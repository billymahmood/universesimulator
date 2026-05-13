"""Tests for sim/polymers.py — Phase-5 polymer analysis (issue #5).

Coverage groups:

  * Topology — chain length, branching, ring counts on hand-built bond
    graphs whose answers a chemist can verify by inspection.
  * Motifs — amino-acid-like (glycine connectivity) and nucleotide-like
    (pyrimidine-ring + phosphate placeholder) graph fingerprints.
  * Aggregate — polymer_stats HUD output: max / mean chain, motif counts.
  * Chemistry-driven growth — a row of cold carbon atoms within bonding
    range should grow into a chain after running ``chemistry.update``.
    This is the integration smoke test for the acceptance criterion
    "carbon-rich warm regions reliably grow chains of length ≥ 10"
    (issue #5).
  * Breakage — hot carbon chains lose bonds and shrink, satisfying the
    "near-monomer in plasma / cold vacuum" half of the acceptance
    criterion (here represented by the hostile high-KE regime, where the
    bond_break_factor stretches push bonds past dissociation).
"""

from __future__ import annotations
import numpy as np
import pytest

from sim.world import World
from sim.elements import ELEMENTS
from sim.particle import Bond
from sim.polymers import (
    Polymer,
    identify_polymers,
    polymer_of,
    polymer_stats,
)
from sim import chemistry


# ---------------------------------------------------------------------------
# Test helpers — mirror tests/test_molecules.py style
# ---------------------------------------------------------------------------

def _place(world, sym: str, pos, vel=None):
    if vel is None:
        vel = np.zeros(3)
    return world.add_particle(
        ELEMENTS[sym], np.array(pos, dtype=float), np.array(vel, dtype=float),
    )


def _bond(world, i: int, j: int, order: int = 1):
    """Append a Bond between i and j and update bond_counts consistently."""
    e_i = ELEMENTS[world.element_of(i).symbol]
    e_j = ELEMENTS[world.element_of(j).symbol]
    r_eq = e_i.covalent_radius + e_j.covalent_radius
    world.bonds.append(Bond(i, j, r_eq, 50.0, 0.8, order=order))
    world.bond_counts[i] += order
    world.bond_counts[j] += order


def _chain(world, sym: str, length: int, spacing: float = 130.0):
    """Place a straight chain of ``length`` atoms of ``sym`` along +x and
    pre-bond consecutive pairs.  Returns the list of indices."""
    indices = [_place(world, sym, [k * spacing, 0.0, 0.0]) for k in range(length)]
    for k in range(length - 1):
        _bond(world, indices[k], indices[k + 1])
    return indices


# ---------------------------------------------------------------------------
# Topology — chain length / branches / rings
# ---------------------------------------------------------------------------

class TestChainLength:
    """Linear backbones — diameter equals atom count (exact for trees)."""

    def test_dimer_chain_is_2(self, cfg):
        w = World(cfg)
        _chain(w, 'C', 2)
        polys = identify_polymers(w)
        assert len(polys) == 1
        assert polys[0].chain_length == 2

    def test_decamer_chain_is_10(self, cfg):
        """The acceptance criterion in issue #5 calls out chain length ≥ 10."""
        w = World(cfg)
        _chain(w, 'C', 10)
        polys = identify_polymers(w)
        assert len(polys) == 1
        assert polys[0].chain_length == 10
        assert polys[0].branches == 0
        assert polys[0].rings == 0

    def test_long_chain_no_branches(self, cfg):
        w = World(cfg)
        _chain(w, 'C', 25)
        polys = identify_polymers(w)
        assert polys[0].chain_length == 25
        assert polys[0].branches == 0


class TestBranches:
    """Vertices of degree ≥ 3 are branching points (Flory 1953 §III)."""

    def test_tee_shape_one_branch(self, cfg):
        # Main chain of 5 carbons + a side carbon attached to the middle.
        #            C
        #            |
        #   C-C-C-C-C
        w = World(cfg)
        backbone = _chain(w, 'C', 5)
        side = _place(w, 'C', [2 * 130.0, 130.0, 0.0])
        _bond(w, backbone[2], side)
        polys = identify_polymers(w)
        assert polys[0].branches == 1
        assert polys[0].rings == 0
        # Longest shortest path: side -> mid -> end (3 atoms via mid + 2 = 5)
        # The diameter is end-of-backbone -> other-end-of-backbone = 5.
        assert polys[0].chain_length == 5

    def test_star_three_branches(self, cfg):
        # Central C with four arms (degree 4): only the centre is degree ≥3.
        w = World(cfg)
        centre = _place(w, 'C', [0, 0, 0])
        arms = [
            _place(w, 'C', [ 130.0,  0.0, 0.0]),
            _place(w, 'C', [-130.0,  0.0, 0.0]),
            _place(w, 'C', [ 0.0,  130.0, 0.0]),
            _place(w, 'C', [ 0.0, -130.0, 0.0]),
        ]
        for a in arms:
            _bond(w, centre, a)
        polys = identify_polymers(w)
        assert polys[0].branches == 1   # only the centre has degree ≥ 3
        assert polys[0].chain_length == 3   # arm -> centre -> arm


class TestRings:
    """Cyclomatic number r = E − V + 1 per connected component."""

    def test_pure_chain_no_rings(self, cfg):
        w = World(cfg)
        _chain(w, 'C', 6)
        assert identify_polymers(w)[0].rings == 0

    def test_six_membered_carbon_ring(self, cfg):
        # Closed hexagon: 6 atoms, 6 bonds → cyclomatic = 1.
        w = World(cfg)
        ring = [_place(w, 'C', [np.cos(t), np.sin(t), 0.0])
                for t in np.linspace(0.0, 2 * np.pi, 7)[:-1]]
        for k in range(6):
            _bond(w, ring[k], ring[(k + 1) % 6])
        p = identify_polymers(w)[0]
        assert p.rings == 1
        assert 6 in p.ring_sizes

    def test_fused_bicyclic_two_rings(self, cfg):
        # Naphthalene-like: two 6-rings sharing one edge → 10 atoms, 11 bonds,
        # cyclomatic = 11 − 10 + 1 = 2.
        w = World(cfg)
        ring_a = [_place(w, 'C', [k, 0.0, 0.0]) for k in range(6)]
        for k in range(5):
            _bond(w, ring_a[k], ring_a[k + 1])
        _bond(w, ring_a[5], ring_a[0])           # close ring A
        # Ring B shares the edge (ring_a[2], ring_a[3])
        ring_b = [_place(w, 'C', [k, 1.0, 0.0]) for k in range(4)]
        _bond(w, ring_a[3], ring_b[0])
        _bond(w, ring_b[0], ring_b[1])
        _bond(w, ring_b[1], ring_b[2])
        _bond(w, ring_b[2], ring_b[3])
        _bond(w, ring_b[3], ring_a[2])           # close ring B
        p = identify_polymers(w)[0]
        assert p.rings == 2


# ---------------------------------------------------------------------------
# Motifs — graph-only subgraph fingerprints
# ---------------------------------------------------------------------------

class TestAminoAcidLikeMotif:
    """Glycine-like connectivity: H₂N–Cα–C(=O)–OH backbone, abstracted.

    The motif fires when the same component contains an N atom *and* a
    C atom that participates in a C=O double bond (the carbonyl / carboxyl
    fingerprint).  Lehninger 6e §3.1.
    """

    def test_glycine_skeleton_is_amino_acid_like(self, cfg):
        # N–Cα–C(=O)–O   plus the C=O double bond.
        w = World(cfg)
        n  = _place(w, 'N', [   0.0, 0.0, 0.0])
        ca = _place(w, 'C', [ 150.0, 0.0, 0.0])
        cb = _place(w, 'C', [ 300.0, 0.0, 0.0])
        o1 = _place(w, 'O', [ 450.0, 0.0, 0.0])   # =O
        o2 = _place(w, 'O', [ 300.0, 150.0, 0.0]) # –OH oxygen
        _bond(w, n,  ca, order=1)
        _bond(w, ca, cb, order=1)
        _bond(w, cb, o1, order=2)                 # carbonyl
        _bond(w, cb, o2, order=1)
        p = identify_polymers(w)[0]
        assert 'amino-acid-like' in p.motifs

    def test_pure_carbonyl_without_nitrogen_is_not_amino_acid(self, cfg):
        # Formaldehyde H₂C=O — no N, must not be flagged.
        w = World(cfg)
        c  = _place(w, 'C', [0.0,    0.0, 0.0])
        o  = _place(w, 'O', [150.0,  0.0, 0.0])
        h1 = _place(w, 'H', [0.0,  150.0, 0.0])
        h2 = _place(w, 'H', [0.0, -150.0, 0.0])
        _bond(w, c, o, order=2)
        _bond(w, c, h1)
        _bond(w, c, h2)
        assert 'amino-acid-like' not in identify_polymers(w)[0].motifs

    def test_amine_alone_is_not_amino_acid(self, cfg):
        # Methylamine CH₃-NH₂ — N present, but no C=O.  Must not fire.
        w = World(cfg)
        n  = _place(w, 'N', [   0.0,    0.0, 0.0])
        c  = _place(w, 'C', [ 150.0,    0.0, 0.0])
        h1 = _place(w, 'H', [   0.0,  150.0, 0.0])
        h2 = _place(w, 'H', [-150.0,    0.0, 0.0])
        _bond(w, n, c)
        _bond(w, n, h1)
        _bond(w, n, h2)
        assert 'amino-acid-like' not in identify_polymers(w)[0].motifs


class TestNucleotideLikeMotif:
    """Pyrimidine / purine / imidazole / pyrrole fingerprint + a P atom.

    The strict rule (sim/polymers.py:_detect_motifs) requires the C+N
    pair to appear *inside* a 5- or 6-membered ring — checking the
    surrounding component would mark every benzene-with-external-amine
    as nucleotide-like, which is wrong.
    """

    def test_pyrimidine_with_phosphate_is_nucleotide_like(self, cfg):
        # 6-ring of alternating C / N + an external P attached to one C.
        w = World(cfg)
        ring = []
        for k in range(6):
            sym = 'C' if k % 2 == 0 else 'N'
            ring.append(_place(w, sym, [k * 100.0, 0.0, 0.0]))
        for k in range(6):
            _bond(w, ring[k], ring[(k + 1) % 6])
        p_idx = _place(w, 'P', [0.0, 200.0, 0.0])
        _bond(w, ring[0], p_idx)
        p = identify_polymers(w)[0]
        assert 'nucleotide-like' in p.motifs
        assert 6 in p.ring_sizes

    def test_imidazole_with_external_phosphate_is_nucleotide_like(self, cfg):
        """5-ring of C / N (imidazole / histidine side-chain scaffold —
        Lehninger 6e §3.2) plus an external P atom → must fire.  Covers
        the 5-membered branch of the motif, distinct from the 6-ring
        pyrimidine case above."""
        w = World(cfg)
        ring = []
        for k, sym in enumerate(['C', 'N', 'C', 'N', 'C']):
            ring.append(_place(w, sym, [k * 100.0, 0.0, 0.0]))
        for k in range(5):
            _bond(w, ring[k], ring[(k + 1) % 5])
        p_idx = _place(w, 'P', [0.0, 200.0, 0.0])
        _bond(w, ring[0], p_idx)
        p = identify_polymers(w)[0]
        assert 'nucleotide-like' in p.motifs
        assert 5 in p.ring_sizes

    def test_ring_without_phosphorus_is_not_nucleotide_like(self, cfg):
        # Pyrimidine-style C/N ring but no P → just an aromatic, not a
        # nucleotide fingerprint.
        w = World(cfg)
        ring = []
        for k in range(6):
            sym = 'C' if k % 2 == 0 else 'N'
            ring.append(_place(w, sym, [k * 100.0, 0.0, 0.0]))
        for k in range(6):
            _bond(w, ring[k], ring[(k + 1) % 6])
        assert 'nucleotide-like' not in identify_polymers(w)[0].motifs

    def test_chain_with_phosphorus_no_ring_is_not_nucleotide_like(self, cfg):
        # Linear chain with a P somewhere — no ring, so not a nucleotide.
        w = World(cfg)
        atoms = [_place(w, 'C', [k * 130.0, 0.0, 0.0]) for k in range(4)]
        atoms.append(_place(w, 'P', [4 * 130.0, 0.0, 0.0]))
        for k in range(4):
            _bond(w, atoms[k], atoms[k + 1])
        assert 'nucleotide-like' not in identify_polymers(w)[0].motifs

    def test_pure_carbon_ring_with_external_amine_phosphate_is_not_nucleotide(self, cfg):
        """Regression: a 6-ring of pure carbon (benzene scaffold) with an
        external –NH–P substituent. The component contains C, N, *and* P,
        but the 6-ring itself is all C — there is no N inside the ring,
        so the motif must NOT fire.

        Earlier surrogate implementations checked the component as a
        whole and incorrectly flagged this as nucleotide-like. The strict
        check (motif requires C+N in the ring vertex set itself) is the
        whole point of the refactor.
        """
        w = World(cfg)
        ring = [_place(w, 'C', [100.0 * np.cos(t), 100.0 * np.sin(t), 0.0])
                for t in np.linspace(0.0, 2 * np.pi, 7)[:-1]]
        for k in range(6):
            _bond(w, ring[k], ring[(k + 1) % 6])
        n_atom = _place(w, 'N', [0.0,  200.0, 0.0])
        p_atom = _place(w, 'P', [0.0,  400.0, 0.0])
        _bond(w, ring[0], n_atom)
        _bond(w, n_atom,  p_atom)
        p = identify_polymers(w)[0]
        assert 'nucleotide-like' not in p.motifs
        # Sanity: the 6-ring is still detected; the motif just doesn't
        # claim it.
        assert 6 in p.ring_sizes


# ---------------------------------------------------------------------------
# Aggregate / HUD stats
# ---------------------------------------------------------------------------

class TestPolymerStats:
    def test_empty_world(self, cfg):
        w = World(cfg)
        s = polymer_stats(w)
        assert s == {
            'n_polymers':   0,
            'max_chain':    0,
            'mean_chain':   0.0,
            'motif_counts': {},
        }

    def test_only_unbonded_atoms(self, cfg):
        w = World(cfg)
        for k in range(5):
            _place(w, 'H', [k * 500.0, 0.0, 0.0])
        s = polymer_stats(w)
        assert s['n_polymers'] == 0

    def test_max_and_mean_chain(self, cfg):
        # Build three components: chains of lengths 3, 5, 10.
        w = World(cfg)
        _chain(w, 'C', 3,  spacing=130.0)
        # Offset the second chain so it doesn't touch the first
        for j in range(5):
            _place(w, 'C', [j * 130.0, 1000.0, 0.0])
        base = w.n
        for k in range(4):
            _bond(w, base - 5 + k, base - 5 + k + 1)
        # Third chain
        for j in range(10):
            _place(w, 'C', [j * 130.0, 2000.0, 0.0])
        base = w.n
        for k in range(9):
            _bond(w, base - 10 + k, base - 10 + k + 1)
        s = polymer_stats(w)
        assert s['n_polymers'] == 3
        assert s['max_chain']  == 10
        assert s['mean_chain'] == pytest.approx((3 + 5 + 10) / 3.0)

    def test_motif_counts_aggregate(self, cfg):
        """Two glycine-like clusters and one neutral chain → motif count 2."""
        w = World(cfg)
        # Two amino-acid skeletons at different y offsets
        for y_off in (0.0, 1000.0):
            n  = _place(w, 'N', [   0.0, y_off,        0.0])
            ca = _place(w, 'C', [ 150.0, y_off,        0.0])
            cb = _place(w, 'C', [ 300.0, y_off,        0.0])
            o1 = _place(w, 'O', [ 450.0, y_off,        0.0])
            o2 = _place(w, 'O', [ 300.0, y_off + 150., 0.0])
            _bond(w, n,  ca); _bond(w, ca, cb)
            _bond(w, cb, o1, order=2); _bond(w, cb, o2)
        # A separate neutral hexane-like backbone (no N or P)
        atoms = [_place(w, 'C', [k * 130.0, 2000.0, 0.0]) for k in range(6)]
        for k in range(5):
            _bond(w, atoms[k], atoms[k + 1])
        s = polymer_stats(w)
        assert s['motif_counts'].get('amino-acid-like') == 2


# ---------------------------------------------------------------------------
# polymer_of
# ---------------------------------------------------------------------------

class TestPolymerOf:
    def test_returns_containing_polymer(self, cfg):
        w = World(cfg)
        indices = _chain(w, 'C', 5)
        p = polymer_of(w, indices[2])
        assert p is not None
        assert p.chain_length == 5
        assert set(p.indices) == set(indices)

    def test_returns_none_for_lone_atom(self, cfg):
        w = World(cfg)
        lone = _place(w, 'C', [0, 0, 0])
        _place(w, 'C', [1000.0, 0.0, 0.0])    # other atom, unrelated
        assert polymer_of(w, lone) is None


# ---------------------------------------------------------------------------
# Chemistry-driven chain growth (acceptance criterion smoke test)
# ---------------------------------------------------------------------------

class TestChainGrowthFromChemistry:
    """A row of cold carbon atoms within bond-formation range must
    spontaneously bond into a chain via ``chemistry.update``.

    Real carbon parameters: bond_dissociation_self = 348 kJ/mol (CRC), χ =
    2.55 (Pauling). The Pauling formula gives D(C-C) = 348 kJ/mol exactly
    (Δχ = 0). With the test cfg's slow velocity threshold of 1000 SU and
    zero initial velocities, every adjacent pair is below the threshold,
    so consecutive bonds form on the first call.
    """

    def test_cold_carbon_row_forms_chain(self, cfg):
        w = World(cfg)
        # Spacing < bond_formation_factor × 2 × rc_C
        # = 1.4 × 154 ≈ 215 pm; we use 130 pm so each pair is well inside.
        n = 12
        for k in range(n):
            _place(w, 'C', [k * 130.0, 0.0, 0.0])
        chemistry.update(w)
        # The chain should reach the full row, length 12.
        polys = identify_polymers(w)
        # All n atoms in one component
        assert len(polys) == 1
        assert polys[0].size == n
        # Chain length must clear the issue-#5 ≥ 10 bar.
        assert polys[0].chain_length >= 10

    def test_hot_carbon_row_does_not_form_chain(self, cfg):
        """Same geometry but with relative velocities above the
        bond-velocity threshold → no bonds. The kinetic regime
        determines bond formation, not a flag."""
        w = World(cfg)
        v_thresh = cfg.chemistry.bond_velocity_threshold
        for k in range(8):
            # Alternating velocities ±1.5×v_thresh in x → relative
            # velocity between neighbours is 3×v_thresh > threshold.
            sign = 1.0 if k % 2 == 0 else -1.0
            _place(w, 'C', [k * 130.0, 0.0, 0.0], vel=[sign * 1.5 * v_thresh, 0.0, 0.0])
        chemistry.update(w)
        polys = identify_polymers(w)
        # No long chain should form — at most diatomic pairs sneak through
        # at the extremes, and even those have neighbours moving too fast.
        max_chain = max((p.chain_length for p in polys), default=0)
        assert max_chain < 5


# ---------------------------------------------------------------------------
# Chain breakage under hostile conditions
# ---------------------------------------------------------------------------

class TestChainBreakageUnderStress:
    """The acceptance criterion requires that polymer length distribution is
    near-monomer in plasma / cold vacuum. We simulate that here by
    pre-building a chain and then letting bond-break_factor × eq_length
    cull bonds when atoms get stretched apart."""

    def test_stretched_chain_loses_bonds(self, cfg):
        w = World(cfg)
        # Build a 6-atom carbon chain at normal spacing, then displace
        # alternating atoms far enough that every bond exceeds
        # bond_break_factor × eq_length and breaks on the next chemistry
        # step.
        atoms = _chain(w, 'C', 6, spacing=130.0)
        rc_c = ELEMENTS['C'].covalent_radius
        break_dist = cfg.chemistry.bond_break_factor * 2 * rc_c * 1.5  # well past
        # Push every other atom out in +y so all 5 bonds stretch
        for idx, a in enumerate(atoms):
            if idx % 2 == 1:
                w.positions[a, 1] = break_dist

        n_bonds_before = len(w.bonds)
        chemistry._break_bonds(w)
        n_bonds_after = len(w.bonds)
        assert n_bonds_after < n_bonds_before
        # Final state should have shorter chains — definitely below 6.
        max_chain = max(
            (p.chain_length for p in identify_polymers(w)),
            default=0,
        )
        assert max_chain < 6


# ---------------------------------------------------------------------------
# Full-simulation acceptance — issue #5 criterion exercised through the
# integrated step loop (gravity, thermal pressure, vdW, chemistry,
# reactions). The earlier ``TestChainGrowthFromChemistry`` exercised
# ``chemistry.update`` in isolation; here we run the same physics the
# simulator runs at runtime, end to end.
# ---------------------------------------------------------------------------

class TestFullSimulationAcceptance:
    """Acceptance criterion from issue #5 verified against ``world.step``:

        "Carbon-rich warm regions reliably grow chains of length ≥ 10."
        "Polymer length distribution shows a long tail in good niches,
         vs near-monomer in plasma/cold vacuum."

    The two tests below cover the two halves: positive (cold carbon
    cluster, full step loop → chain ≥ 10) and negative (hot kinetic
    regime → near-monomeric).
    """

    @staticmethod
    def _make_closed_carbon_world(cfg, count: int, spacing: float,
                                  velocity_fn=None):
        """Build a closed-system World with `count` carbon atoms in a row.

        Closes the system by zeroing the injection rate so the test
        doesn't accumulate stray H/He from cfg.injection. Disables
        fusion (focuses the test on bond chemistry — fusion is
        exercised by tests/test_chemistry.py) and accretion (gravita-
        tional collapse into rocky/stellar bodies is a different
        physical regime; for a "carbon-rich warm region" the relevant
        process is covalent bonding, not coagulation — disabling it
        keeps the AC test from being shadowed by accretion eating the
        atoms before chemistry can bond them).
        """
        cfg.simulation.injection_rate = 0
        cfg.chemistry.fusion_enabled = False
        cfg.accretion.accretion_radius = 0.0
        w = World(cfg)
        for k in range(count):
            vel = velocity_fn(k) if velocity_fn is not None else np.zeros(3)
            w.add_particle(
                ELEMENTS['C'],
                np.array([k * spacing, 0.0, 0.0], dtype=float),
                np.asarray(vel, dtype=float),
            )
        return w

    def test_cold_carbon_cluster_grows_chain_through_full_step_loop(self, cfg):
        """A row of cold carbons subjected to the full physics stack
        (gravity + thermal pressure + vdW + chemistry + reactions) must
        develop a chain of ≥ 10 atoms and keep it stable across
        multiple integration steps.

        With zero initial velocity, ``add_thermal_pressure`` is exactly
        zero on the first step (the force law scales with relative KE,
        sim/thermal.py:67), so the first chemistry.update sees the row
        with bond_velocity_threshold cleanly satisfied and bonds the
        consecutive carbons. The Morse well then keeps them together
        through the velocity-Verlet integration that follows.
        """
        w = self._make_closed_carbon_world(cfg, count=15, spacing=130.0)

        for _ in range(30):
            w.step(0.01)

        polys = identify_polymers(w)
        max_chain = max((p.chain_length for p in polys), default=0)
        assert max_chain >= 10, (
            f'expected chain length >= 10 after 30 full-sim steps, '
            f'got {max_chain}; polymers={[(p.size, p.chain_length) for p in polys]}'
        )

    def test_hot_carbon_cluster_stays_near_monomeric(self, cfg):
        """The negative half of the AC: the same carbon row launched
        with random velocities well above ``bond_velocity_threshold`` in
        every direction must NOT grow a meaningful chain. The relative
        velocity between any pair exceeds the bond-formation kinetic
        gate, so bonds rarely form and break quickly when they do.
        """
        rng = np.random.default_rng(20260512)
        v_thresh = cfg.chemistry.bond_velocity_threshold

        def hot_velocity(_k: int) -> np.ndarray:
            # 3× threshold in a random direction → pair-relative speed
            # well above the bond-formation gate.
            return rng.standard_normal(3) * 3.0 * v_thresh

        w = self._make_closed_carbon_world(
            cfg, count=15, spacing=200.0, velocity_fn=hot_velocity,
        )

        for _ in range(20):
            w.step(0.01)

        polys = identify_polymers(w)
        max_chain = max((p.chain_length for p in polys), default=0)
        # A handful of brief diatomic encounters are physically possible,
        # but no chain that would clear the AC bar.
        assert max_chain < 5, (
            f'hot regime should stay near-monomeric, got chain={max_chain}'
        )

    def test_pre_built_chain_disintegrates_under_high_thermal_kick(self, cfg):
        """A 15-atom pre-built chain injected with enough thermal energy
        that pair-relative speeds exceed the bond_break_factor stretch
        within a few steps loses bonds. Verifies the simulator's chain
        breakage behaviour against the AC's "near-monomer in plasma"
        clause."""
        cfg.simulation.injection_rate = 0
        cfg.chemistry.fusion_enabled = False
        cfg.accretion.accretion_radius = 0.0   # focus on bond dynamics

        w = World(cfg)
        atoms = _chain(w, 'C', 15, spacing=130.0)
        assert identify_polymers(w)[0].chain_length == 15

        # Inject large outward radial velocities — pulls atoms apart and
        # stretches bonds past the dissociation threshold.
        rng = np.random.default_rng(20260513)
        v_thresh = cfg.chemistry.bond_velocity_threshold
        for a in atoms:
            w.velocities[a] = rng.standard_normal(3) * 5.0 * v_thresh

        for _ in range(20):
            w.step(0.01)

        polys = identify_polymers(w)
        max_chain = max((p.chain_length for p in polys), default=0)
        assert max_chain < 15, (
            f'high-kinetic stress should break the chain, got {max_chain}'
        )


# ---------------------------------------------------------------------------
# Dataclass invariants
# ---------------------------------------------------------------------------

class TestPolymerDataclass:
    def test_size_matches_indices(self, cfg):
        w = World(cfg)
        _chain(w, 'C', 4)
        p = identify_polymers(w)[0]
        assert isinstance(p, Polymer)
        assert p.size == len(p.indices) == 4

    def test_indices_are_sorted(self, cfg):
        w = World(cfg)
        _chain(w, 'C', 5)
        p = identify_polymers(w)[0]
        assert list(p.indices) == sorted(p.indices)
