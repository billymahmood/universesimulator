# Universe Simulator

A 3D real-time particle physics simulation that grows a universe from nothing — empty space, atoms injected at cosmic abundances, and **every observable outcome emerging from real physics**. No hardcoded galaxy shapes, no scripted star formation, no fake nebulae painted on a backdrop. Every numerical value cites NIST, AME 2020, CRC Handbook, or a named formula.

**Status:** 11+ physics modules, 360+ tests, ~6 phases worth of emergent behaviour. Real radiation transport. Real chemistry. Real stellar nucleosynthesis ending at iron. See [PHILOSOPHY.md](PHILOSOPHY.md) for the project principle.

```bash
git clone https://github.com/billymahmood/universesimulator.git
cd universesimulator
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python main.py
```

## What this simulator does (all real physics, all emergent)

| Mechanism | What it is | Reference |
|---|---|---|
| **N-body gravity** | Symmetric O(N²/2), Newton's 3rd law, optional numba JIT 30–50× faster | Newton |
| **Primordial structure** | Hubble expansion + density seeds → tidal-torque rotation, no hardcoded spin | Standard cosmology |
| **Atomic data** | NIST IE, AME 2020 isotope masses, CRC bond energies, Cordero covalent radii, CRC polarisabilities | All cited per-row in `sim/elements.py` |
| **Covalent bonding** | Pauling: `D(A-B) = √(D_AA·D_BB) + 96·Δχ²` (kJ/mol) with bond order single/double/triple | Pauling 1932; Atkins |
| **van der Waals** | London dispersion `F = 6·C₆/r⁷` with `C₆ = (3/2)·α·α'·I·I'/(I+I')` from real polarisability data | London 1937 |
| **Heterogeneous catalysis** | Bond formation accelerated near solid surfaces | Sabatier 1911, Langmuir 1918 |
| **Chemical reactions** | Electronegativity-driven single-bond swaps, real Arrhenius kinetics, ΔE released as KE | Pauling + Eyring |
| **Nuclear fusion** | Coulomb barrier `V_C = scale·Z₁Z₂/r_nuc` (Krane), Gamow quantum tunnelling (Clayton), Q-values from AME 2020 masses, products from Z+A search | Burbidge×4 + Clayton |
| **Outgassing** | Hot bodies emit light volatiles at escape velocity (real Jeans escape) | – |
| **Accretion** | Bound non-fusing pairs merge, composition history tracked element-by-element | – |
| **Differentiation** | Body colour = mass-weighted CPK of constituents; heavy core, light atmosphere halo | – |
| **States of matter** | Solid / liquid / gas / plasma from local density + KE + bonding + ionisation | – |
| **Ionisation** | Per-element NIST IE thresholds with hysteresis on recombination | NIST ASD |
| **Stefan-Boltzmann radiation** | Halo size ∝ T², Planck blackbody RGB (Tanner Helland's CIE fit) | Stefan, Planck |
| **Inverse-square illumination** | Gas lit by stars; reflection-nebula colouring from each emitter's Planck colour | Rybicki & Lightman |
| **Beer-Lambert extinction** | Cold gas absorbs light along line of sight (dark nebulae) | Rybicki & Lightman §1 |
| **Doppler + Hubble redshift** | Apparent T = T / (1 + v_radial/c) on rendered colour | Wien's law |
| **Wavelength-dependent reddening** | Per-RGB extinction with Fitzpatrick 1999 ISM weights | Fitzpatrick 1999 |
| **HDR tone mapping** | Reinhard exposure curve so dim sources stay visible | Real CCD response |
| **Explicit photons** | Stochastic thermal emission at rate ∝ T⁴, finite-c straight-line propagation, swept-volume absorption, real radiation pressure | Stefan-Boltzmann + Rybicki §1 |
| **Polymer analysis** | Chain length via BFS diameter, branching count, ring detection from cyclomatic complexity, motif tagging (amino-acid-like, nucleotide-like) on the bond graph | Aho-Hopcroft-Ullman; Diestel §1.9; Flory 1953; Lehninger 6e |

Every halo around every particle is the result of real Stefan-Boltzmann emission scaling, coloured by the real Planck-curve. There is no painted starfield — the void is empty.

## Roadmap

The project is organised into two milestones on GitHub:

- **[Universe Roadmap](https://github.com/billymahmood/universesimulator/milestone/1)** — Phases 1–6 of emergent behaviour (atoms → states of matter → molecules → planets → chemistry → polymers → life), plus cross-cutting concerns (save/load, presets, time tiers).
- **[Advanced Physics (multi-week each)](https://github.com/billymahmood/universesimulator/milestone/2)** — long-horizon items each requiring careful design: explicit electrons, multi-bond reactions, stereochemistry, gravitational lensing, volumetric ray-marching, wavelength tracking, scattering.

Progress so far:

| | |
|---|---|
| Phase 1 — States of matter | ✅ |
| Phase 2 — Molecule recognition | ✅ |
| Phase 3 — Planetary differentiation | ✅ |
| Phase 4 — Chemical reactions | ✅ |
| Phase 5 — Polymers | ✅ |
| Phase 6 — Emergent life | ⬜ blocked on dependency chain below |
| Photorealistic rendering pipeline | 6/8 sub-tasks done (issue #12) |
| De-hardcoding pass | 6 of 11 items done (issue #11) |

### Path to Phase 6 — Emergent life

Phase 6 is the project's capstone (autocatalytic chemistry → template replication → selection pressure). For it to be **truly emergent** rather than hand-coded "Avida-style" artificial life, four prerequisites must land first, in a specific order so we never have to rewrite earlier work:

| Order | Issue | What it is | Effort |
|:-:|---|---|---|
| 1 | [#13 Explicit electrons](https://github.com/billymahmood/universesimulator/issues/13) | Replace implicit `max_bonds` with real electron orbitals. Foundational chemistry rewrite — bond formation, ionisation, photon emission, spectral lines all become electron-driven. | 4–8 weeks |
| 2 | [#14 Multi-bond reactions](https://github.com/billymahmood/universesimulator/issues/14) | Real metabolism breaks/forms several bonds at once (`2H₂ + O₂ → 2H₂O`). Required for autocatalytic network closure. | 3–4 weeks |
| 3 | [#15 Stereochemistry / bond angles](https://github.com/billymahmood/universesimulator/issues/15) | 3D molecular shape — needed for template-based replication. With #13 in place, hybridisation emerges from orbital geometry. | 3–4 weeks |
| 4 | [#23 Lineage tracking](https://github.com/billymahmood/universesimulator/issues/23) | Per-particle `parent_id` so replication events can be traced through generations. Independent of the chemistry stack — can land in parallel. | ~1 week |
| 5 | [#6 Phase 6 — Emergent life](https://github.com/billymahmood/universesimulator/issues/6) | Autocatalysis (6a) → template replication (6b) → selection pressure (6c) → HUD + lineage viz (6d). All four sub-phases emergent from chemistry, no hardcoded "alive" tag. | ~3 weeks |

**Total time to Phase 6:** ~10–17 weeks of prerequisite work (depending on contributor parallelism) + ~3 weeks of life work = ~13–20 weeks.

Why this order matters: doing #14 or #15 *before* #13 would mean rewriting them once electrons land. Doing #6 *before* the prerequisites either forces us to violate [PHILOSOPHY.md](PHILOSOPHY.md) by hardcoding life-like behaviour, or it produces something that doesn't quite work.

If you want to contribute toward Phase 6: pick whichever of #13 / #14 / #15 / #23 hasn't been claimed yet, comment on the issue that you're taking it, branch + PR per [CONTRIBUTING.md](CONTRIBUTING.md).

### Where to start (full priority order for contributors)

All 15 open issues, ranked. Pick from the top tier if it's available; drop down a tier if those are claimed or if the scope/effort doesn't match what you can commit.

**Tier 1 — Critical path to Phase 6 (highest leverage on the project's main goal)**

Doing these in this order means the project's capstone (emergent life) becomes possible. Tagged [`blocks-phase-6`](https://github.com/billymahmood/universesimulator/labels/blocks-phase-6) on GitHub.

| Order | Issue | Effort | Status |
|:-:|---|---|---|
| 1 | [#13 Explicit electrons](https://github.com/billymahmood/universesimulator/issues/13) | 4–8 weeks | Foundational chemistry rewrite. Biggest task on the board. |
| 2 | [#14 Multi-bond reactions](https://github.com/billymahmood/universesimulator/issues/14) | 3–4 weeks | After #13 lands. Unlocks real metabolism. |
| 3 | [#15 Stereochemistry / bond angles](https://github.com/billymahmood/universesimulator/issues/15) | 3–4 weeks | After #14. Needed for template-based replication. |
| 4 | [#23 Lineage tracking](https://github.com/billymahmood/universesimulator/issues/23) | ~1 week | Can land **in parallel** with any of the above. Good first issue. |
| 5 | [#6 Phase 6 — Emergent life](https://github.com/billymahmood/universesimulator/issues/6) | ~3 weeks | After all four prerequisites. |

**Tier 2 — Useful infrastructure (parallel work, easy to pick up)**

Doesn't block Phase 6 but improves the project. Tagged [`effort-small`](https://github.com/billymahmood/universesimulator/labels/effort-small) on GitHub.

| Issue | Effort | Why it's useful |
|---|---|---|
| [#8 Save / load world state](https://github.com/billymahmood/universesimulator/issues/8) | ~1 week | Long Phase 6 runs survive restarts. Good first issue. |
| [#9 Scenario presets](https://github.com/billymahmood/universesimulator/issues/9) | ~1 week | Protoplanetary disk, primordial soup — useful for Phase 6 testing. |
| [#18 Wavelength tracking on photons](https://github.com/billymahmood/universesimulator/issues/18) | 3–4 weeks | Extends the photon module — opens up spectral lines. |
| [#20 Photon-photon interactions](https://github.com/billymahmood/universesimulator/issues/20) | ~1 hour | Mostly a docstring no-op (cross-sections are unobservably small in our regime). Good for someone wanting a tiny first commit. |

**Tier 3 — Beautiful but optional (luxury items)**

Pick if you have a specific demo or visual you want to see.

| Issue | Effort | What it gets you |
|---|---|---|
| [#16 Gravitational lensing](https://github.com/billymahmood/universesimulator/issues/16) | 3–4 weeks | Einstein rings, real GR light bending. |
| [#17 Volumetric ray-marching](https://github.com/billymahmood/universesimulator/issues/17) | 4–8 weeks | Cosmology-quality gas rendering. |
| [#19 Compton + Rayleigh scattering](https://github.com/billymahmood/universesimulator/issues/19) | 3–5 weeks | Explains why the sky is blue; needs #18 first. |
| [#21 Wavelength-dependent cross-sections](https://github.com/billymahmood/universesimulator/issues/21) | 3–5 weeks | Real ISM reddening, Strömgren spheres; needs #18 first. |
| [#10 Time tiers](https://github.com/billymahmood/universesimulator/issues/10) | 1–2 weeks | Fast-forward through quiet periods of the sim. |
| [#12 Photorealistic rendering umbrella](https://github.com/billymahmood/universesimulator/issues/12) | – | Tracks #18-21; 6/8 sub-tasks already done. |

**Good first issues (small, clear scope, low risk)**

If you're new to the codebase, start with one of these — they have well-defined scope and won't require deep architectural changes. Tagged [`good first issue`](https://github.com/billymahmood/universesimulator/labels/good%20first%20issue) on GitHub.

- [#23 Lineage tracking](https://github.com/billymahmood/universesimulator/issues/23) — copies the existing `composition` pattern; pure addition
- [#8 Save / load](https://github.com/billymahmood/universesimulator/issues/8) — file I/O, no physics changes
- [#20 Photon-photon](https://github.com/billymahmood/universesimulator/issues/20) — literally a docstring + reference

**How to claim a task**

1. Comment on the issue saying you're taking it (so others don't duplicate effort)
2. Branch + PR per [CONTRIBUTING.md](CONTRIBUTING.md)
3. The repo owner reviews and approves

## Controls

| Key / Mouse | Action |
|---|---|
| Left-drag / Arrow keys | Orbit camera |
| Right-drag / Scroll | Zoom |
| Page Up / Page Down | Zoom in / out |
| Middle-drag | Pan |
| `Space` | Pause / resume |
| `+` / `−` | Speed up / slow down |
| `F` | Cinematic camera (auto orbit) |
| `P` | Toggle astronomy palette (Planck colour on cores) |
| `R` | Reset camera |
| `Q` / `Escape` | Quit |
| Left-click | Inspect a particle (composition, bonds, molecule, kind) |

The window title HUD shows: particle count, bonds, fusions, reactions, accretions, outgassing, photon counts, state distribution (S/L/G/P), total energy + drift %, velocity-clamp count, simulation time, FPS, time-scale multiplier, top elements present, top molecules formed.

## Configuration

All parameters live in [`settings.yaml`](settings.yaml) and are documented inline. The major sections:

- `simulation` — particle counts, injection rate, time step
- `physics` — gravity, softening, max velocity
- `thermal` — pressure, VdW (London formula), ionisation, plasma factor
- `chemistry` — bond formation/breaking, Coulomb barrier, Gamow tunnelling, Pauling scale, catalysis
- `reactions` — Arrhenius kinetics for bond swaps
- `outgassing` — Jeans-escape volatile release
- `accretion` — coagulation radius
- `cosmology` — Hubble parameter, primordial seeds
- `injection` — cosmic element abundances (H 74%, He 24%, traces)
- `renderer` — illumination, extinction, Doppler/Hubble redshift, exposure
- `photons` — explicit radiation transport
- `physics` calibration scales — see [PHILOSOPHY.md](PHILOSOPHY.md) for which constants are derived vs measured vs calibrated

## Performance

Pure NumPy gravity runs at N=1000 on most modern hardware. For larger sims:

- **N=1k–10k:** install [numba](https://numba.pydata.org/) — gravity auto-detects it and switches to a parallel JIT loop. No code changes needed.
  ```bash
  pip install numba
  ```
  Measured speedup (Apple Silicon):

  | N | NumPy | numba | speedup |
  |---:|---:|---:|---:|
  | 100 | 0.73 ms | 0.11 ms | 6.3× |
  | 500 | 6.43 ms | 0.19 ms | 33.8× |
  | 1000 | 19.6 ms | 0.57 ms | 34.2× |
  | 2000 | 70.4 ms | 1.35 ms | 52.0× |

  Both backends produce numerically identical force fields (verified to 1e-10 in `tests/test_physics_gravity.py`).
- **N=10k+:** even with numba, O(N²) is the bottleneck. A Barnes–Hut tree (O(N log N)) is on the roadmap. Contributions welcome.

## Architecture

```
main.py              — entry point, config loading
settings.yaml        — all tunable parameters (inline-documented)
PHILOSOPHY.md        — emergence-over-hardcoding principle
CONTRIBUTING.md      — workflow, physics rigor, citation rules
sim/
  elements.py        — periodic table (NIST + AME + CRC + Pauling)
  nuclear.py         — Coulomb barrier, Gamow tunnelling, Q-values, product search
  particle.py        — Bond class (Morse potential, bond order)
  spatial.py         — O(1) neighbour-lookup grid (shared)
  physics.py         — gravity (NumPy + optional numba) + Morse bond forces
  thermal.py         — kinetic-pressure analogue + plasma repulsion
  vdw.py             — London dispersion (real polarisability data)
  chemistry.py       — bond formation/breaking, fusion, radiation kicks, catalysis
  reactions.py       — single-bond electronegativity-driven swaps
  accretion.py       — coagulation of bound non-fusing pairs
  outgassing.py      — Jeans-escape volatile release
  photons.py         — explicit radiation transport (emit + propagate + absorb)
  states.py          — solid / liquid / gas / plasma classifier
  molecules.py       — connected-component identification, Hill formulae
  illumination.py    — inverse-square flux from emitters to receivers
  extinction.py      — Beer-Lambert line-of-sight absorption
  polymers.py        — chain length / branching / rings / motif tagging
  diagnostics.py     — total KE, PE, momentum, drift %
  injector.py        — particle birth (position, velocity, element sampling)
  world.py           — state, velocity-Verlet integrator, full step loop
  viewer.py          — vispy 3D renderer + HDR tone mapping + Doppler redshift + click-inspect
tests/               — 360+ tests covering every module above
```

## Contributing

**Pull requests are welcome.** Read [CONTRIBUTING.md](CONTRIBUTING.md) first — it covers:

1. The branch + PR workflow (no direct push to `main`)
2. Physics-rigor rules (every formula needs a citation; every constant needs a source or a calibration label; tests should validate physics, not just code paths)
3. The de-hardcoding philosophy

The repo owner reviews and approves every PR. Issues for proposed work live in the two milestones above; pick one, comment that you're working on it, branch and PR.

If you want to discuss architecture before coding (especially for `multi-week` / `needs-architect` items), open an issue with the `question` label.

## License

MIT — see [LICENSE](LICENSE).
