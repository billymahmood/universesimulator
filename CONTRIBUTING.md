# Contributing to Universe Simulator

Welcome! This project simulates the universe from real physics — no hardcoded outcomes, every formula traceable to a textbook or measurement. Contributions are very welcome, but they have to keep that bar.

Read [PHILOSOPHY.md](PHILOSOPHY.md) before submitting your first PR. It's short.

## Looking for a task to pick up?

The full prioritised list of open work is in the [**"Where to start" section of the README**](../README.md#where-to-start-full-priority-order-for-contributors) — 15 open issues sorted into three tiers plus a "good first issue" pickup list.

**Quick decisions:**
- **First contribution?** Pick a [`good first issue`](https://github.com/billymahmood/universesimulator/labels/good%20first%20issue) — they're small, well-scoped, and low risk.
- **Want maximum impact?** Pick a [`blocks-phase-6`](https://github.com/billymahmood/universesimulator/labels/blocks-phase-6) issue — they unlock the project's capstone (emergent life).
- **Have a weekend?** Filter by [`effort-small`](https://github.com/billymahmood/universesimulator/labels/effort-small) — work that fits in one focused session.

Before starting: comment on the issue saying you're taking it so we don't duplicate effort.

## Workflow

1. **Fork the repo** (or for direct collaborators, clone it).
2. **Create a feature branch**:
   ```bash
   git checkout -b your-feature-name
   ```
3. **Make your changes** — write tests as you go.
4. **Run the full test suite locally** before pushing:
   ```bash
   .venv/bin/python -m pytest tests/ -q
   ```
   Don't open a PR if tests fail.
5. **Open a PR against `main`**.
6. **The repo owner (Billy Mahmood) reviews and approves** all PRs. Direct push to `main` is restricted to the owner.

## Branch protection

The `main` branch requires:
- All status checks pass (tests)
- At least one approval from the repo owner before merge

Don't try to force-push to `main` — the protection rule blocks it. If your PR has the wrong commits, force-push to your feature branch instead.

## Physics rigor rules

This is **not** an art project. We are explicit about what's real physics and what isn't. Every contribution must respect three rules:

### Rule 1: Every formula needs a textbook citation in code comments

Bad:
```python
F = 6 * c6 / r**7  # van der Waals
```

Good:
```python
# London dispersion (F. London, Trans. Faraday Soc. 33, 8 (1937)):
#   U(r) = −C₆/r⁶,   C₆ = (3/2)·α₁·α₂ · I₁·I₂ / (I₁+I₂)
# Force on i toward j: F = dU/dr = 6·C₆/r⁷
F = 6.0 * c6 / r**7
```

If you can't cite a formula, you probably shouldn't be writing it.

### Rule 2: Every numeric value needs either a source or a calibration label

Bad:
```python
'H': Element(..., 436.0, 13.598, 0.667)  # what are these numbers?
```

Good:
```python
'H': Element(
    ..., 
    bond_dissociation_self=436.0,   # kJ/mol, CRC Handbook 95th ed. §9
    ionization_energy=13.598,        # eV, NIST ASD
    polarisability=0.667,            # Å³, CRC Handbook §10
)
```

If a value is a calibration knob (not real-world data), label it as such in a comment:
```python
# This scale exists only because the sim runs in custom time-scaled units;
# with consistent SI units it would vanish. Calibrated so H+H Coulomb barrier
# ≈ 4000 sim energy.
coulomb_barrier_scale: 9600.0
```

### Rule 3: Tests should validate physics, not just code paths

Bad:
```python
def test_function_returns_a_number():
    assert some_physics_fn(...) is not None
```

Good:
```python
def test_helium_harder_to_ionize_than_hydrogen():
    """At matched KE, H ionises but He doesn't because He's first
    ionization energy (24.587 eV) is 1.8× higher than H's (13.598 eV)."""
    # ... place an H and an He at the same KE
    # ... run ionization update
    assert w.ionized[h_idx]
    assert not w.ionized[he_idx]
```

A physicist reading your test should be able to verify the result against a textbook. "It does what the code says it does" is not a useful test.

## Adding new physics

Common pitfalls when adding a new feature:

- **Don't add a flag or boolean** if the behaviour can come from a continuous physical quantity. Bad: `can_react: bool`. Good: compare activation energy to local kinetic energy.
- **Don't add a hand-tuned constant** if a real measurement exists. Bad: `bond_strength_multiplier = 1.5`. Good: `D_ij = pauling_formula(α_i, α_j, I_i, I_j)` using NIST data.
- **Don't paint visual effects** that don't come from particles in the simulator. We had a fake starfield once. It got removed.

When in doubt, ask in an issue before coding.

## Setting up your dev environment

```bash
# Clone (or fork-clone)
git clone https://github.com/billymahmood/universesimulator.git
cd universesimulator

# Python 3.11+ in a virtual environment
python -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate

# Core dependencies
pip install -r requirements.txt

# Optional but recommended for N > 1000 sims:
pip install numba

# Run the simulator
python main.py

# Run tests
python -m pytest tests/ -q
```

The project currently has **350+ tests**. New features should add tests for the physics they introduce.

## Issue triage / labels

When opening an issue, use these labels where applicable:

- `roadmap` — work explicitly on the project's planned phases (Phase 1–6)
- `multi-week` — scope that genuinely doesn't fit in one session; requires planning
- `needs-architect` — needs design discussion before coding starts
- `dehardcoding` — replacing a hardcoded constant/table/flag with derived physics
- `visual-physics` — photorealistic rendering work (illumination, redshift, lensing, etc.)
- `cross-cutting` — touches many modules
- `performance` — speed / memory work
- `phase-1-thermo` through `phase-6-life` — for tasks tied to a specific roadmap phase

The current roadmap milestones:
- **Universe Roadmap** — the main near-term plan
- **Advanced Physics (multi-week each)** — long-horizon items that each require careful planning

## Commit message style

Follow what's in the existing log (`git log --oneline`). Lower-case type prefix (`feat`, `fix`, `docs`, `dehardcoding`), short imperative subject, optional body explaining *why*.

Good commits:
- `feat: bond order — single / double / triple bonds`
- `dehardcoding: fusion gate from Coulomb barrier (#11)`
- `feat(phase-3): planetary differentiation by composition (closes #3)`

## Code review process

The repo owner reviews every PR. Expect questions about:
- Where the physics comes from (cite your sources)
- Whether the tests validate the actual physics
- Whether the code respects the emergence-over-hardcoding principle

PRs that don't engage with these aren't going to merge. PRs that do — even imperfectly — will get good feedback.

## Questions?

Open a GitHub issue with the `question` label or comment on an existing related issue. Discussions on architecture (especially for `needs-architect` items) are welcome.
