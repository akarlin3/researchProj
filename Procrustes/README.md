# Procrustes

*Misspecification-aliasing of a calibrated error bar* — fitting a bi-exponential
IVIM model on non-bi-exponential truth keeps **marginal** coverage but breaks the
**conditional** coverage of the *well-identified* tissue-diffusion map D, distinct
from Gauge's within-model identifiability wall.

- [`procrustes-core/`](procrustes-core/) — the clean-room core: hardened GATE
  B/C/D results and the observable misspecification diagnostic. Ground truth is
  the Lattice DRO (seed-generated, no data files).
- [`procrustes-core/RESULTS.md`](procrustes-core/RESULTS.md) — the hardened,
  load-bearing numbers (GATE B/C/D).
- [`procrustes-core/POSITIONING.md`](procrustes-core/POSITIONING.md) — the novelty
  gate record (vs Gauge, Lei 2018, Barber 2021, Wang–Tamir–Bush 2026, IVIM model
  selection).

Status: **submission-ready but HELD** — GATE B (apples-to-apples Gauge
separation), GATE C (diagnostic reach), and GATE D (robustness envelope) all
PASS, and the manuscript compiles (`procrustes-core/paper/procrustes.pdf`).
Submission is withheld pending Gauge's publication
(`procrustes-core/release_gate.py`, `procrustes-core/SUBMISSION_HOLD`); run
`bash procrustes-core/reproduce.sh` to re-validate everything in one command.
Venue TBC.
