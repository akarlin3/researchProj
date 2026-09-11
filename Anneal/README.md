# Anneal — chimera collapse, aging, and finite-size scaling

*Paper:* **"Chimera Collapse Ages: Topology-Dependent Finite-Size Scaling in
Mean-Field and Ring Oscillator Systems"** (submitted to *Nonlinear Dynamics*).

Anneal studies *chimera death* — the spontaneous collapse of the coexisting
synchronized/desynchronized state in coupled identical oscillators — and asks
which features of that collapse are universal versus topology-dependent. The
work grew out of a real-time additive music synthesizer whose two-population
Sakaguchi–Kuramoto engine drives partial amplitudes (hence the origin repo name
`annealMusic`), which prompted the operational question of what a "chimera
death" detector actually measures. The study runs matched long-duration
survival experiments on two canonical substrates — the two-population mean
field and the nonlocally coupled ring near its existence boundary — under a
pre-registered protocol, and validates against a parameter-free reduced
two-population flow.

## Headline findings

- The standard first-passage death criterion undercounts true lifetimes by
  ~50% (up to 98% of near-synchrony excursions are self-healing "grazes").
- Collapse *ages* at the hazard level in both topologies (Weibull shape
  k > 1 in all 15 pre-registered cells) and arrives via a breath-synchronized
  ratchet.
- Finite-size scaling is topology-dependent: mean-field lifetimes plateau
  around 139 s while near-boundary ring lifetimes *decrease* with N.
- The reduced flow predicts the breath period and per-cycle dynamics but
  underpredicts finite-N lifetimes by a constant 3.2×, leaving the
  prolongation mechanism as a sharply constrained open problem.

## Layout

- [`anneal-hazard/`](anneal-hazard/) — the main pre-registered ring hazard
  experiment (numba engine, survival/hazard library, CP1–CP4 validation
  gates, [`PREREGISTRATION.md`](anneal-hazard/PREREGISTRATION.md),
  [`SUMMARY.md`](anneal-hazard/SUMMARY.md)).
- [`tools/`](tools/) — supporting campaigns:
  [`absorption-recampaign/`](tools/absorption-recampaign/) (graze vs.
  absorption), [`reduced-ode/`](tools/reduced-ode/) (reduced two-population
  flow), [`breath-phase/`](tools/breath-phase/),
  [`manifold-probe/`](tools/manifold-probe/),
  [`noise-test/`](tools/noise-test/) (ruled-out mechanisms),
  [`kinetic-probe/`](tools/kinetic-probe/),
  [`chimera-campaign/`](tools/chimera-campaign/),
  [`beta010-slice/`](tools/beta010-slice/),
  [`transient-tests/`](tools/transient-tests/),
  [`paper-figures/`](tools/paper-figures/), and `unify_checks.py`.
- [`paper/`](paper/) — Springer `sn-jnl` manuscript, critical-review
  supplement, cover letter.
- [`critical_review/`](critical_review/) — analysis scripts backing the
  critical-review supplement (CP1–CPn batches, calibration, figures).
- `paper_figures/` + `*_results/` (`absorption_results/`,
  `campaign_results/`, `kinetic_results/`, `manifold_results/`,
  `noise_results/`, `phase_results/`, `reduced_results/`,
  `transient_results/`) — rendered figures and segregated result caches per
  experiment.

Each subdirectory's own README/`SUMMARY.md`, where present, is authoritative
for that experiment's methodology and status; this file is the entry point.
See the root [`README.md`](../README.md#anneal--chimera-collapse-aging-and-finite-size-scaling)
for how Anneal fits into the rest of the monorepo, and its
[provenance table](../README.md) for import history (Anneal was split from
the `annealMusic` science subtree with full commit history preserved).
