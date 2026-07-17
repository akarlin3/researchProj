# Anneal — chimera collapse, aging, and finite-size scaling

**Paper:** *"Chimera Collapse Ages: Topology-Dependent Finite-Size Scaling in Mean-Field and
Ring Oscillator Systems"* (submitted to *Nonlinear Dynamics*).

Anneal studies *chimera death* — the spontaneous collapse of the coexisting
synchronized/desynchronized state in coupled identical oscillators — and asks which features
of that collapse are universal versus topology-dependent. The work grew out of a real-time
additive music synthesizer whose two-population Sakaguchi–Kuramoto engine drives partial
amplitudes (hence the origin repo name `annealMusic`), which prompted the operational question
of what a "chimera death" detector actually measures. The study runs matched long-duration
survival experiments on two canonical substrates — the two-population mean field and the
nonlocally coupled ring near its existence boundary — under a pre-registered protocol, and
validates against a parameter-free reduced two-population flow.

**Headline findings:** the standard first-passage death criterion undercounts true lifetimes by
~50% (up to 98% of near-synchrony excursions are self-healing "grazes"); collapse *ages* at the
hazard level in both topologies (Weibull shape k > 1 in all 15 pre-registered cells) and arrives
via a breath-synchronized ratchet; but finite-size scaling is topology-dependent (mean-field
lifetimes plateau ~139 s while near-boundary ring lifetimes *decrease* with N). The reduced flow
predicts the breath period and per-cycle dynamics but underpredicts finite-N lifetimes by a
constant 3.2×, leaving the prolongation mechanism as a sharply constrained open problem.

## Layout

- [`anneal-hazard/`](anneal-hazard/) — the main pre-registered ring hazard experiment: numba
  engine, survival/hazard library, CP1–CP4 validation gates, `PREREGISTRATION.md`, `SUMMARY.md`
  (verdict: **STRUCTURED**).
- [`tools/`](tools/) — supporting campaigns: `absorption-recampaign/` (graze vs. absorption),
  `reduced-ode/` (reduced two-population flow), `breath-phase/`, `manifold-probe/`,
  `noise-test/` (ruled-out mechanisms), `chimera-campaign/`, `kinetic-probe/`, `beta010-slice/`,
  `paper-figures/`, `xlabel/`, and `unify_checks.py`.
- [`paper/`](paper/) — Springer `sn-jnl` manuscript sources, `critical_review_supplement.tex`/
  `.pdf`, conversion/reconstruction changelogs, `appendix-drafts.md`.
- [`critical_review/`](critical_review/) — the critical-review response campaign (CP1–CP3
  batches, calibration, figures, `CHANGELOG.md`, `verification/`).
- [`paper_figures/`](paper_figures/) — rendered manuscript figures plus `FIGURES_REPORT.md`.
- `*_results/` (`absorption_results/`, `campaign_results/`, `kinetic_results/`,
  `manifold_results/`, `noise_results/`, `phase_results/`, `reduced_results/`,
  `transient_results/`) — segregated result caches per supporting experiment.
- `death_of_a_chimera*.tex`/`.pdf`, `chimera_collapse_ages_NLD_v6.docx`,
  `anneal_NLD_cover_letter_FINAL_v2.tex` — manuscript drafts and submission cover letter at the
  project root.

## Reproduce

The pre-registered headline result lives in `anneal-hazard/`; see
[`anneal-hazard/README.md`](anneal-hazard/README.md) for its own reproduce steps
(`pip install -r requirements.txt`, then the `src.cp3_ensemble` / `src.cp4_analysis` drivers).
Figure regeneration for the manuscript is driven from the `Anneal/` project root via
`python3 tools/paper-figures/run_all.py` (see [`tools/paper-figures/README.md`](tools/paper-figures/README.md)).
