# Anneal — chimera collapse, aging, and finite-size scaling

*Paper:* **"Chimera Collapse Ages: Topology-Dependent Finite-Size Scaling in Mean-Field
and Ring Oscillator Systems"** (submitted to *Nonlinear Dynamics*).

Anneal studies *chimera death* — the spontaneous collapse of the coexisting
synchronized/desynchronized state in coupled identical oscillators — and asks which
features of that collapse are universal versus topology-dependent. Runs matched
long-duration survival experiments on two canonical substrates (two-population mean
field; nonlocally coupled ring near its existence boundary) under a pre-registered
protocol, validated against a parameter-free reduced two-population flow.

See the root [README.md](../README.md#anneal--chimera-collapse-aging-and-finite-size-scaling)
for the full writeup and headline findings.

## Layout
- `anneal-hazard/` — the main pre-registered ring hazard experiment (numba engine, survival/hazard library, CP1–CP4 validation gates, `PREREGISTRATION.md`). Has its own [README](anneal-hazard/README.md).
- `tools/` — supporting campaigns: `absorption-recampaign/` (graze vs. absorption), `reduced-ode/` (reduced two-population flow), `breath-phase/`, `manifold-probe/`, `noise-test/` (ruled-out mechanisms), `paper-figures/`, and more.
- `paper/` — Springer `sn-jnl` manuscript, critical-review supplement, cover letter.
- `critical_review/` — critical-review supplement materials.
- `paper_figures/` + `*_results/` (`absorption_results/`, `campaign_results/`, `kinetic_results/`, `manifold_results/`, `noise_results/`, `phase_results/`, `reduced_results/`, `transient_results/`) — rendered figures and segregated result caches per experiment.
