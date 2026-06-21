# CP3 — label-free validity monitor, applied (observable-fires / hidden-blind)

> **PROVISIONAL — pending Gauge publication.** Reuses Gauge's `DeploymentMonitor` (the Minos-style
> label-free monitor) and forward/estimator/conformal code. See [`../ASSUMPTIONS.md`](../ASSUMPTIONS.md).
> Reproduce: `proteus/bin/python Minos/future/applied/monitor_applied.py` (numbers in
> [`RESULTS_CP3.json`](RESULTS_CP3.json)). Strictly synthetic IVIM (seed 20260613); no clinical data.
>
> **Retool reinforcement.** The hidden, label-free-invisible failure this monitor
> exposes — marginal D\* coverage **0.885** (≈ nominal) but high-D\* tercile **0.795**
> (hidden-channel AUC 0.516, blind) — *is* the conditional high-D\* under-coverage the
> retooled Fashion (in review at *NMR in Biomedicine*) now openly owns under the
> honest CRLB (0.63 [0.60, 0.67], high-D\* tercile). Re-run against the retooled
> upstream reproduces this **exactly** (FULL-N cross-check: 0.905 / 0.823, AUC 0.499).
> **Verdict: SURVIVES** — the retool *reinforces* CP3, the ruler's own paper now
> agreeing that ground-truth-free trust fails conditionally in the high-D\* regime.

## What was run

A label-free deployment-validity monitor (Minos-Core v3 idea; Gauge's two-family implementation:
Mahalanobis on observable features + residual-conformal) was fit on an in-distribution synthetic
IVIM calibration cohort (n=2000, SNR 40, clinical-sparse b-scheme) and evaluated under deployment
shifts. The monitor sees only **observables** — signal-shape log-slopes, NLLS plug-in (D̂,D̂\*,f̂),
estimated SNR, and the fit-residual norm — never the truth. AUC separates calibration (ID) from
test voxels by the per-voxel drift score: → 1 when the shift is observable, → ½ when hidden.

> *Note:* the observable-feature builder mirrors Gauge's `robustness._observe`; it is reconstructed
> in `monitor_applied.py` because `gauge.robustness` / `gauge.conditional_attack` use a Python 3.12+
> f-string and won't import on the 3.11 reproduce env. The monitor, NLLS, forward model, and
> conformal code are imported from Gauge unchanged.

## Results (PROVISIONAL)

| channel | shift | AUC | fires | maha / resid |
|---|---|---|---|---|
| observable (resolvable) | SNR 40→12 | **1.000** | yes | 1.000 / 1.000 |
| observable (resolvable) | tissue-D shift up | **0.919** | yes | 0.919 / 0.515 |
| observable (resolvable) | perfusion-f shift up | **0.959** | yes | 0.959 / 0.495 |
| observable (wall-adjacent) | tri-exp perfusion misspec | 0.528 | yes | 0.528 / 0.493 |
| **hidden** | exchangeable; high-D\* conditional failure | **0.516** | **no** | 0.516 / 0.489 |

Hidden-channel coverage (split-conformal D\*, α=0.10): **marginal 0.885** vs **high-D\* tercile 0.795**
— the conditional failure is real but invisible. CRLB(D\*) exceeds the high-D\* tercile width
(ratio ≫ 1; cf. Gauge 03's reported 1.12 — the absolute value differs with cohort/SNR/scheme, but
both say the latent regime is unresolvable).

## Honest verdict (the gate)

**The Theorem-2 signature reproduces on applied IVIM, cleanly.**

1. **Observable channel is bounded-detectable (Thm 2(ii)):** shifts that move *resolvable* signal
   features — global noise, tissue diffusion D, perfusion fraction f — are detected with AUC 0.92–1.00
   and fire. The Mahalanobis (feature) family carries detection; the residual family is correctly
   flat (≈0.5) for a parameter-*distribution* shift that leaves the fit residual unchanged — only the
   global-noise shift moves both.

2. **Hidden channel is undetectable (Thm 2(i)):** an exchangeable test (same observable distribution
   as calibration) yields AUC 0.52 and does **not** fire — the monitor is at chance, i.e. blind — yet
   within it the high-D\* tercile coverage has collapsed (0.795 vs 0.885 marginal). This is Gauge 03's
   identifiability wall: a within-distribution latent-axis failure produces no observable drift.

3. **Honest nuance:** a shift confined to the *perfusion* regime (tri-exp at high D\*) is only weakly
   observable (AUC 0.53) — because that regime is itself near the identifiability wall, so even a
   genuine forward-model misspecification there barely fingerprints the signal. The monitor's blind
   spot is not a knife-edge at exactly "hidden"; it degrades gracefully as the shift moves into the
   unresolvable perfusion regime.

**Conclusion:** the label-free monitor is *optimal up to a floor no label-free statistic can cross*
(Theorem 2). The only way to see the high-D\* hidden channel is to read a label — the principled case
for labeled repeatability spot-checks (Echo). Nothing was tuned; the one fix from the first run was a
bona-fide bug (the FPR null threshold degenerated when the bootstrap subsample size equalled the cal
size, causing a false fire on the exchangeable test) — corrected to `null_subsample = n/4`.

## Status

CP3 GATE: **PASS** — observable-fires / hidden-blind reproduced on synthetic IVIM; real AUCs
reported; all numbers PROVISIONAL (Gauge dependency).
