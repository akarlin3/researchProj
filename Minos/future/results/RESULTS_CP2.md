# CP2 — applied decision–calibration gap on Fashion-calibrated IVIM posteriors

> **PROVISIONAL — pending Fashion publication** (retooled, in review at *NMR in
> Biomedicine*). Every number here consumes Fashion's uncertainty generators (in
> review). See [`../ASSUMPTIONS.md`](../ASSUMPTIONS.md). Reproduce:
> `proteus/bin/python Minos/future/applied/gap_applied.py` (raw numbers in
> [`RESULTS_CP2.json`](RESULTS_CP2.json)). Clean synthetic data only (Fashion IVIM simulator,
> pancreatic-anchor priors); no `pancData3`/MSK/clinical data.
>
> **Retool survival (re-run against the retooled upstream).** Minos imports no
> Fashion/Caliper *ruler* code, so the decision-gap headline re-runs unchanged: the
> max coverage-calibration regret stays **3.2 utility units** and the gap (`tau* ≠
> tau_stat`, driven by posterior skewness γ) persists, **concentrated in the
> high-D\* tercile** — the same regime where the retool's honest *conditional* D\*
> under-coverage (0.63 [0.60, 0.67], high-D\* tercile) lives. The dropped marginal
> 0.30/0.67 was never a Minos input. **Verdict: SURVIVES**, now honestly sized.

## What was run

A synthetic IVIM cohort (varying truths from Fashion's pancreatic priors, b = clinical-sparse,
SNR 40) fit with `OGC_AmsterdamUMC_biexp`, with Fashion uncertainty from two generators:
**MCMC** (n=300; Fashion's skew-aware/calibrated recipe — bounded credible posterior) and
**bootstrap** (n=500; resampling SD). For each parameter (D\*, f) and cost asymmetry
λ ∈ {1,2,3,4} we computed, on the real (θ, μ, σ):

- `tau_stat` — scale making the central-90% interval `[μ ± z·τ·σ]` cover at nominal rate;
- `tau*` — scale maximising the realised treat/spare/escalate utility (decision core reused
  read-only from `minos.decision` / `minos.utility`);
- **raw gap** `G_raw = tau* − tau_stat`;
- **standardized gap** `G_std = tau*_cal − 1` after coverage-calibrating first
  (`σ_cal = tau_stat·σ`), the clean test of Theorem 1's `G = (1/6)|z*(λ)|·γ`;
- `γ` = skewness of the standardised error `u = (θ−μ)/σ`.

Fashion's **Laplace/Hessian** σ for D\* was *excluded*: it blows up (σ→0 ⇒ |u|~10⁷) at the
identifiability wall — a numerical instability consistent with Gauge's high-D\* wall, documented
rather than patched.

## Headline numbers (PROVISIONAL)

| gen | param | λ | tau_stat | std(u) | γ | G_raw | G_std | G_theory | regret* |
|---|---|---|---|---|---|---|---|---|---|
| mcmc | D\* | 3 | 2.00 | 1.70 | +0.02 | +2.25 | +1.15 | +0.001 | +0.94 |
| mcmc | f  | 3 | 1.60 | 1.55 | +0.42 | +0.57 | +0.54 | +0.030 | +0.001 |
| boot | D\* | 3 | 1.44 | 2.24 | +2.99 | +3.27 | +2.23 | +0.218 | +2.56 |
| boot | f  | 3 | 1.46 | 1.48 | −0.21 | −1.26 | −0.80 | −0.015 | +0.000 |

\*regret = utility lost by using the coverage-calibrated bar instead of the decision-optimal
one (≥ 0 by construction; D\* costs are large because D\* dominates the decision here).

## Honest verdict (the gate)

1. **The qualitative theory claim HOLDS.** On real Fashion-calibrated IVIM posteriors,
   coverage-calibration is **not** decision-calibration: `tau* ≠ tau_stat`, a real gap exists,
   and using the reported bar for the decision leaves substantial regret (up to ≈ +3.2 utility
   units for D\*). This is the central Minos thesis, confirmed on applied data.

2. **Fashion's reported σ is itself miscalibrated** (mean `tau_stat = 1.62`, mean `std(u) = 1.75`):
   the reported D\*/f uncertainty **under-disperses** the true error. This reproduces Fashion's
   D\* under-coverage finding and extends it to f. The **raw** gap is dominated by this coverage
   miscalibration, not by the skew term.

3. **The quantitative leading-order law `G = (1/6)|z*(λ)|·γ` does NOT transfer.** Even after
   coverage-calibration, for the cells with non-negligible skew (|γ|>0.3) the median
   `|G_std / G_theory| ≈ 17×`. Sign agrees in 10/12 cells, but a 17× magnitude gap means
   sign-agreement is **not** quantitative support. Cause: real IVIM errors are **biased**
   (E[u] ≠ 0; e.g. −0.71 for MCMC D\*) and **over-dispersed** (std(u) > 1) and heavy-tailed —
   deviations the small-skew, zero-mean, unit-variance idealization abstracts away. The measured
   skew is also strongly **generator-dependent** (MCMC D\* γ≈0 vs bootstrap D\* γ≈3), so "the
   posterior skew" is not a single well-defined applied quantity.

**Nothing was tuned.** This is the honest outcome the CP2 gate asks for: the *mechanism* the
theory names (a decision-vs-coverage gap) is real and clinically material on Fashion-calibrated
posteriors; the *clean scaling law* is an idealization that does not survive contact with biased,
over-dispersed IVIM error distributions. CP4 reports the theory as the *idealized* account and
the applied gap as the *empirical, miscalibration-dominated* reality — not as a clean confirmation.

## Degeneracies flagged

- `bootstrap / f / λ=3`: standardized `tau*` hit the grid edge (decision near-insensitive for f
  at this threshold) — flagged, not silently dropped.

## Status

CP2 GATE: **PASS** — runs on Fashion-calibrated posteriors; all numbers PROVISIONAL; theory
comparison reported honestly (qualitative hold, quantitative miss).
