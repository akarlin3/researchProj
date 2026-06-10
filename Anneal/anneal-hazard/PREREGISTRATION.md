# Pre-registration of the CP4 hazard-shape test
Committed BEFORE running the CP3 ensemble / any CP4 fit. Frozen with the config below.

## Frozen design
- Model: nonlocal-ring (top-hat), r=0.15, alpha=pi/2-beta.
- Conditions: beta in {0.110,0.115,0.120,0.125,0.130} x N in {32,64,128} = 15 cells, M=300 each.
- Event: rho_std(t) < eps_std=0.04 held for dt_hold=50; right-censored at T_max=12000.
- Death-time estimator and all survival code: src/survival.py (KM+Greenwood, Nelson-Aalen,
  Epanechnikov hazard, censored exponential/Weibull MLE, profile-likelihood CI on k, LRT).

## Primary endpoint (pre-specified)
Per condition, the Weibull shape parameter k̂ with profile-likelihood 95% CI from the
right-censored MLE. Directional hypothesis fixed in advance: **k > 1** (increasing hazard).

## Decision rule (pre-specified, from the project spec)
A condition is **STRUCTURED** iff ALL three hold:
1. k̂ 95% CI excludes 1, AND
2. LRT 2(ll_Weibull - ll_exp) gives p < 0.05 (chi^2_1), AND
3. kernel-smoothed hazard ĥ(t) (Epanechnikov, with a bandwidth-sensitivity pass) is visibly
   non-constant.
Otherwise **FLAT** (memoryless / constant hazard = the known attractor-escape result).

## Headline claim being tested (pre-committed)
Chimera collapse in this near-boundary regime is a STRUCTURED (non-memoryless) hazard with
k>1 — consistent with a post-homoclinic chaotic-transient collapse that has a characteristic
approach time, NOT a memoryless escape. Pre-registered corroborating predictions:
- The k>1 signal is CONSISTENT across the 15 cells (report all; no cell selection).
- Mechanism: per-run dwell_stat (terminal committed-descent duration, band rho_std<0.10)
  correlates with lifetime / hazard contribution (Spearman, sign pre-specified: longer-lived
  runs are NOT instantaneous escapes -> positive structure).
- Supporting GOF (not gating): ln S(t) vs t nonlinearity (low R^2, runs-test p<0.05) and a
  KS/AD test of lifetimes vs the fitted exponential.

## Multiple comparisons
15 conditions. We do NOT cherry-pick a single significant cell; the claim rests on the
CONSISTENT direction and significance pattern across cells. All 15 (k̂, CI, LRT p) reported.

## Scope / what we are NOT claiming (per reviewer caveat)
This is a hazard-STRUCTURE claim (memoryless vs structured), NOT a critical-exponent claim.
Dropping beta=0.10 and capping T_max=12000 censors the long-lived low-beta tail, and the
sweep (0.110-0.130) is a short lever arm against an uncertain boundary beta_c~0.13-0.14 — too
little to fit a divergence/scaling exponent of lifetime vs (beta_c - beta). The stored
(tau,event) rows permit later recovery of the low-beta regime via Kaplan-Meier / RMST if a
critical-exponent claim is ever desired; that would require re-introducing the dropped corner
with a much larger T_max.

## N-dependence (reported, not gated)
Lifetime DECREASES with N in this regime (near/just inside the boundary). We report how the
hazard SHAPE (k̂, ĥ(t)) varies with N; we do not claim WO-style lifetime-grows-with-N (that
regime is unmeasurably long-lived at feasible compute).
