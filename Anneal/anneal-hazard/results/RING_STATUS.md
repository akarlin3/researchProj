# Ring-model pivot — engine + regime validated (CP0-redo / CP1 for the ring)

## What was built
- `src/ring_model.py` — numpy reference engine (top-hat nonlocal ring, exact O(N) mean field).
- `src/ring_fast.py` — numba RK4 integrator, O(N) circular running-sum mean field.
  Validated vs numpy (max|Δrho_std| ~ 1e-6 short-time; identical to ~1e-13 step-for-step
  before chaotic divergence). ~10x faster.
- `src/ring_detector.py` — death = rho_std(t) (spatial std of local coherence) below eps_std,
  held for dt_hold. Plus precollapse stats and a near-collapse dwell statistic.
- `src/survival.py` — KM+Greenwood, Nelson-Aalen, Epanechnikov hazard, censored
  exponential/Weibull MLE (profile-likelihood CI on k), Weibull-vs-exp LRT, lnS linearity +
  runs test. Self-test recovers k=1.0 (CI incl. 1) for exponential data and k=1.81 (CI excl. 1,
  LRT p~0) for Weibull k=1.8 data.

## A detection bug was found and fixed (important)
The early-stop terminated integration exactly when the hold window completed, truncating the
stored trace so the detector measured the below-eps run as ~48 t.u. (just under dt_hold=50) and
reported CENSORED even when the chimera had clearly collapsed (rho_std -> 0.02). This made
chimeras look far more stable than they are. Fixed: early-stop now evaluates on decimated
samples (same data the detector sees) and stores a buffer past hold completion. After the fix,
seed 20260613 correctly registers death at tau=374.5, and death fractions jumped from ~0 to ~100%.

## Confirmed phenomenology (beta sweep, r=0.15, detection fixed)
Clear chimeras (rho_std ~ 0.15, sync floor ~ 0.018) that genuinely collapse:
- beta=0.06: 0 deaths to N=192 at T=15000 (deep-stable; lifetime >> 15000).
- beta=0.10-0.12: ~100% deaths, clean lifetimes ~200-7000.

Lifetime-vs-N (this regime, fixed detector):
```
beta=0.10  N=24:med14285(50%cens)  N=48:4376  N=96:3908  N=192:1197
beta=0.12  N=32:1472  N=48:1752  N=64:464  N=96:722  N=128:434   (all ~100% died)
```
Lifetime DECREASES with N in the measurable regime (we sit near/just inside the chimera
boundary, the opposite side from WO's deep-stable "lifetime grows with N" regime, which for us
is unmeasurably long). Direction is a finding, not a blocker: the hazard SHAPE vs N is still
fully recoverable, and the criticality axis is beta -> beta_c^- (~0.13-0.14).

## Preview of the scientific answer
The N=48 lifetime CDF shows a clear lag before deaths begin, then a steep rise — sub-exponential
early = a NON-memoryless (structured) hazard signature (Weibull k>1 expected). The full CP4
analysis will quantify this and its dependence on N and beta.

## Proposed operating point (see config.yaml, pending confirm)
r=0.15, beta_sweep {0.09..0.13} (criticality), N {32,64,128}, T_max=8000, eps_std=0.04,
dt_hold=50, M_pilot=30, M_full=300.

Open question flagged for you: lifetime-vs-N is non-monotonic at beta=0.12 (near-boundary
chaos). Acceptable for a hazard-shape study, but if you want monotone "lifetime grows with N"
(more WO-canonical) it needs a deeper-stable/smaller-N hunt with larger T_max (more compute).
