# Anneal — chimera-death hazard experiment

Survival/hazard analysis of chimera-state collapse: is the death of a finite-N chimera a
**memoryless** escape (constant hazard, exponential lifetimes) or a **structured** one
(non-constant hazard), and how does that depend on system size N and proximity to the
chimera's stability boundary?

**Verdict: STRUCTURED.** See [`SUMMARY.md`](SUMMARY.md). Weibull shape k̂ > 1 (95% CI excludes
1) and the exponential is rejected (LRT p < 0.05) in all 15 (β, N) conditions; the structure
strengthens with N and toward criticality β → β_c. The test was pre-registered:
[`PREREGISTRATION.md`](PREREGISTRATION.md).

## Model
The two-population mean-field Sakaguchi–Kuramoto chimera was found to be a *stable attractor*
with no finite-N hazard (see `results/CP1_FINDINGS.md`), so the study uses the canonical
chimera-death system — the **nonlocally-coupled ring** (Wolfrum & Omel'chenko 2011):

    dθ_k/dt = -(1/(2P+1)) Σ_{|j-k|<=P, ring} sin(θ_k - θ_j + α),  α = π/2 - β,  P = round(r·N)

integrated O(N)/step (circular running-sum mean field) with fixed-step RK4 (numba). Death =
collapse of the chimera's spatial structure: rho_std(t) (spatial std of the local order
parameter) falls below ε_std and stays below for dt_hold; right-censored at T_max.

## Layout
- `config.yaml` — resolved configuration.
- `src/` — engine (`ring_model.py` numpy ref, `ring_fast.py` numba), detector
  (`ring_detector.py`), survival/hazard library (`survival.py`), and checkpoint drivers
  (`cp1_validate.py`, `cp2_pilot.py`, `cp3_ensemble.py`, `cp4_analysis.py`) plus regime-finding
  diagnostics.
- `results/` — figures, `ensemble.csv` (4500 per-run rows), fit CSV/JSON, resolved config.
  `results/traces/` (decimated traces for ε re-sweep) is gitignored — regenerate with CP3.

## Reproduce
    pip install -r requirements.txt
    python3 -m src.survival            # self-test of the survival estimators
    python3 -m src.cp3_ensemble --verify   # RNG-independence checks
    python3 -m src.cp3_ensemble        # full ensemble (4500 runs, ~3 min, writes ensemble.csv + traces)
    python3 -m src.cp4_analysis        # survival & hazard analysis -> figures + cp4_fits.csv

All randomness is per-run `np.random.default_rng(seed)` with globally-unique seeds logged in
`ensemble.csv`, so every trajectory is reproducible and the gitignored traces are regenerable.

## Scope
This is a hazard-**structure** claim (memoryless vs structured), not a critical-**exponent**
claim — see the caveats in `SUMMARY.md`.
