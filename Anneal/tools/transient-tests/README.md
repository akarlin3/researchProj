# transient-tests — deterministic-transient tests of the A=0.5 chimera

Tests Avery's hypothesis that, at A=0.5, the finite-N chimera collapse is a
**deterministic transient**: the breathing trajectory spirals outward, each
high-R excursion closer to sync capture. Two falsifiable predictions:

1. **Ratchet** — per-cycle breath maxima `M_k` (k-th local max of `min(R₁,R₂)`)
   increase monotonically within a run.
2. **Collective-IC predictability** — the absorption lifetime `t_abs` is
   (near-)deterministically set by the collective initial condition, far beyond
   the manifold probe's null on D₀ (|ρ|≤0.15).

New code only. All prior result directories are read-only. Peak detection reuses
PR #41's canonical detector (`absorption-recampaign/breath.mjs detectPeaks`)
verbatim; initial conditions are recomputed from the logged seed via
`seedChimera` with **zero integration** (the IC is a pure function of the seed).

## Run

```bash
node   tools/transient-tests/extract.mjs     # ~1.5 min (re-traces 160 A=0.2 runs)
python3 tools/transient-tests/analysis.py    # stats + figure + report
```

Everything is driven by `transient.config.json` and is deterministic (fixed CV
seed, fixed bootstrap RNG).

## What each stage does

`extract.mjs` (the only integration here is the A=0.2 re-trace):
- **CP0** trace coverage per (N,A) → `coverage.json`.
- **CP1 (A=0.5)** `M_k` sequences from the on-disk traces → `cp1_mk_a05.jsonl`.
- **CP1 (A=0.2)** re-traces ≥40 persistent (`abs_censored`) seeds/point over the
  full `t_max` window → `cp1_mk_a02.jsonl`. A **determinism gate** checks the
  re-traced `t_graze`/censor flags against the logged campaign bit-for-bit
  (`determinism_gate.json`).
- **CP2** collective t=0 IC features (`R_incoh0`, `R_sync0`, `Δφ0`, `|Δφ0|`)
  from the seed, plus first-cycle `M₁`/`T_b,1` where traced → `cp2_features.jsonl`.

`analysis.py`:
- **CP1** monotone (lenient ≤1 inversion + strict) ratchet fractions, grand mean
  per-cycle increment ⟨ΔM⟩ with a cluster-over-runs bootstrap CI, ensemble ⟨M_k⟩
  forward and backward (`cp1_ensemble_mk.csv`), and the `log(θ−M_k)` spiral-out
  slope distribution over the sub-θ approach peaks. A=0.2 stationary contrast.
- **CP2** per-N Spearman ρ for each t=0 feature, 5-fold CV R² (linear and
  quadratic+interaction) on t=0 features alone, and on the traced subset the
  ΔR² from adding the first-cycle features — all printed against the D₀ null.
- **CP3** `transient_decider.{png,pdf}` (a: ensemble ⟨M_k⟩ vs k with the A=0.2
  band; b: predicted vs actual log t_abs colored by N) and
  `TRANSIENT_REPORT.md` with explicit pass/fail verdicts.

## Headline verdicts (see `transient_results/TRANSIENT_REPORT.md`)

- **Ratchet: YES.** ⟨ΔM⟩ is positive with a strictly-positive bootstrap CI at
  all four N; the strict (0-inversion) fraction climbs 0.18→1.00 and ⟨ΔM⟩
  0.031→0.091 with N. Spiral-out slopes are negative for 100% of fittable runs.
  The A=0.2 never-absorbers are flat (⟨ΔM⟩≈0, ratchet fraction 0.00) under the
  identical detector — so the A=0.5 drift is not a detector artifact.
- **Predictability: supported relative to the null, partial in absolute terms.**
  `|Δφ0|` predicts log t_abs with Spearman |ρ|≈0.43–0.59 (vs the D₀ null
  |ρ|≤0.08) — a strong, clean, previously-unreported predictor. CV R² of the
  collective t=0 model is ≈0.2–0.44 (null R²≈0.007), i.e. two orders of
  magnitude above the manifold-distance null, but well short of 1, so these
  low-dimensional summaries fix only part of the lifetime. The first-cycle
  features add no consistent ΔR².

## Outputs (`transient_results/`)

| file | contents |
|---|---|
| `TRANSIENT_REPORT.md` | full report with verdicts |
| `transient_decider.{png,pdf}` | the two-panel deciding figure |
| `coverage.json` | CP0 trace coverage |
| `determinism_gate.json` | re-trace vs logged-label check |
| `cp1_mk_a05.jsonl`, `cp1_mk_a02.jsonl` | per-run M_k sequences |
| `cp1_ratchet.csv`, `cp1_a02_contrast.csv`, `cp1_ensemble_mk.csv` | CP1 tables |
| `cp2_features.jsonl`, `cp2_spearman.csv`, `cp2_r2.csv` | CP2 features and tables |
