# Breath-phase clustering + absorption check

Tests Avery's **breath-synchronized-collapse** hypothesis on the merged finite-N
chimera collapse-time campaign (`tools/chimera-campaign/`):

1. **Do collapse times cluster at a specific phase of the breathing cycle?** (CP2)
   If collapse is a per-pass event keyed to the breath's high-R excursion, the
   collapse phase φ_c relative to the preceding R_incoh peak should be non-uniform.
2. **Are the θ-crossings the criterion counts as collapse true absorptions, not
   grazes?** (CP3) Post-crossing min(R₁,R₂) should stay absorbed (recovery ≈ 0).

Everything here is **new** code that imports the shipped-identical integrator
primitives via `../manifold-probe/trace.mjs` → `../chimera-campaign/integrator.mjs`
(verified bit-identical to `src/audio/chimera.ts`). The shipped voice / supervisor /
collapse criterion are untouched; the campaign JSONL is read-only input. Dynamics
only — no audio path.

## Pipeline

```
node tools/breath-phase/cp1_trace.mjs      # CP1: traced re-runs + determinism gate
python3 tools/breath-phase/analysis.py     # CP2 + CP3 + figures + PHASE_REPORT.md
```

`cp1_trace.mjs` re-integrates the 100 lowest-id non-censored campaign seeds per point
(N ∈ {8,16,32,64} × A ∈ {0.5,0.2}) at the campaign stride (sampleStride=0.1, **required**
so the collapse criterion reproduces the logged lifetime bit-for-bit), recording
R_incoh(t) = min(R₁,R₂) until collapse + a 60 s tail. The **determinism gate** asserts
every traced lifetime equals its logged campaign lifetime.

`analysis.py` (needs `numpy scipy matplotlib`):

- **CP2** — breath period T_b = median peak-to-peak interval of smoothed-R_incoh maxima
  over the pre-collapse window _excluding the final cycle_; collapse phase
  φ_c = 2π·(t_collapse − t_prevpeak)/T_b (peak = φ 0); **Rayleigh test** for
  non-uniformity (z = n·R̄², p ≈ e^(−z)(1+(2z−z²)/4n); Mardia & Jupp 2000). Runs not
  completing ≥2 full cycles are excluded and counted as the early-collapse fraction.
- **CP3** — post-crossing min R_incoh over the ≥60 s tail: recovery (dips below 0.80
  sustained ≥W after the crossing), terminal absorption (final 5 s above 0.80), graze
  occupancy, and the sub-W graze rate per run-hour.

## Config & reproducibility

All knobs live in `phase.config.json` (model/criterion mirror the campaign; seed
selection, smoothing, prominence, cycle threshold, recovery band). Figures are
reproducible from config + the committed `phase_results/cp1_traces.jsonl`.

## Outputs (`phase_results/`)

| file                             | contents                                                   |
| -------------------------------- | ---------------------------------------------------------- |
| `cp1_traces.jsonl`               | one compact row per traced run (R_incoh series + lifetime) |
| `cp1_determinism.json`           | determinism-gate verdict + per-point match counts          |
| `cp2_rose.{png,pdf}`             | per-point + pooled rose plots of collapse phase            |
| `cp2_example_trace.png`          | example breath trace with detected peaks + collapse        |
| `cp3_absorption.{png,pdf}`       | post-crossing min-R_incoh distributions per point          |
| `cp2_table.csv`, `cp3_table.csv` | the report tables                                          |
| `phase_metrics.json`             | all computed metrics + verdicts                            |
| `PHASE_REPORT.md`                | the final report                                           |

## Headline result

Breath-locking is **supported wherever testable** (all four A=0.5 points reject
uniformity, R̄ tightening 0.43→0.69 with N; pooled p≈0), but A=0.2 is dominated by
near-total early collapse (≤2 breaths) and cannot be tested. **CP3 refutes the
absorption prediction:** recovery is large (max 98%), so the W=5 s criterion frequently
dates collapse to a long **graze** that re-forms before (usually) re-absorbing — the
campaign's lifetimes are first-passage-to-a-long-graze times, and the criterion should
be hardened before they are read as absorption times. See `phase_results/PHASE_REPORT.md`.
