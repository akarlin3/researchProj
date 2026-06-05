# Absorption-grade re-measurement

Re-measures the finite-N two-population Sakaguchi–Kuramoto chimera collapse-time
campaign under an **absorption-grade** criterion, after PR #41 showed the published
criterion's θ-crossings recover up to 98% of the time. The campaign's published
"lifetimes" are **first-long-graze** times; this module defines a criterion that
only counts a θ-crossing as a collapse if it is **not** followed by recovery within
a verification horizon `T_v`, re-runs the full campaign under it (both labels from
one trace per run), and determines which headline results survive.

New code only. The shipped voice/supervisor/criterion (`src/audio/chimera*.ts`),
the published campaign (`campaign_results/`), and the PR #40/#41 results
(`manifold_results/`, `phase_results/`) are **read-only**. The dynamics core is
the verbatim shipped RK4 (`../chimera-campaign/integrator.mjs`); `t_graze`
reproduces the published campaign **bit-for-bit** (determinism gate).

## Pipeline

```
# CP1 — labeling unit tests (synthetic R(t))
node --test labeling.test.mjs

# CP2 — pilot: t_abs censoring + CP1 (T_v, recThresh) sensitivity
node pilot.mjs

# CP2 — full re-campaign (2100 runs, both labels + T_b + n_grazes), ~90s/4 cores
node sweep.mjs

# CP2 — determinism gate (t_graze == published campaign, bit-for-bit)
node determinism.mjs

# CP3 — phase-clustering subset (R_incoh traces for absorbed runs)
node phase_trace.mjs

# CP4 — supervisor over-trigger replay (shipped detector, no behavior change)
node supervisor_replay.mjs

# CP3 + CP4 + ABSORPTION_REPORT.md (figures, tables, survive/revise/retire)
python3 analysis.py
```

## Files

| File | Role |
| --- | --- |
| `absorption.config.json` | Single committed config (criterion, sensitivity grid, sweeps). |
| `labeling.mjs` | Streaming two-timescale `Labeler` (`t_graze`, `t_abs`) — single source of truth. |
| `labeling.test.mjs` | CP1 unit tests on synthetic R(t). |
| `breath.mjs` | Per-run breath period `T_b` (PR #41 estimator ported to JS; bit-identical cross-check). |
| `tracer.mjs` | Integrates one run, drives the labeler online, computes `T_b`. |
| `sweep.mjs` | Worker-pool full re-campaign → `absorption_campaign.jsonl`. |
| `pilot.mjs` | CP2 censoring + CP1 sensitivity (one integration, re-labeled under every knob). |
| `determinism.mjs` | t_graze == published campaign gate. |
| `phase_trace.mjs` | CP3 subset: R_incoh traces for true-absorption phase. |
| `supervisor_replay.mjs` | CP4 shipped-detector over-trigger measurement. |
| `analysis.py` | CP3 survival/Weibull/geometric/phase/graze + CP4 table + `ABSORPTION_REPORT.md`. |

Outputs land in `absorption_results/`. See `absorption_results/ABSORPTION_REPORT.md`
for the headline survive/revise/retire verdict.

## Absorption criterion (default)

`t_abs` = run-start time of the first `W`-sustained θ-crossing of `R_incoh =
min(R₁,R₂)` that is **not** followed by recovery within `T_v` seconds, where
recovery = `R_incoh < recoveryThreshold` sustained for `≥ recoveryWindowSec` after
the crossing's confirmation. Defaults: θ=0.85, W=5 s, T_v=120 s,
recoveryThreshold=0.80, recoveryWindowSec=5 s. A crossing whose full `T_v` window
is unavailable at `t_max` leaves the run right-censored for `t_abs`.
