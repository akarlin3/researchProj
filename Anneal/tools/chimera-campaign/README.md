# Finite-N chimera collapse-time campaign

A dynamics-only measurement harness for the time-to-collapse distribution of the
shipped two-population Sakaguchi–Kuramoto chimera voice, versus oscillator count
N, with right-censoring. Built for a physics paper (target: _Chaos_/AIP).

**No audio anywhere in this path.** It reuses the shipped integrator's math
(`src/audio/chimera.ts`) — verified bit-identical by a cross-check test — but
runs headless, deterministic, and fast. The shipped voice and its supervisor are
untouched; everything here is new modules.

## Pipeline

| Stage              | File                                 | What it does                                                                                                                                         |
| ------------------ | ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| Integrator         | `integrator.mjs`                     | Self-contained ESM port of the two-population RK4 mean-field core + mulberry32 + canonical chimera seed. Plain `node`, no deps.                      |
| Cross-check        | `integrator.crosscheck.test.ts`      | Asserts the integrator reproduces `src/audio/chimera.ts` bit-for-bit (zero drift).                                                                   |
| Collapse criterion | `collapse.mjs` + `collapse.test.ts`  | Paper-grade collapse rule: `min(R₁,R₂) > θ` sustained for `W`; lifetime = first-crossing time, else right-censored. Unit-tested on synthetic traces. |
| Online runner      | `runner.mjs`                         | Integrates with early-exit at collapse → collapsing runs cost ~lifetime, not `t_max`.                                                                |
| CP1 gate           | `validate.mjs`                       | Reproduces the documented basin (N=8→63%, N=32→92% at A=0.5,β=0.05) under the old collapse definition.                                               |
| CP2 robustness     | `robustness.mjs`                     | (θ,W) grid + scaling-shape robustness → `campaign_results/cp2_robustness.{csv,md}`.                                                                  |
| CP3 sweep          | `sweep.mjs` + `campaign.config.json` | Worker-pool, resumable, append-only JSONL campaign driver with ETA.                                                                                  |
| CP4 dt check       | `cp4_dt.mjs`                         | 40 seeds at dt and dt/4 (N=8,32) → `campaign_results/cp4_dt.jsonl`.                                                                                  |
| CP5 analysis       | `analysis.py`                        | KM + exponential-censored MLE, CIs, τ(N), power-vs-exp fit, figures, tables.                                                                         |

## Run it

```bash
# 1. validate the runner against the documented basin (≈30 s)
node tools/chimera-campaign/validate.mjs

# 2. unit + cross-check tests
npx vitest run tools/chimera-campaign/

# 3. CP2 robustness table
node tools/chimera-campaign/robustness.mjs

# 4. the campaign itself (resumable; re-run to top up). ~8 s on 4 cores.
node tools/chimera-campaign/sweep.mjs
#    smoke first:   node tools/chimera-campaign/sweep.mjs --seeds 5
#    more cores:    node tools/chimera-campaign/sweep.mjs --workers 16

# 5. CP4 timestep data
node tools/chimera-campaign/cp4_dt.mjs

# 6. survival analysis + figures (writes campaign_results/)
pip install --break-system-packages numpy scipy matplotlib pandas
python3 tools/chimera-campaign/analysis.py --cp4
```

The campaign is **resumable**: `sweep.mjs` skips any `(A, β, N, seed, θ, W, dt,
t_max)` already in the output JSONL, so a killed run continues where it left off,
and bumping `seeds` in `campaign.config.json` only runs the new seeds.

## Determinism

Every row is reproducible from `(seed, params, dt, θ, W)` plus the logged
`git_hash` and `runner_version`. The only randomness is the integer-seeded
mulberry32 initial-phase draw; there is no wall-clock coupling. The model is
exactly homogeneous (identical natural frequency, **δω = 0**), so collapse is
driven purely by finite-N fluctuations.

## Model

Two equal populations of `Np` oscillators, identical ω = 0 (rotating frame):

```
dθᵢ^σ/dt = μ·R_σ·sin(Φ_σ − θᵢ^σ − α) + ν·R_σ'·sin(Φ_σ' − θᵢ^σ − α)
```

μ=(1+A)/2, ν=(1−A)/2, α=π/2−β, RK4 at dt=0.05 (the shipped control-rate step).
Canonical chimera seed: pop 1 a tight synchronized cluster, pop 2 incoherent.
