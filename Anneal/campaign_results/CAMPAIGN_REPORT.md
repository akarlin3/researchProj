# Finite-N chimera collapse-time campaign — report

Dynamics-only measurement of the time-to-collapse distribution of the shipped
two-population Sakaguchi–Kuramoto chimera voice versus oscillator count N, with
right-censoring. No audio in the campaign path; the shipped voice and supervisor
are untouched. All numbers are reproducible from `(seed, params, dt, θ, W)` +
`git_hash` + `runner_version`.

Code: `tools/chimera-campaign/`. Artifacts: this directory.

---

## CP0 — Audit recap

| #   | Finding                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Integrator `src/audio/chimera.ts` `chimeraStep` — state = `2·Np` phases, ω=0 mean-field, **RK4**, **dt=0.05** (`DRIFT_DT`), TypeScript.                                                                                                                                                                                                                                                                                                                                |
| 2   | μ=(1+A)/2, ν=(1−A)/2, **α=π/2−β** — confirmed exactly.                                                                                                                                                                                                                                                                                                                                                                                                                 |
| 3   | **δω = 0 exactly; δω/μ = 0.** Identical-ω by construction (no ω term, no per-oscillator frequency). No inharmonicity/detune in the dynamics (detune explicitly zeroed for the chimera path; harmonic partials are an audio-only mapping, excluded here). No per-block timing in the math. The spread-ω heterogeneity lives in a _different_ engine. **The model is exactly homogeneous — collapse is driven purely by finite-N fluctuations, not frequency disorder.** |
| 4   | Canonical seed `seedChimera`: anchor = `rng·2π`; pop 1 = cluster ±0.25 rad; pop 2 = uniform. RNG (mulberry32) consumed anchor → Np jitter → Np incoherent.                                                                                                                                                                                                                                                                                                             |
| 5   | Supervisor collapse detector: `isChimeraAlive` (`max R>0.9 ∧ min R<0.85`); "collapsed" = not-alive sustained ≥ 2.0 s.                                                                                                                                                                                                                                                                                                                                                  |
| 6   | Basin map = `examples/probes/chimera_probe.mjs`; results in `docs/CHIMERA_CHARACTERIZATION.md`. "Out of basin" = fracLive < 0.8.                                                                                                                                                                                                                                                                                                                                       |
| 7   | Headless path already exists (the probe); minimal surface = `seedChimera, chimeraStep, orderParam, isChimeraAlive` + mulberry32.                                                                                                                                                                                                                                                                                                                                       |
| 8   | Fully deterministic from `(seed, params, dt)`; mulberry32; no wall-clock on the campaign path.                                                                                                                                                                                                                                                                                                                                                                         |

**δω/μ = 0.** (Standalone deliverable.)

---

## CP1 — Headless runner + validation gate

`integrator.mjs` reproduces `src/audio/chimera.ts` **bit-for-bit** (cross-check
test: identical seed phases and RK4 trajectory over 4000 steps, maxdiff = 0).

**Gate (old collapse definition, A=0.5, β=0.05, 24 seeds, seed0=5000):**

| N   | measured          | documented | verdict |
| --- | ----------------- | ---------- | ------- |
| 8   | **62.5%** (15/24) | 63%        | PASS    |
| 32  | **91.7%** (22/24) | 92%        | PASS    |

**GATE: PASS** — exact reproduction of the documented basin.

---

## CP2 — Paper-grade criterion + robustness

Criterion: collapse at first t where `min(R₁,R₂) > θ` sustained for `W`; lifetime
= that first-crossing time, else right-censored at `t_max`. Defaults **θ=0.85**
(the supervisor's INCOH_LO, just above the measured ~0.83 healthy-breathing
envelope) and **W=5 s** (≫ the breath's fast phase). Unit-tested on synthetic
R(t) traces (sustained-excursion, spike-rejection, run-reset, censoring).

- **Absolute τ** is criterion-sensitive: across the 3×3 (θ,W) grid it moves
  smoothly/monotonically over ~52% of its mean (stricter θ,W date collapse
  later) — see `cp2_robustness.md`.
- **Scaling shape** is **ROBUST**: at loose/default/strict criteria the plateau
  ratio τ(64)/τ(8) = 1.65 / 1.37 / 1.31 (spread 0.34) — all agree on _weak
  sub-linear growth → plateau_, none near exponential (orders of magnitude) or
  even linear (≈8×).

**Verdict: the scaling conclusion the paper relies on is criterion-independent.**

---

## CP3 — Campaign

Spec in `tools/chimera-campaign/campaign.config.json` (single committed
artifact). Worker-pool, resumable, append-only JSONL.

- Primary: N∈{4,8,16,24,32,48,64}, **A=0.5, β=0.05, 200 seeds/point**.
- Secondary: same N, **A=0.20, β=0.05, 100 seeds/point**.
- **t_max = 2000 s** — from pilots (max observed lifetime ~223 s), so N=4–16
  collapse almost entirely uncensored; early-exit makes the large cap free.
- Output row: `(N, A, beta, dt, seed, lifetime, censored, theta, W, t_max,
git_hash, runner_version, wall_ms)`.

**Compute estimate (measured):** full campaign = **2100 runs in ~8 s on 4 cores**
(~260 runs/s; 2.5M RK4 steps). Trivially runnable locally; no need for Origin.
A 10× denser campaign (2000 seeds/point) would still be ~1–2 min.

---

## CP4 — Timestep sanity check

N=8 and N=32 (A=0.5), 40 seeds at dt=0.05 vs dt/4=0.0125.

| N   | KS p (uncensored) | τ̂(dt)               | τ̂(dt/4)             | verdict |
| --- | ----------------- | ------------------- | ------------------- | ------- |
| 8   | 1.000             | 54.9 s [40.3, 76.9] | 54.9 s [40.3, 76.9] | PASS    |
| 32  | 1.000             | 68.5 s [50.3, 95.8] | 68.5 s [50.3, 95.8] | PASS    |

**Verdict: PASS** — lifetimes are identical to 0.1 s; KM curves overlay, τ̂ CIs
coincide. The shipped dt is fully converged; the collapse statistics are not a
timestep artifact. (`cp4_dt_km.png/pdf`, `cp4_verdict.md`.) dt was **not** switched.

---

## CP5 — Survival analysis (outputs)

`lifelines` would not build here (autograd-gamma wheel failure), so Kaplan–Meier
(Greenwood CIs) and the exponential-censored MLE (exact χ²/Garwood-Poisson CIs)
are implemented directly in `analysis.py`.

**Headline result — τ(N) is sub-exponential, not the ring-topology law.**

| Regime | exp rate c (τ∝e^{cN}) | power p (τ∝Nᵖ) | ΔAIC (pow−exp) | preferred |
| ------ | --------------------- | -------------- | -------------- | --------- |
| A=0.5  | 0.0047 (≈0)           | 0.151          | −2.95          | power-law |
| A=0.2  | −0.0097 (≈0)          | −0.230         | −0.25          | power-law |

At **A=0.5**, τ rises from ~39 s (N=4) to ~72 s (N=16) then **plateaus at ~64–68 s**
through N=64 — τ(64)/τ(4)=1.74 over a 16× N range (exponential would be orders of
magnitude). At **A=0.2**, τ is flat at ~28–30 s. The KM panels show _why_: small-N
survival is memoryless (immediate decay), while large-N survival develops a
**shoulder** (a near-flat plateau then a steeper drop) — finite-N delays collapse
onset but the characteristic lifetime saturates. **The two-population topology
does not exhibit the exponential τ(N) scaling of ring chimeras.**

### Output inventory (`campaign_results/`)

| File                                    | Contents                                               |
| --------------------------------------- | ------------------------------------------------------ |
| `collapse_campaign.jsonl`               | 2100 raw runs (one row each).                          |
| `summary_table.{csv,md}`                | Per-(A,N): KM median + exponential-MLE τ̂ with 95% CIs. |
| `fit_comparison.md`                     | Power-law vs exponential τ(N) AIC comparison.          |
| `tau_vs_N.{png,pdf}`                    | τ(N) semilog-y, both A, χ² CIs.                        |
| `km_survival_primary.{png,pdf}`         | KM curves per N, A=0.5 (one panel per N).              |
| `km_survival_secondary.{png,pdf}`       | KM curves per N, A=0.2.                                |
| `cp2_robustness.{csv,md}`               | (θ,W) grid + scaling-shape robustness.                 |
| `cp4_dt_km.{png,pdf}`, `cp4_verdict.md` | Timestep sanity check.                                 |

---

## Anomalies / deferred

- **A=0.20, N=8 exponential-MLE τ̂ = 109 s is an outlier** (vs KM median 25.9 s and
  neighbours ~28 s). It is inflated by 4/100 right-censored survivors (each adding
  2000 s of exposure to the MLE numerator). This is expected behaviour of the
  exponential MLE under a heavy censored tail and is exactly why **both** estimators
  are reported — the KM median is the robust readout there. Not a bug.
- **Large-N survival is non-exponential** (the KM shoulder ⇒ increasing hazard,
  Weibull-like). The exponential MLE is therefore a summary mean-lifetime, not a
  distributional fit; a Weibull/aging fit is a natural follow-up if the paper wants
  to quantify the shoulder. Deferred.
- **No log-rank test** (no `lifelines`): CP4 uses a two-sample KS on the uncensored
  subsets plus τ̂ CI overlap. Sufficient for the PASS verdict here (p=1.0, identical
  lifetimes); a log-rank could be added if `lifelines` is installed on Origin.
- Censoring is rare at A=0.5 (0 across all N at t_max=2000 s) and only appears at
  small N for A=0.2 — consistent with the sub-exponential, saturating lifetimes.
