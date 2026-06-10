# Kinetic-theory probe harness (anneal v6, exploratory branch)

Verification rig for the one remaining principled candidate for the
N-independent ≈3.2× finite-size prolongation at the operating corner
(A = 0.5, β = 0.05): a **finite-size kinetic-theory correction to the
collective flow** (system-size expansion of the order-parameter dynamics,
Buice–Chow / Hildebrand–Buice–Chow class).

The analytical derivation is a **human task** (CP2). This harness only
integrates and scores whatever the derivation yields, untuned. Three
mechanisms are already excluded under the same bar: additive 1/√N collective
noise (Appendix B), multiplicative breath-locked noise, and off-manifold
Watanabe–Strogatz constants (`paper/revision-data-gated/results_mech.json`).

## What the human must supply (CP2)

The body of `f_corr` in `tools/kinetic-probe/f_corr.py`:

```python
def f_corr(rho1, rho2, psi, N) -> (drift_vec, diff_matrix)
```

- `drift_vec` — ndarray shape `(3,)`: the additional deterministic drift the
  derivation yields on (ρ₁, ρ₂, ψ), in s⁻¹.
- `diff_matrix` — ndarray shape `(3, 3)`: the noise-**amplitude** matrix `B` in

  dX = [rhs₃d(X) + drift(X, N)] dt + B(X, N) dW,

  with `W` a 3-vector of independent standard Wiener processes. The
  Fokker–Planck diffusion matrix is `D = B Bᵀ`; if the derivation yields `D`,
  supply any matrix square root (e.g. Cholesky). Zero matrix = drift-only
  correction.

Also set `F_CORR_META["supplied_by"] = "human"` plus a derivation reference;
`run_probe.py` refuses to run without it (the CP2 gate).

**Conventions (must match):**
- State (ρ₁, ρ₂, ψ) as in `tools/reduced-ode/reduced_core.rhs_3d`:
  μ = (1+A)/2, ν = (1−A)/2, α = π/2 − β, ω = 0; operating corner
  A = 0.5, β = 0.05.
- 1:1 clock — one model time unit = one second; all rates in s⁻¹
  (the corner spiral has σ = +0.01243 s⁻¹, ω = 0.3277 s⁻¹).
- All N-dependence lives inside `f_corr`; the harness sweeps
  N ∈ {8, 16, 32, 64} and never rescales the returned values.
- ρ₁, ρ₂ are kept in [0,1] by reflection (the Appendix B operator); ψ is
  unbounded. If the derivation requires a different boundary treatment near
  ρ → 1 (the homoclinic ghost / capture boundary), say so explicitly — that
  is a physics statement the human must make, not a harness default to lean on.
- **Coefficients are theory-fixed.** No free parameters tuned toward 3.2×.
  If a coefficient is genuinely undetermined, declare its range in
  `F_CORR_META["undetermined_coefficients"]` and provide
  `F_CORR_VARIANTS = [(tag, callable), ...]` spanning it; the runner scores
  every variant and reports the **range** of outcomes without selecting one.
- "The derivation does not close" or "no N-independent term at physical
  order" are valid CP2 answers; they are recorded as a clean negative.

## How a supplied correction is scored

```
python3 tools/kinetic-probe/run_probe.py        # CP3, after CP2 input
```

1. **Anchor gate** — committed per-IC deterministic (DOP853) capture times
   from `reduced_results/reduced_runs.jsonl` must median-match the Appendix B
   anchors {8: 41.65, 16: 42.05, 32: 43.05, 64: 64.0} s exactly; 3 rows/N are
   recomputed with DOP853 and must match bit-identically.
2. **Zero-correction gate** — with `f_corr ≡ 0` the Euler–Maruyama driver
   (dt = 0.01, sampling stride 0.1 s) must reproduce the deterministic capture
   times to < 0.3 % median per N (asserted; same gate as Appendix B).
3. **Sweep** — N ∈ {8,16,32,64} × 200 realizations from the §6.4 seed-mapped
   collective ICs, seeds `9_000_000 + cell*100_000 + N*1000 + j`, t_max = 2000 s.
   Capture = first sustained crossing of min(ρ₁,ρ₂) > θ = 0.85 with no dip
   below 0.8 for 5 s (the finite-N absorption-grade label); censoring at
   t_max is counted and fed to the censored fits.
4. **Score** (`kinetic_results/score.json`) — the four **pre-committed**
   conditions, identical in kind to the bar the excluded mechanisms faced:

   | # | Condition | Threshold |
   |---|-----------|-----------|
   | 1 | prolongation factor: median over N of (median capture)/(per-N deterministic median) | within **[2.9, 3.5]** |
   | 2 | N-independence: CV of the per-N factor across N ∈ {8,16,32,64} | **< 0.15** |
   | 3 | breath-phase locking: Rayleigh p of capture phases (canonical PR#41 breath detector) | **< 0.05 at every N** |
   | 4 | rising-in-cycles hazard: censored-Weibull k_cyc (cycles n = t/T_b), 95 % profile CI | **CI low end > 1 at every N** |

   Verdict: **PASS** = all four; **PARTIAL** = (1) holds but not all four;
   **FAIL** = (1) fails. All four are always reported individually. No
   coefficient is adjusted between runs.

**Context for honest reading (non-binding diagnostics in the score):** under
per-N referencing the *measured* finite-N system itself shows factors
3.44 / 3.16 / 3.03 / 1.98 (CV 0.190), because the reduced ensemble's median
crosses one extra breath cycle at N = 64 (43 → 64 s) while the measured
plateau stays flat. The score therefore also reports the factor against the
pooled deterministic reference (the manuscript's 3.2× convention) and the CV
over N ∈ {8,16,32} — reported, never substituted for conditions (1)–(2).

## Pre-registered secondary criterion (frozen 2026-06-10, before any CP2 input)

The primary four conditions are kept verbatim for apples-to-apples comparison
with the excluded mechanisms — but the measured system itself scores only
PARTIAL under them (cond 2: CV 0.190; cond 3: Rayleigh p = 0.165 at N=32 and
0.099 at N=64). A correction that passed cond 3 at *every* N would in fact be
**wrong about the data**. So the scorer also asks, as a separate frozen
criterion approved by the human before any theory existed: **does the
correction reproduce the measured per-N pattern?**

| # | Criterion | Threshold |
|---|-----------|-----------|
| S1 | per-N factor vs measured factors 3.444 / 3.158 / 3.028 / 1.984 (incl. the N=64 dip) | relative deviation ≤ **20 %** at every N |
| S2 | breath-phase pattern as measured | Rayleigh p < 0.05 at N = 8 **and** 16; R̄(8) > R̄(64) |
| S3 | cycle-hazard shape as measured | k_cyc 95 % CI overlaps the measured CI at every N (2.13 [1.91,2.36], 2.15 [1.92,2.40], 2.11 [1.86,2.37], 2.89 [2.55,3.24]) |

Secondary verdict: `matches_measured_pattern` iff S1 ∧ S2 ∧ S3. Reported
alongside the primary verdict; neither overrides the other. Tolerance
rationale for S1 (fixed in advance): same-mechanism seed-replication of the
Appendix B null reproduced committed factors to 2–26 % (worst in broad cells),
so 20 % bounds sampling noise, while the additive null at physical amplitude
misses the measured factors by ~70 %. Targets load from the committed
`paper/revision-data-gated/results_mech.json` (experiment 2b "actual" rows =
genuine finite-N Eq.-1 runs), never hand-entered. The null self-test asserts
the additive mechanism FAILS this secondary criterion too.

## Null self-test (why you can trust a FAIL/PASS from this rig)

```
python3 tools/kinetic-probe/null_test.py
```

Feeds the harness the **already-excluded** additive 1/√N noise as `f_corr`
(zero drift, B = (c/√N)·I) at the physical amplitude c = 0.05 and the strong
amplitude c = 0.2, and asserts the scorer rejects it, reproducing the
committed Appendix B result: no ~3× prolongation at the physical amplitude
(condition 1 fails), N-dependent prolongation at c = 0.2 (condition 2 fails,
factor(N=8)/factor(N=64) > 1.5). It also asserts the driver is bit-identical
to the Appendix B integrator (`tools/noise-test/em_core.em_run`) for this
mechanism at equal seeds. Outputs `kinetic_results/score_null_c*.json` and
`kinetic_results/null_selftest.json`.

## Files

| File | Role |
|------|------|
| `f_corr.py` | CP2 plug-in (human analytical input; zero placeholder) |
| `f_corr_zero.py` | the zero correction used by the determinism gate |
| `harness.py` | EM driver (rhs₃d + pluggable correction), ICs, gates, sweep |
| `scorer.py` | four-condition scorer → `score.json` (thresholds frozen) |
| `null_test.py` | additive-noise null self-test (must FAIL) |
| `run_probe.py` | CP3 runner for the supplied derivation |
| `kinetic_results/` | runs (jsonl), scores (json) — all seeds fixed |

Everything is regenerable: fixed seeds, committed inputs
(`reduced_results/reduced_runs.jsonl`, `noise_results/noise_results.json`),
no hand-entered numbers besides the frozen thresholds and the read-only
measured context. This branch does not touch the manuscript.
