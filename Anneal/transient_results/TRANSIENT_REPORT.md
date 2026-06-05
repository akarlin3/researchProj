# TRANSIENT_REPORT — deterministic-transient tests of the A=0.5 chimera

_Generated 2026-06-05T01:50:49.279327+00:00 · config `tools/transient-tests/transient.config.json`_

**Hypothesis (Avery's):** at A=0.5 the chimera is a *deterministic transient* — the breathing trajectory spirals outward, each high-R excursion closer to sync capture. Predictions: (1) per-cycle breath maxima M_k ratchet upward within runs; (2) lifetime is (near-)deterministically set by the collective initial condition, far beyond the manifold probe's null on D₀. Nulls reported straight.

Peak detection: PR #41 canonical detector (`breath.mjs detectPeaks`) on the smoothed pre-absorption min(R₁,R₂). ICs recomputed from the logged seed (`seedChimera`, zero integration). Determinism gate (re-trace vs logged labels): graze 160/160, abs-censor 160/160 → **PASS**.

## CP0 — coverage (light audit)

| N | A | trace records | with full R_incoh | absorbed/traced | campaign total | persistent |
|---|---|---|---|---|---|---|
| 8 | 0.2 | 100 | 25 | 25 | 100 | 75 |
| 16 | 0.2 | 100 | 22 | 22 | 100 | 78 |
| 32 | 0.2 | 100 | 29 | 29 | 100 | 71 |
| 64 | 0.2 | 100 | 27 | 27 | 100 | 73 |
| 8 | 0.5 | 100 | 100 | 100 | 200 | 0 |
| 16 | 0.5 | 100 | 100 | 100 | 200 | 0 |
| 32 | 0.5 | 100 | 100 | 100 | 200 | 0 |
| 64 | 0.5 | 100 | 100 | 100 | 200 | 0 |

A=0.5 traces exist only for N∈{8,16,32,64} (the phase subset; N∈{4,24,48} were not trace-dumped, so CP1's per-cycle envelope is reported on those four points). All A=0.5 traced runs absorb. A=0.2 persistent runs carry no on-disk trace and were re-traced (160 runs, 40/point) for the contrast.

## CP1 — the ratchet test (A=0.5)

A run is *testable* for monotonicity if it has ≥3 detected peaks (≥2 increments). It *ratchets* (lenient) if M_k is non-decreasing allowing ≤1 single-step inversion(s); *strict* allows none. ⟨ΔM⟩ is the grand mean per-cycle increment with a cluster (over-runs) bootstrap 95% CI — the robust statistic, immune to the short-sequence leniency of the monotone fraction.

| N | runs | testable | ratchet frac (≤1 inv) | strict (0 inv) | ⟨ΔM⟩ | 95% CI | median cycles |
|---|---|---|---|---|---|---|---|
| 8 | 100 | 83 | 0.65 | 0.18 | +0.0314 | [+0.0263, +0.0365] | 4 |
| 16 | 100 | 73 | 0.84 | 0.53 | +0.0518 | [+0.0459, +0.0577] | 3 |
| 32 | 100 | 63 | 1.00 | 0.89 | +0.0783 | [+0.0732, +0.0836] | 2 |
| 64 | 100 | 80 | 1.00 | 1.00 | +0.0913 | [+0.0884, +0.0941] | 2 |

**Spiral-out fit** — per-run linear fit of log(θ − M_k) vs cycle index k over the sub-θ *approach* peaks (θ=0.85; capturing/grazing peaks at or above θ are dropped, not the run). A clean outward spiral gives a negative slope (exponential approach to the boundary). Runs with <3 sub-θ peaks are counted apart.

| N | runs fit | median slope | frac slope<0 | runs <3 sub-θ peaks |
|---|---|---|---|---|
| 8 | 26 | -0.5479 | 1.00 | 74 |
| 16 | 24 | -0.4973 | 1.00 | 76 |
| 32 | 30 | -0.5883 | 1.00 | 70 |
| 64 | 33 | -0.7671 | 1.00 | 67 |

## CP1 — A=0.2 contrast (the never-absorbers)

Same M_k extraction over the full t_max window on re-traced persistent runs. Under the hypothesis these should be **bounded/stationary — no ratchet**: mean increment ≈ 0 and the per-run slope of M_k vs k ≈ 0.

| N | runs | median cycles | ⟨ΔM⟩ | 95% CI | median per-run slope (M_k~k) | ratchet frac |
|---|---|---|---|---|---|---|
| 8 | 40 | 227 | -0.00067 | [-0.00079, -0.00054] | -1.51e-04 | 0.00 |
| 16 | 40 | 202 | -0.00078 | [-0.00085, -0.00071] | -2.19e-04 | 0.00 |
| 32 | 40 | 181 | -0.00097 | [-0.00112, -0.00084] | -2.67e-04 | 0.00 |
| 64 | 40 | 79 | -0.00191 | [-0.00254, -0.00147] | -1.04e-03 | 0.00 |

**Artifact control.** The same detector and pre-event windowing produce ⟨ΔM⟩≈0 here (persistent runs, hundreds of cycles) but a strongly positive ⟨ΔM⟩ at A=0.5. The A=0.5 upward drift is therefore not an artifact of the peak detector or of slicing on the pre-absorption prefix — it is a property of the absorbing dynamics.

## CP2 — collective-IC predictability of log t_abs

Target = log t_abs. Per-N: Spearman ρ for each single t=0 feature; 5-fold CV R² of OLS on the t=0 collective features alone (linear and quadratic+interactions); and — on the traced subset where the first-cycle features exist — CV R² adding M₁ and T_b,1. D₀ null = the manifold probe's |ρ| (read-only benchmark).

### Single-feature Spearman ρ (t=0 collective features) vs the D₀ null

| N | n | ρ(R_incoh0) | ρ(R_sync0) | ρ(\|Δφ0\|) | D₀ null \|ρ\| |
|---|---|---|---|---|---|
| 4 | 200 | -0.109 | -0.185 | +0.032 | 0.081 |
| 8 | 200 | -0.135 | -0.099 | -0.445 | 0.041 |
| 16 | 200 | -0.010 | -0.003 | -0.428 | 0.039 |
| 24 | 200 | -0.051 | -0.063 | -0.507 | 0.072 |
| 32 | 200 | +0.061 | +0.106 | -0.582 | 0.049 |
| 48 | 200 | -0.086 | +0.056 | -0.594 | 0.072 |
| 64 | 200 | -0.058 | -0.011 | -0.538 | 0.024 |

### CV R² — collective IC predicting log t_abs

| N | n (all) | R² t=0 linear | R² t=0 quad | n (traced) | R² t=0 (traced) | R² t=0+1st-cycle | ΔR² |
|---|---|---|---|---|---|---|---|
| 4 | 200 | 0.017 | -0.053 | — | — | — | — |
| 8 | 200 | 0.195 | 0.184 | 96 | 0.145 | 0.255 | +0.109 |
| 16 | 200 | 0.149 | 0.237 | 96 | 0.103 | 0.487 | +0.384 |
| 24 | 200 | 0.212 | 0.223 | — | — | — | — |
| 32 | 200 | 0.296 | 0.340 | 87 | 0.229 | 0.356 | +0.127 |
| 48 | 200 | 0.299 | 0.438 | — | — | — | — |
| 64 | 200 | 0.200 | 0.242 | 97 | 0.105 | -0.008 | -0.114 |

R² is out-of-fold (5-fold). The deterministic ceiling: the dynamics is fully deterministic given the *complete* IC (all 2N phases) plus constants; these collective summaries are a lossy projection of that IC, so residual unpredictability mixes genuine constants-influence with information thrown away by the projection — R² here is a **lower bound** on collective-IC determinism, not the deterministic ceiling itself.

## CP3 — deciding figure

`transient_decider.png` / `.pdf`. Pooled out-of-fold R² (t=0 collective IC, all N) = 0.207.

## Verdicts

- **Prediction (1): per-cycle ratchet — YES (ratchet supported).** A=0.5 mean ratcheting fraction 0.87; mean per-cycle increment ⟨ΔM⟩ over N = +0.0632, with a strictly-positive bootstrap 95% CI at **4/4** N-points. Both the increment and the ratchet fraction grow with N (the spiral tightens for larger populations). The spiral-out fit corroborates this where there are ≥3 sub-θ approach peaks (see table); short, fast captures at large N leave few sub-θ peaks to fit.
- **A=0.2 contrast:** ⟨ΔM⟩ = -0.00108 over a median of 192 cycles, ratchet fraction 0.00 at every N — **stationary, no ratchet**, exactly as the hypothesis predicts for the never-absorbers. Same detector and windowing as A=0.5, so the A=0.5 drift is not a detector artifact.
- **Prediction (2): collective-IC predictability.** The single feature **|Δφ₀|** (initial inter-population phase difference) carries it: Spearman |ρ| median 0.52 across N≥8 (monotone, negative — larger initial phase separation ⇒ shorter life), versus the manifold-probe D₀ null |ρ|≤0.081. CV R² of the t=0 collective model: linear ⟨0.196⟩, quadratic ⟨0.230⟩ (best 0.438) — i.e. **R² ≈ 0.23 versus the null's R² ≈ 0.007**, two orders of magnitude above the manifold-distance null.
- **First-cycle features (M₁, T_b,1):** ⟨ΔR²⟩ = +0.127 on the traced subset, but inconsistent across N (range -0.114…+0.384); they help at some N and not others. No clean 'set up within the first breath' signal beyond what t=0 already provides.
- **Predictability verdict:** collective IC (driven by |Δφ₀|) predicts log t_abs **far beyond the D₀ null** — the transient is collective-IC-organized — but the absolute R² (~0.2–0.44) means these low-dimensional summaries fix only part of the lifetime; the rest lives in the finer IC structure these projections discard. Supported in the ≫-null sense the prompt's discriminator asks for, not as full determinism from collective summaries alone.

### Overall
**Prediction (1) [ratchet]: YES (ratchet supported).** **Prediction (2) [collective-IC predictability]: SUPPORTED relative to the D₀ null, partial in absolute terms** — |Δφ₀| is a strong, clean, previously-unreported predictor of lifetime. The A=0.2 stationary contrast holds and rules out a detector artifact. Refuted sub-claims: the first-cycle features add no consistent predictive power, and at large N captures are too fast (≤2–3 cycles) to exhibit a long resolvable spiral — the ratchet there is real in ⟨ΔM⟩ but compressed into few cycles. Nulls reported straight.

