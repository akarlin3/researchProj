# RESULTS — Minos-Core

Every number below was printed by `python experiments/run_all.py` in this session
(seed `20240517`, `n_voxels = 1_000_000`). Re-running from the clean seed reproduces them
exactly. No number here was hand-entered from anywhere but that run's stdout.

**Config.** `t1=0.0, t2=2.0, k_under=2.0, k_over=1.0, s=0.5`; prior weights
`(0.35, 0.30, 0.35)`, means `(-1.0, 1.0, 3.0)`, std `1.0`; shift `alpha=0.5, beta=5.0`;
gate `q=0.995` (threshold `g*=2.5758`); `delta_test=1.5`.

## GATE 1 — decision core
| quantity | value |
|---|---|
| `EU_point` | `-0.059371` |
| `EU_posterior` | `-0.055138` |
| `EU_oracle` | `+0.000000` |
| EVPI-analog (`EU_oracle - EU_posterior`) | `+0.055138` |
| EVPI-analog, degenerate `s=1e-4` | `2.428e-09` |

Ordering `EU_oracle ≥ EU_posterior ≥ EU_point` holds, and the EVPI-analog collapses to
≈0 (`2.4e-9`) as the posterior approaches a point mass at `theta_true`. The oracle utility
is exactly zero by construction, so the EVPI-analog equals the posterior policy's expected
regret.

## GATE 2 — Value of Calibration
| quantity | value |
|---|---|
| `VoC(tau=1)` | `+0.000e+00` |
| `VoC(tau=0.5)` (2× over-confident) | `+0.001169` |
| `VoC(tau=2.0)` (2× under-confident) | `+0.004361` |
| `argmin VoC` | `tau = 1.00` |
| value of using the error bar (`EU_post(1) - EU_point`) | `+0.004233` |
| central-interval coverage @ `tau=1` (nominal 0.5 / 0.8 / 0.9) | `0.4997 / 0.8001 / 0.9003` |
| ECE @ `tau=1.0` / `0.6` / `1.5` | `0.00015 / 0.21483 / 0.14004` |

`VoC` is exactly 0 at `tau=1`, strictly positive away from it, and minimised at `tau=1`.
Under-confidence costs more decision utility than equal-factor over-confidence
(`0.004361 > 0.001169`) because over-hedging escalates without bound while under-hedging
saturates at the point policy. The calibrated error bar is itself worth `+0.004233` of utility
over ignoring it. Coverage at `tau=1` matches nominal to 3 decimals and ECE is ~`1.5e-4`;
mis-scaling the error bar drives ECE up by 2–3 orders of magnitude.

## GATE 3 — trust-gate
| quantity | value |
|---|---|
| `VoTG(delta=0)` | `-0.006471` |
| `VoTG(delta=1.5)` | `+0.121806` |
| detection `AUC(delta=0)` | `0.4995` |
| detection `AUC(delta=1.5)` | `0.8555` |
| posterior regret @ shift | `+2.278442` |
| gated regret @ shift | `+2.156635` |

Without shift the gate is near-neutral (`VoTG ≈ -0.0065`, the 0.5% false-positive cost) and
its detector is at chance (`AUC ≈ 0.50`). Under the deployment shift the over-confident,
downward-biased posterior incurs large regret (`2.278`); detecting the shift and escalating
the flagged voxels recovers `+0.1218` of utility (gated regret `2.157 < 2.278`), with detection
`AUC = 0.8555` (matching the analytic `Φ(delta/√2)`).

## Figures (vector PDF, `figures/`)
- **`fig_a_regret_vs_tau.pdf`** — (a) decision regret of the posterior policy vs calibration
  quality `tau`, against the point-policy baseline; minimised at `tau=1`.
- **`fig_b_voc_evpi.pdf`** — (b) `VoC(tau)` (zero and minimal at `tau=1`, growing both ways)
  with the EVPI-analog (`0.0551`) drawn as the irreducible-regret reference.
- **`fig_c_gate_roc_votg.pdf`** — (c) trust-gate ROC at `delta=1.5` (AUC `0.856`, operating
  point at `q=0.995`) and `VoTG(delta)` showing the near-zero hold at `delta=0`, a small
  break-even dip, then growth.
- **`fig_d_utility_bars.pdf`** — (d) expected utility of `{point, posterior, gated, oracle}`
  in-distribution and under shift; the gated bar beats the posterior bar only under shift.

## What the toy demonstrates
Calibration of the per-voxel error bar has **positive decision value in its own right**: at
fixed point estimate, scaling the reported uncertainty away from its calibrated width strictly
lowers expected utility (`VoC > 0`), and the calibrated bar is worth measurable utility over
ignoring uncertainty entirely. Separately, when a deployment shift makes the reported
uncertainty *untrustworthy* (over-confident and biased), a trust-gate that detects the shift
and acts conservatively **recovers utility** (`VoTG > 0` past a break-even shift). Both
quantities price the error bar — its width (VoC) and its trustworthiness (VoTG) — not the point
prediction, which is what distinguishes Minos from decision-curve net-benefit and from
population EVPI/EVPPI (see `POSITIONING.md`).
