# projSentinel — CP0 separation gate: RESULTS

**Verdict: 🔴 RED.** On the mandated Matrix substrate the regret-targeted decision-stop
does **not** robustly halt earlier than a WATCH-style coverage-changepoint alarm. The
pre-registered wedge is refuted. *The separation was the paper — it is not there, so
there is no paper.*

Reproduce: `bash reproduce.sh` (or `python experiments/run_cp0.py`). All numbers below
are from `GLOBAL_SEED = 20260623`, proteus env (numpy 2.3.5), with both stopping rules
calibrated to the **same anytime false-alarm budget** `δ = 0.05`.

## The fair test
The contribution rule and the baseline are placed on equal footing — both **sequential
on the accumulated sequence**, both anytime-`δ` calibrated — so the only difference is
*what each accumulates evidence about*:

| rule | statistic accumulated | calibration | stops? |
|---|---|---|---|
| **regret-stop** (ours) | CUSUM on the regret-targeted monitor `M_k` (stakes-weighted near the decision threshold) | anytime-`δ` running-max quantile (MC on no-drift courses) | yes |
| **WATCH** (baseline) | conformal test martingale on the whole nonconformity-score stream | Ville `c = 1/δ` | yes (alarm) |
| **ACI / conformal-PID** (baseline) | online level `α_t` | — | **never** (recalibrate forever) |

`gap ≡ t_watch − t_regret`. The wedge needs `gap > 0` (regret-stop **earlier**) with a
bootstrap CI excluding 0, ACI never stopping while holding coverage, in the named regime
(coverage held + decision value dead).

## 1. Instrument control (synthetic, dense decision band) — the harness is not blind
| case | t_regret | t_watch | gap | 95% CI | holds | dead | separated |
|---|---|---|---|---|---|---|---|
| dense band, default drift | 2 | 3 | **+1** | [0, +2] | ✓ | ✓ | **False** |

Even an idealized patient with dense mass at the threshold yields a gap whose CI touches
0 once the regret-stop is a *fair* sequential test. (A naive **per-session** regret rule
spuriously "separates" here — that result was an artifact of comparing a single-shot
threshold against a sequential martingale, not of regret-targeting.)

## 2. Matrix verdict (mandated substrate)
**Multi-patient, default regime (band 0.8, rate 0.01):**

| patient | t_regret | t_watch | gap | 95% CI | separated |
|---|---|---|---|---|---|
| 20260623 | 14 | 9 | **−5** | [−1, +12] | False |
| 20260624 | 10 | 13 | +3 | [+3, +13] | True |
| 20260625 | 7 | 10 | +3 | [−2, +11] | False |
| 20260626 | 10 | 10 | 0 | [−2, +10] | False |

The gap **sign flips across patients** (−5 … +3); bootstrap CIs are wide and mostly
straddle 0. Only **1 / 4** patients separates.

**Drift concentration × magnitude sweep (patient 20260623), separation rate 1 / 20 cells.**
Concentrated drift (band 0.3–0.5) → neither monitor fires (`None`/`None`): the
near-threshold density signal is too sparse to accumulate. Strong drift → both fire the
same session (`gap → 0`). The lone separating cell (band 2.0, rate 0.005: gap +2,
CI [0, +6]) sits at the **diffuse** end — the opposite of the regret-targeting
mechanism — and does not survive as a regime.

## 3. Why the wedge dies (the kill mechanism)
A WATCH-style conformal **test martingale** detects the very nonconformity-score-
distribution shift that *drives* decision-value collapse, on a timescale comparable to
(often faster than) the regret CUSUM. The premise that coverage-validity monitoring is
**blind** to decision-value death is false against a sequential martingale: the
near-threshold reported drift that `M` targets necessarily perturbs the score stream
WATCH monitors, so the two fire near-simultaneously. ACI behaves exactly as published
(holds coverage by widening, never stops) — but that only re-confirms it is not the
relevant competitor; WATCH is, and WATCH ties or wins.

## Pre-registered REFUTE — status
- **R-ACI** (ACI reproduces the halt): n/a as a stop — ACI never stops; it holds coverage by widening as designed.
- **R-WATCH** (coverage alarm reproduces/beats the regret-stop): **FIRES** — WATCH ties or beats the regret-stop in 19/20 sweep cells and 3/4 patients.
- **R-REGIME** (named regime absent): the regime is reachable, but the *timing separation within it* is not.

→ **RED.** Reported honestly; not rescued by cherry-picking the favorable patient/corner
or by reverting to the unfair per-session rule.
