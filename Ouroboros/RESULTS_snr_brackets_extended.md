# Extended Weak-Form SNR Brackets — below the 20 dB floor (Checkpoint 1)

This report extends the weak-form Grünwald–Letnikov order-recovery sweep **below the
20 dB floor** of `RESULTS_snr_brackets.md`, down to **−5 dB**, on the same 5 dB grid
with **500 noise realizations per cell** and the **identical metric/criterion**
(train on noisy weak features at candidate α, score R² against the clean reference
derivative at the true α; success = exact match `|α̂−α| < 1e-5`; bracket at the
**≥95%** success rate).

It does **not** overwrite the original `RESULTS_snr_brackets.md`. Raw data:
`data/cp1_extended_weak.json`, `data/cp1_brackets.json`.

## 1. Headline result

The original report listed the **α=0.5 weak-form bracket as "N/A"** — but only
because the sweep floored at 20 dB, where the success rate was already 0.994 (no
failing cell existed to bracket against). Extending below the floor reveals the
**true** weak-form bracket for α=0.5, and it is the **same** as for 0.7 and 0.9:

| True $\alpha_t$ | Weak-Form bracket [fail, succeed] (≥95%) | Was (original report) |
| :---: | :---: | :---: |
| 0.5 | **[15 dB, 20 dB]** | N/A (floor artifact) |
| 0.7 | **[15 dB, 20 dB]** | [15 dB, 20 dB] |
| 0.9 | **[15 dB, 20 dB]** | [15 dB, 20 dB] |

**Weak-form helps all three orders, to essentially the same ~20 dB threshold.**
The prior basis claim that weak-form "does nothing for α=0.5" is refuted.

## 2. Full low-SNR success curves (500 realizations/cell)

Columns: success rate; mean error over all trials; **scatter of the error among the
trials that did *not* exactly snap** — the diagnostic that separates genuine
recovery (scatter → 0) from chance grid-snaps (scatter large, ~one grid step).

### True $\alpha_t = 0.5$ (weak-form)
| SNR (dB) | Success | Mean err | Failed-trial scatter |
| :---: | :---: | :---: | :---: |
| −5 | 0.000 | 0.2406 ± 0.0617 | 0.2406 ± 0.0617 |
| 0 | 0.004 | 0.1732 ± 0.0529 | 0.1739 ± 0.0519 |
| 5 | 0.142 | 0.1092 ± 0.0632 | 0.1273 ± 0.0485 |
| 10 | 0.316 | 0.0856 ± 0.0761 | 0.1251 ± 0.0593 |
| 15 | 0.790 | 0.0348 ± 0.0799 | 0.1657 ± 0.0934 |
| 20 | 0.994 | 0.0010 ± 0.0148 | 0.1667 ± 0.0943 |

### True $\alpha_t = 0.7$ (weak-form)
| SNR (dB) | Success | Mean err | Failed-trial scatter |
| :---: | :---: | :---: | :---: |
| −5 | 0.070 | 0.2888 ± 0.1208 | 0.3105 ± 0.0946 |
| 0 | 0.340 | 0.1154 ± 0.1106 | 0.1748 ± 0.0901 |
| 5 | 0.272 | 0.1074 ± 0.0827 | 0.1475 ± 0.0590 |
| 10 | 0.332 | 0.0760 ± 0.0618 | 0.1138 ± 0.0378 |
| 15 | 0.880 | 0.0120 ± 0.0325 | 0.1000 ± 0.0000 |
| 20 | 1.000 | 0.0000 ± 0.0000 | — |

### True $\alpha_t = 0.9$ (weak-form)
| SNR (dB) | Success | Mean err | Failed-trial scatter |
| :---: | :---: | :---: | :---: |
| −5 | 0.382 | 0.1608 ± 0.2168 | 0.2602 ± 0.2240 |
| 0 | 0.482 | 0.0662 ± 0.0946 | 0.1278 ± 0.0971 |
| 5 | 0.276 | 0.1094 ± 0.0913 | 0.1511 ± 0.0722 |
| 10 | 0.494 | 0.0506 ± 0.0500 | 0.1000 ± 0.0000 |
| 15 | 0.944 | 0.0056 ± 0.0230 | 0.1000 ± 0.0000 |
| 20 | 1.000 | 0.0000 ± 0.0000 | — |

## 3. The 0–10 dB wiggle is chance grid-snapping, not recovery — CONFIRMED

The prompt asked whether the apparent low-SNR "successes" (e.g. α=0.9 weak showing
0.48 / 0.28 / 0.49 at 0 / 5 / 10 dB) reflect real recovery or chance grid-snaps.
The extended data **confirms chance grid-snaps**:

1. **Non-monotonic in SNR.** A genuine recovery curve rises monotonically with SNR.
   The 0–10 dB region instead wiggles (α=0.9: 0.482 → 0.276 → 0.494; α=0.7: 0.340 →
   0.272 → 0.332). Real, monotone recovery only sets in at 15 → 20 dB.
2. **Large scatter among the non-snapped trials.** Genuine recovery would show the
   error collapsing toward 0 for *all* trials as SNR rises. Instead, at 10 dB the
   failed-trial scatter is still ≈0.10–0.13 (a full candidate grid step), i.e. the
   selector is choosing essentially at random across neighbouring orders and lands
   on the true one ~⅓–½ of the time by coincidence.
3. **α-ordering of the snap rate matches A(α), not signal.** At low SNR the snap
   rate is *higher* for α=0.9 (≈0.48) than α=0.5 (≈0.00–0.32): the higher-order
   target, with larger derivative amplitude, scatters its argmax more broadly over
   the upper candidates and thus snaps onto 0.9 more often — a sampling effect, not
   recovery.

**Conclusion:** the only real, monotone recovery threshold for weak-form GL is the
**[15 dB, 20 dB]** bracket, common to all three orders. Everything reported below
15 dB in either sweep is coincidence-level and must not be read as recovery.
