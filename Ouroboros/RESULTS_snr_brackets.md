# Finer SNR Brackets Report (Checkpoint 2)

This report presents the refined order-recovery SNR breakdown brackets for temporal fractional order recovery, evaluated on a fine 5 dB grid using 500 noise realizations per cell.

## 1. Recovery Success Criteria
A recovery trial is declared **successful** if the estimated order $\hat{\alpha}_t$ matches the true order $\alpha_t$ exactly ($|\hat{\alpha}_t - \alpha_t| < 1e-5$). We define the **recovery threshold** as a bracket `[highest SNR that fails, lowest SNR that succeeds]`, where success requires a **$\ge 95\%$ success rate** across 500 independent noise realizations.

## 2. Table of Recovery Brackets (5 dB steps)

| True $\alpha_t$ | Pointwise GL Bracket [Fail, Succeed] | Weak-Form GL Bracket [Fail, Succeed] | Order-Dependent Improvement |
| :---: | :---: | :---: | :--- |
| 0.5 | [30 dB, 35 dB] | [15 dB, 20 dB] | Weak-form improves recovery by 15 dB |
| 0.7 | [35 dB, 40 dB] | [15 dB, 20 dB] | Weak-form improves recovery by 20 dB |
| 0.9 | [40 dB, 45 dB] | [15 dB, 20 dB] | Weak-form improves recovery by 25 dB |

> **Note on the $\alpha_t = 0.5$ weak-form bracket.** The original sweep floored this
> cell at 20 dB (0.994 success), so no failing SNR existed to bracket against and the
> cell was previously marked *N/A*. The Checkpoint-1 below-floor extension
> (`ouroboros_cp1_extend.py`; raw data `data/cp1_extended_weak.json`, report
> `RESULTS_snr_brackets_extended.md`) re-ran the identical weak-form pipeline
> (same selector `select_temporal_order_weak_fast`, candidate grid $\{0.2,\dots,1.0\}$
> at 0.1 spacing, 500 realizations, exact-match criterion $|\hat\alpha-\alpha|<1e\text{-}5$,
> $\ge 95\%$ bracket) down to $-5$ dB. It measures a genuine fail edge at **15 dB
> (0.790 success, 105/500 fail $<95\%$)** and **20 dB (0.994)**, yielding the measured
> bracket **[15 dB, 20 dB]** — the same as $\alpha_t = 0.7$ and $0.9$. The apparent
> low-SNR "successes" below 15 dB are non-monotonic in SNR with failed-trial error
> scatter of roughly one candidate-grid step ($\approx 0.1$), i.e. chance grid-snaps,
> **not** genuine recovery. The below-floor rows are reproduced in the
> $\alpha_t = 0.5$ Weak-Form GL section below.

## 3. Full Sweep Results Data

### True $\alpha_t = 0.5$

#### Pointwise GL:
- SNR = 20 dB: Success Rate = 0.000, Error = 0.2000 ± 0.0000
- SNR = 25 dB: Success Rate = 0.000, Error = 0.2000 ± 0.0000
- SNR = 30 dB: Success Rate = 0.000, Error = 0.1000 ± 0.0000
- SNR = 35 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 40 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 45 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 50 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 55 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 60 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000

#### Weak-Form GL:
<!-- Below-floor rows (-5 to 15 dB) from the Checkpoint-1 extension (data/cp1_extended_weak.json); 20 dB and above from the original fine sweep. -->
- SNR = -5 dB: Success Rate = 0.000, Error = 0.2406 ± 0.0617  (chance grid-snap regime; failed-trial scatter 0.2406 ± 0.0617)
- SNR = 0 dB: Success Rate = 0.004, Error = 0.1732 ± 0.0529  (chance grid-snap regime; failed-trial scatter 0.1739 ± 0.0519)
- SNR = 5 dB: Success Rate = 0.142, Error = 0.1092 ± 0.0632  (chance grid-snap regime; failed-trial scatter 0.1273 ± 0.0485)
- SNR = 10 dB: Success Rate = 0.316, Error = 0.0856 ± 0.0761  (chance grid-snap regime; failed-trial scatter 0.1251 ± 0.0593)
- SNR = 15 dB: Success Rate = 0.790, Error = 0.0348 ± 0.0799  (fail edge, <95%; failed-trial scatter 0.1657 ± 0.0934)
- SNR = 20 dB: Success Rate = 0.994, Error = 0.0010 ± 0.0148
- SNR = 25 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 30 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 35 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 40 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 45 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 50 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 55 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 60 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000

### True $\alpha_t = 0.7$

#### Pointwise GL:
- SNR = 20 dB: Success Rate = 0.000, Error = 0.3000 ± 0.0000
- SNR = 25 dB: Success Rate = 0.000, Error = 0.2000 ± 0.0000
- SNR = 30 dB: Success Rate = 0.000, Error = 0.1000 ± 0.0000
- SNR = 35 dB: Success Rate = 0.000, Error = 0.1000 ± 0.0000
- SNR = 40 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 45 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 50 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 55 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 60 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000

#### Weak-Form GL:
- SNR = 0 dB: Success Rate = 0.340, Error = 0.1154 ± 0.1106
- SNR = 5 dB: Success Rate = 0.272, Error = 0.1074 ± 0.0827
- SNR = 10 dB: Success Rate = 0.332, Error = 0.0760 ± 0.0618
- SNR = 15 dB: Success Rate = 0.880, Error = 0.0120 ± 0.0325
- SNR = 20 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 25 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 30 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 35 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 40 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000

### True $\alpha_t = 0.9$

#### Pointwise GL:
- SNR = 40 dB: Success Rate = 0.588, Error = 0.0412 ± 0.0492
- SNR = 45 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 50 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 55 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 60 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 65 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 70 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 75 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 80 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000

#### Weak-Form GL:
- SNR = 0 dB: Success Rate = 0.482, Error = 0.0662 ± 0.0946
- SNR = 5 dB: Success Rate = 0.276, Error = 0.1094 ± 0.0913
- SNR = 10 dB: Success Rate = 0.494, Error = 0.0506 ± 0.0500
- SNR = 15 dB: Success Rate = 0.944, Error = 0.0056 ± 0.0230
- SNR = 20 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 25 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000
- SNR = 30 dB: Success Rate = 1.000, Error = 0.0000 ± 0.0000

* **Plot Citation**: ![Mitigation Comparison Plot](file:///Users/averykarlin/projOuroboros/figures/mitigation_comparison.png)
