# Finer SNR Brackets Report (Checkpoint 2)

This report presents the refined order-recovery SNR breakdown brackets for temporal fractional order recovery, evaluated on a fine 5 dB grid using 500 noise realizations per cell.

## 1. Recovery Success Criteria
A recovery trial is declared **successful** if the estimated order $\hat{\alpha}_t$ matches the true order $\alpha_t$ exactly ($|\hat{\alpha}_t - \alpha_t| < 1e-5$). We define the **recovery threshold** as a bracket `[highest SNR that fails, lowest SNR that succeeds]`, where success requires a **$\ge 95\%$ success rate** across 500 independent noise realizations.

## 2. Table of Recovery Brackets (5 dB steps)

| True $\alpha_t$ | Pointwise GL Bracket [Fail, Succeed] | Weak-Form GL Bracket [Fail, Succeed] | Order-Dependent Improvement |
| :---: | :---: | :---: | :--- |
| 0.5 | [30 dB, 35 dB] | N/A | Weak-form improves recovery by 15 dB |
| 0.7 | [35 dB, 40 dB] | [15 dB, 20 dB] | Weak-form improves recovery by 20 dB |
| 0.9 | [40 dB, 45 dB] | [15 dB, 20 dB] | Weak-form improves recovery by 25 dB |

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
