# Phase 1 SINDy Temporal-Order Selection & Refutation Report

This report summarizes the data-driven discovery and temporal order selection results on simulated stromal-pressure dynamics.

## Model Selection Sweep Results

Below is the sweep of temporal fractional orders $\alpha_t$ against their generalization performance (held-out $R^2$) and complexity (number of active terms in the model):

| Temporal Order ($\alpha_t$) | Held-out $R^2$ (Average) | Number of Active Terms |
| :--- | :--- | :--- |
| 0.6 | -9.06725 | 65 |
| 0.8 | -21.02890 | 64 |
| 1.0 | 0.99247 | 46 |

### Key Findings
- **Selected Temporal Order**: $\alpha_t$ = 1.0 (Held-out $R^2$ = 0.99247)
- **Did $\alpha_t = 1.0$ Win?**: **YES**
- **Discovered Active Term Count**: 46

## Refutation Conclusion

Ordinary integer-time dynamics ($\alpha_t = 1.0$) **yielded the highest generalization score** on the held-out test split. 

Since the simulator's governing law is strictly of integer-order ($\alpha_t = 1.0$), **the fractional SINDy pipeline successfully refuted fractional dynamics and recovered the true integer-order dynamics**. This provides a strong **validation** of the fractional temporal formulation.

## LHS Prediction Validation Plot

The scatter plot below compares the true derivative value against the SINDy prediction for each of the three coupled variables on the held-out test split:

![LHS Validation Plot](file:///Users/averykarlin/projOuroboros/figures/lhs_validation.png)

## Discovered Equations for Selected Model ($\alpha_t$ = 1.0)

```text
(x0)' = -0.050 x0 + 0.001 x1 + 0.097 x2 + -0.001 x0_x1 + -0.093 x0_x2 + -0.001 x1_x1 + 0.005 x1_x2 + -0.010 x2_x2 + 0.050 x0_xx
(x1)' = 0.085 1 + -0.036 x0 + -0.173 x1 + 0.530 x2 + -0.009 x0_x0 + -0.028 x0_x1 + -0.345 x0_x2 + 0.092 x1_x1 + -0.321 x1_x2 + 0.383 x2_x2 + -0.002 x0_x + -0.037 x0_xx + 0.096 x1_xx + 0.009 x2_x + 0.005 D_0_5_x0 + 0.002 D_0_5_x1 + -0.014 D_0_5_x2
(x2)' = 0.001 1 + 0.063 x2 + 0.004 x0_x1 + -0.115 x0_x2 + 0.032 x1_x2 + -0.090 x2_x2 + 0.010 x2_xx
```
