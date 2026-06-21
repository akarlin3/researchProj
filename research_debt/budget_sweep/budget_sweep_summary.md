# Simulation-budget sweep — re-derived summary (HC4/CS3, BANKED)

Re-derived from `cp1_budget_sweep.csv` (prior run product). D\* claimed SD-to-CRLB ratio (median) by SNR, and the overall D\* below-floor fraction (mean over SNR), per training budget.

| budget | SNR10 | SNR20 | SNR50 | SNR100 | D* below-floor |
|---|---|---|---|---|---|
|  50,000 | 0.077 | 0.156 | 0.362 | 0.740 | 0.695 |
| 100,000 | 0.080 | 0.160 | 0.361 | 0.704 | 0.696 |
| 250,000 | 0.079 | 0.159 | 0.380 | 0.716 | 0.699 |
| 500,000 | 0.084 | 0.160 | 0.376 | 0.671 | 0.695 |
| 1,000,000 | 0.081 | 0.160 | 0.368 | 0.692 | 0.695 |

- SNR10 claimed-ratio span across 20× budget range: [0.077, 0.084]
- Overall D\* below-floor span: [0.695, 0.699]

**Verdict:** BUDGET-INVARIANT: the below-floor D* overconfidence is flat across a 20-fold budget range -> NOT simulation starvation; supports the information-theoretic reading.
