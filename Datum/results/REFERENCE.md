# Datum reference numbers (PROVISIONAL)

> **PROVISIONAL.** These numbers are scored on Fashion's calibration ruler, which is *in review at MRM*. They are **not** final reference values and must not be cited as such until the ruler locks. Regenerate with `python revalidate.py --full`.

- Ruler: **Fashion calibration ruler v0.1.0 @ f078802** (in review at MRM (R2 revision) -- NOT finalized)
- Ruler implementation: `caliper.metrics` (read-only)
- Substrate: **Gauge synthetic cohort** (seed 20260613), converted to Caliper `(D, f, D*)` convention
- Task: `datum-ivim-calibration` v1; nominal central interval = 0.90; bootstrap CIs: 1000 resamples, 95%
- Honest gate: numbers are reported as run; **no tuning**, no cherry-picking. All baselines reused from Caliper.

## Headline -- D\* marginal coverage gap (nominal 0.90), with 95% bootstrap CI

| baseline | paradigm | D\* coverage | D\* coverage gap [95% CI] |
|---|---|---|---|
| `nlls_gaussian` | parametric | 0.830 | -0.070 [-0.083, -0.057] |
| `nlls_split_conformal` | conformal | 0.880 | -0.020 [-0.032, -0.009] |
| `reference_segmented` | segmented | 0.206 | -0.694 [-0.707, -0.679] |
| `reference_cqr` | conformal | 0.896 | -0.004 [-0.014, +0.007] |
| `reference_mondrian_cqr` | mondrian | 0.907 | +0.007 [-0.003, +0.016] |
| `maf_raw` | flow | 0.854 | -0.046 [-0.059, -0.034] |
| `maf_cqr` | flow | 0.901 | +0.001 [-0.010, +0.011] |

## The identifiability wall -- high-D\* tercile coverage

| baseline | high-D\* coverage [95% CI] | high-D\* mean width |
|---|---|---|
| `nlls_gaussian` | 0.827 [0.804, 0.849] | 190 |
| `nlls_split_conformal` | 0.771 [0.745, 0.797] | 119 |
| `reference_segmented` | 0.057 [0.044, 0.071] | 19.7 |
| `reference_cqr` | 0.915 [0.897, 0.932] | 195 |
| `reference_mondrian_cqr` | 0.915 [0.897, 0.932] | 195 |
| `maf_raw` | 0.851 [0.829, 0.872] | 54.5 |
| `maf_cqr` | 0.884 [0.863, 0.905] | 59.8 |

## External validation -- OSIPI DRO (independent synthetic phantom)

Analytic, b-flexible baselines only (the trained MAF is excluded: the OSIPI DRO is 7-b and its true D\* is shifted out of the Gauge prior). This tests whether the calibration story survives on a phantom we did not generate.

| baseline | D\* coverage | D\* coverage gap [95% CI] |
|---|---|---|
| `nlls_gaussian` | 0.767 | -0.133 [-0.150, -0.118] |
| `nlls_split_conformal` | 0.892 | -0.008 [-0.020, +0.004] |
| `reference_segmented` | 0.022 | -0.878 [-0.884, -0.872] |
| `reference_cqr` | 0.912 | +0.012 [+0.002, +0.024] |
| `reference_mondrian_cqr` | 0.903 | +0.003 [-0.009, +0.014] |

Full long-form numbers (all params x strata x substrate): [`reference_numbers.csv`](reference_numbers.csv).
