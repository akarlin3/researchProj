# Datum reference numbers (PROVISIONAL)

> **PROVISIONAL.** These numbers are scored on Fashion's calibration ruler, which is *in review at NMR in Biomedicine* (retooled, boundary-railing-first; a scoped, ground-truth-only secondary reported under the honest CRLB). They are **not** final reference values and must not be cited as such until the ruler locks. Regenerate with `python revalidate.py --full`.

- Ruler: **Fashion calibration ruler v0.1.0 @ f078802** (in review at NMR in Biomedicine (retooled, boundary-railing-first) -- NOT finalized)
- Ruler implementation: `caliper.metrics` (read-only)
- Substrate: **Lattice IVIM DRO** (seed 20260619), converted to Caliper `(D, f, D*)` convention
- Task: `datum-ivim-calibration` v2; nominal central interval = 0.90; bootstrap CIs: 1000 resamples, 95%
- Honest gate: numbers are reported as run; **no tuning**, no cherry-picking. All baselines reused from Caliper.

## Headline -- D\* marginal coverage gap (nominal 0.90), with 95% bootstrap CI

| baseline | paradigm | D\* coverage | D\* coverage gap [95% CI] |
|---|---|---|---|
| `nlls_gaussian` | parametric | 0.811 | -0.089 [-0.103, -0.076] |
| `nlls_split_conformal` | conformal | 0.905 | +0.005 [-0.006, +0.016] |
| `reference_segmented` | segmented | 0.251 | -0.649 [-0.663, -0.632] |
| `reference_cqr` | conformal | 0.893 | -0.007 [-0.018, +0.004] |
| `reference_mondrian_cqr` | mondrian | 0.897 | -0.003 [-0.014, +0.008] |
| `maf_raw` | flow | 0.862 | -0.038 [-0.050, -0.025] |
| `maf_cqr` | flow | 0.899 | -0.001 [-0.012, +0.011] |

## The identifiability wall -- high-D\* tercile coverage

| baseline | high-D\* coverage [95% CI] | high-D\* mean width |
|---|---|---|
| `nlls_gaussian` | 0.807 [0.782, 0.833] | 164 |
| `nlls_split_conformal` | 0.798 [0.773, 0.823] | 124 |
| `reference_segmented` | 0.078 [0.062, 0.096] | 19.7 |
| `reference_cqr` | 0.889 [0.870, 0.908] | 184 |
| `reference_mondrian_cqr` | 0.906 [0.888, 0.923] | 188 |
| `maf_raw` | 0.887 [0.867, 0.907] | 61.3 |
| `maf_cqr` | 0.908 [0.890, 0.927] | 64.6 |

## External validation -- OSIPI DRO (independent synthetic phantom)

Analytic, b-flexible baselines only (the trained MAF is excluded: the OSIPI DRO is 7-b and its true D\* is shifted out of the Gauge prior). This tests whether the calibration story survives on a phantom we did not generate.

| baseline | D\* coverage | D\* coverage gap [95% CI] |
|---|---|---|
| `nlls_gaussian` | 0.770 | -0.130 [-0.146, -0.114] |
| `nlls_split_conformal` | 0.898 | -0.002 [-0.014, +0.009] |
| `reference_segmented` | 0.020 | -0.880 [-0.886, -0.875] |
| `reference_cqr` | 0.912 | +0.012 [+0.000, +0.024] |
| `reference_mondrian_cqr` | 0.923 | +0.023 [+0.012, +0.033] |

Full long-form numbers (all params x strata x substrate): [`reference_numbers.csv`](reference_numbers.csv).
