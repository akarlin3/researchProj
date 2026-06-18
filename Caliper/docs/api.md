# API reference

The public surface, by module. Full docstrings live in the source; this is the
map. The central data contract everywhere is the quantile array

```
q_pred   : (n, n_params, n_levels)   # from any estimator's predict_quantiles
q_levels : (n_levels,)               # ascending, in (0, 1)
y_true   : (n, n_params)
```

## `caliper.metrics` — the ruler (numpy only)

| Object | Summary |
|---|---|
| `score_quantiles(y_true, q_pred, q_levels, alpha=0.1, ...)` | per-parameter scorecard → `list[ParamScore]` |
| `ParamScore` | dataclass: `coverage`, `coverage_gap`, `ece`, `sharpness`, `mean_pinball`, `mean_interval_score`, `conditional`, `nominal` |
| `format_scorecard(scores, title)` | render a scorecard as a fixed-width table |
| `empirical_coverage(y, lo, hi)` | fraction of `y` inside `[lo, hi]` |
| `ece_quantile(y, q_pred, q_levels)` | expected calibration error over levels |
| `sharpness(lo, hi)` | mean interval width |
| `pinball_loss(y, q_pred, q)` / `interval_score(lo, hi, y, alpha)` | proper scoring primitives |
| `central_interval(q_pred, q_levels, alpha)` | extract a central `(1-alpha)` interval |
| `tercile_groups(x)` | assign each `x` to a tercile `{0,1,2}` |
| `conditional_coverage(y, lo, hi, groups)` | coverage within each group label |

## `caliper.forward` — synthetic IVIM (numpy only)

| Object | Summary |
|---|---|
| `synthetic_cohort(n=4000, snr=50, noise="rician", seed=0)` | reproducible `Cohort` (PHI-free) |
| `Cohort` | dataclass: `params (n,3)`, `bvalues`, `signals_clean`, `signals`, `snr`, `noise` |
| `ivim_signal(bvalues, D, f, Dstar, s0=1.0)` | bi-exponential forward model |
| `sample_params(n, rng)` | draw `(D, f, D*)` from physiological priors |
| `DEFAULT_BVALUES`, `PARAM_NAMES` | the b-schedule and `("D","f","Dstar")` |

## `caliper.conformal` — coverage correction (numpy only)

| Object | Summary |
|---|---|
| `SplitConformalQuantile(q_levels)` | marginal CQR; `.calibrate(q_cal, y_cal).apply(q_test)` |
| `SplitConformalResidual(alpha)` | point ± Q split-conformal; `.calibrate(point_cal, y_cal).apply(point_test) → (lo, hi)` |
| `MondrianConformalQuantile(q_levels)` | group-conditional CQR; `.calibrate(q_cal, y_cal, groups).apply(q_test, groups)` |
| `conditional_coverage_by_strata(y, lo, hi, strata)` | `{stratum → StratumCoverage(n, coverage, mean_width)}` |
| `format_strata_table(per_method, ...)` | render coverage + width per stratum |
| `conformity_scores(q_lo, q_hi, y)` / `conformal_offset(scores, alpha)` | CQR primitives |

Every wrapper also exposes `calibrate_apply(...)` (calibrate on one split,
correct another, in one call).

## `caliper.estimator_reference` — over-confident segmented fit (numpy only)

| Object | Summary |
|---|---|
| `ReferenceIVIMEstimator(bvalues=..., sigma_D=..., sigma_f=..., sigma_Dstar=...)` | closed-form segmented IVIM fit with deliberately narrow reported quantiles |
| `.predict_point(signals) → (n, 3)` | segmented least-squares point estimate |
| `.predict_quantiles(signals, q_levels) → (n, 3, L)` | over-confident Gaussian quantiles about the point |

## `caliper.estimator_maf` — MAF posterior (requires `[estimator]`: torch)

| Object | Summary |
|---|---|
| `MAFPosterior(n_bvalues, epochs=60, seed=0, ...)` | conditional masked-autoregressive-flow posterior over `(D, f, D*)` |
| `.fit(signals, params)` | train by maximum likelihood (seeded) |
| `.predict_quantiles(signals, q_levels) → (n, 3, L)` | posterior quantiles |
| `.posterior_samples(signals) → (n, n_posterior, 3)` | natural-unit posterior draws |

The module imports without torch; constructing `MAFPosterior` without torch
raises a clear `ImportError`.

## `caliper.benchmark` — evaluation harness (numpy only; MAF used iff torch present)

| Object | Summary |
|---|---|
| `run_grid(estimators=None, calibrations=CALIBRATIONS, snrs=..., seeds=..., ...)` | sweep the grid → tidy long-form `list[dict]` |
| `write_csv(rows, path)` | write rows to `results/benchmark.csv` |
| `summarize(rows, param="Dstar")` | render the headline patterns as text |
| `check_reproducible(**grid_kwargs)` | assert fixed seeds reproduce the table |
| `default_estimators()` / `torch_available()` | registry + torch detection |
| `CALIBRATIONS`, `CSV_COLUMNS`, `STRATUM_NAMES` | grid + schema constants |

CLI: `python -m caliper.benchmark [--quick] [--check] [--out PATH]`.
