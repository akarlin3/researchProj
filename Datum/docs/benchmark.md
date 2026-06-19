# The Datum benchmark task

Datum scores one thing: **is an IVIM method's reported uncertainty honest?** A
method that reports an error bar (or quantiles) on the IVIM parameters `(D, D*, f)`
is judged by whether those intervals actually cover the truth at their nominal
level -- and how wide they have to be to do so.

> **PROVISIONAL.** Every Datum number is scored on **Fashion's calibration ruler**,
> which is *in review at MRM*. Numbers are not final until the ruler locks (see
> [`../ASSUMPTIONS.md`](../ASSUMPTIONS.md)). Regenerate with `python -m datum.run`.

## The fixed task (`datum/task.py::TASK_V1`)

| field | value |
|---|---|
| substrate | Lattice IVIM DRO (`lattice.make_cohort`, seed `20260619`); Gauge cohort as bootstrap; OSIPI DRO for external validation |
| splits | train 3000 / cal 2000 / test 3000 |
| signal | bi-exponential IVIM, 22-b scheme (0–800 s/mm²), Rician noise |
| prediction | quantiles `(n, 3, L)` for `(D, D*, f)` at `quantile_levels` (13 levels) |
| headline interval | central 0.90 (`alpha = 0.10`) — a Fashion nominal level |
| metrics (the ruler) | coverage, coverage-gap, ECE, sharpness, pinball, interval-score |
| conditioning | D\* terciles (the identifiability-wall regime) |
| CIs | 95% nonparametric bootstrap (1000 resamples) on the load-bearing numbers |

The task is **frozen and versioned**: changing the substrate, levels, or metric set
is a new task version and triggers re-validation of every reference number.

## Metrics — Fashion's ruler, via Caliper

Scoring goes through `caliper.metrics.score_quantiles` (read-only), which packages
Fashion's coverage / ECE / sharpness recipe model-agnostically:

- **coverage(L)** — fraction of voxels whose nominal-L interval contains the truth.
  Calibrated ⇔ coverage ≈ L. **coverage-gap** = coverage − nominal (the headline).
- **ECE** — mean |coverage(L) − L| across levels (0 = perfect quantile calibration).
- **sharpness / interval-score / pinball** — width and proper-scoring penalties;
  calibration is cheap with huge intervals, so width is always reported alongside.

All scored **marginally and per-D\* tercile**. The high-D\* tercile is where the
identifiability wall lives.

## The substrate (Lattice primary; Gauge bootstrap; OSIPI external)

The primary substrate is the **Lattice IVIM DRO** (`lattice.make_cohort`, read-only)
— the intended substrate, now built; task **v2** swapped it in for the **v1 Gauge
bootstrap** (`gauge.cohort.generate_cohort`), which is kept runnable as a
cross-check. Lattice's convention matches Gauge's exactly, so `datum.convert` is
unchanged. The **OSIPI DRO** (an independent synthetic phantom, DOI
`10.5281/zenodo.14605039`) is the external-validation substrate via
`python -m datum.osipi_fetch` (download-on-demand, provenance-tracked, raw arrays
git-ignored). Everything is synthetic and PHI-free by construction.

## Convention

Internally the curated baselines run in Caliper's `(D, f, D*)` / 1e-3 convention
(the substrate cohort is converted once, proven exact in `tests/test_convert.py`). The
**public submission interface** speaks the IVIM-natural physical convention
`(D, D*, f)` — submitters never touch the internal convention, and coverage / gap /
ECE are scale-invariant so the comparison is sound.

## Reference numbers

See [`../results/REFERENCE.md`](../results/REFERENCE.md) and
[`../results/reference_numbers.csv`](../results/reference_numbers.csv). The curated
panel reproduces the IVIM calibration story: parametric/segmented/learned error bars
under-cover D\*; conformal restores marginal coverage; the high-D\* wall persists for
marginal conformal and is recovered by CQR/Mondrian only by inflating width.
