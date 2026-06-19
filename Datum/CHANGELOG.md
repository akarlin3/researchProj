# Changelog

All notable changes to Datum are recorded here. Datum follows the monorepo's
provisional-until-the-ruler-locks discipline; see `ASSUMPTIONS.md`.

## [0.1.0] — 2026-06-19 — Substrate swap: Lattice is now primary (task v2)

Lattice (the intended substrate) was merged into the monorepo, so Datum honors its
build-order plan and swaps it in as the primary substrate.

### Changed
- `task.py`: added `TASK_V2` (substrate `lattice`, seed `20260619`) and
  `CURRENT_TASK = TASK_V2`; `TASK_V1` (Gauge) kept frozen for provenance. The
  runner / submission interface / example now default to `CURRENT_TASK`.
- `substrate.py`: implemented `lattice()` (real adapter over `lattice.make_cohort`,
  train/cal/test across an SNR grid for comparability); `gauge_cohort()` demoted to
  the bootstrap substrate (still runnable). `run.py` dispatches on `task.substrate`.
- `manifest.py`: `SUBSTRATE['primary']` is now Lattice (v0.1.0 @ `cbabffe`, seed
  `20260619`); Gauge moved to `SUBSTRATE['bootstrap']`; `planned`/"NOT BUILT"
  removed. `_paths.ensure_deps(names=...)` resolves `lattice` on demand.
- Regenerated reference numbers on Lattice (+ OSIPI external); docs/README/ASSUMPTIONS
  updated. The calibration story is substrate-invariant (NLLS-Gaussian D\* gap −0.089,
  segmented −0.649, MAF −0.038; conformal restores; wall persists for split-conformal,
  CQR/Mondrian recover by inflating width). 28 tests pass.

## [0.1.0] — 2026-06-19 — CP3: submission interface + docs + Casali framing

### Added
- `datum/submit.py` — the standardized submission/scoring interface: `load_task`
  (materialise the task, test truth held out) + `score_submission` (score a new
  method's quantiles on the ruler, with bootstrap CIs, per-D\* tercile coverage, and
  a `vs_reference` ranking against the baselines). Public convention is physical
  `(D, D*, f)`; every result PROVISIONAL-stamped.
- `examples/submit_demo.py` — runnable worked example submitting Gauge's quantile
  regressor (a non-baseline method) and printing its scorecard + ranking.
- `docs/` — `benchmark.md` (the task), `submitting.md` (the contract + example),
  `ruler-as-standard.md` (the Casali differentiation), `release.md` (citable-release
  path, documented, NOT executed). CP3 tests (`test_submit`).

### Notes
- Resource, not a novel-result claim. The Casali framing is positioning: Fashion's
  ruler is a *standard*; Datum is the benchmark that scores any method (including a
  Casali-style one) against it. Numbers remain PROVISIONAL until the ruler locks.

## [0.1.0] — 2026-06-19 — CP2: benchmark build (PROVISIONAL reference numbers)

### Added
- `run.py` — the benchmark runner: drives the curated panel through the ruler on
  the substrate and writes long-form reference numbers + a report. `ci.py` —
  nonparametric bootstrap CIs (none existed in Caliper/Gauge) on the load-bearing
  numbers (marginal coverage-gap & ECE, per-D\* tercile coverage). `convert.py` —
  the Gauge↔Caliper convention adapter (proven exact: forward-model match to 0.0).
- `osipi_fetch.py` — on-demand OSIPI DRO fetch reusing Gauge's pinned DOI/URL/MD5 +
  verifier (Gauge unmodified); `substrate.osipi_dro()` now loads + splits the DRO.
- `results/reference_numbers.csv` + `results/REFERENCE.md` — the PROVISIONAL
  reference numbers (Gauge cohort, all baselines; OSIPI DRO external validation for
  the analytic baselines). `results/osipi_provenance.json` (tracked; raw DRO is not).
- `baselines.py` rewritten as (estimator × calibration) cells wired to Caliper's
  NLLS / segmented-reference / MAF estimators + conformal wrappers. CP2 tests
  (`test_convert`, `test_run`, `test_osipi`).

### Result (PROVISIONAL — scored on Fashion's in-review ruler; no tuning)
- Raw error bars under-cover D\* (NLLS-Gaussian gap −0.07, segmented −0.69, MAF
  −0.05); conformal restores marginal coverage (split −0.02, CQR −0.004, MAF+CQR
  +0.001, Mondrian +0.007). The high-D\* identifiability wall persists for marginal
  split-conformal (0.77) and is recovered by CQR/Mondrian only by inflating width.

## [0.1.0] — 2026-06-19 — CP1: subrepo + README + manifest

### Added
- Datum package scaffolding: `task.py` (frozen `TASK_V1` spec), `substrate.py`
  (read-only Gauge cohort + OSIPI DRO + Lattice swap-in stub), `ruler.py`
  (read-only adapter over `caliper.metrics`), `baselines.py` (curated 7-method
  panel registry), `manifest.py` + `provisional.py` (assumption pins + PROVISIONAL
  stamping), `_paths.py` (sibling-path bootstrap for the read-only deps).
- `ASSUMPTIONS.md` pinning Fashion's ruler (v0.1.0 @ `f078802`, in review at MRM)
  and the substrate (Gauge cohort seed 20260613; OSIPI DRO DOI
  10.5281/zenodo.14605039; Lattice not built).
- `revalidate.py` one-command re-validation; CP1 import/manifest/task tests.

### Notes
- Embedded into the ResearchProj monorepo mirroring Minos: own clean synthetic-only
  history merged with `--allow-unrelated-histories`.
- No reference numbers yet — those are the CP2 deliverable and are PROVISIONAL until
  Fashion's ruler locks.
