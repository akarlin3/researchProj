# Changelog

All notable changes to Caliper are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

No versioned release has been tagged: a citable software release (JOSS / DOI) is
deferred and gated on a separate publication (see [ROADMAP.md](ROADMAP.md)).

## [Unreleased]

### Added
- `caliper.benchmark` — a reproducible evaluation harness sweeping
  `{estimator} × {raw, split, CQR, Mondrian} × {SNR} × {seed}`, scoring each
  cell with the ruler plus D\*-tercile conditional coverage and mean width. Emits
  a tidy long-form `results/benchmark.csv`; fixed seeds reproduce it exactly.
- `examples/benchmark_report.py` — regenerates the benchmark figures *solely*
  from `results/benchmark.csv` (calibration-vs-SNR, coverage-by-D\*-tercile,
  Mondrian width cost).
- `examples/ruler_demo.py` — a numpy-only quickstart for the model-agnostic
  calibration ruler.
- `docs/` — plain-Markdown index and API reference.
- `CHANGELOG.md`, `CONTRIBUTING.md`; packaging metadata (classifiers, project
  URLs); CI estimator (torch) job alongside the numpy matrix.

### Fixed
- `caliper.estimator_maf` now imports cleanly without torch (the `@torch.no_grad()`
  decorator was evaluated at class-definition time and raised `NameError`; it is
  now an in-method `with torch.no_grad():` context). Constructing `MAFPosterior`
  without torch still raises a clear `ImportError`.

### Notes
- Initial toolkit (model-agnostic ruler `caliper.metrics`, synthetic IVIM
  `caliper.forward`, conformal wrappers `caliper.conformal`, the over-confident
  `estimator_reference`, and the MAF posterior `estimator_maf`) predates this
  changelog.
