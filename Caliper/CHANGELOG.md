# Changelog

All notable changes to Caliper are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

No versioned release has been tagged: a citable software release (JOSS / DOI) is
deferred and gated on a separate publication (see [ROADMAP.md](ROADMAP.md)).

## [Unreleased]

### Added
- **Optional, publication-gated reproduction & citation feature (OFF by default).**
  - `caliper.repro_gauge` + `examples/gauge_repro.py` — a numpy-only synthetic
    reproduction of the Gauge conformal-coverage result (marginal CQR restores
    pooled D\* coverage; the high-D\* tercile stays under-covered — the
    identifiability wall; Mondrian buys it back only by inflating width). Reuses
    `caliper.conformal` + `caliper.metrics`; adds no new method. Sits beside the
    existing Fashion reproduction (`caliper.baselines` / `examples/fashion_repro.py`).
  - `caliper.publication` — the single source of truth for the two associated
    manuscripts' (pre-publication) status. A paper is `published` **iff** it has a
    real `paper_doi` (`None` by default), so `publication_enabled()` ships `False`
    and nothing renders as published until a real DOI is filled in. Real Zenodo
    *software* code-archive DOIs are recorded separately and never flip the gate.
  - `CITATION.cff`, `docs/citing.md`, `docs/gauge_reproduction.md` — software +
    pre-publication manuscript citations (`@unpublished` while DOI-less) and the
    Gauge claim → synthetic-result map.
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
