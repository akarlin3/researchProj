# Contributing to Caliper

Thanks for your interest. Caliper is a small, deliberately-scoped toolkit; the
bar for contributions is *correctness, reproducibility, and honest reporting*.

## Scope

Caliper ships only the **un-gated** calibration tooling: the model-agnostic
ruler, synthetic IVIM generator, conformal wrappers, two estimators, and the
evaluation harness. Value-of-information scoring, decision-gap analysis, a
deployment validity monitor, and any citable JOSS/DOI release are **deferred and
gated** on a separate publication — please do not add them here (see
[ROADMAP.md](ROADMAP.md)). Safe, in-scope extensions are listed there too.

**Synthetic data only.** No clinical or external datasets enter this repository.

## Development setup

```bash
pip install -e ".[estimator,dev]"   # core + torch MAF + pytest/ruff/matplotlib
```

The core (`metrics`, `forward`, `conformal`, `benchmark`, `estimator_reference`)
is numpy-only; `estimator_maf` needs torch. Tests that require torch auto-skip
when it is absent.

## Before opening a PR

```bash
ruff check .        # lint (line length 100; see pyproject.toml)
pytest -q           # 50 passed with torch; 46 passed + 1 skipped without it
```

Also confirm the examples still run from a clean state:

```bash
python examples/ruler_demo.py
python examples/conformal_demo.py
python -m caliper.benchmark          # writes results/benchmark.csv
python examples/benchmark_report.py  # regenerates examples/figures/
python examples/demo.py              # MAF path (needs torch)
```

## Ground rules

- **Every number in docs comes from a real, named, fixed-seed run.** No
  aspirational or hand-edited metrics. If you change a result, re-run the script
  that produces it and update the figure/table together.
- **Determinism.** New code paths take an explicit seed and reproduce exactly
  (`caliper.benchmark.check_reproducible` is the pattern). Avoid `Date.now`-style
  nondeterminism.
- **Keep the core numpy-only.** Heavy dependencies (torch, etc.) belong behind an
  optional extra and a lazy import, never in `caliper.metrics`/`forward`/
  `conformal`.
- **Tests accompany behavior.** Add or extend a pytest under `tests/` for any new
  public function; prefer asserting invariants over snapshotting noisy numbers.

## License

By contributing you agree your contributions are licensed under the project's
[MIT license](LICENSE).
