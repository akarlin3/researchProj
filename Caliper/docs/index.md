# Caliper docs

Caliper is a small, model-agnostic **uncertainty-quantification calibration
toolkit** for IVIM diffusion MRI fitting. It does three things:

1. **Score** the calibration of any estimator's predicted quantiles
   (`caliper.metrics` — the ruler).
2. **Correct** miscalibration with split-conformal / CQR / Mondrian wrappers
   (`caliper.conformal`).
3. **Generate** PHI-free synthetic IVIM cohorts to exercise the above
   (`caliper.forward`), plus two estimators (`estimator_reference`,
   `estimator_maf`) under one `predict_quantiles` contract.

> **Scope.** This is the un-gated calibration tooling, released for public use
> but **not** as a citable software release. The value-of-information scoring,
> deployment validity monitor, and a citable JOSS/DOI release are deferred and
> gated on a separate publication — see [`../ROADMAP.md`](../ROADMAP.md).

## Install

```bash
pip install -e .                 # numpy-only core (ruler + forward + conformal + reference)
pip install -e ".[estimator]"    # + the MAF posterior (torch)
pip install -e ".[dev]"          # + pytest, ruff
```

## Examples (each: one command, fixed seed)

| Script | What it shows | Needs torch |
|---|---|---|
| `examples/ruler_demo.py` | the model-agnostic ruler on a toy estimator | no |
| `examples/conformal_demo.py` | raw → CQR → Mondrian, the D\* conditional result | no |
| `examples/demo.py` | MAF posterior → split-conformal under SNR shift | **yes** |
| `python -m caliper.benchmark` | the full eval grid → `results/benchmark.csv` | no |
| `examples/benchmark_report.py` | regenerate figures from the CSV | no |

## Reference

- [API reference](api.md) — the public surface of each module.
- [`../README.md`](../README.md) — quickstart, the conformal result, and the
  benchmark summary, with every number traced to the script that produced it.
