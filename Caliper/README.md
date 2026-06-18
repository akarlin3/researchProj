# Caliper

**An IVIM uncertainty-quantification calibration toolkit.**

Caliper is a small, reviewer-oriented Python package for *measuring and
correcting* the calibration of uncertainty estimates in intravoxel incoherent
motion (IVIM) diffusion MRI. It ships three composable pieces:

1. **`caliper.metrics`** — a **model-agnostic calibration ruler** (numpy-only):
   coverage, quantile ECE, sharpness, pinball/interval score, and
   group-conditional coverage. It knows nothing about IVIM or any particular
   estimator — feed it true values and predicted quantiles.
2. **`caliper.estimator_maf`** — a conditional **masked-autoregressive-flow**
   posterior over `(D, f, D*)` given the multi-b signal decay (optional torch
   extra). **`caliper.estimator_reference`** is a torch-free stand-in: an
   over-confident segmented-fit IVIM estimator exposing the same
   `predict_quantiles` contract, so the calibration story runs on numpy alone.
3. **`caliper.conformal`** — **split-conformal**, **CQR**, and **Mondrian
   (group-conditional) CQR** wrappers that coverage-correct *any* estimator
   exposing `predict_quantiles`, plus `conditional_coverage_by_strata` for
   reading coverage *and width* per stratum.

All data is **synthetic and PHI-free**, generated in-repo with fixed seeds
(`caliper.forward`). There are no clinical-data dependencies.

> Every number in this README is produced by a fixed-seed script in this repo:
> the MAF results by `examples/demo.py` (torch), and the **Conformal
> calibration** numbers below by `examples/conformal_demo.py` (numpy only).
> Re-run them to reproduce exactly.

---

## Install

```bash
# core ruler + forward model + conformal wrapper (numpy only)
pip install -e .

# add the MAF estimator (pulls in torch)
pip install -e ".[estimator]"
```

Python 3.10–3.12. The core (`metrics`, `forward`, `conformal`) is numpy-only;
torch is required only for `estimator_maf`.

---

## Quickstart

```bash
python examples/demo.py
```

This runs the full pipeline — synthetic IVIM → MAF posterior → split-conformal →
calibration scorecard — under a realistic **deployment shift** (the flow is
trained for high-SNR fitting but evaluated at lower SNR, where model-based IVIM
UQ is known to be over-confident).

### The model-agnostic ruler API

```python
import numpy as np
from caliper import metrics as M

# y_true:   (n, n_params)
# q_pred:   (n, n_params, n_levels)   <- from any estimator's predict_quantiles
# q_levels: (n_levels,)               ascending in (0, 1)
scores = M.score_quantiles(y_true, q_pred, q_levels, alpha=0.10,
                           param_names=["D", "f", "Dstar"],
                           conditioning=y_true)      # tercile-conditional probe
print(M.format_scorecard(scores))
# each ParamScore exposes: coverage, coverage_gap, ece, sharpness,
# mean_pinball, mean_interval_score, conditional (per-tercile coverage)
```

### Conformal wrapper over any estimator

```python
from caliper.conformal import SplitConformalQuantile

cq = SplitConformalQuantile(q_levels).calibrate(q_cal, y_cal)
q_corrected = cq.apply(q_test)   # coverage-corrected quantiles, same shape
```

---

## Results (from `examples/demo.py`)

Nominal central coverage **0.900** (90% intervals, α = 0.10). Held-out synthetic
test set, deployment shift train-SNR 60 → test-SNR 25.

### Raw MAF is over-confident (the known model-based UQ result)

| param | coverage | gap | ECE | sharpness |
|-------|---------:|------:|------:|----------:|
| D     | 0.528 | −0.372 | 0.139 | 0.281 |
| f     | 0.550 | −0.350 | 0.133 | 0.072 |
| D\*   | 0.555 | −0.345 | 0.124 | 26.86 |

The raw posterior intervals are far too tight: ~53–56% empirical coverage
against a 90% target. This is expected and is reported honestly — **not** tuned.

### Split-conformal restores **marginal** coverage

| param | raw coverage | conformal coverage | raw \|gap\| | conformal \|gap\| |
|-------|-------------:|-------------------:|------------:|------------------:|
| D     | 0.528 | 0.890 | 0.372 | **0.010** |
| f     | 0.550 | 0.876 | 0.350 | **0.024** |
| D\*   | 0.555 | 0.903 | 0.345 | **0.003** |

Marginal coverage is restored to within ≤0.024 of nominal for every parameter.

### Honest caveat: **conditional** coverage is *not* restored

Conformal applies a single marginal offset, so it cannot fix coverage that
varies across the parameter range. Post-conformal conditional coverage by
true-D\* tercile:

```
   Dstar  g0(low)=0.972  g1(mid)=0.929  g2(high)=0.810
```

The high-D\* tercile still under-covers (0.810 vs 0.900) while low-D\*
over-covers (0.972). This is the **irreducible identifiability limit** of IVIM
`D*` — pseudo-diffusion is weakly constrained by the signal — and it is a
property of the data, not a bug in the wrapper. Caliper's job is to *measure*
this faithfully, which it does.

---

## Conformal calibration (torch-free — `examples/conformal_demo.py`)

The same calibration story runs **without torch**, using the over-confident
segmented-fit `ReferenceIVIMEstimator` in place of the MAF. This is the
`conformal_demo.py` path and the source of the numbers in this section.

### API

```python
import numpy as np
from caliper import metrics as M, conformal as C
from caliper.estimator_reference import ReferenceIVIMEstimator
from caliper.forward import synthetic_cohort, PARAM_NAMES

levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
cal, test = synthetic_cohort(n=4000, snr=40, seed=1), synthetic_cohort(n=9000, snr=40, seed=2)
est = ReferenceIVIMEstimator()
q_cal, q_test = est.predict_quantiles(cal.signals, levels), est.predict_quantiles(test.signals, levels)

# (a) split conformal on the absolute residual of a point predictor
sc = C.SplitConformalResidual(alpha=0.10).calibrate(est.predict_point(cal.signals), cal.params)
lo, hi = sc.apply(est.predict_point(test.signals))            # point ± Q

# (b) CQR — input-adaptive, restores MARGINAL coverage
cq = C.SplitConformalQuantile(levels).calibrate(q_cal, cal.params)
q_cqr = cq.apply(q_test)

# (c) Mondrian / group-conditional CQR — restores CONDITIONAL coverage per group
groups_cal = M.tercile_groups(cal.params[:, 2])              # D* terciles
strata     = M.tercile_groups(test.params[:, 2])
mq = C.MondrianConformalQuantile(levels).calibrate(q_cal, cal.params, groups_cal)
q_mond = mq.apply(q_test, strata)

# read coverage AND mean width per stratum
lo, hi = M.central_interval(q_cqr[:, 2, :], levels, 0.10)
C.conditional_coverage_by_strata(test.params[:, 2], lo, hi, strata)
```

### CQR restores **marginal** coverage (nominal 0.900)

| param | raw coverage | CQR coverage | raw \|gap\| | CQR \|gap\| |
|-------|-------------:|-------------:|------------:|------------:|
| D     | 0.676 | 0.902 | 0.224 | **0.002** |
| f     | 0.435 | 0.901 | 0.465 | **0.001** |
| D\*   | 0.359 | 0.903 | 0.541 | **0.003** |

The raw reference estimator is over-confident on every parameter (reported
quantiles too narrow); CQR restores marginal coverage to within **≤0.003** of
nominal. (`SplitConformalResidual` does the same from the point estimate alone.)

### …but **conditional** coverage is not — the D\* tercile result

Coverage **and mean interval width** of the 90% `D*` interval, stratified by
true-D\* tercile:

| method | low-D\* cov | width | mid-D\* cov | width | high-D\* cov | width |
|--------|------:|------:|------:|------:|------:|------:|
| raw          | 0.655 | 19.7 | 0.359 | 19.7 | 0.062 | 19.7 |
| marginal CQR | 0.951 | 215  | 0.875 | 215  | 0.882 | 215  |
| Mondrian CQR | 0.893 | 58.7 | 0.909 | 261  | 0.902 | 227  |

- **Marginal CQR** restores *pooled* D\* coverage (0.903) but applies one global
  width everywhere, so the well-identified **low-D\* tercile over-covers
  (0.951)** while the poorly-identified **high-D\* tercile under-covers (0.882)**.
  Conditional coverage is not delivered.
- **Mondrian CQR** restores per-tercile coverage (0.893 / 0.909 / 0.902) **only
  by inflating width**: high-D\* intervals are **3.87×** the low-D\* width.

Conformal guarantees marginal coverage unconditionally; conditional coverage
costs sharpness, and at high `D*` — the identifiability wall — the trade is
steep. The gap is the finding, reported as-is, not tuned away.

---

## What's in the box

```
caliper/
  metrics.py             # numpy-only calibration ruler (the canonical core)
  forward.py             # bi-exponential IVIM model + synthetic cohorts
  estimator_reference.py # over-confident segmented-fit IVIM estimator  [numpy]
  estimator_maf.py       # conditional MAF posterior over (D, f, D*)    [torch]
  conformal.py           # split-conformal / CQR / Mondrian + strata diagnostics
examples/
  demo.py                # MAF end-to-end pipeline (fixed seeds)        [torch]
  conformal_demo.py      # conformal + D* tercile result (fixed seeds)  [numpy]
tests/                   # pytest: metrics, forward, conformal, reference (numpy)
                         #         + estimator_maf (auto-skips without torch)
```

Run the tests:

```bash
pip install -e ".[dev]"
pytest -q          # 38 tests (1 MAF test auto-skips without torch)
```

## License

MIT — see [LICENSE](LICENSE).

## Roadmap

See [ROADMAP.md](ROADMAP.md). Value-of-information, decision-gap, and
deployment validity-monitor functionality are **deliberately deferred** and not
implemented here.
