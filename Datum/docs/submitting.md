# Submitting a method to Datum

Datum scores **any** IVIM uncertainty method on the same frozen task, so you can see
how your method's calibration compares to the reference baselines. You implement one
function; Datum does the rest.

> All reported numbers are **PROVISIONAL** until Fashion's ruler locks.

## The contract

Predict central quantiles for the held-out test signals, in the physical
`(D, D*, f)` convention (mm²/s), at the task's quantile levels:

```
q_test : np.ndarray, shape (n_test, 3, L)
    q_test[i, p, j] = predicted q_levels[j] quantile of parameter p for voxel i,
    with p = 0 -> D, 1 -> D*, 2 -> f, and quantiles ascending along axis 2.
```

You get the calibration split (with ground truth) and the test signals (truth held
out). Anything that turns signals into quantiles is fair game: parametric fits,
conformal wrappers, learned posteriors, your own method.

## Minimal example

```python
import numpy as np
from datum.submit import load_task, score_submission

td = load_task()                      # td.b, td.q_levels, td.alpha,
                                      # td.cal_signals, td.cal_params (truth),
                                      # td.test_signals  (test truth is held out)

def my_method(td):
    q = td.empty_prediction()         # (n_test, 3, L)
    # ... fill q with your predicted quantiles for (D, D*, f) ...
    return np.sort(q, axis=2)         # quantiles must be non-decreasing

result = score_submission("my-method", my_method(td))
print(result.summary())
```

`score_submission` returns a `SubmissionResult` with, per parameter:
marginal **coverage**, **coverage-gap** (with 95% bootstrap CI), **ECE**,
**sharpness**; **per-D\*-tercile** coverage and width; and a **`vs_reference`**
ranking on |D\* coverage gap| against the curated baselines. Every field is
PROVISIONAL-stamped.

## A complete, runnable worked example

[`examples/submit_demo.py`](../examples/submit_demo.py) submits Gauge's
gradient-boosted **quantile regressor** (a method that is *not* one of Datum's
baselines) and prints its scorecard and ranking:

```
python Datum/examples/submit_demo.py          # full task
python Datum/examples/submit_demo.py --quick  # tiny smoke
```

## Rules of the road

- **Use only `cal` for calibration.** The test ground truth is held out by the
  interface; don't reconstruct it.
- **Report width with coverage.** Coverage is trivial to inflate; Datum always
  reports sharpness/width next to it (and the high-D\* width is the honest cost of
  beating the identifiability wall).
- **Numbers are PROVISIONAL.** Do not cite a Datum score as final until the ruler
  locks; re-run `python -m datum.run` after any change to the pins.
