#!/usr/bin/env python
"""Worked example: submit a *new* method to Datum and get it scored.

The method here -- Gauge's gradient-boosted **conditional quantile regressor** -- is
deliberately NOT one of Datum's curated baselines, to show the submission flow for
an external method. It is trained on the task's calibration split and predicts test
quantiles in the physical (D, D*, f) convention; Datum scores it on Fashion's ruler
and reports where it lands relative to the reference baselines.

Run (from the monorepo root; works with no install / no PYTHONPATH):
    python Datum/examples/submit_demo.py          # full task
    python Datum/examples/submit_demo.py --quick  # tiny cohort smoke
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

# Put the Datum package root on sys.path so this runs as a plain script
# (`python Datum/examples/submit_demo.py`) without an install or PYTHONPATH.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from datum import _paths  # noqa: E402

_paths.ensure_deps()

from gauge.estimators import IVIMQuantileRegressor  # noqa: E402  (a reused, non-baseline method)

from datum.submit import load_task, score_submission  # noqa: E402
from datum.task import TASK_V1  # noqa: E402


def gauge_qr_submission(td):
    """Train Gauge's quantile regressor on the cal split; predict test quantiles.

    Returns (n_test, 3, L) for (D, D*, f) at td.q_levels.
    """
    qr = IVIMQuantileRegressor(list(td.q_levels))
    qr.fit(td.cal_signals, td.cal_params)            # cal_params is (D, D*, f) physical
    n, L = td.n_test, td.q_levels.size
    q = np.empty((n, 3, L), dtype=float)
    for p in range(3):
        for j, ql in enumerate(td.q_levels):
            q[:, p, j] = qr.predict_quantile(td.test_signals, p, float(ql))
    return np.sort(q, axis=2)                         # enforce non-crossing quantiles


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="tiny cohort smoke")
    args = ap.parse_args(argv)

    task = TASK_V1
    if args.quick:
        task = replace(TASK_V1, n_train=100, n_cal=300, n_test=300, n_bootstrap=100)

    print("Loading Datum task (test truth held out)...")
    td = load_task(task=task)
    print(f"  b-values: {td.b.size}; cal: {td.cal_signals.shape}; "
          f"test: {td.test_signals.shape}; levels: {td.q_levels.size}")

    print("Training + predicting with Gauge quantile regressor (a non-baseline method)...")
    q_test = gauge_qr_submission(td)

    result = score_submission("gauge-quantile-regressor", q_test, task=task)
    print()
    print(result.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
