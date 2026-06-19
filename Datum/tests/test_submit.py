"""CP3 gate: the submission/scoring interface works and is PROVISIONAL-flagged."""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from datum import submit
from datum.task import TASK_V1

TASK = replace(TASK_V1, n_train=100, n_cal=300, n_test=300, n_bootstrap=50)


def test_load_task_shapes_and_holds_out_truth():
    td = submit.load_task(task=TASK)
    assert td.cal_signals.shape == (300, td.b.size)
    assert td.cal_params.shape == (300, 3)            # cal truth provided
    assert td.test_signals.shape == (300, td.b.size)
    assert td.q_levels.size == len(TASK.quantile_levels)
    assert not hasattr(td, "test_params")             # test truth held out
    assert td.empty_prediction().shape == (300, 3, td.q_levels.size)


def _wide_submission(td, half=5.0):
    """A deliberately conservative submission: wide symmetric intervals everywhere."""
    # crude point = column means of cal truth, broadcast; wide quantile spread.
    point = td.cal_params.mean(axis=0)                # (3,)
    L = td.q_levels.size
    q = np.empty((td.n_test, 3, L))
    z = (np.asarray(td.q_levels) - 0.5) * 2 * half
    for p in range(3):
        q[:, p, :] = point[p] + z * abs(point[p])
    return np.sort(q, axis=2)


def test_score_submission_runs_and_is_provisional():
    td = submit.load_task(task=TASK)
    res = submit.score_submission("wide", _wide_submission(td), task=TASK)
    assert res.provisional is True
    assert "Fashion" in res.ruler
    assert set(res.per_param) == {"D", "Dstar", "f"}
    for p in res.per_param:
        lo, hi = res.ci[p]["coverage_gap"]
        assert lo <= res.per_param[p]["coverage_gap"] <= hi + 1e-9
    assert ("Dstar", "dstar_high") in res.by_tercile
    assert isinstance(res.summary(), str) and "PROVISIONAL" in res.summary()


def test_score_submission_validates_shape():
    td = submit.load_task(task=TASK)
    with pytest.raises(ValueError):
        submit.score_submission("bad", np.zeros((td.n_test, 3, 2)), task=TASK)


def test_vs_reference_present_when_csv_exists():
    """If reference numbers have been generated, ranking context is populated."""
    td = submit.load_task(task=TASK)
    res = submit.score_submission("wide", _wide_submission(td), task=TASK)
    # reference_numbers.csv is committed, so vs_reference should be non-empty.
    assert "submission_dstar_gap" in res.vs_reference
