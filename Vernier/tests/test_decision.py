"""Experiment B decision-value (PROVISIONAL, Minos lens): runs, deterministic, monotone.

Skipped automatically if scipy (a Minos dependency) is unavailable.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")  # Minos's utility module imports scipy.stats

from vernier import _paths
from vernier.schemes import by_name, EFFICIENCY_FRONTIER


def _cfg_and_params(seed=0, n=3000):
    _paths.add_caliper()
    from caliper.forward import sample_params
    from vernier.decision import make_decision_config

    params = sample_params(n, np.random.default_rng(seed))
    return make_decision_config(params[:, 2]), params


def test_decision_value_runs_and_is_deterministic():
    from vernier.decision import decision_value

    cfg, _ = _cfg_and_params()
    s = by_name(EFFICIENCY_FRONTIER, "clinical-11")
    a = decision_value(s, n=2000, snr=33.0, seed=0, cfg=cfg)
    b = decision_value(s, n=2000, snr=33.0, seed=0, cfg=cfg)
    assert a.mean_utility == b.mean_utility
    assert a.width_dstar == b.width_dstar
    assert a.mean_utility <= 0.0  # utility is a non-positive cost


def test_more_scan_time_sharpens_and_improves_decisions():
    from vernier.decision import decision_value

    cfg, _ = _cfg_and_params()
    sparse = decision_value(by_name(EFFICIENCY_FRONTIER, "sparse-7"), n=4000, snr=33.0, seed=0, cfg=cfg)
    dense = decision_value(by_name(EFFICIENCY_FRONTIER, "dense-22"), n=4000, snr=33.0, seed=0, cfg=cfg)
    # denser sampling -> sharper corrected D* and better (less negative) decisions
    assert dense.width_dstar < sparse.width_dstar
    assert dense.mean_utility > sparse.mean_utility
