"""CP1 gate: the frozen task spec is well-formed and substrate wiring resolves."""
from __future__ import annotations

import numpy as np

from datum.baselines import BASELINES, panel_keys
from datum.task import QUANTILE_LEVELS, TASK_V1


def test_task_spec_well_formed():
    assert TASK_V1.version == "v1"
    assert TASK_V1.substrate == "gauge_cohort"
    assert TASK_V1.seed == 20260613
    assert 0.0 < TASK_V1.alpha < 1.0
    assert len(TASK_V1.param_names) == 3
    assert tuple(sorted(QUANTILE_LEVELS)) == tuple(QUANTILE_LEVELS)


def test_task_central_level_is_a_fashion_nominal_level():
    central = round(1.0 - TASK_V1.alpha, 3)
    assert central in TASK_V1.nominal_levels  # 0.90 is one of Fashion's LEVELS


def test_baseline_registry_covers_task():
    assert set(TASK_V1.baselines) == set(panel_keys())
    assert len(BASELINES) >= 6
    # spans paradigms from under-coverer to conditional fix
    paradigms = {b.paradigm for b in BASELINES.values()}
    assert {"parametric", "conformal", "mondrian"} <= paradigms


def test_substrate_registry_has_named_substrate():
    from datum.substrate import SUBSTRATES
    assert TASK_V1.substrate in SUBSTRATES
    assert {"gauge_cohort", "osipi_dro", "lattice"} <= set(SUBSTRATES)


def test_gauge_cohort_smoke():
    """Tiny end-to-end wiring check: substrate -> right shapes (no scoring yet)."""
    from datum.substrate import gauge_cohort
    sub = gauge_cohort(n_train=8, n_cal=8, n_test=8)
    assert set(sub.signals) == {"train", "cal", "test"}
    assert sub.params["test"].shape == (8, 3)          # (D, D*, f)
    assert sub.signals["test"].shape[0] == 8
    assert np.all(np.isfinite(sub.signals["test"]))
    assert sub.provenance["seed"] == 20260613
