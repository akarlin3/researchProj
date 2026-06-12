"""Calibration diagnostics on the toy: central-interval coverage and ECE."""
from __future__ import annotations

import numpy as np

from minos.config import MinosConfig
from minos.diagnostics import central_interval_coverage, ece
from minos.generative import make_population
from minos.seeding import make_rng

CFG = MinosConfig(n_voxels=400_000)
LEVELS = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95])


def _pop(cfg=CFG):
    return make_population(cfg, make_rng(cfg.seed))


def test_calibrated_coverage_matches_nominal():
    base = _pop()
    for lvl in (0.5, 0.8, 0.9):
        cov = central_interval_coverage(base, CFG, lvl, tau=1.0)
        assert abs(cov - lvl) < 0.005


def test_overconfident_undercovers_and_underconfident_overcovers():
    base = _pop()
    lvl = 0.9
    assert central_interval_coverage(base, CFG, lvl, tau=0.6) < lvl - 0.02
    assert central_interval_coverage(base, CFG, lvl, tau=1.5) > lvl + 0.02


def test_ece_minimal_when_calibrated():
    base = _pop()
    ece_cal = ece(base, CFG, LEVELS, tau=1.0)
    ece_over = ece(base, CFG, LEVELS, tau=0.6)
    ece_under = ece(base, CFG, LEVELS, tau=1.5)
    assert ece_cal < 0.005
    assert ece_over > ece_cal
    assert ece_under > ece_cal


def test_shift_degrades_calibration_even_at_tau_one():
    base = _pop()
    ece_id = ece(base, CFG, LEVELS, tau=1.0, delta=0.0, shift=True)
    ece_shift = ece(base, CFG, LEVELS, tau=1.0, delta=1.5, shift=True)
    assert ece_shift > ece_id + 0.05
