"""Value-of-information core: policy expected utilities, EVPI-analog, VoC.

Encodes GATE 1 (EU ordering, degenerate EVPI) and GATE 2 (VoC limits)."""
from __future__ import annotations

import numpy as np

from minos.config import MinosConfig
from minos.generative import make_population
from minos.seeding import make_rng
from minos.voi import (
    evpi,
    expected_utility,
    posterior_eu_curve,
    value_of_error_bar,
    voc,
)

CFG = MinosConfig(n_voxels=200_000)


def _pop(cfg=CFG):
    return make_population(cfg, make_rng(cfg.seed))


# ---- GATE 1 ---------------------------------------------------------------------
def test_eu_ordering_oracle_posterior_point():
    base = _pop()
    eu_point = expected_utility("point", base, CFG)
    eu_post = expected_utility("posterior", base, CFG, tau=1.0)
    eu_oracle = expected_utility("oracle", base, CFG)
    assert eu_oracle >= eu_post - 1e-9
    assert eu_post >= eu_point - 1e-9
    assert np.isclose(eu_oracle, 0.0, atol=1e-9)  # oracle loss is identically zero


def test_evpi_nonnegative():
    base = _pop()
    assert evpi(base, CFG, tau=1.0) >= -1e-9


def test_evpi_vanishes_as_posterior_to_point_mass():
    cfg = MinosConfig(n_voxels=100_000, s=1e-4)
    base = make_population(cfg, make_rng(cfg.seed))
    assert evpi(base, cfg, tau=1.0) < 1e-2


def test_value_of_using_error_bar_is_positive_when_calibrated():
    base = _pop()
    assert value_of_error_bar(base, CFG) > 0.0


# ---- GATE 2 ---------------------------------------------------------------------
def test_voc_zero_at_tau_one():
    base = _pop()
    assert np.isclose(voc(base, CFG, 1.0), 0.0, atol=1e-9)


CFG_BIG = MinosConfig(n_voxels=1_000_000)  # VoC near tau=1 is small; resolve it


def test_voc_minimised_at_calibration_and_grows_away_from_it():
    base = make_population(CFG_BIG, make_rng(CFG_BIG.seed))
    taus = np.round(np.arange(0.5, 2.0001, 0.05), 3)
    curve = np.array([voc(base, CFG_BIG, t) for t in taus])
    # VoC >= 0 everywhere (within MC tolerance) and is exactly 0 at tau=1.
    assert curve.min() >= -5e-5
    one = int(np.argmin(np.abs(taus - 1.0)))
    assert np.isclose(curve[one], 0.0, atol=1e-9)
    # The calibration point is the optimum: argmin sits at tau=1.
    assert int(np.argmin(curve)) == one
    # Strictly increasing as tau moves away from 1 on each side, where the signal
    # is well above the MC noise floor (|tau-1| >= 0.1).
    offs = [0.1, 0.2, 0.3, 0.4, 0.5]
    left = [voc(base, CFG_BIG, round(1 - d, 3)) for d in offs]
    right = [voc(base, CFG_BIG, round(1 + d, 3)) for d in offs]
    assert np.all(np.diff(left) > 0)
    assert np.all(np.diff(right) > 0)


def test_posterior_eu_curve_peaks_at_calibration():
    base = make_population(CFG_BIG, make_rng(CFG_BIG.seed))
    taus = np.round(np.arange(0.5, 2.0001, 0.05), 3)
    eus = posterior_eu_curve(base, CFG_BIG, taus)
    one = int(np.argmin(np.abs(taus - 1.0)))
    assert int(np.argmax(eus)) == one
