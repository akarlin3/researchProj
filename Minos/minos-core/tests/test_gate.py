"""Trust-gate: signal, threshold, VoTG, and shift-detection AUC.

Encodes GATE 3: VoTG(0) ~ 0, VoTG(delta) > 0 under shift, gated regret < posterior
regret, and detection AUC > 0.5 + margin."""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from minos.config import MinosConfig
from minos.gate import (
    detection_auc,
    expected_utility_gated,
    gate_signal,
    gate_threshold,
    votg,
)
from minos.generative import make_population
from minos.seeding import make_rng
from minos.voi import expected_utility

CFG = MinosConfig(n_voxels=600_000)
DELTA_TEST = 1.5


def _pop(cfg=CFG):
    return make_population(cfg, make_rng(cfg.seed))


def test_gate_signal_standardizes_with_training_stats():
    cfg = MinosConfig(w_train_mean=0.0, w_train_std=1.0)
    w = np.array([-1.0, 0.0, 2.0])
    assert np.allclose(gate_signal(w, cfg), w)


def test_gate_threshold_is_training_quantile():
    assert np.isclose(gate_threshold(CFG), norm.ppf(CFG.q_gate))


# ---- GATE 3 ---------------------------------------------------------------------
def test_votg_approximately_zero_without_shift():
    base = _pop()
    assert abs(votg(base, CFG, delta=0.0)) < 0.02


def test_votg_positive_under_shift():
    base = _pop()
    assert votg(base, CFG, delta=DELTA_TEST) > 0.05


def test_gated_regret_below_posterior_regret_under_shift():
    base = _pop()
    # regret = EU(oracle) - EU(policy) = -EU(policy), since oracle utility is 0.
    eu_post = expected_utility("posterior", base, CFG, delta=DELTA_TEST, shift=True)
    eu_gated = expected_utility_gated(base, CFG, delta=DELTA_TEST, shift=True)
    assert (-eu_gated) < (-eu_post)


def test_detection_auc_is_chance_without_shift():
    base = _pop()
    assert abs(detection_auc(base, CFG, 0.0) - 0.5) < 0.01


def test_detection_auc_above_chance_under_shift():
    base = _pop()
    assert detection_auc(base, CFG, DELTA_TEST) > 0.7


def test_detection_auc_matches_gaussian_theory():
    base = _pop()
    for d in (0.5, 1.0, 1.5):
        assert abs(detection_auc(base, CFG, d) - norm.cdf(d / np.sqrt(2))) < 0.01
