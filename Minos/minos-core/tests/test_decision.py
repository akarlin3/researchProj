"""Bayes action: argmax expected utility under the reported posterior."""
from __future__ import annotations

import numpy as np

from minos.config import MinosConfig
from minos.decision import bayes_action
from minos.utility import Action

CFG = MinosConfig()


def test_point_action_correct_deep_in_each_region():
    mu = np.array([CFG.t1 - 2.0, 0.5 * (CFG.t1 + CFG.t2), CFG.t2 + 2.0])
    got = bayes_action(mu, 0.0, CFG)
    assert list(got) == [Action.SPARE, Action.TREAT, Action.ESCALATE]


def test_posterior_hedges_toward_costlier_error_near_lower_threshold():
    # Just below t1 the point rule spares; with a calibrated error bar the
    # asymmetric cost makes treating the better hedge.
    mu = CFG.t1 - 0.05
    assert bayes_action(mu, 0.0, CFG) == Action.SPARE
    assert bayes_action(mu, CFG.s, CFG) == Action.TREAT


def test_wider_posterior_pushes_toward_escalation():
    # At a fixed mu between the thresholds, growing the reported spread eventually
    # makes escalation (guarding the heavy under-treatment slope) optimal.
    mu = 0.5 * (CFG.t1 + CFG.t2)
    assert bayes_action(mu, 0.05, CFG) == Action.TREAT
    assert bayes_action(mu, 5.0, CFG) == Action.ESCALATE


def test_bayes_action_is_vectorised():
    mu = np.linspace(-5, 7, 50)
    got = bayes_action(mu, CFG.s, CFG)
    assert got.shape == mu.shape
    assert set(np.unique(got)).issubset({int(a) for a in Action})
