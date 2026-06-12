"""Utility primitives: U(a, theta), the expected-positive-part EPP, and EU(a | q)."""
from __future__ import annotations

import numpy as np

from minos.config import MinosConfig
from minos.seeding import make_rng
from minos.utility import Action, epp, expected_utility_under_q, utility

CFG = MinosConfig()


def test_epp_zero_sigma_is_relu():
    m = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
    assert np.allclose(epp(m, 0.0), np.maximum(0.0, m))


def test_epp_zero_mean_is_sigma_over_root_2pi():
    sigma = 0.7
    assert np.isclose(epp(0.0, sigma), sigma / np.sqrt(2 * np.pi))


def test_epp_matches_monte_carlo():
    rng = make_rng(1)
    m, sigma = 0.3, 0.9
    draws = m + sigma * rng.standard_normal(2_000_000)
    mc = np.mean(np.maximum(0.0, draws))
    assert abs(epp(m, sigma) - mc) < 2e-3


def test_utility_correct_action_is_zero_per_region():
    # spare correct below t1; treat correct in (t1,t2); escalate correct above t2.
    assert np.isclose(utility(Action.SPARE, CFG.t1 - 1.0, CFG), 0.0)
    assert np.isclose(utility(Action.TREAT, 0.5 * (CFG.t1 + CFG.t2), CFG), 0.0)
    assert np.isclose(utility(Action.ESCALATE, CFG.t2 + 1.0, CFG), 0.0)


def test_oracle_utility_identically_zero():
    # max_a U(a, theta) == 0 for every theta (best action attains zero loss).
    theta = np.linspace(-6, 8, 1401)
    best = np.max(
        np.stack([utility(a, theta, CFG) for a in Action], axis=0), axis=0
    )
    assert np.allclose(best, 0.0, atol=1e-12)


def test_under_treatment_costs_more_than_over_treatment():
    # symmetric distance d past each threshold: sparing a severe case (under) must
    # cost strictly more than escalating a mild case (over).
    d = 0.8
    under = -utility(Action.SPARE, CFG.t1 + d, CFG)      # should-treat, spared
    over = -utility(Action.ESCALATE, CFG.t2 - d, CFG)    # should-treat, escalated
    assert under > over > 0


def test_expected_utility_reduces_to_pointwise_when_sigma_zero():
    mu = np.linspace(-4, 6, 200)
    for a in Action:
        assert np.allclose(
            expected_utility_under_q(a, mu, 0.0, CFG), utility(a, mu, CFG)
        )


def test_expected_utility_matches_monte_carlo():
    rng = make_rng(2)
    mu, sigma = 0.4, 0.6
    theta = mu + sigma * rng.standard_normal(2_000_000)
    for a in Action:
        mc = np.mean(utility(a, theta, CFG))
        assert abs(expected_utility_under_q(a, mu, sigma, CFG) - mc) < 3e-3
