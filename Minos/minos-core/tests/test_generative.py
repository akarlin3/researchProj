"""Latent field generator + measurement model + common-random-number realisation."""
from __future__ import annotations

import numpy as np

from minos.config import MinosConfig
from minos.generative import make_population, realise
from minos.seeding import make_rng

CFG = MinosConfig(n_voxels=200_000)


def test_population_spans_all_three_regions():
    base = make_population(CFG, make_rng(0))
    theta = base.theta
    assert (theta < CFG.t1).mean() > 0.1
    assert ((theta >= CFG.t1) & (theta < CFG.t2)).mean() > 0.1
    assert (theta >= CFG.t2).mean() > 0.1


def test_population_mean_matches_mixture_mean():
    base = make_population(CFG, make_rng(0))
    expected = float(np.dot(CFG.mix_weights, CFG.mix_means))
    assert abs(base.theta.mean() - expected) < 0.02


def test_population_is_deterministic_given_seed():
    a = make_population(CFG, make_rng(0)).theta
    b = make_population(CFG, make_rng(0)).theta
    assert np.array_equal(a, b)


def test_realise_in_distribution_is_unbiased_with_intrinsic_spread():
    base = make_population(CFG, make_rng(0))
    mu, w = realise(base, CFG, delta=0.0, shift=False)
    assert abs((mu - base.theta).mean()) < 0.01
    assert abs((mu - base.theta).std() - CFG.s) < 0.01
    assert abs(w.mean()) < 0.01


def test_realise_under_shift_biases_down_and_moves_feature():
    base = make_population(CFG, make_rng(0))
    delta = 1.0
    mu, w = realise(base, CFG, delta=delta, shift=True)
    # downward bias of magnitude beta*s*delta, inflated spread, feature mean ~ delta.
    assert np.isclose((mu - base.theta).mean(), -CFG.beta * CFG.s * delta, atol=0.02)
    assert (mu - base.theta).std() > CFG.s * (1 + CFG.alpha * delta) - 0.05
    assert np.isclose(w.mean(), delta, atol=0.02)


def test_realise_crn_only_changes_with_delta_not_theta():
    base = make_population(CFG, make_rng(0))
    mu0, _ = realise(base, CFG, delta=0.0, shift=False)
    mu1, _ = realise(base, CFG, delta=0.5, shift=True)
    # theta is fixed; only the error transform changes -> residual is a deterministic
    # function of the same base normal draws (no re-sampling of theta).
    assert np.array_equal(base.theta, make_population(CFG, make_rng(0)).theta)
    assert not np.array_equal(mu0, mu1)
