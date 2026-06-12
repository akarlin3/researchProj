"""Value-of-information core: policy expected utilities, the EVPI-analog,
the value of using the error bar, and the Value of Calibration (VoC).

Math: DESIGN.md Sections 4 and 6. All expectations use the CRN base draws so that
differences (VoC) are low-variance. The gated policy and VoTG live in :mod:`gate`.
"""
from __future__ import annotations

import numpy as np

from .config import MinosConfig
from .decision import bayes_action
from .generative import BaseDraws, realise
from .utility import Action, utility

# Policies handled here (the gated policy is in gate.py).
_NONGATED = ("point", "posterior", "oracle")


def realised_utility(actions: np.ndarray, theta: np.ndarray, cfg: MinosConfig) -> np.ndarray:
    """Per-voxel ``U(action_i, theta_i)`` for an integer action array."""
    cols = np.stack([utility(a, theta, cfg) for a in Action], axis=0)  # (3, n)
    return np.take_along_axis(cols, actions[None, :], axis=0)[0]


def policy_actions(
    policy: str, base: BaseDraws, cfg: MinosConfig, *, tau: float = 1.0,
    delta: float = 0.0, shift=False,
) -> np.ndarray:
    """Action chosen by ``policy`` for every voxel (non-gated policies only)."""
    if policy == "oracle":
        return np.asarray(bayes_action(base.theta, 0.0, cfg))
    mu, _ = realise(base, cfg, delta=delta, shift=shift)
    if policy == "point":
        return np.asarray(bayes_action(mu, 0.0, cfg))
    if policy == "posterior":
        return np.asarray(bayes_action(mu, tau * cfg.s, cfg))
    raise ValueError(f"unknown policy {policy!r} (gated policy lives in gate.py)")


def expected_utility(
    policy: str, base: BaseDraws, cfg: MinosConfig, *, tau: float = 1.0,
    delta: float = 0.0, shift=False,
) -> float:
    """Monte-Carlo expected utility ``E[U(a_policy, theta)]`` of a policy."""
    actions = policy_actions(policy, base, cfg, tau=tau, delta=delta, shift=shift)
    return float(np.mean(realised_utility(actions, base.theta, cfg)))


def evpi(base: BaseDraws, cfg: MinosConfig, *, tau: float = 1.0,
         delta: float = 0.0, shift=False) -> float:
    """EVPI-analog ``= EU(oracle) - EU(posterior)`` (expected posterior regret)."""
    eu_oracle = expected_utility("oracle", base, cfg)
    eu_post = expected_utility("posterior", base, cfg, tau=tau, delta=delta, shift=shift)
    return eu_oracle - eu_post


def value_of_error_bar(base: BaseDraws, cfg: MinosConfig) -> float:
    """``EU(posterior | tau=1) - EU(point)`` at ``delta=0`` (calibrated)."""
    eu_post = expected_utility("posterior", base, cfg, tau=1.0)
    eu_point = expected_utility("point", base, cfg)
    return eu_post - eu_point


def posterior_eu_curve(base: BaseDraws, cfg: MinosConfig, taus, *,
                       delta: float = 0.0, shift=False) -> np.ndarray:
    """``EU(posterior | tau)`` over a grid of ``tau`` (CRN: one realisation)."""
    mu, _ = realise(base, cfg, delta=delta, shift=shift)
    out = np.empty(len(taus))
    for i, tau in enumerate(taus):
        actions = np.asarray(bayes_action(mu, tau * cfg.s, cfg))
        out[i] = float(np.mean(realised_utility(actions, base.theta, cfg)))
    return out


def voc(base: BaseDraws, cfg: MinosConfig, tau: float) -> float:
    """Value of Calibration ``VoC(tau) = EU(posterior|tau=1) - EU(posterior|tau)``."""
    eu1 = expected_utility("posterior", base, cfg, tau=1.0)
    eu_tau = expected_utility("posterior", base, cfg, tau=tau)
    return eu1 - eu_tau
