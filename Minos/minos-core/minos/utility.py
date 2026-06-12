"""Action set, decision utility ``U(a, theta)``, and the closed-form expected
utility ``EU(a | q)`` for a Gaussian reported posterior ``q = N(mu, sigma^2)``.

Math: DESIGN.md Sections 1 and 4.
"""
from __future__ import annotations

from enum import IntEnum

import numpy as np
from scipy.stats import norm

from .config import MinosConfig


class Action(IntEnum):
    SPARE = 0
    TREAT = 1
    ESCALATE = 2


def _relu(x):
    return np.maximum(0.0, x)


def epp(m, sigma):
    """Expected positive part ``E[max(0, Y)]`` for ``Y ~ N(m, sigma^2)``.

    ``EPP(m, sigma) = m*Phi(m/sigma) + sigma*phi(m/sigma)`` for ``sigma > 0``;
    the ``sigma -> 0`` limit is ``relu(m)``. Vectorised over ``m`` (and ``sigma``).
    """
    m = np.asarray(m, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    m, sigma = np.broadcast_arrays(m, sigma)
    # Safe divisor avoids 0/0 warnings; the sigma==0 branch falls back to relu(m).
    safe = np.where(sigma > 0, sigma, 1.0)
    z = m / safe
    gaussian = m * norm.cdf(z) + sigma * norm.pdf(z)
    out = np.where(sigma > 0, gaussian, _relu(m))
    return out if out.ndim else float(out)


def utility(action: Action, theta, cfg: MinosConfig):
    """``U(a, theta)`` — piecewise-linear asymmetric mismatch utility (<= 0)."""
    theta = np.asarray(theta, dtype=float)
    if action == Action.SPARE:
        u = -cfg.k_under * _relu(theta - cfg.t1)
    elif action == Action.TREAT:
        u = -cfg.k_over * _relu(cfg.t1 - theta) - cfg.k_under * _relu(theta - cfg.t2)
    elif action == Action.ESCALATE:
        u = -cfg.k_over * _relu(cfg.t2 - theta)
    else:  # pragma: no cover - exhaustive
        raise ValueError(action)
    return u if u.ndim else u.item()


def expected_utility_under_q(action: Action, mu, sigma, cfg: MinosConfig):
    """``EU(a | q)`` for ``q = N(mu, sigma^2)`` in closed form via :func:`epp`."""
    if action == Action.SPARE:
        u = -cfg.k_under * epp(mu - cfg.t1, sigma)
    elif action == Action.TREAT:
        u = -cfg.k_over * epp(cfg.t1 - np.asarray(mu, float), sigma) - cfg.k_under * epp(
            np.asarray(mu, float) - cfg.t2, sigma
        )
    elif action == Action.ESCALATE:
        u = -cfg.k_over * epp(cfg.t2 - np.asarray(mu, float), sigma)
    else:  # pragma: no cover - exhaustive
        raise ValueError(action)
    return u
