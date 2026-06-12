"""Bayes action under a Gaussian reported posterior ``q = N(mu, sigma^2)``.

``a*(q) = argmax_a EU(a | q)``, evaluated in closed form (no inner Monte Carlo).
"""
from __future__ import annotations

import numpy as np

from .config import MinosConfig
from .utility import Action, expected_utility_under_q


def bayes_action(mu, sigma, cfg: MinosConfig):
    """Return the utility-maximising action(s) for ``q = N(mu, sigma^2)``.

    ``mu`` may be scalar or array; ``sigma`` is broadcast against it. Returns an
    int array (``Action`` values) for array ``mu``, or a scalar ``int`` otherwise.
    Ties resolve to the lower-index (more conservative-to-spare) action via argmax.
    """
    mu = np.asarray(mu, dtype=float)
    eus = np.stack(
        [expected_utility_under_q(a, mu, sigma, cfg) for a in Action], axis=0
    )
    idx = np.argmax(eus, axis=0)
    return idx if idx.ndim else int(idx)
