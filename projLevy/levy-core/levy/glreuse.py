"""Grunwald-Letnikov fractional-derivative operators -- REUSED READ-ONLY from Ouroboros.

PROVENANCE
----------
Copied verbatim (clean-room, read-only reuse) from:
    Ouroboros/ouroboros_fractional_sindy.py  (gl_weights, gl_derivative_time,
                                              gl_derivative_space; AGPLv3, fully synthetic)

These are *not* Levy's contribution. They are kept here so the clean-room subrepo is
self-contained and so the noise-amplification law A(alpha) = dt^{-2 alpha} ||w(alpha)||^2
(reused from Ouroboros/ouroboros_noise_analysis.py) can be reproduced as a cross-check
against the b-indexed Fisher information that *is* Levy's contribution.

Levy's net-new layer (forward.py / fisher.py / identifiability.py / wall.py) does NOT
depend on these for the stretched-exponential lead lane; they are wired in for the
Phase-3 joint (alpha, beta) work and for the A(alpha) sanity anchor.
"""
from __future__ import annotations

import numpy as np


def gl_weights(order, n):
    """Grunwald-Letnikov weights w_k(order), k=0..n-1, via the standard recurrence.

    REUSED from Ouroboros/ouroboros_fractional_sindy.py:7 (read-only).
    """
    w = np.zeros(n)
    w[0] = 1.0
    for m in range(1, n):
        w[m] = w[m - 1] * (1.0 - (order + 1.0) / m)
    return w


def noise_amplification(alpha, nt, dt):
    """A(alpha) = dt^{-2 alpha} * sum_k w_k(alpha)^2  -- the GL noise-amplification law.

    REUSED (re-expressed) from Ouroboros/ouroboros_noise_analysis.py:20-45
    (``compute_analytic_factors``). Kept as a sanity anchor: it is the *time-domain
    derivative-operator* noise penalty, conceptually distinct from the b-indexed
    signal-decay Fisher information Levy derives. Used only as a cross-reference.
    """
    w = gl_weights(alpha, nt)
    return (dt ** (-2.0 * alpha)) * float(np.sum(w ** 2))
