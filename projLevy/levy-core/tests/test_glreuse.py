"""Reused Grunwald-Letnikov layer (read-only from Ouroboros) -- provenance + A(alpha) anchor."""
import numpy as np

from levy import glreuse


def test_gl_weights_recurrence():
    w = glreuse.gl_weights(0.5, 6)
    assert w[0] == 1.0
    # w[m] = w[m-1]*(1-(alpha+1)/m)
    for m in range(1, 6):
        assert np.isclose(w[m], w[m - 1] * (1 - 1.5 / m))


def test_noise_amplification_increases_with_alpha():
    # A(alpha) = dt^{-2 alpha} ||w||^2 grows with alpha (Ouroboros result, cross-check anchor)
    a_lo = glreuse.noise_amplification(0.3, nt=200, dt=0.01)
    a_hi = glreuse.noise_amplification(0.9, nt=200, dt=0.01)
    assert a_hi > a_lo
