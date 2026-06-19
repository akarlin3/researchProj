"""Physical-property gates on the forward generators."""

import numpy as np
import pytest

from lattice import generators as G
from lattice.cohort import DEFAULT_BVALUES as B


def test_signal_at_b0_is_one():
    assert G.ivim_biexp(0.0, 1.5e-3, 30e-3, 0.2) == pytest.approx(1.0)


def test_f_zero_is_monoexponential():
    s = G.ivim_biexp(B, 1.5e-3, 30e-3, 0.0)
    assert np.allclose(s, np.exp(-B * 1.5e-3), rtol=1e-12)


def test_dstar_equals_d_is_monoexponential_for_any_f():
    for f in (0.05, 0.2, 0.4):
        s = G.ivim_biexp(B, 1.2e-3, 1.2e-3, f)
        assert np.allclose(s, np.exp(-B * 1.2e-3), rtol=1e-12)


def test_signal_strictly_decreasing_in_b():
    s = G.ivim_biexp(B, 1.5e-3, 40e-3, 0.25)
    assert np.all(np.diff(s) < 0)


def test_high_b_logslope_recovers_D():
    D = 1.3e-3
    s = G.ivim_biexp(B, D, 50e-3, 0.2)
    hi = B >= 300
    slope = np.polyfit(B[hi], np.log(s[hi]), 1)[0]
    assert -slope == pytest.approx(D, rel=1e-3)


def test_all_families_b0_is_one():
    assert G.ivim_dispersion_gamma(0.0, 1.5e-3, 30e-3, 4.0, 0.2) == pytest.approx(1.0)
    assert G.ivim_dispersion_lognormal(np.array([0.0]), 1.5e-3, 30e-3, 0.5, 0.2)[0] == pytest.approx(1.0)
    assert G.ivim_stretched(0.0, 1.5e-3, 30e-3, 0.2, 0.7) == pytest.approx(1.0)
    assert G.ivim_triexp(0.0, 1.5e-3, 30e-3, 0.2, 120e-3, 0.3) == pytest.approx(1.0)


def test_rician_noise_floor_on_zero_signal():
    rng = np.random.default_rng(0)
    zeros = np.zeros(200_000)
    snr = 50.0
    mag = G.add_rician_noise(zeros, snr, rng)
    sigma = 1.0 / snr
    assert mag.mean() == pytest.approx(sigma * np.sqrt(np.pi / 2), rel=2e-2)


def test_vectorised_shape():
    n, nb = 17, len(B)
    D = np.full(n, 1.5e-3)
    Dstar = np.full(n, 30e-3)
    f = np.full(n, 0.2)
    assert G.ivim_biexp(B, D, Dstar, f).shape == (n, nb)
    assert G.ivim_dispersion_lognormal(B, D, Dstar, np.full(n, 0.5), f).shape == (n, nb)
    assert G.ivim_triexp(B, D, Dstar, f, 4 * Dstar, np.full(n, 0.3)).shape == (n, nb)
