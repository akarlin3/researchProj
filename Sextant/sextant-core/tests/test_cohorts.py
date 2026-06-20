import numpy as np

from sextant.cohorts import load_array_cohort
from sextant.railing import analyze_cohort


def _synthetic_volume():
    b = np.array([0.0, 50.0, 500.0, 800.0])
    H = W = 24
    vol = np.abs(np.random.default_rng(0).normal(2.0, 1.0, size=(H, W, 2, 4)))  # air
    D, Dstar, f, S0 = 1.2e-3, 0.03, 0.15, 500.0
    sig = S0 * ((1 - f) * np.exp(-b * D) + f * np.exp(-b * Dstar))
    vol[8:16, 8:16, :, :] = sig   # central body block (clean IVIM)
    return vol, b


def test_array_cohort_normalises_and_finds_body():
    vol, b = _synthetic_volume()
    coh = load_array_cohort("synthetic", vol, b)
    assert list(coh.bvals) == [0.0, 50.0, 500.0, 800.0]
    assert np.allclose(coh.fit_signals[:, 0], 1.0)            # normalised by b=0
    assert coh.n_high_snr >= 120                              # the 8x8x2 body block


def test_array_cohort_sorts_unsorted_bvals():
    vol, b = _synthetic_volume()
    perm = [2, 0, 3, 1]
    coh = load_array_cohort("perm", vol[..., perm], b[perm])
    assert list(coh.bvals) == [0.0, 50.0, 500.0, 800.0]


def test_clean_signals_barely_rail():
    vol, b = _synthetic_volume()
    coh = load_array_cohort("synthetic", vol, b)
    r = analyze_cohort(coh, bounds="tight")
    assert r.frac_railed < 0.1                                # clean interior D* -> little railing
