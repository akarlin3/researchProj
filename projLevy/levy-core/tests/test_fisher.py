"""Fisher / CRLB layer -- Levy's net-new contribution."""
import numpy as np

from levy import fisher, forward


def test_rician_fim_approaches_gaussian_at_high_snr():
    b = forward.StretchedExp().theta  # placeholder, replaced below
    b = np.array([0.0, 500.0, 1000.0, 2000.0, 3000.0])
    theta = np.array([1.0, 1.5e-3, 0.75])
    fim_g = fisher.fisher_matrix(b, theta, snr=500.0, model="gaussian")
    fim_r = fisher.fisher_matrix(b, theta, snr=500.0, model="rician")
    assert np.allclose(fim_g, fim_r, rtol=2e-2)


def test_crlb_positive_definite_at_reasonable_snr():
    r = fisher.crlb(np.array([0.0, 500, 1000, 2000, 3000.0]), np.array([1.0, 1.5e-3, 0.75]), snr=50)
    assert np.all(np.diag(r.cov) > 0)
    assert r.se_alpha > 0


def test_crlb_alpha_decreases_with_snr():
    b = np.array([0.0, 500, 1000, 2000, 3000.0])
    theta = np.array([1.0, 1.5e-3, 0.75])
    cv = [fisher.crlb(b, theta, snr).cv_alpha for snr in (10, 30, 60, 120)]
    assert all(cv[i] > cv[i + 1] for i in range(len(cv) - 1)), cv


def test_more_bvalues_tighten_alpha():
    theta = np.array([1.0, 1.5e-3, 0.75])
    few = fisher.crlb(np.array([0.0, 1000, 2000.0]), theta, snr=40).cv_alpha
    many = fisher.crlb(np.linspace(0, 3000, 16), theta, snr=40).cv_alpha
    assert many < few


def test_rician_crlb_looser_than_gaussian_at_low_snr():
    b = np.linspace(0, 3000, 12)
    theta = np.array([1.0, 1.5e-3, 0.75])
    cv_g = fisher.crlb(b, theta, snr=8, model="gaussian").cv_alpha
    cv_r = fisher.crlb(b, theta, snr=8, model="rician").cv_alpha
    assert cv_r > cv_g  # Rician information is lower -> bound is looser
