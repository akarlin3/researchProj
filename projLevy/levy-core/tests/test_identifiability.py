"""Empirical identifiability: MLE recovery, profile-likelihood CI, bootstrap.

CRLB is the *bound*; these tests confirm the estimator actually approaches it at high SNR
(so the wall, when it appears at low SNR, is a real information limit, not an estimator bug).
"""
import numpy as np

from levy import forward, identifiability, noise, seeding


def test_mle_recovers_truth_at_high_snr():
    rng = seeding.make_rng(7)
    truth = forward.StretchedExp(S0=1.0, D=1.5e-3, alpha=0.75)
    b = np.linspace(0, 3000, 16)
    sigma = noise.sigma_from_snr(truth.S0, snr=200.0)
    nu = forward.signal(b, truth.theta)
    # average several high-SNR realizations
    a_hats = []
    for _ in range(8):
        M = noise.rician_sample(nu, sigma, rng)
        theta_hat, _, _ = identifiability.mle(b, M, sigma, truth.theta)
        a_hats.append(theta_hat[forward.IDX["alpha"]])
    assert abs(np.mean(a_hats) - truth.alpha) < 0.05


def test_profile_ci_coverage_and_tightness_at_high_snr():
    # A 95% CI from one realization need NOT contain truth; test COVERAGE across realizations.
    rng = seeding.make_rng(3)
    truth = forward.StretchedExp()
    b = np.linspace(0, 3000, 16)
    sigma = noise.sigma_from_snr(truth.S0, snr=150.0)
    nu = forward.signal(b, truth.theta)
    contained, widths = 0, []
    K = 16
    for _ in range(K):
        M = noise.rician_sample(nu, sigma, rng)
        ci = identifiability.profile_ci_alpha(b, M, sigma, truth.theta)
        assert np.isfinite(ci.lo) and np.isfinite(ci.hi)
        widths.append(ci.width)
        if ci.lo <= truth.alpha <= ci.hi:
            contained += 1
    assert contained >= 12, f"coverage too low: {contained}/{K}"  # ~95% nominal, allow slack
    assert np.median(widths) < 0.2  # tight at high SNR


def test_bootstrap_se_grows_as_snr_drops():
    truth = forward.StretchedExp()
    b = np.linspace(0, 3000, 12)
    hi = identifiability.parametric_bootstrap(truth, b, snr=80, n_boot=60, rng=seeding.make_rng(1))
    lo = identifiability.parametric_bootstrap(truth, b, snr=15, n_boot=60, rng=seeding.make_rng(2))
    assert lo.se > hi.se
