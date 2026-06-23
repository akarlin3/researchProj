"""Forward model + Jacobian correctness."""
import numpy as np

from levy import forward


def test_alpha_one_is_monoexponential():
    b = np.array([0.0, 500.0, 1000.0, 2000.0])
    theta = np.array([1.0, 1.5e-3, 1.0])
    got = forward.signal(b, theta)
    expected = np.exp(-b * 1.5e-3)  # S0=1
    assert np.allclose(got, expected, rtol=1e-12)


def test_b0_returns_S0():
    theta = np.array([2.3, 1.5e-3, 0.6])
    assert np.isclose(forward.signal([0.0], theta)[0], 2.3)


def test_signal_monotone_decreasing_in_b():
    b = np.linspace(0, 3000, 50)
    s = forward.signal(b, np.array([1.0, 1.5e-3, 0.75]))
    assert np.all(np.diff(s) <= 1e-12)


def test_jacobian_matches_finite_difference():
    b = np.array([0.0, 100.0, 500.0, 1000.0, 2000.0, 3000.0])
    theta = np.array([1.0, 1.5e-3, 0.75])
    J = forward.jacobian(b, theta)
    # central finite differences, per-parameter step scaled to the parameter magnitude
    steps = np.array([1e-6, 1e-9, 1e-6])
    for k in range(3):
        tp = theta.copy(); tp[k] += steps[k]
        tm = theta.copy(); tm[k] -= steps[k]
        fd = (forward.signal(b, tp) - forward.signal(b, tm)) / (2 * steps[k])
        assert np.allclose(J[:, k], fd, rtol=1e-4, atol=1e-7), f"param {k}"


def test_jacobian_b0_only_informs_S0():
    J = forward.jacobian([0.0], np.array([1.0, 1.5e-3, 0.75]))
    assert np.isclose(J[0, forward.IDX["S0"]], 1.0)
    assert np.isclose(J[0, forward.IDX["D"]], 0.0)
    assert np.isclose(J[0, forward.IDX["alpha"]], 0.0)
