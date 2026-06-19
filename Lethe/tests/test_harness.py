"""Unit tests for echo_repeat.harness and the invivo plug-in fit."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from echo_repeat import harness, invivo  # noqa: E402


def test_harness_self_test_all_pass():
    # the locked CP1 method self-test must pass end-to-end
    res = harness.self_test(n=20000, level=0.10, seed=0)
    assert res["ALL_PASS"], res


def test_bias_invariance_explicit():
    no_bias = harness.simulate(20000, bias=0.0, seed=5)
    big_bias = harness.simulate(20000, bias=1e-2, seed=5)  # same seed, only bias differs
    c0 = harness.st.test_retest_coverage(no_bias.est_a, no_bias.est_b,
                                         no_bias.lo_a, no_bias.hi_a, no_bias.lo_b, no_bias.hi_b)
    cb = harness.st.test_retest_coverage(big_bias.est_a, big_bias.est_b,
                                         big_bias.lo_a, big_bias.hi_a, big_bias.lo_b, big_bias.hi_b)
    assert abs(c0 - cb) < 1e-9          # coverage literally identical: bias cancels in Delta


def test_segmented_fit_recovers_truth():
    b = invivo.ACRIN_BVALS
    D, Dstar, f = 1.2e-3, 30e-3, 0.10
    clean = invivo.ivim_signal(b, D, Dstar, f, S0=1.0)
    eD, eDs, ef = invivo.segmented_fit(clean, b)
    assert abs(eD - D) < 2e-4
    assert abs(ef - f) < 0.05


def test_synthetic_cohort_shapes():
    truth, est = invivo.synthetic_cohort(n=300, snr=12.0, seed=0)
    assert truth.shape[1] == 3 and est.shape == truth.shape
    assert np.all(np.isfinite(est))


def test_build_deployer_offsets_positive():
    dep = invivo.build_deployer(level=0.10, n_cal=500, seed=0)
    assert dep.offsets.shape == (3,)
    assert np.all(dep.offsets > 0)
    lo, hi = dep.apply(np.array([1.2e-3, 30e-3, 0.1]))
    assert np.all(hi > lo)
