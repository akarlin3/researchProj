"""Unit tests for echo_repeat.statistic -- the scientific core."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from echo_repeat import statistic as st  # noqa: E402


def test_norm_roundtrip():
    for p in (0.025, 0.1, 0.5, 0.9, 0.975):
        assert abs(st.norm_cdf(st.norm_ppf(p)) - p) < 1e-6


def test_analytic_repeat_coverage_below_nominal():
    # a perfectly measurement-scaled 90% interval covers the *repeat* at ~0.755, not 0.90
    c = st.analytic_repeat_coverage(0.10, scale=1.0, model_to_meas_ratio=0.0)
    assert 0.74 < c < 0.77
    # non-repeating model spread pushes coverage up toward 1
    assert st.analytic_repeat_coverage(0.10, model_to_meas_ratio=2.0) > c
    # narrower width pushes coverage down
    assert st.analytic_repeat_coverage(0.10, scale=0.5) < c


def test_coverage_basic():
    est_a = np.array([1.0, 2.0, 3.0])
    est_b = np.array([1.05, 2.5, 3.0])
    lo_a, hi_a = est_a - 0.2, est_a + 0.2
    # B in A: |1.05-1|=.05 in; |2.5-2|=.5 out; |3-3|=0 in -> 2/3 (one direction)
    c = st.test_retest_coverage(est_a, est_b, lo_a, hi_a, symmetrize=False)
    assert abs(c - 2.0 / 3.0) < 1e-9


def test_spearman_invariant_to_monotone_rescale():
    rng = np.random.default_rng(0)
    w = rng.random(50) + 0.1
    d = w * 2.0 + rng.random(50) * 0.01
    r1 = st.spearman(w, d)
    r2 = st.spearman(5.0 * w, d)          # pure positive rescale -> identical rank corr
    assert abs(r1 - r2) < 1e-12


def test_coverage_changes_under_rescale_while_spearman_fixed():
    # the distinctness-from-Gauge property, at unit level
    rng = np.random.default_rng(1)
    n = 5000
    truth = np.full(n, 1.0)
    sig = 0.1
    ea = truth + rng.normal(0, sig, n)
    eb = truth + rng.normal(0, sig, n)
    hw = 1.645 * sig
    c_full = st.test_retest_coverage(ea, eb, ea - hw, ea + hw, eb - hw, eb + hw)
    c_half = st.test_retest_coverage(ea, eb, ea - hw / 2, ea + hw / 2, eb - hw / 2, eb + hw / 2)
    assert c_full - c_half > 0.1        # scale matters
    width_full = np.full(n, 2 * hw)
    width_half = np.full(n, hw)
    dabs = np.abs(eb - ea)
    # constant widths -> spearman is nan/degenerate; use varying widths instead
    wv = hw * (1 + rng.random(n))
    assert abs(st.spearman(wv, dabs) - st.spearman(2 * wv, dabs)) < 1e-12


def test_bca_ci_brackets_point():
    rng = np.random.default_rng(2)
    x = rng.normal(0.8, 0.05, 200)
    point, lo, hi = st.bca_ci(lambda idx: float(np.mean(x[idx])), len(x),
                              n_boot=1500, seed=3)
    assert lo < point < hi
    assert abs(point - x.mean()) < 1e-9
