"""Tests for the interval score (proper scoring rule for central intervals)."""
import numpy as np
import pytest

from gauge.conformal import interval_score


def test_interval_score_inside_is_just_width():
    # y inside [l, u] -> no penalty, score == width.
    assert interval_score(0.0, 2.0, 1.0, alpha=0.1) == pytest.approx(2.0)


def test_interval_score_penalizes_miss_above():
    # y above u: width + (2/alpha)*(y-u).
    assert interval_score(0.0, 2.0, 3.0, alpha=0.1) == pytest.approx(2.0 + 20 * 1.0)


def test_interval_score_penalizes_miss_below_symmetrically():
    assert interval_score(0.0, 2.0, -1.0, alpha=0.1) == pytest.approx(2.0 + 20 * 1.0)


def test_interval_score_vectorized():
    lo = np.array([0.0, 0.0, 0.0])
    hi = np.array([2.0, 2.0, 2.0])
    y = np.array([1.0, 3.0, -1.0])
    got = interval_score(lo, hi, y, alpha=0.1)
    np.testing.assert_allclose(got, [2.0, 22.0, 22.0])


def test_interval_score_smaller_alpha_penalizes_misses_more():
    s_small = interval_score(0.0, 1.0, 2.0, alpha=0.05)
    s_big = interval_score(0.0, 1.0, 2.0, alpha=0.20)
    assert s_small > s_big
