import numpy as np

from sextant.bootstrap import bootstrap_fraction


def test_degenerate_all_true():
    ci = bootstrap_fraction(np.ones(200, bool), n_boot=500)
    assert ci.point == 1.0 and ci.lo == 1.0 and ci.hi == 1.0


def test_ci_brackets_point_and_is_reproducible():
    rng = np.random.default_rng(0)
    x = rng.random(2000) < 0.5
    a = bootstrap_fraction(x, n_boot=2000)
    b = bootstrap_fraction(x, n_boot=2000)
    assert a.lo <= a.point <= a.hi
    assert (a.lo, a.hi) == (b.lo, b.hi)          # seeded -> identical
    assert abs(a.point - 0.5) < 0.05


def test_ci_width_shrinks_with_n():
    rng = np.random.default_rng(1)
    small = bootstrap_fraction(rng.random(100) < 0.5, n_boot=1000)
    large = bootstrap_fraction(rng.random(10000) < 0.5, n_boot=1000)
    assert (large.hi - large.lo) < (small.hi - small.lo)
