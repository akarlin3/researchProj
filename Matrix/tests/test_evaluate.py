"""Evaluation helpers: bootstrap CI and AUROC."""
import numpy as np

from matrix.evaluate import bootstrap_ci, auroc


def test_bootstrap_ci_brackets_point_and_is_seed_stable():
    x = np.concatenate([np.zeros(50), np.ones(50)])
    p1 = bootstrap_ci(x, seed=0)
    p2 = bootstrap_ci(x, seed=0)
    assert p1 == p2                              # seeded -> reproducible
    point, lo, hi = p1
    assert abs(point - 0.5) < 1e-9
    assert lo < point < hi


def test_bootstrap_ci_constant_sample_zero_width():
    point, lo, hi = bootstrap_ci(np.zeros(20), seed=3)
    assert point == 0.0 and lo == 0.0 and hi == 0.0


def test_auroc_perfect_and_chance():
    labels = np.array([False] * 10 + [True] * 10)
    perfect = np.arange(20)                      # higher score == positive
    assert auroc(perfect, labels) == 1.0
    flipped = np.arange(20)[::-1]
    assert auroc(flipped, labels) == 0.0


def test_auroc_undefined_single_class():
    assert np.isnan(auroc(np.arange(5), np.zeros(5, bool)))
