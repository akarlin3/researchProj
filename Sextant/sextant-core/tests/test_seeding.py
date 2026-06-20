import numpy as np

from sextant.seeding import GLOBAL_SEED, make_rng


def test_global_seed_matches_monorepo_convention():
    assert GLOBAL_SEED == 20260613


def test_make_rng_is_deterministic():
    a = make_rng().integers(0, 1_000_000, size=50)
    b = make_rng().integers(0, 1_000_000, size=50)
    assert np.array_equal(a, b)


def test_make_rng_respects_explicit_seed():
    assert not np.array_equal(make_rng(1).random(10), make_rng(2).random(10))
