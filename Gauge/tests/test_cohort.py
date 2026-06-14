"""Tests for the labeled synthetic IVIM cohort generator."""
import numpy as np

from gauge.cohort import generate_cohort, D_RANGE, DSTAR_RANGE, F_RANGE
from gauge.forward import DEFAULT_B_VALUES


def test_cohort_is_deterministic_under_seed():
    a = generate_cohort(200, 100, 100, snr_grid=(20, 50), seed=42)
    b = generate_cohort(200, 100, 100, snr_grid=(20, 50), seed=42)
    np.testing.assert_array_equal(a.signals["train"], b.signals["train"])
    np.testing.assert_array_equal(a.params["test"], b.params["test"])
    np.testing.assert_array_equal(a.snr["cal"], b.snr["cal"])


def test_cohort_shapes_and_split_sizes():
    c = generate_cohort(300, 120, 150, snr_grid=(10, 30, 60), seed=1)
    nb = DEFAULT_B_VALUES.size
    assert c.signals["train"].shape == (300, nb)
    assert c.signals["cal"].shape == (120, nb)
    assert c.signals["test"].shape == (150, nb)
    assert c.params["train"].shape == (300, 3)
    assert c.snr["test"].shape == (150,)


def test_params_within_published_ranges():
    c = generate_cohort(2000, 10, 10, snr_grid=(20,), seed=7)
    P = c.params["train"]
    assert P[:, 0].min() >= D_RANGE[0] and P[:, 0].max() <= D_RANGE[1]
    assert P[:, 1].min() >= DSTAR_RANGE[0] and P[:, 1].max() <= DSTAR_RANGE[1]
    assert P[:, 2].min() >= F_RANGE[0] and P[:, 2].max() <= F_RANGE[1]


def test_snr_drawn_from_grid():
    grid = (10, 25, 80)
    c = generate_cohort(1000, 10, 10, snr_grid=grid, seed=3)
    assert set(np.unique(c.snr["train"]).astype(int)).issubset(set(grid))


def test_b0_signal_near_one_at_high_snr():
    # With S0=1 and high SNR, the b=0 column should sit just above 1
    # (Rician floor) but close to it.
    c = generate_cohort(500, 10, 10, snr_grid=(100,), seed=5)
    b0 = c.signals["train"][:, 0]
    assert 0.97 < b0.mean() < 1.03
