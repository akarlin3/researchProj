"""The secondary ruler is reused read-only from Fashion and is scoped to data
with ground truth (it cannot be applied to the real abdomen)."""
import numpy as np

from sextant import ruler


def test_ruler_requires_ground_truth():
    # Documents, in code, why railing (no truth needed) is the stronger primary.
    assert ruler.requires_ground_truth() is True


def test_well_calibrated_gaussian_has_near_nominal_coverage():
    rng = np.random.default_rng(0)
    truth = 0.03
    sigma = 0.01
    estimates = truth + rng.normal(0, sigma, size=20000)
    cov = ruler.coverage(estimates, truth, sigma)
    # empirical coverage should track the nominal levels closely
    for level, emp in cov.items():
        assert abs(emp - level) < 0.03
    assert ruler.ece(cov) < 0.02


def test_levels_exposed():
    assert 0.95 in list(ruler.levels())
