"""Continuity gates: every family reduces to bi-exponential at its limit.

Three families reduce *exactly* (residual == 0); gamma reduces asymptotically.
"""

import numpy as np
import pytest

from lattice.cohort import continuity_residual


@pytest.mark.parametrize("family", ["dispersion_lognormal", "stretched", "triexp"])
def test_exact_reduction_to_biexp(family):
    resid = continuity_residual(family, n=2000)
    assert resid == 0.0, f"{family} continuity residual {resid:.3e} != 0"


def test_gamma_asymptotic_reduction():
    resid = continuity_residual("dispersion_gamma", n=2000)
    assert resid < 1e-6, f"gamma k->inf residual {resid:.3e} not < 1e-6"


def test_biexp_identity():
    assert continuity_residual("biexp", n=500) == 0.0


def test_gamma_monotone_in_k():
    # Larger k -> smaller deviation from bi-exp.
    from lattice.cohort import sample_params, DEFAULT_BVALUES, DEFAULT_SEED, _forward
    rng = np.random.default_rng(DEFAULT_SEED)
    p = sample_params(800, rng)
    base = _forward("biexp", DEFAULT_BVALUES, p, {})
    devs = []
    for k in (2.0, 10.0, 100.0, 1e4):
        d = _forward("dispersion_gamma", DEFAULT_BVALUES, p, {"k": k})
        devs.append(np.max(np.abs(d - base)))
    assert all(devs[i] > devs[i + 1] for i in range(len(devs) - 1))
