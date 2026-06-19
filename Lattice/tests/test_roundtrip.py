"""Clean-signal round-trip gate: bi-exp truth is recovered within tolerance.

Requires scipy (the 'selfcheck' extra); skipped if unavailable.
"""

import numpy as np
import pytest

pytest.importorskip("scipy")

from lattice import make_cohort
from lattice.selfcheck import fit_biexp_nlls, clean_roundtrip_error
from lattice.cohort import DEFAULT_BVALUES as B


CASES = [
    (0.8e-3, 20e-3, 0.10),
    (1.2e-3, 40e-3, 0.20),
    (1.6e-3, 60e-3, 0.30),
    (2.0e-3, 80e-3, 0.15),
    (1.0e-3, 100e-3, 0.35),
]


@pytest.mark.parametrize("D,Dstar,f", CASES)
def test_clean_recovery_returns_truth(D, Dstar, f):
    from lattice.generators import ivim_biexp
    s = ivim_biexp(B, D, Dstar, f)
    D_hat, Dstar_hat, f_hat = fit_biexp_nlls(B, s)
    assert D_hat == pytest.approx(D, rel=1e-2)
    assert Dstar_hat == pytest.approx(Dstar, rel=1e-2)
    assert f_hat == pytest.approx(f, rel=1e-2)


def test_clean_recovery_grid_max_error_small():
    c = make_cohort("biexp", n=60, noise="none", seed=11)
    report = clean_roundtrip_error(c)
    assert report["max_rel_overall"] < 1e-2, report
