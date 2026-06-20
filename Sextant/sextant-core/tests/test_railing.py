import numpy as np
import pytest

from conftest import extract_dir, osipi_available
from sextant.cohorts import load_nifti_cohort
from sextant.railing import analyze_cohort, fit_dstar, rail_mask


def test_rail_mask_lower_upper_and_nan():
    d = np.array([0.001, 0.20, 0.05, np.nan])  # lower-rail, upper-rail, interior, nan
    m, lo, hi = rail_mask(d, lower=0.0033, upper=0.1485)
    assert list(m) == [True, True, False, False]
    assert list(lo) == [True, False, False, False]
    assert list(hi) == [False, True, False, False]


def test_fit_dstar_shape_and_finiteness():
    R_b = [0, 10, 20, 30, 50, 75, 100, 150, 400, 600]
    sigs = np.array([
        0.85 * np.exp(-np.array(R_b) * 1.2e-3) + 0.15 * np.exp(-np.array(R_b) * 0.03)
        for _ in range(5)])
    out = fit_dstar(sigs)
    assert out.shape == (5,)
    assert np.isfinite(out).all()


@pytest.mark.data
def test_homogeneous_cohort_reproduces_fashion():
    if not osipi_available():
        pytest.skip("OSIPI data not downloaded; run scripts/fetch_osipi.py")
    D = extract_dir()
    coh = load_nifti_cohort("abdomen_homogeneous", f"{D}/abdomen.nii.gz",
                            f"{D}/mask_abdomen_homogeneous.nii.gz", f"{D}/abdomen.bval")
    assert coh.n_high_snr == 1618           # Fashion's reported high-SNR count
    r = analyze_cohort(coh, bounds="tight")
    assert abs(r.frac_railed - 0.547) < 0.005   # Fashion's reported 54.7%
    # wide bounds must still rail a substantial minority (not a tight-bounds artefact)
    rw = analyze_cohort(coh, bounds="wide")
    assert rw.frac_railed > 0.20
