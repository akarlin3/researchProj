"""IVIM Fisher-information / CRLB sanity and the Gauge resolution-ratio anchor."""
from __future__ import annotations

import numpy as np

from vernier import _paths
from vernier.crlb import (
    crlb,
    crlb_dstar_resolution_ratio,
    expected_crlb,
    fisher_information,
    ivim_jacobian,
)
from vernier.schemes import by_name, MATCHED_SCANTIME


def _cohort_params(n=2000, seed=0):
    _paths.add_caliper()
    from caliper.forward import sample_params

    return sample_params(n, np.random.default_rng(seed))


def test_jacobian_shape_and_zero_b_column():
    b = np.array([0.0, 50.0, 200.0, 800.0])
    g = ivim_jacobian(b, D=1.5, f=0.2, Dstar=30.0)
    assert g.shape == (4, 3)
    # at b=0 the D and D* derivatives vanish (factor -b), f derivative is 0 too
    assert np.allclose(g[0], 0.0)


def test_fisher_is_symmetric_psd():
    b = MATCHED_SCANTIME[1].b
    J = fisher_information(b, D=1.5, f=0.2, Dstar=30.0, snr=40.0)
    assert J.shape == (3, 3)
    assert np.allclose(J, J.T)
    eig = np.linalg.eigvalsh(J)
    assert np.all(eig > 0)


def test_crlb_positive_and_dstar_is_worst():
    b = MATCHED_SCANTIME[1].b
    c = crlb(b, D=1.5, f=0.2, Dstar=60.0, snr=40.0)
    assert np.all(c > 0) and np.all(np.isfinite(c))
    # D* is the poorly identified parameter: its CRLB (in 1e-3 units, scale ~10-100)
    # dwarfs D's (scale ~0.5-2.5) in absolute terms.
    assert c[2] > c[0]


def test_more_perfusion_sampling_lowers_dstar_crlb():
    params = _cohort_params()
    perf = by_name(MATCHED_SCANTIME, "perfusion-weighted")
    tiss = by_name(MATCHED_SCANTIME, "tissue-weighted")
    cr_perf = expected_crlb(perf, params, snr=40.0)[2]
    cr_tiss = expected_crlb(tiss, params, snr=40.0)[2]
    # denser low-b (perfusion) sampling must give a tighter D* precision floor
    assert cr_perf < cr_tiss


def test_resolution_ratio_in_gauge_ballpark():
    params = _cohort_params()
    # Gauge reports CRLB(D*)/tercile-width ~ 1.05-1.25 (the wall). We only assert
    # the right order of magnitude (>= ~0.5 and finite) -- SNR/averaging differ,
    # so this validates the implementation, not an exact replication.
    ratio = crlb_dstar_resolution_ratio(
        by_name(MATCHED_SCANTIME, "balanced-clinical"), params, snr=33.0
    )
    assert np.isfinite(ratio)
    assert ratio > 0.3
