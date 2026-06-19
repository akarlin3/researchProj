"""The feasibility gate runs, is deterministic, and returns a pre-registered verdict."""
from __future__ import annotations

import numpy as np

from vernier import _paths
from vernier.feasibility import run_gate, select_matched_crlb
from vernier.schemes import CANDIDATE_POOL


def _params(n=3000, seed=0):
    _paths.add_caliper()
    from caliper.forward import sample_params

    return sample_params(n, np.random.default_rng(seed))


def test_select_matched_crlb_returns_matched_set():
    params = _params()
    matched = select_matched_crlb(CANDIDATE_POOL, params, snr=33.0, tol=0.10, min_keep=3)
    # selection is against the POOL median (the same basis the function uses)
    pool_med = np.median([_crlb(s, params) for s in CANDIDATE_POOL])
    cr = np.array([_crlb(s, params) for s in matched])
    assert np.all(np.abs(cr - pool_med) / pool_med <= 0.10 + 1e-9)
    assert len(matched) >= 3
    # all matched schemes share scan-time (matched-scan-time pool is all 11 b)
    assert len({s.n_b for s in matched}) == 1


def _crlb(scheme, params):
    from vernier import crlb

    return float(crlb.expected_crlb(scheme, params, 33.0)[2])


def test_gate_runs_and_is_deterministic():
    params = _params()
    matched = select_matched_crlb(CANDIDATE_POOL, params, snr=33.0)
    g1 = run_gate(matched, n=1500, snr=33.0, seed=0, n_boot=200)
    g2 = run_gate(matched, n=1500, snr=33.0, seed=0, n_boot=200)
    assert g1.verdict in ("PASS", "FAIL")
    assert g1.verdict == g2.verdict
    assert np.isclose(g1.delta_sharp, g2.delta_sharp)
    assert np.isclose(g1.delta_cond, g2.delta_cond)
    # CIs are ordered and gaps are non-negative (range statistics)
    assert g1.delta_sharp >= 0 and g1.delta_cond >= 0
    assert g1.delta_sharp_ci[0] <= g1.delta_sharp_ci[1]


def test_marginal_coverage_is_restored_by_conformal():
    # sanity: conformal restores ~nominal marginal D* coverage for every scheme
    params = _params()
    matched = select_matched_crlb(CANDIDATE_POOL, params, snr=33.0)
    g = run_gate(matched, n=4000, snr=33.0, seed=0, n_boot=50)
    for s in g.schemes:
        assert abs(s.cov_dstar - 0.90) < 0.05, f"{s.name} marginal cov {s.cov_dstar}"
        # raw (pre-conformal) coverage is far below nominal (over-confident)
        assert s.cov_dstar_raw < 0.6
