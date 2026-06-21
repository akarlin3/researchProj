"""CP2 self-consistency gates: the rebuild must recover truth on clean signals.

These are the "runs, self-consistent" gate for CP2: clean-signal round-trip recovers
truth (forward + NLLS), the analytic Jacobian matches finite differences, and the
independent ruler computes coverage exactly on a constructed case.
"""
from __future__ import annotations

import numpy as np
import pytest

from gnomon import forward as F
from gnomon import metrics as M
from gnomon import cohort


B = np.array([0, 10, 20, 40, 80, 200, 400, 800], dtype=float)  # clinical-sparse 8b


def test_forward_b0_is_s0():
    s = F.ivim(B, D=1.2e-3, Dstar=30e-3, f=0.2, s0=1.0)
    assert np.isclose(s[0], 1.0)  # at b=0 both exponentials are 1 -> S0


def test_analytic_jacobian_matches_finite_difference():
    p = F.to_fit(D=1.2e-3, Dstar=25e-3, f=0.2, s0=1.0)  # scaled (S0,D3,f,Ds3)
    J = F.ivim_jac_scaled(p, B)
    eps = 1e-6
    Jfd = np.empty_like(J)
    for k in range(4):
        pp = p.copy(); pm = p.copy()
        pp[k] += eps; pm[k] -= eps
        Jfd[:, k] = (F.ivim_design(pp, B) - F.ivim_design(pm, B)) / (2 * eps)
    assert np.allclose(J, Jfd, atol=1e-6, rtol=1e-4)


def test_clean_round_trip_recovers_truth():
    from gnomon.nlls import NLLSEstimator
    truths = np.array([[1.0e-3, 20e-3, 0.15],
                       [1.5e-3, 40e-3, 0.25],
                       [0.8e-3, 60e-3, 0.10]])
    clean = F.ivim(B, truths[:, 0], truths[:, 1], truths[:, 2], s0=1.0)
    est = NLLSEstimator(B)
    fit = est.fit(clean)  # noise-free
    # D and f recover tightly; Dstar a touch looser (it is the hard one) but close.
    assert np.allclose(fit.params[:, 0], truths[:, 0], rtol=1e-3, atol=1e-5), "D"
    assert np.allclose(fit.params[:, 2], truths[:, 2], rtol=2e-2, atol=1e-2), "f"
    assert np.allclose(fit.params[:, 1], truths[:, 1], rtol=5e-2, atol=5e-3), "Dstar"
    assert not fit.dstar_railed.any()  # clean signals should not rail


def test_continuity_no_perfusion_is_monoexponential():
    # f = 0 -> the fast compartment vanishes -> pure mono-exponential exp(-bD).
    s = F.ivim(B, D=1.3e-3, Dstar=40e-3, f=0.0, s0=1.0)
    assert np.allclose(s, np.exp(-B * 1.3e-3), atol=1e-12)


def test_cohort_uses_lattice_substrate():
    co = cohort.make_headline_cohort(n_noise=5)
    assert co.substrate == "lattice"
    assert co.n == 3 * 3 * 5          # 3 truths x 3 SNR x 5 noise
    assert co.signals.shape == (co.n, len(B))
    # b0 column ~ 1 after the S0=1 forward + noise (within a few sigma)
    assert abs(co.signals_clean[:, 0].mean() - 1.0) < 1e-9


def test_driver_runs_end_to_end_and_emits_verdict(tmp_path):
    # Fast path: no flow, no real-data fetch; tiny MCMC. Exercises the full
    # synthetic pipeline + verdict machinery and writes results JSON.
    from gnomon import reproduce
    res = reproduce.run(n_noise=20, run_flow=False, run_real=False,
                        mcmc={"burn": 300, "keep": 400, "thin": 1}, verbose=False,
                        out_dir=tmp_path)
    assert res["verdict"] in ("REPRODUCES", "PARTIAL", "DOES NOT REPRODUCE")
    assert (tmp_path / "reproduction.json").exists()
    # The shape-correct quantile interval should recover near-nominal D* coverage
    # (loose bound here: this is a tiny noisy smoke run, not the CP3 gate).
    assert res["synthetic"]["T3c_mcmc_quantile_Dstar"]["coverage"] > 0.70


def test_ruler_coverage_is_exact_on_constructed_case():
    # 4 points, central 0.5 interval [lo,hi]; 3 of 4 truths inside -> coverage 0.75.
    q_levels = np.array([0.25, 0.5, 0.75])
    y = np.array([[0.0], [0.5], [2.0], [0.4]])           # (n,1)
    qp = np.zeros((4, 1, 3))
    qp[:, 0, 0] = -1.0  # lo (0.25)
    qp[:, 0, 1] = 0.0
    qp[:, 0, 2] = 1.0   # hi (0.75)
    scores = M.score_quantiles(y, qp, q_levels, alpha=0.5, param_names=["x"])
    assert scores[0].coverage == pytest.approx(0.75)
    assert scores[0].nominal == pytest.approx(0.5)
