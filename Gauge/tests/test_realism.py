"""Tests for the realism cohort (realistic joint prior + measurement nuisance).

Fast tests only: the nuisance forward physics, the byte-identical continuity to the
existing bi-exponential cohort, the realistic-prior shape (correlated / skewed), the
fixed-vs-per-distribution tercile policy, and determinism of the draws. The full
Arm-A / Arm-B pipeline is exercised by the seeded report + multiseed sweep, not in
unit tests.
"""
import numpy as np
import pytest

from gauge.cohort import generate_cohort, DSTAR_RANGE
from gauge.conditional_attack import _regime_from_true
from gauge.forward import (ivim_signal, add_rician_noise, add_ncchi_noise,
                           partial_volume_mix, apply_motion_phase,
                           DEFAULT_B_VALUES)
import gauge.realism as R

B = DEFAULT_B_VALUES


# --------------------------------------------------------------------------- #
# Continuity to the existing bi-exponential cohort (the built-in gate).
# --------------------------------------------------------------------------- #
def test_ncchi_L1_equals_rician_byte_identical():
    """Non-central chi at L=1 reduces to Rician byte-for-byte (same rng order)."""
    S = ivim_signal(B[None, :], np.array([[1.2e-3]]), np.array([[30e-3]]),
                    np.array([[0.2]]))
    snr = np.array([[30.0]])
    r1 = np.random.default_rng(5)
    r2 = np.random.default_rng(5)
    assert np.array_equal(add_rician_noise(S, snr, r1),
                          add_ncchi_noise(S, snr, r2, L=1))


def test_ncchi_more_coils_lifts_floor():
    """More coils (L>1) raises the magnitude noise floor (non-Rician)."""
    S = np.zeros((4000, 1))                       # pure-noise voxels: see the floor
    snr = np.full((4000, 1), 20.0)
    r1 = np.random.default_rng(0)
    r8 = np.random.default_rng(0)
    m1 = add_ncchi_noise(S, snr, r1, L=1).mean()
    m8 = add_ncchi_noise(S, snr, r8, L=8).mean()
    assert m8 > m1


def test_partial_volume_zero_is_identity():
    """frac=0 partial volume returns the input unchanged."""
    S = ivim_signal(B[None, :], np.array([[1.2e-3]]), np.array([[30e-3]]),
                    np.array([[0.2]]))
    assert np.array_equal(partial_volume_mix(S, B, 0.0), S)


def test_partial_volume_reshapes_signal_keeps_b0():
    """Free-fluid mixing reshapes the decay (lifts low-b, lowers high-b) but leaves
    b=0 unchanged -- a genuine partial-volume confound on the bi-exp fit."""
    S = ivim_signal(B[None, :], np.array([[1.5e-3]]), np.array([[30e-3]]),
                    np.array([[0.2]]))
    mix = partial_volume_mix(S, B, 0.3)
    assert mix[0, 0] == pytest.approx(S[0, 0])    # b=0 untouched (both = 1.0)
    assert np.max(np.abs(mix - S)) > 1e-3         # the decay shape is contaminated
    assert mix[0, -1] < S[0, -1]                  # fast free water lowers high-b tail


def test_motion_phase_zero_is_identity():
    """phase=0 leaves the (real) signal unchanged."""
    S = ivim_signal(B, 1.2e-3, 30e-3, 0.2)
    assert np.allclose(apply_motion_phase(S, 0.0).real, S)
    assert np.allclose(apply_motion_phase(S, 0.0).imag, 0.0)


def test_nuisance_eta0_equals_rician_byte_identical():
    """nuisance_signal(eta=0) == add_rician_noise(ivim_signal) byte-for-byte."""
    params = np.array([[1.2e-3, 30e-3, 0.2], [0.8e-3, 60e-3, 0.1]])
    snr = np.array([30.0, 20.0])
    clean = ivim_signal(B[None, :], params[:, 0:1], params[:, 1:2], params[:, 2:3])
    r1 = np.random.default_rng(7)
    r2 = np.random.default_rng(7)
    ref = add_rician_noise(clean, snr[:, None], r1)
    got = R.nuisance_signal(B, params, snr, r2, 0.0)
    assert np.array_equal(ref, got)


def test_continuity_gate_uniform_eta0_reproduces_cohort():
    """Uniform prior + zero nuisance reproduces generate_cohort to numerical zero."""
    cont = R.continuity_gate(seed=20260613)
    assert cont["uniform_eta0_max_abs_err"] == 0.0
    assert cont["clean_eta0_max_abs_err"] == 0.0


def test_build_cohort_uniform_eta0_byte_identical_to_generate_cohort():
    """The full uniform+eta=0 build is byte-identical, split by split."""
    sizes = (R.N_TRAIN_A, R.N_CAL_A, R.N_TEST_A)
    coh = generate_cohort(*sizes, seed=20260613)
    sig, par, sn = R._build_cohort(20260613, "uniform", 0.0, sizes, B)
    for s in ("train", "cal", "test"):
        assert np.array_equal(sig[s], coh.signals[s])
        assert np.array_equal(par[s], coh.params[s])
        assert np.array_equal(sn[s], coh.snr[s])


# --------------------------------------------------------------------------- #
# Realistic prior shape: correlated + skewed, D* CV in the published band.
# --------------------------------------------------------------------------- #
def test_published_prior_deterministic():
    """The published prior is a pure function of the rng seed."""
    p1, s1 = R.sample_published_prior(500, np.random.default_rng(20260613))
    p2, s2 = R.sample_published_prior(500, np.random.default_rng(20260613))
    assert np.array_equal(p1, p2) and np.array_equal(s1, s2)


def test_published_prior_skew_and_cv():
    """D* carries a within-cohort CV in the published 50-110% band and is skewed."""
    p, snr = R.sample_published_prior(40000, np.random.default_rng(0))
    D, Dstar, f = p[:, 0], p[:, 1], p[:, 2]
    cv_dstar = Dstar.std() / Dstar.mean()
    assert 0.5 <= cv_dstar <= 1.1                     # within-subject D* CV band
    # log-normal D* is right-skewed: mean > median.
    assert Dstar.mean() > np.median(Dstar)
    # f stays a valid fraction.
    assert f.min() > 0.0 and f.max() < 1.0
    # SNR is realistic and clipped.
    assert snr.min() >= R._PUB["snr_clip"][0] - 1e-9
    assert snr.max() <= R._PUB["snr_clip"][1] + 1e-9


def test_published_prior_is_correlated():
    """The joint is genuinely correlated (D*,f positive; D,f negative)."""
    p, _ = R.sample_published_prior(40000, np.random.default_rng(1))
    D, Dstar, f = p[:, 0], p[:, 1], p[:, 2]
    assert np.corrcoef(Dstar, f)[0, 1] > 0.15        # perfusion params co-vary +
    assert np.corrcoef(D, f)[0, 1] < -0.05           # restricted tissue -> lower f


def test_uniform_source_matches_cohort_ranges():
    """The 'uniform' source draws within the published physiological ranges."""
    p, snr = R._draw_prior(5000, np.random.default_rng(0), "uniform")
    assert p[:, 1].min() >= DSTAR_RANGE[0] and p[:, 1].max() <= DSTAR_RANGE[1]


# --------------------------------------------------------------------------- #
# Nuisance magnitude monotonicity.
# --------------------------------------------------------------------------- #
def test_nuisance_departs_more_with_eta():
    """Larger eta drives the noisy signal further from the clean bi-exp (in mean)."""
    rng_p = np.random.default_rng(3)
    params, snr = R._draw_prior(3000, rng_p, "uniform")
    clean = ivim_signal(B[None, :], params[:, 0:1], params[:, 1:2], params[:, 2:3])
    prev = -1.0
    for eta in (0.0, 0.4, 0.8):
        sig = R.nuisance_signal(B, params, snr, np.random.default_rng(9), eta)
        dep = float(np.mean(np.abs(sig - clean)))
        assert dep >= prev                            # monotone non-decreasing
        prev = dep


# --------------------------------------------------------------------------- #
# Tercile policy: fixed vs per-distribution boundaries (item 4).
# --------------------------------------------------------------------------- #
def test_fixed_edges_are_uniform_quantiles():
    """Fixed edges are the uniform DSTAR_RANGE 1/3, 2/3 quantiles."""
    lo, hi = DSTAR_RANGE
    assert R.FIXED_DSTAR_EDGES == (pytest.approx(lo + (hi - lo) / 3),
                                   pytest.approx(lo + 2 * (hi - lo) / 3))


def test_regime_fixed_vs_perdist_differ_on_skewed_prior():
    """Fixed boundaries hold the physical regime; per-dist boundaries move with it."""
    p, _ = R.sample_published_prior(20000, np.random.default_rng(2))
    dstar = p[:, 1]
    _, edges_perd = _regime_from_true(dstar)                       # quantile edges
    _, edges_fixed = _regime_from_true(dstar, edges=R.FIXED_DSTAR_EDGES)
    assert np.allclose(edges_fixed, R.FIXED_DSTAR_EDGES)
    # the skewed prior's upper quantile edge sits below the fixed physical edge.
    assert edges_perd[-1] < R.FIXED_DSTAR_EDGES[-1]
    # fixed hi-bin is rarer than the per-distribution hi tercile (~1/3 by design).
    prev_fixed = float(np.mean(dstar >= R.FIXED_DSTAR_EDGES[-1]))
    assert prev_fixed < 1.0 / 3.0


# --------------------------------------------------------------------------- #
# NF1 calibration-size control: up-sample only the hi-D* CALIBRATION count, at
# fixed (realistic) within-bin distribution. Fast mechanism tests only.
# --------------------------------------------------------------------------- #
def test_upsample_hi_dstar_only_fixed_bin_and_deterministic():
    """The up-sampler returns exactly n_target realistic voxels, all in the fixed
    hi-D* bin, and is a pure function of the rng."""
    n = 600
    p1, s1 = R._draw_hi_dstar_realistic(n, np.random.default_rng(909))
    p2, s2 = R._draw_hi_dstar_realistic(n, np.random.default_rng(909))
    assert p1.shape == (n, 3) and s1.shape == (n,)
    assert np.array_equal(p1, p2) and np.array_equal(s1, s2)        # deterministic
    assert np.all(p1[:, 1] >= R.FIXED_DSTAR_EDGES[-1])             # all hi-D*
    assert R._draw_hi_dstar_realistic(0, np.random.default_rng(0))[0].shape == (0, 3)


def test_upsampled_hi_matches_native_realistic_not_uniform():
    """Up-sampled hi-D* voxels share the NATIVE realistic hi-D* law, not the uniform
    one: the (D, D*, f) marginals are KS-closer to a realistic hi-D* reference than to
    a uniform hi-D* reference (the fidelity gate's discriminator)."""
    edge = R.FIXED_DSTAR_EDGES[-1]
    up, _ = R._draw_hi_dstar_realistic(3000, np.random.default_rng(909))
    ref_real, _ = R._draw_hi_dstar_realistic(12000, np.random.default_rng(123))
    # uniform hi-D* reference.
    up_u, have, acc = None, 0, []
    rng = np.random.default_rng(456)
    while have < 12000:
        p, _ = R._draw_prior(8000, rng, "uniform")
        m = p[:, 1] >= edge
        acc.append(p[m]); have += int(m.sum())
    ref_unif = np.concatenate(acc)[:12000]
    from scipy.stats import ks_2samp
    for j in range(3):
        d_real = ks_2samp(up[:, j], ref_real[:, j]).statistic
        d_unif = ks_2samp(up[:, j], ref_unif[:, j]).statistic
        assert d_real < d_unif                          # closer to realistic


def test_fidelity_gate_passes_on_realistic_upsample():
    """The fidelity gate PASSES for genuine realistic-generator hi-D* draws and FAILS
    if the up-sampled voxels are (illegitimately) uniform hi-D* voxels."""
    up_real, _ = R._draw_hi_dstar_realistic(2000, np.random.default_rng(909))
    assert R._fidelity_gate(up_real, seed=20260613)["passes"] is True
    # illegitimate substitute: uniform hi-D* voxels must NOT pass.
    edge = R.FIXED_DSTAR_EDGES[-1]
    rng, acc, have = np.random.default_rng(99), [], 0
    while have < 2000:
        p, _ = R._draw_prior(8000, rng, "uniform")
        m = p[:, 1] >= edge
        acc.append(p[m]); have += int(m.sum())
    up_unif = np.concatenate(acc)[:2000]
    assert R._fidelity_gate(up_unif, seed=20260613)["passes"] is False


def test_control_verdict_branches():
    """The control verdict thresholds map gap-closed fraction -> A / PARTIAL / B."""
    # native 0.50, uniform anchor 0.80 -> gap 0.30.
    assert R._control_verdict(0.50, 0.80, 0.80)["branch"] == "A"      # 100% closed
    assert R._control_verdict(0.50, 0.515, 0.80)["branch"] == "B"     # 5% closed
    assert R._control_verdict(0.50, 0.62, 0.80)["branch"] == "PARTIAL"  # 40% closed
