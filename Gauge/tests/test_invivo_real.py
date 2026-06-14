"""Tests for the ADDITIVE real-data in-vivo path (no network; tiny synthetic input).

These check the mechanics of the real-data code path without touching the
synthetic seed-20260613 pipeline:
  * the no-nibabel 4D loader masks + S0-normalizes + subsamples correctly;
  * the hybrid run applies the AS-DEPLOYED monitor (synthetic 22-value calibration)
    to real-scheme features and re-fits the CQR band at the real b-scheme;
  * the run makes NO coverage claim (has_gt is False, no phantom_coverage key).
A genuinely out-of-calibration 4-b "real-like" source is an exchangeability break,
so the as-deployed Mahalanobis monitor family fires.
"""
import numpy as np

from gauge.invivo import (load_dwi_npy4d, run_real_hybrid, synthetic_stand_in,
                          _find_retest_pairs)

REAL_B = np.array([0.0, 100.0, 600.0, 800.0])      # the ACRIN-6698 scheme


def _fake_4d(tmp_path, nx=12, ny=12, nz=4):
    """A small (X,Y,Z,4) volume: bright foreground block + dark air border."""
    rng = np.random.default_rng(0)
    vol = np.zeros((nx, ny, nz, 4), float)
    s0 = 800.0
    # mono-exponential-ish decay across the 4 b-values for the foreground block
    decay = np.exp(-REAL_B * 1.0e-3)
    vol[2:-2, 2:-2, :, :] = s0 * decay[None, None, None, :]
    vol += rng.normal(0, 5.0, vol.shape) ** 0  # keep deterministic shape; no-op add
    vp = tmp_path / "signals_4d.npy"
    bp = tmp_path / "bvals.txt"
    np.save(vp, vol.astype(np.float32))
    np.savetxt(bp, REAL_B, fmt="%g")
    return str(vp), str(bp)


def test_npy4d_loader_masks_and_normalizes(tmp_path):
    vp, bp = _fake_4d(tmp_path)
    signals, b = load_dwi_npy4d(vp, bp, max_voxels=10_000, seed=0)
    assert np.allclose(b, REAL_B)
    assert signals.shape[1] == 4
    # b=0 foreground mask rejects the dark air border (8x8 block kept of 12x12x4)
    assert signals.shape[0] == 8 * 8 * 4
    # every kept voxel is normalized to S0=1 at b=0
    assert np.allclose(signals[:, 0], 1.0, atol=1e-6)


def test_npy4d_loader_subsamples(tmp_path):
    vp, bp = _fake_4d(tmp_path)
    signals, _ = load_dwi_npy4d(vp, bp, max_voxels=50, seed=1)
    assert signals.shape[0] == 50


def test_hybrid_runs_no_coverage_claim_and_monitor_fires():
    # an out-of-calibration 4-b "real-like" source (broadened D*/f prior, low SNR)
    sig22, _, _ = synthetic_stand_in(n=300)
    # resample onto the real 4-b scheme by re-simulating at REAL_B from the same
    # out-of-calibration draw is overkill here; instead use a direct 4-b draw:
    from gauge.robustness import _draw_params, _draw_snr, _simulate
    rng = np.random.default_rng(7)
    params = _draw_params(300, rng, prior={"Dstar": (10e-3, 130e-3), "f": (0.05, 0.45)})
    snr = _draw_snr(300, rng, (8.0, 12.0, 18.0))
    sig = _simulate(REAL_B, params, snr, rng, model="triexp", triexp=(4.0, 0.20))
    sig = np.clip(sig / sig[:, :1], 0.0, None)             # S0-normalize

    out, (lo, hi, widths, theta) = run_real_hybrid(
        sig, REAL_B, "test 4-b real-like", seed=20260613, n_train=300, n_cal=300)

    # NO coverage claim: no ground truth is used or reported
    assert out["has_gt"] is False
    assert "phantom_coverage" not in out
    # qualitative band widths are positive and finite
    assert np.all(widths > 0) and np.all(np.isfinite(out["dstar_width_q"]))
    assert out["dstar_widths"].shape == (300,)
    # the as-deployed b-independent (Mahalanobis) family fires on the transfer
    assert out["monitor_maha_fires"] is True
    assert out["monitor_maha_auc"] > 0.6
    assert out["monitor_fires"] is True
    # real b-scheme is recorded for provenance
    assert out["b_real"] == [0.0, 100.0, 600.0, 800.0]


def test_loader_tumor_mask_rejects_degenerate_voxels(tmp_path):
    # a 3x3x1x4 volume where one masked voxel is pure air (all-zero -> S0=0)
    vol = np.zeros((3, 3, 1, 4), float)
    vol[0, 0, 0, :] = 800.0 * np.exp(-REAL_B * 1e-3)   # one good voxel
    # vol[1,1,0,:] stays all-zero (degenerate)
    mask = np.zeros((3, 3, 1), bool)
    mask[0, 0, 0] = True
    mask[1, 1, 0] = True                                # includes the air voxel
    np.save(tmp_path / "signals_4d.npy", vol.astype(np.float32))
    np.savetxt(tmp_path / "bvals.txt", REAL_B, fmt="%g")
    np.save(tmp_path / "tumor_mask.npy", mask)
    sig, _ = load_dwi_npy4d(str(tmp_path / "signals_4d.npy"),
                            str(tmp_path / "bvals.txt"),
                            mask_path=str(tmp_path / "tumor_mask.npy"))
    assert sig.shape[0] == 1                            # degenerate voxel dropped
    assert np.isfinite(sig).all()


def test_find_retest_pairs(tmp_path):
    root = tmp_path / "invivo"
    for patient, n_exams in [("P1", 2), ("P2", 2), ("P3", 1)]:
        for k in range(n_exams):
            d = root / patient / f"study{k}"
            d.mkdir(parents=True)
            np.save(d / "signals_4d.npy", np.ones((2, 2, 1, 4), np.float32))
            np.save(d / "tumor_mask.npy", np.ones((2, 2, 1), bool))
    pairs = _find_retest_pairs(str(root))
    patients = {p for p, _ in pairs}
    assert patients == {"P1", "P2"}                    # P3 (1 exam) excluded
    assert all(len(exs) == 2 for _, exs in pairs)
