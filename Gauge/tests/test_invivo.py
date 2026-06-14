"""Tests for the CP2 in-vivo demonstration pipeline.

The in-vivo demo makes NO coverage claims (no ground truth in vivo); these tests
check the pipeline mechanics on tiny inputs: S0 normalization, the generic
loader, and that the deployed pipeline produces positive interval widths and a
monitor that fires on an out-of-calibration source (the synthetic->deployment
exchangeability break that forbids transferring the guarantee).
"""
import numpy as np

from gauge.forward import DEFAULT_B_VALUES
from gauge.invivo import (_normalize_to_s0, load_signals_npy, synthetic_stand_in,
                          deployed_calibration, run_demo)


def test_normalize_to_s0_sets_b0_to_one():
    b = np.array([0.0, 50.0, 200.0, 800.0])
    sig = np.array([[2.0, 1.6, 1.0, 0.4], [10.0, 8.0, 5.0, 2.0]])
    norm = _normalize_to_s0(sig, b)
    assert np.allclose(norm[:, 0], 1.0)
    assert np.all(norm >= 0)


def test_load_signals_npy_roundtrip(tmp_path):
    b = np.array([0.0, 100.0, 800.0])
    sig = np.array([[1.0, 0.6, 0.2], [2.0, 1.2, 0.4]])
    sp = tmp_path / "sig.npy"
    bp = tmp_path / "b.txt"
    np.save(sp, sig)
    np.savetxt(bp, b)
    loaded, bb = load_signals_npy(str(sp), str(bp))
    assert loaded.shape == (2, 3)
    assert np.allclose(bb, b)
    assert np.allclose(loaded[:, 0], 1.0)               # normalized to S0


def test_stand_in_is_labeled_and_shaped():
    sig, b, params = synthetic_stand_in(n=50)
    assert sig.shape == (50, len(DEFAULT_B_VALUES))
    assert params.shape == (50, 3)
    assert np.all(np.isfinite(sig))


def test_demo_pipeline_runs_and_monitor_fires_on_transfer():
    cal = deployed_calibration(n_train=300, n_cal=300)
    sig, b, params = synthetic_stand_in(n=250)
    r = run_demo(sig, b, cal, "test stand-in", true_params=params)
    # positive, finite interval widths
    assert np.all(r["width_median"] > 0)
    assert r["dstar_widths"].shape == (250,)
    assert np.all(np.isfinite(r["dstar_width_q"]))
    # the out-of-calibration stand-in is an observable break -> monitor fires
    assert r["monitor_fires"]
    assert r["monitor_auc"] > 0.6
    # phantom coverage is reported (separately) and is a 3-vector
    assert r["phantom_coverage"].shape == (3,)
