"""The segmented IVIM fit: recovers the QoI (f) and inflates sigma_f under low SNR."""
import numpy as np

from matrix import MatrixConfig, Twin, TUMOR
from matrix.fit import fit_scan


def test_fit_recovers_f_with_small_bias():
    cfg = MatrixConfig()
    twin = Twin.build(cfg)
    mu, sig = fit_scan(twin.scan(np.random.default_rng(1)), cfg)
    good = ~twin.lowsnr
    bias = np.mean(mu["f"][good] - twin.f[good])
    assert abs(bias) < 0.03                       # f is well-recovered at good SNR


def test_fit_sigma_f_inflates_in_lowsnr_zone():
    cfg = MatrixConfig()
    twin = Twin.build(cfg)
    _, sig = fit_scan(twin.scan(np.random.default_rng(2)), cfg)
    assert sig["f"][twin.lowsnr].mean() > 1.5 * sig["f"][~twin.lowsnr].mean()


def test_fit_outputs_all_params_finite():
    cfg = MatrixConfig()
    twin = Twin.build(cfg)
    mu, sig = fit_scan(twin.scan(np.random.default_rng(3)), cfg)
    for k in ("D", "Dstar", "f"):
        assert np.all(np.isfinite(mu[k])) and np.all(np.isfinite(sig[k]))
        assert mu[k].shape == (cfg.n_voxels,)
