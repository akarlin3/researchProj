"""Twin reproducibility, geometry, and the synthetic forward model."""
import numpy as np

from matrix import MatrixConfig, Twin, TUMOR, OAR, NORMAL
from matrix.forward import ivim_signal, simulate_scan


def test_twin_bit_reproducible_from_seed():
    cfg = MatrixConfig()
    a, b = Twin.build(cfg), Twin.build(cfg)
    for fld in ("D", "Dstar", "f", "labels", "highdstar", "lowsnr", "snr_map", "dose"):
        assert np.array_equal(getattr(a, fld), getattr(b, fld))


def test_twin_differs_by_seed():
    cfg = MatrixConfig()
    a = Twin.build(cfg)
    c = Twin.build(cfg.with_seed(cfg.seed + 1))
    assert not np.array_equal(a.f, c.f)


def test_twin_has_all_tissue_classes_and_untrust_overlaps_tumor():
    twin = Twin.build(MatrixConfig())
    for lab in (NORMAL, TUMOR, OAR):
        assert np.any(twin.labels == lab)
    # the untrustworthy (low-SNR) zone must clip the tumour, else the trust gate
    # has no real TREAT decision to suppress.
    assert np.any((twin.labels == TUMOR) & twin.lowsnr)
    assert np.any((twin.labels == TUMOR) & ~twin.lowsnr)


def test_lowsnr_zone_has_lower_snr():
    twin = Twin.build(MatrixConfig())
    assert twin.snr_map[twin.lowsnr].max() < twin.snr_map[~twin.lowsnr].min()


def test_ivim_signal_shape_and_monotonicity():
    bvals = np.array([0.0, 50.0, 200.0, 800.0])
    s = ivim_signal(bvals, D=[1e-3], Dstar=[20e-3], f=[0.2])
    assert s.shape == (4, 1)
    assert np.all(np.diff(s[:, 0]) < 0)          # signal decays with b
    assert np.isclose(s[0, 0], 1.0)              # S(b=0) == S0


def test_simulate_scan_shape_and_per_voxel_snr():
    cfg = MatrixConfig()
    rng = np.random.default_rng(0)
    snr_map = np.array([cfg.snr, cfg.snr_low])
    scan = simulate_scan(cfg.bvals, [1e-3, 1e-3], [20e-3, 20e-3], [0.2, 0.2],
                         snr_map, cfg.n_noise, rng)
    assert scan.shape == (len(cfg.bvals), 2, cfg.n_noise)
    # the low-SNR voxel must be noisier at b0
    assert scan[0, 1].std() > scan[0, 0].std()


def test_apply_plan_devascularises_boosted_voxels():
    twin = Twin.build(MatrixConfig())
    f_before = twin.f.copy()
    new_dose = twin.dose.copy()
    new_dose[0] = twin.cfg.dose_baseline + twin.cfg.dose_boost   # full boost
    twin.apply_plan(new_dose)
    assert twin.f[0] < f_before[0]               # boosted voxel loses perfusion
    assert np.isclose(twin.f[1], f_before[1])    # untouched voxel unchanged
