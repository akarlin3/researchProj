"""
Tier-2 test suite.

Fast pure-numeric and estimator-free tests run by default; one small uGUIDE
integration test (a few epochs) exercises the Arm-2 path. No network or external
data is required -- the real-data ingestion is tested via synthetic NIfTI volumes.
"""
import os

import numpy as np
import torch
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

from sibyl.forward_model.ivim import (
    ivim_biexponential, adc_monoexp_fit, ACRIN_B_SCHEME, DEFAULT_B_SCHEME, ACRIN_TO_DENSE_IDX,
)
from sibyl.data.synthetic import sample_prior
from sibyl.data.acrin_reference import generate_acrin_id_reference, normalize_by_b0, signal_features
from sibyl.detectors.signal_space import SignalSpaceDetector
from sibyl.data.imputation import segmented_ivim_fit, impute_dense
from sibyl.data.units import synthetic_id_units, synthetic_invivo_units, concat_units, UnitTable
from sibyl.data.acrin_ingest import reconcile_b0, roi_mean_per_b, assemble_4b_vector
from sibyl.metrics.eval import adc_repeatability, spearman_coupling, partial_spearman, degeneracy_fpr


# ----------------------------- forward model / ADC -----------------------------
def test_adc_monoexp_recovers_known_adc():
    adc_true = np.array([1.2e-3, 0.8e-3, 1.6e-3, 2.0e-3])
    S = torch.exp(-ACRIN_B_SCHEME.view(1, -1) * torch.tensor(adc_true).view(-1, 1))
    _, adc_hat = adc_monoexp_fit(S, ACRIN_B_SCHEME)
    assert np.allclose(adc_hat.numpy(), adc_true, rtol=1e-4)


def test_adc_fit_handles_single_vector():
    S = torch.exp(-ACRIN_B_SCHEME * 1.0e-3)
    _, adc = adc_monoexp_fit(S, ACRIN_B_SCHEME)
    assert adc.shape == (1,) and np.isclose(adc.item(), 1.0e-3, rtol=1e-4)


# ----------------------------- ingestion pure helpers --------------------------
def test_reconcile_b0_averages():
    b0 = np.array([1.0, 1.2, 0.8])
    assert np.isclose(reconcile_b0(b0), 1.0)


def test_roi_mean_per_b():
    vol = np.zeros((2, 3, 3))
    vol[0] = 5.0
    vol[1, 1, 1] = 9.0
    mask = np.zeros((3, 3), dtype=bool)
    mask[1, 1] = True
    out = roi_mean_per_b(vol, mask)
    assert np.allclose(out, [5.0, 9.0])


def test_assemble_4b_reconciles_multiple_b0():
    per_b = {0: [1.0, 1.2, 0.8], 100: 0.9, 600: 0.5, 800: 0.4}
    vec = assemble_4b_vector(per_b)
    assert np.allclose(vec, [1.0, 0.9, 0.5, 0.4])


# ----------------------------- imputation --------------------------------------
def test_imputation_recovers_dense_and_splices_real():
    theta = sample_prior(500, seed=7)
    f, D, Ds = theta[:, 0], theta[:, 1], theta[:, 2]
    dense_true = ivim_biexponential(DEFAULT_B_SCHEME, f, D, Ds).numpy()
    sig4 = ivim_biexponential(ACRIN_B_SCHEME, f, D, Ds).numpy()
    sig4n = normalize_by_b0(sig4)

    dense_imp = impute_dense(sig4n)
    # Real ACRIN positions are spliced exactly.
    assert np.allclose(dense_imp[:, list(ACRIN_TO_DENSE_IDX)], sig4n, atol=1e-9)
    # Imputed positions are close to the true dense signal on clean data.
    imp_idx = [j for j in range(10) if j not in ACRIN_TO_DENSE_IDX]
    rel = np.abs(dense_imp[:, imp_idx] - dense_true[:, imp_idx]) / np.clip(dense_true[:, imp_idx], 1e-6, None)
    assert np.median(rel) < 0.02


def test_segmented_fit_recovers_f_and_D_on_clean():
    theta = sample_prior(500, seed=3)
    f, D, Ds = theta[:, 0], theta[:, 1], theta[:, 2]
    sig4 = ivim_biexponential(ACRIN_B_SCHEME, f, D, Ds).numpy()
    fh, Dh, Dsh = segmented_ivim_fit(normalize_by_b0(sig4))
    assert np.median(np.abs(fh - f.numpy())) < 0.02
    assert np.median(np.abs(Dh - D.numpy())) < 1e-4
    assert np.all(np.isfinite(Dsh))  # D* is allowed to be poor but must stay finite


# ----------------------------- Arm-1 detector ----------------------------------
def test_signal_detector_separates_id_from_shifted():
    _, _, feat_ref = generate_acrin_id_reference(4000, snr=40.0, seed=1)
    _, _, feat_cal = generate_acrin_id_reference(1500, snr=40.0, seed=2)
    det = SignalSpaceDetector(method="knn", k=15).fit(feat_ref, calib_features=feat_cal)

    _, _, feat_id = generate_acrin_id_reference(500, snr=40.0, seed=3)        # in-distribution
    _, _, feat_ood = generate_acrin_id_reference(500, snr=8.0, seed=4)        # low-SNR shift
    s_id, s_ood = det.score(feat_id), det.score(feat_ood)
    assert s_ood.mean() > s_id.mean()
    # conformal p-values for ID are roughly uniform -> mean near 0.5.
    p_id = det.conformal_pvalue(feat_id)
    assert 0.35 < p_id.mean() < 0.65


# ----------------------------- trust + coupling --------------------------------
def test_adc_repeatability_zero_when_identical():
    adc = np.array([1.0e-3, 1.5e-3, 0.8e-3])
    rep = adc_repeatability(adc, adc.copy())
    assert rep["wCV"] == 0.0 and rep["RC"] == 0.0
    assert np.allclose(rep["abs_diff"], 0.0)


def test_headline_coupling_sign_and_controls():
    """Estimator-free end-to-end: OOD score positively tracks ADC unrepeatability."""
    det = SignalSpaceDetector(method="knn", k=15)
    _, _, fr = generate_acrin_id_reference(6000, snr=40.0, seed=1)
    _, _, fc = generate_acrin_id_reference(2000, snr=40.0, seed=2)
    det.fit(fr, calib_features=fc)

    units = concat_units(synthetic_id_units(300, snr=40.0, seed=10),
                         synthetic_invivo_units(300, seed=11))
    feats = signal_features(normalize_by_b0(units.sig_test))
    scores = det.score(feats)

    _, adc_t = adc_monoexp_fit(torch.tensor(units.sig_test), ACRIN_B_SCHEME)
    _, adc_r = adc_monoexp_fit(torch.tensor(units.sig_retest), ACRIN_B_SCHEME)
    rep = adc_repeatability(adc_t.numpy(), adc_r.numpy())

    coup = spearman_coupling(scores, rep["cov"])
    assert coup["rho"] > 0.3 and coup["p"] < 1e-6

    mean_adc = 0.5 * (adc_t.numpy() + adc_r.numpy())
    pc = partial_spearman(scores, rep["cov"], mean_adc)
    assert pc["rho_partial"] > 0.25  # survives controlling for ADC level

    # Degeneracy control: low-f ID units not over-flagged.
    idm = units.is_synth_id
    deg = degeneracy_fpr(scores[idm], units.theta[idm][:, 0] < np.percentile(units.theta[idm][:, 0], 20))
    assert deg["degenerate_fpr"] < 0.15


# ----------------------------- UnitTable round-trip ----------------------------
def test_unit_table_save_load(tmp_path):
    u = synthetic_invivo_units(20, seed=0)
    p = tmp_path / "u.npz"
    u.save(p)
    v = UnitTable.load(p)
    assert np.allclose(u.sig_test, v.sig_test)
    assert np.array_equal(u.qa_pass, v.qa_pass)


# ----------------------------- NIfTI ingestion round-trip ----------------------
def test_nifti_ingestion_roundtrip(tmp_path):
    sitk = pytest.importorskip("SimpleITK")
    from sibyl.data.acrin_ingest import roi_4b_from_nifti

    # Build a synthetic 4D DWI with TWO b0 volumes (tests reconciliation) + 100/600/800.
    Z, Y, X = 4, 5, 6
    bvals = np.array([0, 0, 100, 600, 800])
    true_roi = np.array([1.0, 0.9, 0.55, 0.30])  # [b0(after avg), b100, b600, b800]
    vols = []
    # two b0 volumes that average to 1.0 inside the ROI
    vols.append(np.full((Z, Y, X), 0.8))
    vols.append(np.full((Z, Y, X), 1.2))
    for s in true_roi[1:]:
        vols.append(np.full((Z, Y, X), s))
    dwi = np.stack(vols, axis=0).astype(np.float32)  # [5, Z, Y, X]
    mask = np.zeros((Z, Y, X), dtype=np.uint8)
    mask[1:3, 1:3, 1:3] = 1

    dwi_p = tmp_path / "dwi.nii.gz"
    bval_p = tmp_path / "dwi.bval"
    mask_p = tmp_path / "mask.nii.gz"
    sitk.WriteImage(sitk.GetImageFromArray(dwi), str(dwi_p))
    sitk.WriteImage(sitk.GetImageFromArray(mask), str(mask_p))
    np.savetxt(bval_p, bvals[None], fmt="%d")

    vec = roi_4b_from_nifti(dwi_p, bval_p, mask_p)
    assert np.allclose(vec, true_roi, atol=1e-5)


# ----------------------------- uGUIDE Arm-2 integration (small) ----------------
@pytest.mark.slow
def test_arm2_uguide_integration(tmp_path):
    """Train a tiny dense estimator and confirm Arm-2 produces finite family scores."""
    from sibyl.experiments.tier2 import (
        train_dense_estimator, fit_dense_family_detectors, arm2_scores,
    )
    cfg = train_dense_estimator(tmp_path / "dense", epochs=3, n_train=1500, seed=42)
    det1, det2 = fit_dense_family_detectors(cfg, n_calib=800, seed=43)
    u = synthetic_invivo_units(16, seed=5)
    out = arm2_scores(cfg, det1, det2, u.sig_test)
    assert out["family1"].shape == (16,) and out["family2"].shape == (16,)
    assert np.all(np.isfinite(out["family1"])) and np.all(np.isfinite(out["family2"]))
