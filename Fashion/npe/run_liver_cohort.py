"""
run_liver_cohort.py
===================
In-vivo COHORT validation (Checkpoint IV) of the amortized NPE's held-out-b
calibration, on an open multi-subject liver IVIM dataset.

Dataset: "The role of IVIM parameters in characterization of hypovascular liver
lesions" (Dryad/Zenodo record 4408313, CC0; 77 patients; GE Signa 1.5T; b-values
{0, 10, 25, 50, 200, 500, 800} s/mm^2). DICOM -> NIfTI via dcm2niix; each b-value
is its own 2-volume series [b_X, b_0], so every diffusion image carries its own
b=0 reference.

Why this dataset is the right cohort test
-----------------------------------------
The trained NPE is amortized over the clinical-sparse scheme
{0,10,20,30,50,100,200,400,600,800,1000}. The liver cohort's b-scheme differs, so
this is an OFF-SCHEME, multi-subject test of the main-text claim that amortized
calibration does not transfer across acquisitions (Figure 4C). NLLS, which refits
per voxel, is the calibrated reference. We therefore expect: NPE under-covers
held-out b across the cohort; NLLS stays near nominal.

Method (mirrors npe/run_f_realdata.py exactly; helpers imported from it)
-----------------------------------------------------------------------
Per patient:
  * Assemble a per-voxel normalized attenuation curve over b = {0,10,25,50,200,
    500,800}: each b_X series is normalized by its OWN in-series b=0 volume
    (robust to inter-series scaling); b=0 -> 1.0 by construction.
  * Foreground/liver-parenchyma voxels: b0 above a robust threshold (exclude air)
    and monotone-plausible, finite, positive signal across all b. A fixed-seed
    random subsample (default 300) is taken for tractable per-voxel NLLS.
  * Per-patient noise sigma from a background (air) ROI of the highest-b trace
    image; per-voxel SNR = S0 / sigma, clipped to [8,100] (as in run_g). The
    images are registered traces ("RTr"), so SNR_avg = SNR (no sqrt(N_dir) factor,
    unlike the 6-direction brain).
  * Held-out-b posterior-predictive coverage: fit on a b-subset, propagate
    posterior samples through the forward model, score the fraction of held-out
    signals covered. fit = {0,25,200,800}; held-out = {10,50,500}.

Outputs:
  npe/liver_cohort_coverage.csv   per-patient coverage (NPE + NLLS) at each level
  npe/liver_cohort_summary.json   cohort aggregate stats
"""
from __future__ import annotations

import os
import re
import sys
import csv
import json
import glob
import argparse
import numpy as np
import torch
import nibabel as nib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Identical methodology: import the exact helpers used by the single-subject check.
from run_f_realdata import (
    fit_biexp_nlls, biexp_signal, add_rician_noise, sample_nlls_asymptotic,
)
from npe_prior import get_processed_prior, invert_theta
from train_npe import pack_x, SNRWrapperEmbedding

ALL_BVALS = np.array([0, 10, 25, 50, 200, 500, 800], float)
BVALS_FIT = np.array([0, 25, 200, 800], float)
BVALS_HELDOUT = np.array([10, 50, 500], float)
NOMINAL_LEVELS = [0.50, 0.68, 0.80, 0.90, 0.95]


def _b_from_name(fname):
    m = re.search(r"B=(\d+)", os.path.basename(fname))
    return float(m.group(1)) if m else None


def load_patient_curve(pat_nifti_dir):
    """Return (norm_signal[H,W,Z,7], s0[H,W,Z]) over ALL_BVALS, or (None,None).

    Each b-value is a 2-volume series [b_X, b_0]. The b-value is taken from the
    series filename (the per-series .bval is occasionally degenerate, e.g. [0,0]).
    Within a series the b=0 reference is the higher-mean volume. All required
    b-series must share the same spatial shape (else cross-series voxel
    correspondence is undefined and the patient is skipped).
    """
    series = {}
    for f in sorted(glob.glob(os.path.join(pat_nifti_dir, "*.nii.gz"))):
        if "ADC" in os.path.basename(f):
            continue  # scanner-derived ADC maps are not raw DWI
        b = _b_from_name(f)
        if b is None or b <= 0:
            continue
        img = nib.load(f).get_fdata()
        if img.ndim != 4 or img.shape[3] != 2:
            continue
        m0, m1 = float(np.nanmean(img[..., 0])), float(np.nanmean(img[..., 1]))
        b0_idx = 0 if m0 >= m1 else 1     # b=0 has the higher mean signal
        bx_idx = 1 - b0_idx
        series[b] = (img[..., bx_idx], img[..., b0_idx])
    needed = [b for b in ALL_BVALS if b > 0]
    if any(b not in series for b in needed):
        return None, None
    shapes = {series[b][0].shape for b in needed}
    if len(shapes) != 1:
        return None, None                  # mismatched geometry across b-series
    shape = series[needed[0]][0].shape
    s0 = np.mean([series[b][1] for b in needed], axis=0)
    norm = np.zeros(shape + (len(ALL_BVALS),), float)
    for k, b in enumerate(ALL_BVALS):
        if b == 0:
            norm[..., k] = 1.0
        else:
            bx_vol, b0_vol = series[b]
            with np.errstate(divide="ignore", invalid="ignore"):
                norm[..., k] = np.where(b0_vol > 0, bx_vol / b0_vol, np.nan)
    return norm, s0


def select_voxels(norm, s0, n_voxels, rng):
    """Foreground liver voxels with valid attenuation; subsample; return signals+S0."""
    s0_thresh = 0.15 * np.nanpercentile(s0, 99)
    finite = np.isfinite(norm).all(axis=-1)
    # plausible attenuation: all normalized values in (0, 1.5]
    plausible = finite & (norm[..., 1:] > 0).all(axis=-1) & (norm[..., 1:] <= 1.5).all(axis=-1)
    fg = (s0 > s0_thresh) & plausible
    coords = np.argwhere(fg)
    if len(coords) == 0:
        return None, None
    if len(coords) > n_voxels:
        coords = coords[rng.choice(len(coords), size=n_voxels, replace=False)]
    sig = np.array([norm[c[0], c[1], c[2], :] for c in coords])     # (n,7)
    s0v = np.array([s0[c[0], c[1], c[2]] for c in coords])
    return sig, s0v


def estimate_snr(sig, s0v):
    """Per-patient SNR (on the normalized signal) from the median NLLS-fit residual.

    Backgrounds are zero-filled, so air-ROI noise is unrecoverable. Instead we fit
    the biexponential to each voxel's full 7-b normalized curve and take the MEDIAN
    per-voxel residual RMS (dof = N_b - 3) as the normalized-signal noise level
    sigma. The median is robust to the minority of poorly-fitting voxels (vessels,
    lesions, inter-series motion); SNR = 1/sigma, clipped to [8,100], applied to all
    voxels (one acquisition -> one noise level on the [0,1]-normalized curve).
    """
    n = sig.shape[0]
    idx = np.arange(n) if n <= 120 else np.linspace(0, n - 1, 120).astype(int)
    dofc = len(ALL_BVALS) - 3
    rms = []
    for i in idx:
        theta = fit_biexp_nlls(ALL_BVALS, sig[i])
        resid = biexp_signal(ALL_BVALS, *theta) - sig[i]
        rms.append(np.sqrt(np.sum(resid ** 2) / dofc))
    sigma = float(np.median(rms))
    if not np.isfinite(sigma) or sigma <= 0:
        return np.full(n, 30.0), None
    snr_scalar = float(np.clip(1.0 / sigma, 8.0, 100.0))
    return np.full(n, snr_scalar), sigma


def coverage_for_patient(sig, snr, posterior, set_size_setter, prior_bounds,
                         log_dstar, n_samples, rng):
    """Held-out-b posterior-predictive coverage (NPE + NLLS) for one patient."""
    idx_fit = [int(np.where(ALL_BVALS == b)[0][0]) for b in BVALS_FIT]
    idx_held = [int(np.where(ALL_BVALS == b)[0][0]) for b in BVALS_HELDOUT]
    fit_signals = sig[:, idx_fit]            # (n, 4)
    heldout_signals = sig[:, idx_held]       # (n, 3)
    n = fit_signals.shape[0]

    # NLLS on fit subset
    nlls_fits = np.array([fit_biexp_nlls(BVALS_FIT, fit_signals[i]) for i in range(n)])
    nlls_samples = np.zeros((n_samples, n, 3))
    for i in range(n):
        nlls_samples[:, i, :] = sample_nlls_asymptotic(
            nlls_fits[i], BVALS_FIT, snr[i], num_samples=n_samples, bounds=prior_bounds, rng=rng)

    # NPE posterior on the fit set (off-scheme b-values via the set embedding)
    bvals_fit_tiled = np.tile(BVALS_FIT[None, :], (n, 1))
    obs = torch.as_tensor(np.stack([bvals_fit_tiled, fit_signals], axis=-1), dtype=torch.float32)
    snr_ctx = torch.as_tensor(np.log10(snr)[:, None], dtype=torch.float32)
    x = pack_x(obs, snr_ctx, "set")
    set_size_setter(len(BVALS_FIT))
    with torch.no_grad():
        s_npe = posterior.sample_batched((n_samples,), x=x, reject_outside_prior=False)
    s_npe = invert_theta(s_npe, log_dstar=log_dstar).cpu().numpy()
    s_npe = np.clip(s_npe, prior_bounds[0], prior_bounds[1])

    # Predict held-out signals; RTr trace image -> snr_avg = snr (no sqrt(N_dir)).
    # Vectorized over (samples, voxels, held-out b); identical math to the per-voxel
    # biexp_signal + add_rician_noise used in run_f_realdata.
    bh = BVALS_HELDOUT[None, None, :]                  # (1,1,3)
    sigma = (1.0 / snr)[None, :, None]                 # (1,n,1), S0=1

    def predict(theta):
        D, Ds, f = theta[..., 0:1], theta[..., 1:2], theta[..., 2:3]   # (ns,n,1)
        clean = f * np.exp(-bh * Ds) + (1.0 - f) * np.exp(-bh * D)     # (ns,n,3)
        n_re = rng.normal(0.0, 1.0, size=clean.shape) * sigma
        n_im = rng.normal(0.0, 1.0, size=clean.shape) * sigma
        return np.sqrt((clean + n_re) ** 2 + n_im ** 2)

    pred_npe = predict(s_npe)
    pred_nlls = predict(nlls_samples)

    rows = {}
    for lvl in NOMINAL_LEVELS:
        a = 1.0 - lvl
        lo, hi = a / 2.0 * 100, (1 - a / 2.0) * 100
        cov_npe = ((heldout_signals >= np.percentile(pred_npe, lo, axis=0)) &
                   (heldout_signals <= np.percentile(pred_npe, hi, axis=0))).mean()
        cov_nlls = ((heldout_signals >= np.percentile(pred_nlls, lo, axis=0)) &
                    (heldout_signals <= np.percentile(pred_nlls, hi, axis=0))).mean()
        rows[lvl] = (float(cov_npe), float(cov_nlls))
    return rows


def main():
    ap = argparse.ArgumentParser(description="Liver-cohort in-vivo held-out-b coverage (NPE vs NLLS).")
    ap.add_argument("--model", default="npe/npe_posterior_setB.pt")
    ap.add_argument("--nifti-root", default="download/liver_cohort/nifti")
    ap.add_argument("--n-voxels", type=int, default=300)
    ap.add_argument("--n-samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-csv", default="npe/liver_cohort_coverage.csv")
    ap.add_argument("--out-json", default="npe/liver_cohort_summary.json")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    np.random.seed(args.seed); torch.manual_seed(args.seed)

    model_path = args.model
    if not os.path.exists(model_path):
        alt = "npe/" + os.path.basename(model_path)
        model_path = alt if os.path.exists(alt) else model_path
    sys.modules['__main__'].SNRWrapperEmbedding = SNRWrapperEmbedding
    print(f"Loading NPE from {model_path}...")
    posterior = torch.load(model_path, map_location="cpu", weights_only=False)
    log_dstar = bool(posterior.prior.support.base_constraint.lower_bound[1] < 0)
    get_processed_prior(device="cpu", log_dstar=log_dstar)
    prior_bounds = ([0.2e-3, 3.0e-3, 0.0], [3.0e-3, 0.15, 0.5])

    estimator = posterior.posterior_estimator
    std_mod = estimator.net._embedding_net[0]
    orig_mean, orig_std = std_mod._mean.clone(), std_mod._std.clone()

    def set_size_setter(K):
        estimator._condition_shape = torch.Size([K, 3])
        std_mod._mean = orig_mean[:K, :]
        std_mod._std = orig_std[:K, :]

    pat_dirs = sorted(glob.glob(os.path.join(args.nifti_root, "pat*")),
                      key=lambda p: int("".join(ch for ch in os.path.basename(p) if ch.isdigit()) or 0))
    print(f"Found {len(pat_dirs)} patient NIfTI dirs.")

    per_patient = []
    for pdir in pat_dirs:
        pid = os.path.basename(pdir)
        norm, s0 = load_patient_curve(pdir)
        if norm is None:
            print(f"  {pid}: SKIP (missing b-values)"); continue
        sig, s0v = select_voxels(norm, s0, args.n_voxels, rng)
        if sig is None or len(sig) < 20:
            print(f"  {pid}: SKIP (too few foreground voxels)"); continue
        snr, sigma_norm = estimate_snr(sig, s0v)
        sigma_str = f"{sigma_norm:.3f}" if sigma_norm else "FALLBACK"
        try:
            rows = coverage_for_patient(sig, snr, posterior, set_size_setter, prior_bounds,
                                        log_dstar, args.n_samples, rng)
        except Exception as e:
            print(f"  {pid}: SKIP (NPE sampling failed: {type(e).__name__})"); continue
        rec = {"patient": pid, "n_voxels": int(len(sig)), "snr_median": float(np.median(snr))}
        for lvl in NOMINAL_LEVELS:
            rec[f"npe_cov_{lvl}"] = rows[lvl][0]
            rec[f"nlls_cov_{lvl}"] = rows[lvl][1]
        per_patient.append(rec)
        print(f"  {pid}: n={rec['n_voxels']} sigma={sigma_str} snr~{rec['snr_median']:.0f} "
              f"NPE@.95={rec['npe_cov_0.95']:.3f} NLLS@.95={rec['nlls_cov_0.95']:.3f}")

    if not per_patient:
        print("No patients processed."); return

    # Write per-patient CSV
    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    fields = (["patient", "n_voxels", "snr_median"] +
              [f"npe_cov_{l}" for l in NOMINAL_LEVELS] + [f"nlls_cov_{l}" for l in NOMINAL_LEVELS])
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(per_patient)

    # Cohort aggregate
    npe95 = np.array([r["npe_cov_0.95"] for r in per_patient])
    nlls95 = np.array([r["nlls_cov_0.95"] for r in per_patient])
    summary = {
        "dataset": "Dryad/Zenodo 4408313 hypovascular liver lesions (GE 1.5T)",
        "n_patients": len(per_patient),
        "bvals_fit": BVALS_FIT.tolist(), "bvals_heldout": BVALS_HELDOUT.tolist(),
        "npe_cov_0.95_mean": float(npe95.mean()), "npe_cov_0.95_median": float(np.median(npe95)),
        "npe_cov_0.95_iqr": [float(np.percentile(npe95, 25)), float(np.percentile(npe95, 75))],
        "nlls_cov_0.95_mean": float(nlls95.mean()), "nlls_cov_0.95_median": float(np.median(nlls95)),
        "nlls_cov_0.95_iqr": [float(np.percentile(nlls95, 25)), float(np.percentile(nlls95, 75))],
        "frac_patients_npe_undercovers_nlls": float((npe95 < nlls95).mean()),
        "per_level_mean": {
            str(l): {"npe": float(np.mean([r[f"npe_cov_{l}"] for r in per_patient])),
                     "nlls": float(np.mean([r[f"nlls_cov_{l}"] for r in per_patient]))}
            for l in NOMINAL_LEVELS},
    }
    with open(args.out_json, "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 70)
    print(f"LIVER COHORT (N={summary['n_patients']} patients) held-out-b coverage @ nominal 0.95")
    print(f"  NPE  mean {summary['npe_cov_0.95_mean']:.3f}  median {summary['npe_cov_0.95_median']:.3f}  IQR {summary['npe_cov_0.95_iqr']}")
    print(f"  NLLS mean {summary['nlls_cov_0.95_mean']:.3f}  median {summary['nlls_cov_0.95_median']:.3f}  IQR {summary['nlls_cov_0.95_iqr']}")
    print(f"  Patients where NPE under-covers NLLS: {summary['frac_patients_npe_undercovers_nlls']*100:.0f}%")
    print("=" * 70)
    print(f"Wrote {args.out_csv} and {args.out_json}")


if __name__ == "__main__":
    main()
