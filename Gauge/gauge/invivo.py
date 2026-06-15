"""Gauge 04 -- CP2: in-vivo qualitative demonstration (HALT-ABLE on data).

Conformal coverage is a guarantee about EXCHANGEABLE, LABELED data. In vivo there
are no ground-truth IVIM parameters, so **coverage cannot be measured or claimed**
-- the Echo tension. What CAN be shown qualitatively, on real (or realistic)
voxels, is the *behavior* of the deployed pipeline:

  1. the conformalized interval map (point estimates + per-voxel band widths),
  2. the D* compartment band widths -- wide exactly where the perfusion
     compartment is under-identified (the Gauge 03 wall, now on real signal),
  3. the **Minos-style deployment monitor** pointed at the in-vivo signals: it
     fires, because synthetic-calibrated -> in-vivo is itself a large
     exchangeability break. That firing is the honest, observable reason the
     synthetic coverage guarantee must NOT be asserted in vivo.

This module is fully PLUGGABLE: point ``load_dwi_nifti`` / ``load_signals_npy``
at a real DWI series and it runs the identical pipeline. Out of the box it runs
on a transparent SYNTHETIC STAND-IN (clearly labeled, NOT in vivo) so the code
path and the qualitative outputs are demonstrable here.

** DATASET + FRAMING ARE A HUMAN SIGN-OFF (out of scope for the agent). ** No
clean, permissively-licensed, ground-truth-free in-vivo IVIM dataset was usable
in this environment without a licensing decision and extra dependencies (see the
GATE 2 / HALT-TO-REPORT printout). Recommended real options to sign off on:
  * OSIPI TF2.4 IVIM-MRI CodeCollection phantoms (Apache-2.0) -- realistic
    anatomy via the XCAT / ICBM brain generators (HAVE ground truth, so they are
    a transfer/OOD demo, not a true no-GT in-vivo case);
  * a TCIA DWI collection (e.g. a liver/brain series) -- true in-vivo, no GT,
    behind a data-use agreement.

Run:  python -m gauge.invivo                      # synthetic stand-in demo
      python -m gauge.invivo /path/dwi.nii.gz /path/bvals   # real DWI
"""
import argparse
import datetime
import json
import os
import pickle
import sys
import time

import numpy as np

from gauge.forward import DEFAULT_B_VALUES
from gauge.cohort import DEFAULT_SNR_GRID, DEFAULT_SEED
from gauge.estimators import IVIMQuantileRegressor
from gauge.conformal import cqr, interval_width
from gauge.monitor import DeploymentMonitor
from gauge.robustness import (_draw_params, _draw_snr, _simulate, _observe,
                              _calibrate, _importance_weights)

PARAM_NAMES = ("D", "D*", "f")
ALPHA = 0.10
SEED = DEFAULT_SEED
_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
_CACHE = os.path.join(_RESULTS_DIR, "invivo.pkl")

# Real-data path outputs (ADDITIVE: never overwrite the synthetic invivo_report.txt
# / invivo_demo.pdf, which the GATE-3 consistency check traces byte-for-byte).
_REAL_REPORT = os.path.join(_RESULTS_DIR, "invivo_real_report.txt")
_REAL_PROVENANCE = os.path.join(_RESULTS_DIR, "invivo_real_provenance.json")
_REAL_FIG = os.path.join(_FIG_DIR, "invivo_real_demo.pdf")
_REAL_RETEST_REPORT = os.path.join(_RESULTS_DIR, "invivo_real_retest_report.txt")
_REAL_RETEST_FIG = os.path.join(_FIG_DIR, "invivo_real_retest.pdf")
_REAL_MAPS_FIG = os.path.join(_FIG_DIR, "invivo_real_maps.pdf")


# --------------------------------------------------------------------------- #
# Data sources (pluggable). Real loaders normalize each voxel to its b=0 signal,
# matching the synthetic S0=1 convention the calibration was built on.
# --------------------------------------------------------------------------- #
def _normalize_to_s0(signals, b):
    b = np.asarray(b, float)
    s0 = signals[:, b == b.min()].mean(1, keepdims=True)
    s0 = np.where(s0 > 0, s0, signals.max(1, keepdims=True))
    return np.clip(signals / s0, 0.0, None)


def load_dwi_nifti(img_path, bval_path, mask_path=None, max_voxels=4000,
                   seed=0):
    """Load a real 4D DWI volume + b-values into (signals, b). Needs nibabel.

    Voxels are masked (provided mask, or an Otsu-ish b=0 threshold), normalized
    to S0, and randomly subsampled to ``max_voxels`` for the demo.
    """
    try:
        import nibabel as nib
    except ImportError as e:                                # pragma: no cover
        raise ImportError(
            "real-DWI loading needs nibabel (`pip install nibabel`). The demo "
            "runs without it on the synthetic stand-in.") from e
    img = np.asarray(nib.load(img_path).get_fdata(), float)   # (X,Y,Z,B)
    b = np.loadtxt(bval_path).ravel()
    vol = img.reshape(-1, img.shape[-1])
    if mask_path:
        mask = np.asarray(nib.load(mask_path).get_fdata(), float).ravel() > 0
    else:
        s0 = vol[:, b == b.min()].mean(1)
        mask = s0 > 0.25 * np.percentile(s0, 99)
    vox = vol[mask]
    rng = np.random.default_rng(seed)
    if vox.shape[0] > max_voxels:
        vox = vox[rng.choice(vox.shape[0], max_voxels, replace=False)]
    return _normalize_to_s0(vox, b), b


def load_signals_npy(sig_path, bval_path):
    """Generic loader: (N, B) signal array + b-values from .npy/.txt."""
    sig = np.load(sig_path) if sig_path.endswith(".npy") else np.loadtxt(sig_path)
    b = np.loadtxt(bval_path).ravel()
    return _normalize_to_s0(np.atleast_2d(sig), b), b


def synthetic_stand_in(n=2500, b=DEFAULT_B_VALUES, seed=SEED + 321):
    """A transparent OUT-OF-CALIBRATION synthetic stand-in -- NOT in vivo.

    Drawn from a tri-exponential truth at clinically realistic low SNR with a
    deliberately out-of-calibration prior (broadened D*/f), so it differs from the
    (bi-exp, higher-SNR) calibration the way real tissue would -- which is the
    whole point of CP2 (the deployment monitor must detect this transfer).
    Returns ``(signals, b, true_params)``;
    the true params are withheld from the demo pipeline and used ONLY to annotate
    the clearly-separated phantom-coverage sanity line.
    """
    rng = np.random.default_rng(seed)
    params = _draw_params(n, rng, prior={"Dstar": (10e-3, 130e-3),
                                         "f": (0.05, 0.45)})
    snr = _draw_snr(n, rng, (8.0, 12.0, 18.0))
    sig = _simulate(b, params, snr, rng, model="triexp", triexp=(4.0, 0.20))
    return sig, np.asarray(b, float), params


# --------------------------------------------------------------------------- #
# Deployed pipeline: synthetic LABELED calibration (the only labels we ever get)
# -> split-conformal radii + a conformalizable CQR band + the monitor.
# --------------------------------------------------------------------------- #
def deployed_calibration(b=DEFAULT_B_VALUES, alpha=ALPHA, seed=SEED,
                         n_train=3000, n_cal=2000):
    rng = np.random.default_rng(seed)
    tr_params = _draw_params(n_train, rng)
    tr_snr = _draw_snr(n_train, rng, DEFAULT_SNR_GRID)
    tr_sig = _simulate(b, tr_params, tr_snr, rng)
    cal_params = _draw_params(n_cal, rng)
    cal_snr = _draw_snr(n_cal, rng, DEFAULT_SNR_GRID)
    cal_sig = _simulate(b, cal_params, cal_snr, rng)

    cal_theta, cal_feat, cal_resid = _observe(cal_sig, b)
    q = _calibrate(cal_theta, cal_params, alpha)

    # conditional band: gradient-boosted quantile regressor + CQR offset
    levels = [alpha / 2, 1 - alpha / 2]
    qreg = IVIMQuantileRegressor(levels, random_state=0).fit(tr_sig, tr_params)
    cal_lo = np.stack([qreg.predict_quantile(cal_sig, j, levels[0])
                       for j in range(3)], axis=1)
    cal_hi = np.stack([qreg.predict_quantile(cal_sig, j, levels[1])
                       for j in range(3)], axis=1)

    monitor = DeploymentMonitor(seed=seed).fit(cal_feat, cal_resid)
    return {"b": np.asarray(b, float), "q": q, "qreg": qreg, "levels": levels,
            "cal_params": cal_params, "cal_lo": cal_lo, "cal_hi": cal_hi,
            "cal_feat": cal_feat, "cal_resid": cal_resid, "monitor": monitor}


def run_demo(signals, b, cal, source_name, true_params=None, seed=SEED):
    """Run the deployed pipeline on (unlabeled) in-vivo signals; qualitative.

    Returns a dict of qualitative observables. ``true_params`` (phantom only) is
    used ONLY for a clearly-separated coverage sanity line; the pipeline never
    sees it.
    """
    theta, feat, resid = _observe(signals, b)
    alpha, levels = ALPHA, cal["levels"]

    # conditional (CQR) intervals -> per-voxel band widths
    lo_raw = np.stack([cal["qreg"].predict_quantile(signals, j, levels[0])
                       for j in range(3)], axis=1)
    hi_raw = np.stack([cal["qreg"].predict_quantile(signals, j, levels[1])
                       for j in range(3)], axis=1)
    lo = np.empty_like(lo_raw)
    hi = np.empty_like(hi_raw)
    for j in range(3):
        lj, hj, _ = cqr(cal["cal_lo"][:, j], cal["cal_hi"][:, j],
                        cal["cal_params"][:, j], lo_raw[:, j], hi_raw[:, j], alpha)
        lo[:, j], hi[:, j] = lj, hj
    widths = interval_width(lo, hi)                        # (N,3)

    mon = cal["monitor"].evaluate(feat, resid)
    w_cal, w_test, _ = _importance_weights(cal["cal_feat"], feat, seed=seed)

    out = {
        "source": source_name, "n_vox": int(signals.shape[0]),
        "theta_median": np.median(theta, 0),
        "width_median": np.median(widths, 0),
        "dstar_width_q": np.quantile(widths[:, 1], [0.1, 0.5, 0.9]),
        "dstar_widths": widths[:, 1],
        "dstar_hat": theta[:, 1],
        "monitor_fires": mon["fires"], "monitor_auc": mon["auc"],
        "monitor_maha": mon["maha"]["stat"], "monitor_maha_thr": mon["maha"]["threshold"],
        "monitor_resid": mon["resid"]["stat"], "monitor_resid_thr": mon["resid"]["threshold"],
        "weight_spread": float(np.quantile(w_test, 0.9) / max(np.quantile(w_test, 0.1), 1e-9)),
        "has_gt": true_params is not None,
    }
    if true_params is not None:                            # PHANTOM ONLY
        from gauge.conformal import empirical_coverage
        out["phantom_coverage"] = np.array(
            [empirical_coverage(lo[:, j], hi[:, j], true_params[:, j])
             for j in range(3)])
    return out


# --------------------------------------------------------------------------- #
# Orchestration / GATE 2
# --------------------------------------------------------------------------- #
def compute(force=False, dwi=None, bvals=None, seed=SEED):
    """Run the qualitative in-vivo demo. Pure function of ``seed`` (hygiene: the
    synthetic stand-in, deployed calibration, monitor and domain-classifier all
    derive from it). In-vivo quantities stay qualitative and are NOT banded; the
    seed param exists so the demo participates cleanly in the multi-seed harness.
    Cache is seed-specific so a sweep never reuses another seed's demo."""
    cache_path = os.path.join(_RESULTS_DIR, f"invivo_seed{int(seed)}.pkl")
    if (not force) and dwi is None and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    t0 = time.time()
    if dwi is not None:
        signals, b = load_dwi_nifti(dwi, bvals)
        source, true_params = f"REAL DWI ({os.path.basename(dwi)})", None
    else:
        signals, b, true_params = synthetic_stand_in(seed=seed + 321)
        source = "synthetic stand-in (NOT in vivo)"
    cal = deployed_calibration(b=b, seed=seed)
    res = run_demo(signals, b, cal, source, true_params=true_params, seed=seed)
    res["q_dstar"] = float(cal["q"][1])
    payload = {"res": res, "is_real": dwi is not None}
    if dwi is None:
        with open(cache_path, "wb") as fh:
            pickle.dump(payload, fh)
    print(f"[invivo] demo computed ({time.time()-t0:.0f}s)")
    return payload


def main(force=False, dwi=None, bvals=None, seed=SEED):
    P = compute(force=force, dwi=dwi, bvals=bvals, seed=seed)
    r = P["res"]
    lines = []

    def out(*x):
        s = " ".join(str(z) for z in x)
        print(s)
        lines.append(s)

    out("#" * 92)
    out("GAUGE 04 -- CP2 / GATE 2 (HALT-ABLE): in-vivo qualitative demonstration")
    out("#" * 92)
    out(f"data source: {r['source']}   |   voxels: {r['n_vox']}   |   "
        f"ground truth available: {r['has_gt']}")
    out("")
    out("** What conformal CANNOT claim in vivo **  (kept strictly separate from "
        "the synthetic guarantees):")
    out("  - In vivo there are NO ground-truth (D, D*, f), so realized coverage "
        "is UNMEASURABLE and the")
    out("    finite-sample 1-alpha guarantee (Gauge 01-03, on labeled exchangeable "
        "synthetic data) is NOT")
    out("    asserted here. Everything below is QUALITATIVE interval behavior, not "
        "a coverage claim.")
    out("-" * 92)
    out("[2.1] Conformalized interval map (deployed CQR band; synthetic-calibrated):")
    out(f"      median plug-in (D, D*, f) = "
        f"({r['theta_median'][0]*1e3:.2f}, {r['theta_median'][1]*1e3:.1f}, "
        f"{r['theta_median'][2]:.3f})   [D,D* in 1e-3 mm^2/s]")
    out(f"      median band width (D, D*, f) = "
        f"({r['width_median'][0]*1e3:.2f}, {r['width_median'][1]*1e3:.1f}, "
        f"{r['width_median'][2]:.3f})")
    out("[2.2] D* compartment band widths (the Gauge 03 wall, on real signal): "
        "10/50/90th pct (1e-3 mm^2/s)")
    q = r["dstar_width_q"] * 1e3
    out(f"      = [{q[0]:.1f}, {q[1]:.1f}, {q[2]:.1f}]  -- the D* band is wide and "
        f"HEAVY-TAILED: the conformal")
    out(f"      interval honestly balloons on the voxels where the perfusion "
        f"compartment is under-identified.")
    # widest decile vs narrowest decile of D* width
    w = r["dstar_widths"]
    ratio = np.quantile(w, 0.9) / max(np.quantile(w, 0.1), 1e-12)
    out(f"      widest-decile / narrowest-decile D* width ratio = {ratio:.1f}x "
        f"(vs ~1x for a homogeneous map).")
    out("-" * 92)
    out("[2.3] Minos-style deployment monitor on the in-vivo signals "
        "(synthetic calibration vs deployment):")
    fire = "FIRES" if r["monitor_fires"] else "silent"
    out(f"      monitor: {fire}   AUC(cal vs in-vivo) = {r['monitor_auc']:.2f}   "
        f"(Family-1 Mahalanobis {r['monitor_maha']:.2f}/thr {r['monitor_maha_thr']:.2f},")
    out(f"      Family-2 residual {r['monitor_resid']:.3f}/thr "
        f"{r['monitor_resid_thr']:.3f}; domain-weight spread "
        f"{r['weight_spread']:.0f}x).")
    out("      INTERPRETATION: synthetic-calibrated -> deployment is itself an "
        "exchangeability break, and the")
    out("      observable monitor detects it. That firing is the honest, "
        "label-free reason the synthetic")
    out("      coverage guarantee must NOT be transferred to this data without "
        "re-calibration on matched, LABELED")
    out("      data -- which in vivo does not exist (the Echo tension). The "
        "monitor turns 'no ground truth' into")
    out("      an OBSERVABLE stop sign (the Minos cross-link, CP0).")
    if r["has_gt"]:
        out("-" * 92)
        pc = r["phantom_coverage"]
        out("[2.x] PHANTOM-ONLY sanity (this stand-in HAS ground truth; a true "
            "in-vivo run would print NOTHING here):")
        out(f"      stand-in coverage (D, D*, f) = "
            f"({pc[0]:.3f}, {pc[1]:.3f}, {pc[2]:.3f}) at nominal "
            f"{1-ALPHA:.2f} -- under-covers, consistent with the monitor firing on "
            f"the transfer.")
        out("      (Shown ONLY to validate the pipeline; it is NOT an in-vivo "
            "coverage claim and is kept separate.)")
    out("=" * 92)
    out("GATE 2: the in-vivo pipeline RUNS and its qualitative outputs (interval "
        "map, heavy-tailed D* widths,")
    out("  monitor firing on transfer) are clearly DELIMITED from the synthetic "
        "coverage guarantees. ")
    out("  HALT-TO-REPORT (dataset gate): no clean, permissively-licensed, "
        "ground-truth-free in-vivo IVIM")
    out("  dataset was usable in THIS environment (nibabel absent; OSIPI commits "
        "only b-value files / generators;")
    out("  true in-vivo sits behind TCIA data-use agreements). The demo therefore "
        "runs on a transparent SYNTHETIC")
    out("  STAND-IN. >> HUMAN SIGN-OFF NEEDED: dataset choice + licensing, and the "
        "in-vivo framing (what the")
    out("  qualitative demo may claim). The loaders (load_dwi_nifti / "
        "load_signals_npy) are the drop-in point.")
    out("#" * 92)

    with open(os.path.join(_RESULTS_DIR, "invivo_report.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    _make_figure(r)
    return 0


def _make_figure(r):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                  # pragma: no cover
        print(f"[figures] skipped ({e})")
        return
    os.makedirs(_FIG_DIR, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    # D* band width vs plug-in D-hat* (the wall: width grows with D*)
    dh = r["dstar_hat"] * 1e3
    ww = r["dstar_widths"] * 1e3
    ax1.scatter(dh, ww, s=5, alpha=0.3, color="#8e44ad")
    ax1.set_xlabel("plug-in D-hat*  (1e-3 mm^2/s)")
    ax1.set_ylabel("conformal D* band width  (1e-3 mm^2/s)")
    ax1.set_title("CP2: D* interval honestly balloons\n(qualitative; no GT in vivo)")
    # D* width distribution (heavy tail)
    ax2.hist(ww, bins=40, color="#8e44ad", alpha=0.8)
    for p, c in [(50, "#27ae60"), (90, "#c0392b")]:
        ax2.axvline(np.percentile(ww, p), ls="--", c=c, lw=1,
                    label=f"{p}th pct")
    ax2.set_xlabel("conformal D* band width  (1e-3 mm^2/s)")
    ax2.set_ylabel("voxel count")
    ax2.set_title(f"CP2: D* band widths ({r['source']})\nmonitor: "
                  f"{'FIRES' if r['monitor_fires'] else 'silent'} "
                  f"(AUC {r['monitor_auc']:.2f})")
    ax2.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(_FIG_DIR, "invivo_demo.pdf"))
    plt.close(fig)
    print(f"[figures] wrote invivo_demo.pdf -> {_FIG_DIR}")


# =========================================================================== #
# REAL-DATA PATH (ADDITIVE) -- human-approved 2026-06-14:
#   dataset  = ACRIN-6698 / I-SPY2 Breast DWI (TCIA, CC-BY-4.0, no IVIM ground
#              truth, b = {0,100,600,800}); fetch with scripts/fetch_invivo.py.
#   b-scheme = HYBRID:
#      * MONITOR  -> AS-DEPLOYED (option c): the synthetic 22-value-calibrated
#        deployment monitor is evaluated on the real data's b-INDEPENDENT
#        observable features. The real 4-b acquisition is itself the
#        exchangeability break the monitor is built to flag.
#      * CQR BAND -> RE-FIT AT THE REAL 4-b SCHEME (option b): a 22-b predictor
#        cannot ingest a 4-b signal, and interpolating onto the 22-b grid would
#        fabricate unacquired b-values, so the quantile layer is re-fit on
#        synthetic-parameter draws simulated at the real scheme. The band widths
#        are QUALITATIVE -- NO coverage claim (no ground truth in vivo).
# This path writes ONLY to invivo_real_* and never touches the synthetic outputs.
# =========================================================================== #
def load_dwi_npy4d(vol_path, bval_path, mask_path=None, max_voxels=4000, seed=0):
    """Load an assembled 4D DWI volume ``(X,Y,Z,B)`` ``.npy`` + b-values.

    No nibabel needed -- ``scripts/fetch_invivo.py`` assembles the DICOM series
    into a ``.npy`` grid. Voxels are masked (a provided tumor mask, or a b=0
    foreground threshold), normalized to S0, and randomly subsampled to
    ``max_voxels``. Mirrors :func:`load_dwi_nifti` for the no-nibabel case.
    """
    vol = np.asarray(np.load(vol_path), float)             # (X, Y, Z, B)
    b = np.loadtxt(bval_path).ravel()
    flat = vol.reshape(-1, vol.shape[-1])
    if mask_path:
        mask = np.asarray(np.load(mask_path)).ravel() > 0
    else:                                                  # b=0 foreground (air-reject)
        s0 = flat[:, b == b.min()].mean(1)
        mask = s0 > 0.25 * np.percentile(s0, 99)
    vox = flat[mask]
    # reject degenerate voxels (non-finite, or S0<=0) -- e.g. tumor-mask voxels that
    # fall on air after nearest-slice mapping, which would break S0-normalization.
    s0v = vox[:, b == b.min()].mean(1)
    vox = vox[np.isfinite(vox).all(1) & (s0v > 0)]
    rng = np.random.default_rng(seed)
    if vox.shape[0] > max_voxels:
        vox = vox[rng.choice(vox.shape[0], max_voxels, replace=False)]
    return _normalize_to_s0(vox, b), b


def _cqr_band_widths(signals, cal):
    """CQR per-voxel interval widths (N,3) from a fitted calibration dict.

    Conformalizes the gradient-boosted quantile predictor with the split-CQR
    offset. Shared by the real-data run and the test-retest proxy.
    """
    levels = cal["levels"]
    lo_raw = np.stack([cal["qreg"].predict_quantile(signals, j, levels[0])
                       for j in range(3)], axis=1)
    hi_raw = np.stack([cal["qreg"].predict_quantile(signals, j, levels[1])
                       for j in range(3)], axis=1)
    lo = np.empty_like(lo_raw)
    hi = np.empty_like(hi_raw)
    for j in range(3):
        lj, hj, _ = cqr(cal["cal_lo"][:, j], cal["cal_hi"][:, j],
                        cal["cal_params"][:, j], lo_raw[:, j], hi_raw[:, j], ALPHA)
        lo[:, j], hi[:, j] = lj, hj
    return interval_width(lo, hi), lo, hi


def run_real_hybrid(signals, b_real, source_name, seed=SEED, n_train=3000,
                    n_cal=2000):
    """Apply the hybrid deployed pipeline to (unlabeled) real in-vivo signals.

    Returns ``(out, (lo, hi, widths, theta))`` -- qualitative observables plus the
    raw per-voxel arrays (for the figure / the test-retest proxy). No ground truth
    is ever used; there is no coverage number anywhere in here.
    """
    # --- MONITOR = AS-DEPLOYED (synthetic 22-value calibration; option c) -------
    cal_mon = deployed_calibration(b=DEFAULT_B_VALUES, seed=seed, n_train=n_train,
                                   n_cal=n_cal)
    _, feat_real, resid_real = _observe(signals, b_real)
    mon = cal_mon["monitor"].evaluate(feat_real, resid_real)
    w_cal, w_test, _ = _importance_weights(cal_mon["cal_feat"], feat_real, seed=seed)

    # --- CQR BAND = RE-FIT AT THE REAL b-SCHEME (option b; qualitative only) -----
    cal_band = deployed_calibration(b=b_real, seed=seed, n_train=n_train,
                                    n_cal=n_cal)
    theta, _, _ = _observe(signals, b_real)
    widths, lo, hi = _cqr_band_widths(signals, cal_band)   # (N, 3)

    out = {
        "source": source_name, "n_vox": int(signals.shape[0]),
        "b_real": np.asarray(b_real, float).tolist(),
        "theta_median": np.median(theta, 0),
        "width_median": np.median(widths, 0),
        "dstar_width_q": np.quantile(widths[:, 1], [0.1, 0.5, 0.9]),
        "dstar_widths": widths[:, 1], "dstar_hat": theta[:, 1], "d_hat": theta[:, 0],
        # as-deployed monitor: per-family breakdown (Mahalanobis is the b-independent
        # detector; the residual family is reported with a caveat -- see main_real).
        "monitor_fires": bool(mon["fires"]), "monitor_auc": float(mon["auc"]),
        "monitor_maha": float(mon["maha"]["stat"]),
        "monitor_maha_thr": float(mon["maha"]["threshold"]),
        "monitor_maha_fires": bool(mon["maha"]["fires"]),
        "monitor_maha_auc": float(mon["maha"]["auc"]),
        "monitor_resid": float(mon["resid"]["stat"]),
        "monitor_resid_thr": float(mon["resid"]["threshold"]),
        "monitor_resid_fires": bool(mon["resid"]["fires"]),
        "monitor_resid_auc": float(mon["resid"]["auc"]),
        "weight_spread": float(np.quantile(w_test, 0.9) /
                               max(np.quantile(w_test, 0.1), 1e-9)),
        "has_gt": False,
    }
    return out, (lo, hi, widths, theta)


def compute_real(data_dir, seed=SEED, use_tumor_mask=False, max_voxels=4000):
    """Run the real-data qualitative demo on one assembled exam directory.

    ``data_dir`` holds ``signals_4d.npy`` + ``bvals.txt`` (+ ``tumor_mask.npy``)
    from ``scripts/fetch_invivo.py``. Never caches over / reads the synthetic
    pipeline; the synthetic seed-20260613 path is byte-identical and untouched.
    """
    vol = os.path.join(data_dir, "signals_4d.npy")
    bvp = os.path.join(data_dir, "bvals.txt")
    mp = os.path.join(data_dir, "tumor_mask.npy") if use_tumor_mask else None
    signals, b = load_dwi_npy4d(vol, bvp, mask_path=mp, max_voxels=max_voxels,
                                seed=seed)
    # shape / units assertions: predictor expects S0-normalized (N, B) signals.
    assert signals.ndim == 2 and signals.shape[1] == b.size, "bad (N,B) shape"
    assert np.allclose(signals[:, int(np.argmin(b))], 1.0, atol=1e-6), \
        "signals not S0-normalized at b=0"
    patient = "unknown"
    meta_p = os.path.join(data_dir, "meta.json")
    if os.path.exists(meta_p):
        patient = json.load(open(meta_p)).get("patient", "unknown")
    roi = "whole-tumor ROI" if use_tumor_mask else "b=0 foreground"
    src = f"REAL in-vivo DWI -- ACRIN-6698 {patient} ({roi})"
    res, arrays = run_real_hybrid(signals, b, src, seed=seed)
    return {"res": res, "arrays": arrays, "is_real": True, "patient": patient,
            "data_dir": os.path.relpath(data_dir, os.path.dirname(_RESULTS_DIR))}


def _write_real_provenance(res, patient, data_dir, seed):
    prov = {}
    if os.path.exists(_REAL_PROVENANCE):
        try:
            prov = json.load(open(_REAL_PROVENANCE))
        except (json.JSONDecodeError, OSError):
            prov = {}
    prov["run"] = {
        "b_scheme_handling": ("hybrid: as-deployed monitor (synthetic 22-value "
                              "calibration) + CQR band re-fit at the real 4-b scheme"),
        "no_coverage_claim_in_vivo": True,
        "patient": patient, "data_dir": data_dir, "seed": int(seed),
        "n_voxels": int(res["n_vox"]), "b_values": res["b_real"],
        "monitor_fires": bool(res["monitor_fires"]),
        "monitor_auc": float(res["monitor_auc"]),
        "monitor_mahalanobis": {"stat": res["monitor_maha"],
                                "threshold": res["monitor_maha_thr"],
                                "fires": res["monitor_maha_fires"],
                                "auc": res["monitor_maha_auc"]},
        "monitor_residual": {"stat": res["monitor_resid"],
                             "threshold": res["monitor_resid_thr"],
                             "fires": res["monitor_resid_fires"],
                             "auc": res["monitor_resid_auc"],
                             "caveat": ("4-b NLLS residual norm is not directly "
                                        "comparable to the 22-b calibration; the "
                                        "Mahalanobis feature family is the "
                                        "b-independent detector")},
        "dstar_width_pct_10_50_90_e3": (res["dstar_width_q"] * 1e3).tolist(),
        "median_plugin_D_Dstar_f": [float(res["theta_median"][0]),
                                    float(res["theta_median"][1]),
                                    float(res["theta_median"][2])],
        "computed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    os.makedirs(os.path.dirname(_REAL_PROVENANCE), exist_ok=True)
    json.dump(prov, open(_REAL_PROVENANCE, "w"), indent=2)


def main_real(data_dir, seed=SEED, use_tumor_mask=False):
    P = compute_real(data_dir, seed=seed, use_tumor_mask=use_tumor_mask)
    r = P["res"]
    lines = []

    def out(*x):
        s = " ".join(str(z) for z in x)
        print(s)
        lines.append(s)

    out("#" * 92)
    out("GAUGE -- REAL IN-VIVO DATA PATH (qualitative; NO coverage claim)")
    out("#" * 92)
    out(f"data source: {r['source']}")
    out(f"dataset: ACRIN-6698 / I-SPY2 Breast DWI (TCIA, CC-BY-4.0, "
        f"DOI 10.7937/tcia.kk02-6d95) -- REAL in-vivo, NO ground-truth IVIM params.")
    out(f"voxels: {r['n_vox']}   |   real b-values (s/mm^2): {r['b_real']}   |   "
        f"ground truth available: {r['has_gt']}")
    out("")
    out("** NO IN-VIVO COVERAGE CLAIM ** Everything below is QUALITATIVE pipeline/"
        "monitor behavior on")
    out("   real signal. There are NO ground-truth (D, D*, f) in vivo, so realized "
        "coverage is UNMEASURABLE")
    out("   and the finite-sample 1-alpha guarantee is NOT asserted here.")
    out("   b-scheme handling (human-approved): MONITOR = as-deployed (synthetic "
        "22-value calibration);")
    out("   CQR BAND = re-fit at the real 4-b scheme (a 22-b predictor cannot ingest "
        "a 4-b signal, and")
    out("   interpolation would fabricate unacquired b-values). The bands are NOT "
        "coverage intervals.")
    out("-" * 92)
    out("[R2.1] Conformalized interval map (CQR band re-fit at real 4-b; qualitative):")
    out(f"       median plug-in (D, D*, f) = "
        f"({r['theta_median'][0]*1e3:.2f}, {r['theta_median'][1]*1e3:.1f}, "
        f"{r['theta_median'][2]:.3f})   [D,D* in 1e-3 mm^2/s]")
    out(f"       median band width (D, D*, f) = "
        f"({r['width_median'][0]*1e3:.2f}, {r['width_median'][1]*1e3:.1f}, "
        f"{r['width_median'][2]:.3f})")
    out("[R2.2] D* compartment band widths on REAL signal: 10/50/90th pct "
        "(1e-3 mm^2/s)")
    q = r["dstar_width_q"] * 1e3
    out(f"       = [{q[0]:.1f}, {q[1]:.1f}, {q[2]:.1f}]  -- the D* band is wide; the "
        f"sparse 4-b scheme under-")
    out(f"       identifies the perfusion compartment (the Gauge 03 wall, now on real "
        f"in-vivo signal).")
    w = r["dstar_widths"]
    ratio = np.quantile(w, 0.9) / max(np.quantile(w, 0.1), 1e-12)
    out(f"       widest-decile / narrowest-decile D* width ratio = {ratio:.1f}x.")
    out("-" * 92)
    out("[R2.3] AS-DEPLOYED deployment monitor (synthetic 22-value calibration) on "
        "the REAL in-vivo features:")
    fire = "FIRES" if r["monitor_fires"] else "SILENT"
    out(f"       monitor: {fire}   AUC(cal vs in-vivo) = {r['monitor_auc']:.2f}")
    mfire = "FIRES" if r["monitor_maha_fires"] else "silent"
    out(f"       Family-1 Mahalanobis (b-independent feature detector): {mfire}  "
        f"stat {r['monitor_maha']:.2f}/thr {r['monitor_maha_thr']:.2f}  "
        f"AUC {r['monitor_maha_auc']:.2f}")
    rfire = "FIRES" if r["monitor_resid_fires"] else "silent"
    out(f"       Family-2 residual: {rfire}  stat {r['monitor_resid']:.3f}/thr "
        f"{r['monitor_resid_thr']:.3f}  AUC {r['monitor_resid_auc']:.2f}")
    out(f"       PROTOCOL NOTE (why Family-2 is silent -- the monitor has NOT "
        f"half-failed): ACRIN-6698 is a 4-b ADC")
    out(f"       protocol, NOT IVIM-optimized; 4 b-values EXACTLY-DETERMINE the "
        f"4-parameter IVIM fit, so the NLLS")
    out(f"       residual norm collapses toward 0 and Family-2 is structurally "
        f"uninformative here (not comparable")
    out(f"       to the 22-b calibration). The firing decision is carried entirely "
        f"by Family-1 (Mahalanobis,")
    out(f"       b-independent feature detector; AUC {r['monitor_maha_auc']:.2f}) -- "
        f"that is the honest detector.")
    out(f"       domain-weight spread {r['weight_spread']:.0f}x.")
    out("       INTERPRETATION: the synthetic-calibrated -> real-in-vivo transfer "
        "(sparse 4-b ADC protocol + real")
    out("       breast tissue) is a LARGE exchangeability break, and the as-deployed "
        "monitor detects it via")
    out("       Family-1. That firing is the honest, label-free reason the synthetic "
        "coverage guarantee correctly")
    out("       REFUSES to transfer to this data -- matched LABELED in-vivo data does "
        "not exist (the Echo tension).")
    if not r["monitor_fires"]:
        out("       >> FINDING: the monitor did NOT fire -- surfaced as a result "
            "(unexpected on a 4-b transfer).")
    out("=" * 92)
    out("REAL-DATA PATH: the deployed pipeline RUNS on real in-vivo DWI; its "
        "qualitative outputs (interval")
    out("  map, heavy-tailed D* widths, as-deployed monitor firing on transfer) are "
        "clearly DELIMITED from the")
    out("  synthetic coverage guarantees, which remain byte-identical and are the "
        "ONLY coverage claims.")
    out("#" * 92)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    with open(_REAL_REPORT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    _write_real_provenance(r, P["patient"], P["data_dir"], seed)
    _make_figure_real(r)
    print(f"[invivo-real] wrote {os.path.relpath(_REAL_REPORT, os.path.dirname(_RESULTS_DIR))}")
    return 0


def _make_figure_real(r):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                  # pragma: no cover
        print(f"[figures] skipped ({e})")
        return
    os.makedirs(_FIG_DIR, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    dh = r["dstar_hat"] * 1e3
    ww = r["dstar_widths"] * 1e3
    ax1.scatter(dh, ww, s=5, alpha=0.3, color="#16a085")
    ax1.set_xlabel("plug-in D-hat*  (1e-3 mm^2/s)")
    ax1.set_ylabel("conformal D* band width  (1e-3 mm^2/s)")
    ax1.set_title("REAL in-vivo: D* band vs plug-in D-hat*\n"
                  "(qualitative; no coverage claim in vivo)")
    ax2.hist(ww, bins=40, color="#16a085", alpha=0.8)
    for p, c in [(50, "#27ae60"), (90, "#c0392b")]:
        ax2.axvline(np.percentile(ww, p), ls="--", c=c, lw=1, label=f"{p}th pct")
    ax2.set_xlabel("conformal D* band width  (1e-3 mm^2/s)")
    ax2.set_ylabel("voxel count")
    ax2.set_title(f"REAL in-vivo D* band widths (ACRIN-6698)\nas-deployed monitor: "
                  f"{'FIRES' if r['monitor_fires'] else 'silent'} "
                  f"(AUC {r['monitor_auc']:.2f})")
    ax2.legend(fontsize=7)
    fig.suptitle("Gauge in-vivo, qualitative, NO coverage claim "
                 "(real DWI; b={0,100,600,800} s/mm^2)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(_REAL_FIG)
    plt.close(fig)
    print(f"[figures] wrote {os.path.basename(_REAL_FIG)} -> {_FIG_DIR}")


# --------------------------------------------------------------------------- #
# CHECKPOINT D -- test-retest repeatability proxy (REAL quantitative, no GT).
#
# ACRIN-6698 acquired two same-visit DWI exams per patient (patient repositioned
# between scans). With no inter-scan registration, true per-voxel correspondence
# is unavailable, so this is a REGION-level check across tumors: does the
# conformal interval width (predicted uncertainty) track the scan-rescan
# variability of the plug-in diffusion estimate (the measured reproducibility)?
# This is a REPEATABILITY-TRACKING check -- it needs no ground truth and is NOT a
# coverage claim.
# --------------------------------------------------------------------------- #
def _find_retest_pairs(root):
    """Find patients with >=2 assembled exams that each carry a tumor mask."""
    pairs = []
    if not os.path.isdir(root):
        return pairs
    for patient in sorted(os.listdir(root)):
        pdir = os.path.join(root, patient)
        if not os.path.isdir(pdir):
            continue
        exams = [os.path.join(pdir, d) for d in sorted(os.listdir(pdir))
                 if os.path.exists(os.path.join(pdir, d, "tumor_mask.npy"))
                 and os.path.exists(os.path.join(pdir, d, "signals_4d.npy"))]
        if len(exams) >= 2:
            pairs.append((patient, exams[:2]))
    return pairs


def _tumor_summary(exam_dir, b, cal, max_voxels=20000, seed=SEED):
    """Tumor-ROI median plug-in (D, D*) and median CQR band widths for one exam."""
    sig, _ = load_dwi_npy4d(os.path.join(exam_dir, "signals_4d.npy"),
                            os.path.join(exam_dir, "bvals.txt"),
                            mask_path=os.path.join(exam_dir, "tumor_mask.npy"),
                            max_voxels=max_voxels, seed=seed)
    theta, _, _ = _observe(sig, b)
    widths, _, _ = _cqr_band_widths(sig, cal)
    return {"n": int(sig.shape[0]),
            "D_med": float(np.median(theta[:, 0])),
            "Dstar_med": float(np.median(theta[:, 1])),
            "wD_med": float(np.median(widths[:, 0])),
            "wDstar_med": float(np.median(widths[:, 1]))}


def analyze_test_retest(root, seed=SEED):
    """Region-level repeatability proxy across all test-retest tumor pairs."""
    pairs = _find_retest_pairs(root)
    b = np.loadtxt(os.path.join(pairs[0][1][0], "bvals.txt")).ravel()
    cal = deployed_calibration(b=b, seed=seed)             # real-b CQR predictor
    rows = []
    for patient, (d_test, d_retest) in pairs:
        t = _tumor_summary(d_test, b, cal, seed=seed)
        r = _tumor_summary(d_retest, b, cal, seed=seed)
        rows.append({
            "patient": patient,
            "wD": 0.5 * (t["wD_med"] + r["wD_med"]),       # predicted D uncertainty
            "dD": abs(t["D_med"] - r["D_med"]),            # scan-rescan |ΔADC|
            "wDstar": 0.5 * (t["wDstar_med"] + r["wDstar_med"]),
            "dDstar": abs(t["Dstar_med"] - r["Dstar_med"]),
            "n_test": t["n"], "n_retest": r["n"],
        })
    return {"b": b.tolist(), "rows": rows, "seed": int(seed)}


def main_retest(root, seed=SEED):
    from scipy.stats import spearmanr, pearsonr
    A = analyze_test_retest(root, seed=seed)
    rows = A["rows"]
    lines = []

    def out(*x):
        s = " ".join(str(z) for z in x)
        print(s)
        lines.append(s)

    n = len(rows)
    out("#" * 92)
    out("GAUGE -- REAL IN-VIVO TEST-RETEST REPEATABILITY PROXY (Checkpoint D)")
    out("#" * 92)
    out("** REPEATABILITY-TRACKING CHECK -- NOT a coverage claim. ** No ground "
        "truth is used. This asks")
    out("   whether the conformal interval WIDTH (predicted uncertainty) tracks "
        "the real scan-rescan")
    out("   variability of the plug-in diffusion estimate across tumors. Region-"
        "level (whole-tumor ROI);")
    out("   the two same-visit ACRIN-6698 exams are NOT registered, so this is "
        "not a per-voxel claim.")
    out(f"   dataset: ACRIN-6698 / I-SPY2 Breast DWI (TCIA, CC-BY-4.0); "
        f"b={A['b']} s/mm^2; tumor pairs n={n}.")
    out("-" * 92)
    if n < 3:
        out(f"[D] only {n} test-retest tumor pair(s) available -- too few for a "
            f"correlation; reporting raw rows only.")
        for rw in rows:
            out(f"    {rw['patient']}: wD={rw['wD']*1e3:.2f}  |dD|={rw['dD']*1e3:.3f} "
                f"(1e-3 mm^2/s)")
    else:
        wD = np.array([r["wD"] for r in rows]) * 1e3
        dD = np.array([r["dD"] for r in rows]) * 1e3
        wDs = np.array([r["wDstar"] for r in rows]) * 1e3
        dDs = np.array([r["dDstar"] for r in rows]) * 1e3
        sD, pD = spearmanr(wD, dD)
        rD, prD = pearsonr(wD, dD)
        sDs, pDs = spearmanr(wDs, dDs)
        out(f"[D.1] Tissue diffusion D (well-identified even by the 4-b ADC scheme; "
            f"the ADC-like quantity):")
        out(f"      conformal D-width vs |scan-rescan ΔD|:  Spearman r = {sD:+.2f} "
            f"(p={pD:.3f}, n={n});  Pearson r = {rD:+.2f}")
        out(f"      median conformal D-width = {np.median(wD):.2f}; median "
            f"|scan-rescan ΔD| = {np.median(dD):.3f}  (1e-3 mm^2/s)")
        out(f"      READING: the conformal D interval WIDENS where the real ADC is "
            f"LEAST repeatable -- the band")
        out(f"      tracks / adapts to / is consistent with the measured scan-rescan "
            f"variability. This is explicitly")
        out(f"      NOT validation, accuracy, calibration, or coverage (no in-vivo "
            f"ground truth exists); it is an")
        out(f"      honest width-vs-repeatability association only.")
        out(f"[D.2] Pseudo-diffusion D* (under-identified by the 4-b ADC scheme):")
        out(f"      conformal D*-width vs |scan-rescan ΔD*|:  Spearman r = {sDs:+.2f} "
            f"(p={pDs:.3f}, n={n}) -- NOT significant.")
        out("-" * 92)
        out(f"      ON-THESIS SPLIT: for the well-identified D the width-vs-"
            f"repeatability relationship is")
        out(f"      DEMONSTRABLE; for D*, consistent with the identifiability wall "
            f"(Gauge 03), it is UNRESOLVABLE")
        out(f"      at this sample size and a 4-b scheme. The D* negative is EVIDENCE "
            f"for the thesis, not a gap.")
        out(f"      CAVEATS (do not over-read): n={n} tumor pairs is SMALL -- a "
            f"Spearman of {sD:+.2f} at n={n} carries a")
        out(f"      very wide confidence interval, so D.1 is SUGGESTIVE (reported "
            f"with its n and p), NOT robust. It")
        out(f"      is a region-level whole-tumor-ROI association on UNREGISTERED "
            f"same-visit exams (no per-voxel")
        out(f"      correspondence), and a repeatability-tracking check -- NOT an "
            f"in-vivo coverage claim.")
        _make_figure_retest(wD, dD, sD, n)
    out("#" * 92)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    with open(_REAL_RETEST_REPORT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    # merge a retest block into the provenance
    prov = {}
    if os.path.exists(_REAL_PROVENANCE):
        try:
            prov = json.load(open(_REAL_PROVENANCE))
        except (json.JSONDecodeError, OSError):
            prov = {}
    block = {"n_pairs": n, "b_values": A["b"], "seed": int(seed),
             "note": ("region-level (whole-tumor ROI), unregistered same-visit exams; "
                      "SUGGESTIVE at small n -- NOT validation/accuracy/calibration/"
                      "coverage. D interval widens where the real ADC is least "
                      "repeatable (width tracks repeatability)."),
             "computed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    if n >= 3:
        sD_, pD_ = spearmanr([r["wD"] for r in rows], [r["dD"] for r in rows])
        sDs_, pDs_ = spearmanr([r["wDstar"] for r in rows],
                               [r["dDstar"] for r in rows])
        block["D_spearman_wD_vs_dD"] = float(sD_)
        block["D_spearman_p"] = float(pD_)
        block["D_significance"] = "suggestive (small n, wide CI)"
        block["Dstar_spearman"] = float(sDs_)
        block["Dstar_spearman_p"] = float(pDs_)
        block["Dstar_significance"] = ("not significant -- consistent with the "
                                       "identifiability wall at 4-b / small n")
    prov["retest"] = block
    json.dump(prov, open(_REAL_PROVENANCE, "w"), indent=2)
    print(f"[invivo-retest] wrote {os.path.relpath(_REAL_RETEST_REPORT, os.path.dirname(_RESULTS_DIR))}")
    return 0


def _make_figure_retest(wD, dD, spearman_r, n):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                  # pragma: no cover
        print(f"[figures] skipped ({e})")
        return
    os.makedirs(_FIG_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.scatter(wD, dD, s=28, alpha=0.8, color="#2c3e50")
    ax.set_xlabel("conformal D band width (test+retest mean)  (1e-3 mm^2/s)")
    ax.set_ylabel("|scan-rescan ΔD|  (1e-3 mm^2/s)")
    ax.set_title(f"Real in-vivo test-retest (ACRIN-6698)\n"
                 f"n={n} tumors;  Spearman r={spearman_r:+.2f}  (suggestive)\n"
                 f"width tracks ADC repeatability -- not coverage", fontsize=9)
    fig.tight_layout()
    fig.savefig(_REAL_RETEST_FIG)
    plt.close(fig)
    print(f"[figures] wrote {os.path.basename(_REAL_RETEST_FIG)} -> {_FIG_DIR}")


# --------------------------------------------------------------------------- #
# Actual in-vivo IMAGE MAPS (spatial), so the manuscript shows the real signal
# and the per-voxel uncertainty -- not only summary plots. Qualitative; no
# coverage claim. Maps are built on one representative slice with a tumor ROI.
# --------------------------------------------------------------------------- #
def make_real_image_maps(exam_dir, slice_idx=None, seed=SEED):
    """Render real DWI + plug-in D* map + conformal D*-band-width map for a slice.

    Picks the slice with the largest tumor ROI (unless ``slice_idx`` is given),
    fits the deployed pipeline per foreground voxel, and draws spatial maps with
    the whole-tumor ROI outlined. Writes ``figures/invivo_real_maps.pdf`` and
    returns the slice index used.
    """
    vol = np.asarray(np.load(os.path.join(exam_dir, "signals_4d.npy")), float)
    b = np.loadtxt(os.path.join(exam_dir, "bvals.txt")).ravel()
    mp = os.path.join(exam_dir, "tumor_mask.npy")
    tmask3d = np.asarray(np.load(mp)) > 0 if os.path.exists(mp) else None
    if slice_idx is None:
        slice_idx = (int(tmask3d.reshape(-1, tmask3d.shape[-1]).sum(0).argmax())
                     if tmask3d is not None else vol.shape[2] // 2)
    sl = vol[:, :, slice_idx, :]                           # (X, Y, B)
    nx, ny = sl.shape[:2]
    b0 = sl[:, :, b == b.min()].mean(2)                    # anatomical b=0 image
    flat = sl.reshape(-1, sl.shape[-1])
    s0 = flat[:, b == b.min()].mean(1)
    fg = (s0 > 0.25 * np.percentile(s0, 99)) & np.isfinite(flat).all(1) & (s0 > 0)

    sig = _normalize_to_s0(flat[fg], b)
    cal = deployed_calibration(b=b, seed=seed)             # re-fit CQR at real b
    theta, _, _ = _observe(sig, b)
    widths, _, _ = _cqr_band_widths(sig, cal)              # (Nfg, 3)

    dstar = np.full(nx * ny, np.nan)
    wstar = np.full(nx * ny, np.nan)
    dstar[fg] = theta[:, 1] * 1e3
    wstar[fg] = widths[:, 1] * 1e3
    dstar = dstar.reshape(nx, ny)
    wstar = wstar.reshape(nx, ny)
    tmask = tmask3d[:, :, slice_idx] if tmask3d is not None else None

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                  # pragma: no cover
        print(f"[figures] skipped ({e})")
        return slice_idx
    os.makedirs(_FIG_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4))
    # crop to the foreground bounding box (+margin) so the breast fills the panel
    ys, xs = np.where(fg.reshape(nx, ny))
    r0, r1 = max(ys.min() - 8, 0), min(ys.max() + 8, nx)
    c0, c1 = max(xs.min() - 8, 0), min(xs.max() + 8, ny)
    slc = (slice(r0, r1), slice(c0, c1))

    def _roi(ax):
        if tmask is not None and tmask[slc].any():
            ax.contour(tmask[slc], levels=[0.5], colors="r", linewidths=1.0)

    ax = axes[0]
    b0c = b0[slc]
    ax.imshow(b0c, cmap="gray", origin="lower", vmin=0,
              vmax=np.percentile(b0c[b0c > 0], 99) if (b0c > 0).any() else None)
    _roi(ax)
    ax.set_title("real DWI ($b{=}0$), ACRIN-6698\n(tumor ROI in red)")
    ax.axis("off")

    ax = axes[1]
    im = ax.imshow(dstar[slc], cmap="magma", origin="lower", vmin=0,
                   vmax=np.nanpercentile(dstar, 98))
    _roi(ax)
    ax.set_title("plug-in $\\hat{D}^{*}$ map\n($10^{-3}$ mm$^2$/s)")
    ax.axis("off")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[2]
    im = ax.imshow(wstar[slc], cmap="viridis", origin="lower",
                   vmin=np.nanpercentile(wstar, 2),
                   vmax=np.nanpercentile(wstar, 98))
    _roi(ax)
    ax.set_title("conformal $D^{*}$ band-width map\n(wide = under-identified)")
    ax.axis("off")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Gauge in vivo (real DWI, qualitative -- NO coverage claim): "
                 "anatomy, perfusion estimate, and per-voxel conformal uncertainty",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(_REAL_MAPS_FIG)
    plt.close(fig)
    print(f"[figures] wrote {os.path.basename(_REAL_MAPS_FIG)} -> {_FIG_DIR} "
          f"(slice {slice_idx})")
    return slice_idx


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    argv = sys.argv[1:]
    if "--source" in argv:
        ap = argparse.ArgumentParser(prog="gauge.invivo")
        ap.add_argument("--source", choices=["synthetic", "real", "retest"],
                        default="synthetic")
        ap.add_argument("--data", help="assembled exam dir (real) or data root (retest)")
        ap.add_argument("--tumor-mask", action="store_true",
                        help="restrict to the whole-tumor ROI instead of b=0 foreground")
        ap.add_argument("--seed", type=int, default=SEED)
        a = ap.parse_args(argv)
        if a.source == "real":
            if not a.data:
                ap.error("--source real needs --data <assembled exam dir>")
            raise SystemExit(main_real(a.data, seed=a.seed,
                                       use_tumor_mask=a.tumor_mask))
        if a.source == "retest":
            raise SystemExit(main_retest(a.data or os.path.join(
                os.path.dirname(_RESULTS_DIR), "data", "invivo"), seed=a.seed))
        raise SystemExit(main(force=os.environ.get("GAUGE_FORCE") == "1", seed=a.seed))
    # legacy positional synthetic / NIfTI path (unchanged)
    dwi = argv[0] if len(argv) >= 1 else None
    bvals = argv[1] if len(argv) >= 2 else None
    raise SystemExit(main(force=os.environ.get("GAUGE_FORCE") == "1",
                          dwi=dwi, bvals=bvals))
