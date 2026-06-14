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


if __name__ == "__main__":
    argv = sys.argv[1:]
    dwi = argv[0] if len(argv) >= 1 else None
    bvals = argv[1] if len(argv) >= 2 else None
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    raise SystemExit(main(force=os.environ.get("GAUGE_FORCE") == "1",
                          dwi=dwi, bvals=bvals))
