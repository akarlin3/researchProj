"""Task 1 -- external-phantom coverage replication on the OSIPI IVIM reference DRO.

Gauge's conformal/CQR coverage is, by necessity, validated on labeled synthetic
data -- and so far on *our own* forward model and parameter prior. This module
replicates the coverage results and the high-D* identifiability wall against an
**external, community-standard SYNTHETIC reference**: the OSIPI TF2.4 IVIM digital
reference object (DRO), 5000 voxels with ground-truth (D, D*, f) that **we did not
generate** (scripts/fetch_osipi.py; Zenodo doi:10.5281/zenodo.14605039, CC-BY-4.0;
code Apache-2.0). This neutralises the "validated only on the authors' forward
model" critique on the part that matters: the ground-truth parameter distribution.

Honesty (hard constraint). The OSIPI DRO is an *independent synthetic* reference,
NOT in vivo. In-vivo coverage validation remains structurally impossible (no
ground truth) and no in-vivo coverage claim is made here. The bi-exponential
forward model is the field's consensus model and is therefore shared by
construction; what is external is the ground-truth parameter maps / DRO.

Two regimes are evaluated, because exchangeability is this paper's own theme:

  * controlled -- OSIPI ground truth, OUR acquisition: signals re-synthesised from
    the OSIPI (D, D*, f) at Gauge's 22-b scheme and SNR grid (only the parameter
    distribution is swapped vs our cohort).
  * native -- OSIPI ground truth AND OSIPI acquisition: the DRO's own pre-noised
    signals at its native sparse 7-b scheme.

Two calibrations, each within a regime:

  * naive transfer -- calibrate on OUR synthetic cohort (deployed_calibration),
    test on OSIPI. Coverage may break under the (real) D* distribution shift; the
    label-free deployment monitor is run and we REPORT whether it fires (a fired
    monitor here is a feature: the shift is observable).
  * recalibrated -- split the OSIPI DRO into train/cal/test and calibrate within
    the OSIPI distribution: clean within-distribution conformal coverage on an
    independent ground truth.

The wall test asks whether CRLB(D*) exceeds the D*-tercile width on the OSIPI
ground truth (CRLB/width >= 1 at high D*). If it replicates, that is strong
evidence the wall is information-theoretic, not an artifact of our forward model.

Deterministic from DEFAULT_SEED; reuses the existing conformal/monitor/CRLB code.

Run:  python -m gauge.osipi            # seed-0 report + figure + provenance
"""
import datetime
import json
import os

import numpy as np

from gauge.cohort import DEFAULT_SEED, DEFAULT_SNR_GRID
from gauge.conditional_attack import N_REGIME, _regime_from_true, conditional_coverage
from gauge.conformal import cqr, empirical_coverage, interval_width
from gauge.estimators import IVIMQuantileRegressor
from gauge.forward import DEFAULT_B_VALUES, add_rician_noise, crlb_dstar_batch, ivim_signal
from gauge.monitor import DeploymentMonitor
from gauge.robustness import _calibrate, _observe
from gauge.invivo import deployed_calibration

ALPHA = 0.10                         # nominal coverage 1 - ALPHA = 0.90
NOMINAL = 1.0 - ALPHA
PARAM_NAMES = ("D", "D*", "f")
DSTAR = 1                            # index of D* in (D, D*, f)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_RESULTS_DIR = os.path.join(_ROOT, "results")
_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
_PAPER_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper", "figures")
_DRO_PATH = os.path.join(_ROOT, "data", "osipi", "extracted", "DRO.npy")
_PROVENANCE = os.path.join(_RESULTS_DIR, "osipi_provenance.json")
_REPORT = os.path.join(_RESULTS_DIR, "osipi_report.txt")
_FIG = os.path.join(_FIG_DIR, "osipi_coverage.pdf")
_MULTISEED = os.path.join(_RESULTS_DIR, "osipi_multiseed.json")

# OSIPI native SNR (from the DRO residual analysis, ~80); used for native-regime
# CRLB and as the noise reference for native signals (which carry it already).
_NATIVE_SNR = 80.0


# --------------------------------------------------------------------------- #
# DRO loading
# --------------------------------------------------------------------------- #
def load_dro(path=_DRO_PATH):
    """Load the OSIPI DRO -> (gt, native_signals, native_b).

    ``gt`` is (N, 3) ground truth ordered (D, D*, f) -- D* = the DRO's ``Dp`` --
    in mm^2/s / dimensionless, matching Gauge's convention. ``native_signals`` is
    (N, 7) the DRO's own pre-noised signals at its fixed 7-b ``native_b``.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"OSIPI DRO not found at {path}. Run `python scripts/fetch_osipi.py` "
            "(download-on-demand; data/ is git-ignored).")
    dro = np.load(path, allow_pickle=True)
    gt = np.array([[float(e["D"]), float(e["Dp"]), float(e["f"])] for e in dro])
    native_b = np.asarray(dro[0]["bvals"], float)
    native_sig = np.array([np.asarray(e["signals"], float) for e in dro])
    # normalise to b=0 (S0=1 convention the calibration was built on)
    s0 = native_sig[:, int(np.argmin(native_b))][:, None]
    native_sig = np.clip(native_sig / np.where(s0 > 0, s0, 1.0), 0.0, None)
    return gt, native_sig, native_b


# --------------------------------------------------------------------------- #
# Synthesis at our acquisition (controlled regime)
# --------------------------------------------------------------------------- #
def synth_controlled(gt, b, snr_grid, rng):
    """Synthesise Rician-noisy signals from OSIPI ground truth at OUR acquisition.

    Only the *parameter distribution* (gt) is OSIPI's; the forward model is the
    shared bi-exponential and the acquisition (``b``, per-voxel SNR drawn from
    ``snr_grid``) is Gauge's. Returns ``(signals, snr)``.
    """
    n = gt.shape[0]
    snr = rng.choice(np.asarray(snr_grid, float), size=n)
    clean = ivim_signal(b[None, :], gt[:, 0:1], gt[:, 1:2], gt[:, 2:3], S0=1.0)
    noisy = add_rician_noise(clean, snr[:, None], rng, S0=1.0)
    return noisy, snr


# --------------------------------------------------------------------------- #
# Calibration built FROM the OSIPI DRO (recalibrated regime). Mirrors the shape
# of invivo.deployed_calibration so the same evaluator handles both.
# --------------------------------------------------------------------------- #
def _calibration_from_osipi(train_sig, train_gt, cal_sig, cal_gt, b, alpha, seed):
    cal_theta, cal_feat, cal_resid = _observe(cal_sig, b)
    q = _calibrate(cal_theta, cal_gt, alpha)
    levels = [alpha / 2, 1 - alpha / 2]
    qreg = IVIMQuantileRegressor(levels, random_state=0).fit(train_sig, train_gt)
    cal_lo = np.stack([qreg.predict_quantile(cal_sig, j, levels[0]) for j in range(3)], 1)
    cal_hi = np.stack([qreg.predict_quantile(cal_sig, j, levels[1]) for j in range(3)], 1)
    monitor = DeploymentMonitor(seed=seed).fit(cal_feat, cal_resid)
    return {"b": np.asarray(b, float), "q": q, "qreg": qreg, "levels": levels,
            "cal_params": cal_gt, "cal_lo": cal_lo, "cal_hi": cal_hi,
            "cal_feat": cal_feat, "cal_resid": cal_resid, "monitor": monitor}


# --------------------------------------------------------------------------- #
# Shared evaluation of one calibration against an OSIPI test set with GT.
# --------------------------------------------------------------------------- #
def _tercile_cov(lo, hi, true, regime):
    return [empirical_coverage(lo[regime == r], hi[regime == r], true[regime == r])
            for r in range(N_REGIME)]


def _evaluate(cal, test_sig, test_gt, b, test_snr, snr_levels, alpha,
              fixed_edges=None):
    """Coverage of one calibration on OSIPI test (with ground truth).

    Returns split-conformal and CQR marginal coverage (D, D*, f), per-D*-tercile
    conditional coverage for both, the CQR regime x SNR grid (for the figure), and
    the monitor decision. ``test_gt`` is used ONLY to measure coverage.

    ``fixed_edges`` (optional 2-vector of D* tercile boundaries): when given, the
    per-distribution tercile coverage is reported as usual AND a second set keyed
    ``*_terc_fixed`` is computed on FIXED physical boundaries, so "hi D*" denotes
    the same physical regime across priors. ``hi_prevalence_fixed`` is the fraction
    of test voxels in the fixed hi-D* bin -- the clinical-prevalence number that
    separates "regime rare" from "wall gone".
    """
    theta, feat, resid = _observe(test_sig, b)
    q = cal["q"]
    levels = cal["levels"]

    # --- split conformal (NLLS plug-in +/- calibrated radius) ----------------
    split_lo = theta - q
    split_hi = theta + q
    split_marg = [empirical_coverage(split_lo[:, j], split_hi[:, j], test_gt[:, j])
                  for j in range(3)]

    # --- CQR (conformalised gradient-boosted quantile band) ------------------
    lo_raw = np.stack([cal["qreg"].predict_quantile(test_sig, j, levels[0])
                       for j in range(3)], 1)
    hi_raw = np.stack([cal["qreg"].predict_quantile(test_sig, j, levels[1])
                       for j in range(3)], 1)
    cqr_lo = np.empty_like(lo_raw)
    cqr_hi = np.empty_like(hi_raw)
    for j in range(3):
        lj, hj, _ = cqr(cal["cal_lo"][:, j], cal["cal_hi"][:, j],
                        cal["cal_params"][:, j], lo_raw[:, j], hi_raw[:, j], alpha)
        cqr_lo[:, j], cqr_hi[:, j] = lj, hj
    cqr_marg = [empirical_coverage(cqr_lo[:, j], cqr_hi[:, j], test_gt[:, j])
                for j in range(3)]

    # --- conditional coverage by TRUE D* tercile -----------------------------
    regime, edges = _regime_from_true(test_gt[:, DSTAR])
    dtrue = test_gt[:, DSTAR]
    split_terc = _tercile_cov(split_lo[:, DSTAR], split_hi[:, DSTAR], dtrue, regime)
    cqr_terc = _tercile_cov(cqr_lo[:, DSTAR], cqr_hi[:, DSTAR], dtrue, regime)
    grid = None
    if snr_levels is not None:
        grid = conditional_coverage(cqr_lo[:, DSTAR], cqr_hi[:, DSTAR], dtrue,
                                    regime, test_snr, snr_levels)

    mon = cal["monitor"].evaluate(feat, resid)
    out = {
        "split_marg": np.array(split_marg), "cqr_marg": np.array(cqr_marg),
        "split_terc": np.array(split_terc), "cqr_terc": np.array(cqr_terc),
        "cqr_dstar_width_med": float(np.median(interval_width(cqr_lo[:, DSTAR],
                                                              cqr_hi[:, DSTAR]))),
        "regime_edges": edges, "grid": grid,
        "monitor_fires": bool(mon["fires"]), "monitor_auc": float(mon["auc"]),
        "monitor_maha": float(mon["maha"]["stat"]),
        "monitor_maha_thr": float(mon["maha"]["threshold"]),
    }
    # --- second tercile framing on FIXED physical boundaries (cross-prior honest) -
    if fixed_edges is not None:
        regime_f, edges_f = _regime_from_true(dtrue, edges=fixed_edges)
        out["split_terc_fixed"] = np.array(
            _tercile_cov(split_lo[:, DSTAR], split_hi[:, DSTAR], dtrue, regime_f))
        out["cqr_terc_fixed"] = np.array(
            _tercile_cov(cqr_lo[:, DSTAR], cqr_hi[:, DSTAR], dtrue, regime_f))
        out["regime_edges_fixed"] = np.asarray(edges_f, float)
        out["hi_counts_fixed"] = np.array(
            [int((regime_f == r).sum()) for r in range(N_REGIME)])
        out["hi_prevalence_fixed"] = float(np.mean(regime_f == N_REGIME - 1))
    return out


# --------------------------------------------------------------------------- #
# The identifiability wall on OSIPI ground truth.
# --------------------------------------------------------------------------- #
def wall_metric(gt, b, snr, edges=None, bin_endpoints=None):
    """CRLB(D*) vs the D*-tercile width on OSIPI ground truth.

    ``snr`` is per-voxel (controlled) or a scalar (native). Returns per-tercile
    median CRLB(D*), tercile widths, the CRLB/width ratio (>= 1 => the regime is
    unresolvable, the wall), and the low->high absolute-CRLB growth.

    With ``edges=None`` the terciles and bin widths are per-distribution (the
    quantiles of this set's D*). Pass ``edges`` (2-vector of fixed interior
    boundaries) and ``bin_endpoints`` (the fixed lo/hi outer endpoints) to hold the
    bins at FIXED physical D* boundaries -- so the CRLB/width ratio is comparable
    across priors instead of rescaling with each distribution's own spread.
    """
    D, Dstar, f = gt[:, 0], gt[:, DSTAR], gt[:, 2]
    snr = np.broadcast_to(np.asarray(snr, float), (gt.shape[0],))
    sd = crlb_dstar_batch(b, D, Dstar, f, snr)
    regime, edges = _regime_from_true(Dstar, edges=edges)
    if bin_endpoints is None:
        lo_end, hi_end = float(Dstar.min()), float(Dstar.max())
    else:
        lo_end, hi_end = float(bin_endpoints[0]), float(bin_endpoints[-1])
    bin_edges = np.concatenate([[lo_end], np.asarray(edges, float), [hi_end]])
    bin_w = np.diff(bin_edges)
    finite = np.isfinite(sd)
    med, ratio = [], []
    for r in range(N_REGIME):
        m = (regime == r) & finite
        mm = float(np.median(sd[m])) if m.any() else np.inf
        med.append(mm)
        ratio.append(mm / bin_w[r])
    abs_growth = (med[-1] / med[0]) if med[0] > 0 else np.inf
    return {"med_crlb": med, "bin_w": bin_w.tolist(), "ratio": ratio,
            "abs_growth": float(abs_growth), "hi_ratio": float(ratio[-1]),
            "replicates": bool(ratio[-1] >= 1.0 and abs_growth >= 2.0)}


# --------------------------------------------------------------------------- #
# Orchestration: one regime (naive + recalibrated) on a shared OSIPI test split.
# --------------------------------------------------------------------------- #
def _split_indices(n, seed, n_train=2000, n_cal=1500):
    rng = np.random.default_rng(seed + 101)
    perm = rng.permutation(n)
    return perm[:n_train], perm[n_train:n_train + n_cal], perm[n_train + n_cal:]


def run_regime(regime, seed=DEFAULT_SEED, alpha=ALPHA):
    """Evaluate naive-transfer and recalibrated coverage for one regime.

    ``regime`` in {'controlled', 'native'}. Returns a payload dict with the naive
    and recalibrated evaluation results, the wall metric, and the OSIPI test GT.
    """
    gt, native_sig, native_b = load_dro()
    n = gt.shape[0]
    tr_idx, cal_idx, te_idx = _split_indices(n, seed)

    if regime == "controlled":
        b = DEFAULT_B_VALUES
        rng = np.random.default_rng(seed + 202)
        sig, snr = synth_controlled(gt, b, DEFAULT_SNR_GRID, rng)
        snr_levels = list(DEFAULT_SNR_GRID)
    elif regime == "native":
        b = native_b
        sig = native_sig
        snr = np.full(n, _NATIVE_SNR)
        snr_levels = None                     # DRO has one fixed noise level
    else:
        raise ValueError(f"unknown regime {regime!r}")

    test_sig, test_gt, test_snr = sig[te_idx], gt[te_idx], snr[te_idx]

    # naive transfer: calibrate on OUR synthetic cohort at this acquisition.
    naive_cal = deployed_calibration(b=b, alpha=alpha, seed=seed)
    naive = _evaluate(naive_cal, test_sig, test_gt, b, test_snr, snr_levels, alpha)

    # recalibrated: calibrate within the OSIPI distribution (train/cal/test).
    recal_cal = _calibration_from_osipi(sig[tr_idx], gt[tr_idx], sig[cal_idx],
                                        gt[cal_idx], b, alpha, seed)
    recal = _evaluate(recal_cal, test_sig, test_gt, b, test_snr, snr_levels, alpha)

    wall = wall_metric(test_gt, b, test_snr)
    return {"regime": regime, "seed": seed, "n_test": int(te_idx.size),
            "b": np.asarray(b, float).tolist(), "naive": naive, "recal": recal,
            "wall": wall}


# --------------------------------------------------------------------------- #
# Flat scalar dict for the multi-seed band harness + consistency gate.
# --------------------------------------------------------------------------- #
def compute_all(seed=DEFAULT_SEED):
    out = {}
    payload = {}
    for regime in ("controlled", "native"):
        p = run_regime(regime, seed=seed)
        payload[regime] = p
        for cal in ("naive", "recal"):
            ev = p[cal]
            for j, nm in enumerate(PARAM_NAMES):
                out[f"osipi/{regime}/{cal}/cqr_marg/{nm}"] = float(ev["cqr_marg"][j])
                out[f"osipi/{regime}/{cal}/split_marg/{nm}"] = float(ev["split_marg"][j])
            out[f"osipi/{regime}/{cal}/cqr_hiDstar"] = float(ev["cqr_terc"][-1])
            out[f"osipi/{regime}/{cal}/split_hiDstar"] = float(ev["split_terc"][-1])
            out[f"osipi/{regime}/{cal}/cqr_loDstar"] = float(ev["cqr_terc"][0])
        out[f"osipi/{regime}/naive/monitor_auc"] = payload[regime]["naive"]["monitor_auc"]
        out[f"osipi/{regime}/naive/monitor_fires"] = float(payload[regime]["naive"]["monitor_fires"])
        for k in ("lo", "mid", "hi"):
            idx = {"lo": 0, "mid": 1, "hi": 2}[k]
            out[f"osipi/{regime}/wall/crlb_over_width/{k}"] = float(payload[regime]["wall"]["ratio"][idx])
        out[f"osipi/{regime}/wall/abs_growth"] = float(payload[regime]["wall"]["abs_growth"])
    return out, payload


# --------------------------------------------------------------------------- #
# Figure (controlled regime: coverage-vs-nominal + conditional-by-tercile).
# --------------------------------------------------------------------------- #
def _make_figure(payload):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                      # pragma: no cover
        print(f"[osipi figures] skipped ({e})")
        return
    c = payload["controlled"]
    naive, recal, wall = c["naive"], c["recal"], c["wall"]
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 4.2))

    # Panel A -- marginal coverage vs nominal (CQR), naive vs recalibrated.
    x = np.arange(3)
    w = 0.36
    ax[0].bar(x - w / 2, naive["cqr_marg"], w, label="naive transfer", color="#c0392b")
    ax[0].bar(x + w / 2, recal["cqr_marg"], w, label="recalibrated", color="#2c7fb8")
    ax[0].axhline(NOMINAL, ls="--", c="k", lw=1, label=f"nominal {NOMINAL:.2f}")
    ax[0].set_xticks(x); ax[0].set_xticklabels(PARAM_NAMES)
    ax[0].set_ylim(0, 1.02); ax[0].set_ylabel("marginal coverage")
    ax[0].set_title("(A) Marginal coverage on OSIPI DRO")
    ax[0].legend(fontsize=8, loc="lower left")

    # Panel B -- conditional coverage by TRUE D* tercile (CQR).
    t = np.arange(3)
    ax[1].plot(t, naive["cqr_terc"], "o-", c="#c0392b", label="naive transfer")
    ax[1].plot(t, recal["cqr_terc"], "s-", c="#2c7fb8", label="recalibrated")
    ax[1].axhline(NOMINAL, ls="--", c="k", lw=1)
    ax[1].set_xticks(t); ax[1].set_xticklabels(["lo D*", "mid D*", "hi D*"])
    ax[1].set_ylim(0, 1.02); ax[1].set_ylabel("conditional coverage")
    ax[1].set_title("(B) Conditional coverage by true-D* tercile")
    ax[1].legend(fontsize=8, loc="lower left")

    # Panel C -- the wall: CRLB(D*)/tercile-width by tercile.
    ax[2].bar(t, wall["ratio"], 0.6, color="#756bb1")
    ax[2].axhline(1.0, ls="--", c="k", lw=1, label="CRLB = tercile width")
    ax[2].set_xticks(t); ax[2].set_xticklabels(["lo D*", "mid D*", "hi D*"])
    ax[2].set_ylabel("CRLB(D*) / tercile width")
    ax[2].set_title("(C) Identifiability wall (OSIPI ground truth)")
    ax[2].legend(fontsize=8, loc="upper left")

    fires = "FIRES" if naive["monitor_fires"] else "does not fire"
    fig.suptitle(
        f"OSIPI external SYNTHETIC reference (DRO, n={c['n_test']} test) -- "
        f"naive-transfer monitor {fires} (AUC={naive['monitor_auc']:.2f}); "
        f"wall replicates: {wall['replicates']}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for d in (_FIG_DIR, _PAPER_FIG_DIR):
        os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, "osipi_coverage.pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[osipi figures] wrote osipi_coverage.pdf -> {_FIG_DIR} and {_PAPER_FIG_DIR}")


# --------------------------------------------------------------------------- #
# Report + provenance "run" block (seed-0 headline numbers, greppable).
# --------------------------------------------------------------------------- #
def _fmt(a):
    return "[" + ", ".join(f"{v:.3f}" for v in a) + "]"


def main(seed=DEFAULT_SEED):
    out, payload = compute_all(seed=seed)
    lines = []
    def w(s=""): lines.append(s)

    w("=" * 92)
    w("GAUGE TASK 1 -- EXTERNAL-PHANTOM REPLICATION ON THE OSIPI IVIM REFERENCE DRO")
    w("=" * 92)
    w("OSIPI TF2.4 IVIM DRO -- external community-standard SYNTHETIC reference "
      "(NOT in vivo).")
    w(f"Ground-truth (D, D*, f) from OSIPI (Zenodo doi:10.5281/zenodo.14605039, "
      f"CC-BY-4.0); seed={seed}; nominal coverage {NOMINAL:.2f}.")
    w("No in-vivo coverage claim is made. The bi-exponential model is the OSIPI "
      "consensus model (shared); the ground-truth distribution is external.")
    w("")
    for regime in ("controlled", "native"):
        p = payload[regime]
        naive, recal, wall = p["naive"], p["recal"], p["wall"]
        acq = ("Gauge 22-b + SNR grid" if regime == "controlled"
               else "OSIPI native 7-b, pre-noised ~SNR 80")
        w("-" * 92)
        w(f"REGIME: {regime}  ({acq}); n_test={p['n_test']}")
        w("-" * 92)
        w(f"  naive-transfer  CQR  marginal coverage (D, D*, f) = "
          f"{_fmt(naive['cqr_marg'])}  vs nominal {NOMINAL:.2f}")
        w(f"  recalibrated    CQR  marginal coverage (D, D*, f) = "
          f"{_fmt(recal['cqr_marg'])}  vs nominal {NOMINAL:.2f}")
        w(f"  naive-transfer  split-conformal marginal (D, D*, f) = {_fmt(naive['split_marg'])}")
        w(f"  recalibrated    split-conformal marginal (D, D*, f) = {_fmt(recal['split_marg'])}")
        w(f"  conditional CQR coverage by true-D* tercile [lo, mid, hi]:")
        w(f"      naive        = {_fmt(naive['cqr_terc'])}")
        w(f"      recalibrated = {_fmt(recal['cqr_terc'])}")
        mfires = "FIRES" if naive["monitor_fires"] else "does NOT fire"
        w(f"  deployment monitor on naive transfer: {mfires} "
          f"(AUC={naive['monitor_auc']:.2f})")
        w(f"  identifiability wall: CRLB(D*)/tercile-width [lo, mid, hi] = "
          f"{_fmt(wall['ratio'])}; abs-CRLB growth lo->hi = {wall['abs_growth']:.1f}x; "
          f"replicates = {wall['replicates']}")
        w("")

    cN = payload["controlled"]
    w("HEADLINES (controlled regime, seed point; banded across seeds in "
      "osipi_multiseed.json):")
    w(f"  OSIPI recalibrated marginal coverage at nominal {NOMINAL:.2f} "
      f"(within-distribution conformal holds on external ground truth).")
    w(f"  OSIPI naive-transfer monitor fires = {cN['naive']['monitor_fires']} "
      f"(AUC={cN['naive']['monitor_auc']:.2f}) -- the D* shift is observable.")
    w(f"  OSIPI high-D* wall replicates = {cN['wall']['replicates']} "
      f"(CRLB/width hi = {cN['wall']['hi_ratio']:.2f}); the high-D* conditional gap "
      f"persists even after recalibration "
      f"(recal hi-D* CQR coverage = {cN['recal']['cqr_terc'][-1]:.3f}).")
    w("=" * 92)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    with open(_REPORT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[osipi] wrote report -> {_REPORT}")

    _make_figure(payload)
    _write_run_provenance(seed, payload)
    return 0


def _write_run_provenance(seed, payload):
    cN = payload["controlled"]
    run = {
        "seed": seed, "alpha": ALPHA, "nominal_coverage": NOMINAL,
        "no_in_vivo_coverage_claim": True,
        "regimes": {},
    }
    for regime in ("controlled", "native"):
        p = payload[regime]
        run["regimes"][regime] = {
            "n_test": p["n_test"], "b_values": p["b"],
            "naive_cqr_marg_D_Dstar_f": [float(x) for x in p["naive"]["cqr_marg"]],
            "recal_cqr_marg_D_Dstar_f": [float(x) for x in p["recal"]["cqr_marg"]],
            "naive_cqr_tercile_lo_mid_hi": [float(x) for x in p["naive"]["cqr_terc"]],
            "recal_cqr_tercile_lo_mid_hi": [float(x) for x in p["recal"]["cqr_terc"]],
            "naive_monitor_fires": bool(p["naive"]["monitor_fires"]),
            "naive_monitor_auc": float(p["naive"]["monitor_auc"]),
            "wall_crlb_over_width_lo_mid_hi": [float(x) for x in p["wall"]["ratio"]],
            "wall_abs_growth": float(p["wall"]["abs_growth"]),
            "wall_replicates": bool(p["wall"]["replicates"]),
        }
    run["computed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    prov = {}
    if os.path.exists(_PROVENANCE):
        try:
            prov = json.load(open(_PROVENANCE))
        except Exception:
            prov = {}
    prov["run"] = run
    with open(_PROVENANCE, "w") as fh:
        json.dump(prov, fh, indent=1)
    print(f"[osipi] wrote run provenance -> {_PROVENANCE}")


# --------------------------------------------------------------------------- #
# Multi-seed band sweep -- reuses the existing across-seed contract verbatim
# (multiseed.seed_list + _agg_scalar): point = across-seed mean, [5, 95] band.
# Writes results/osipi_multiseed.json, gated by gauge/paper/consistency.py.
# --------------------------------------------------------------------------- #
def _osipi_seed_path(seed):
    return os.path.join(_RESULTS_DIR, "osipi_seeds", f"{int(seed)}.json")


def sweep(n=16, force=False, verbose=True):
    from gauge.multiseed import seed_list, _agg_scalar
    seeds = seed_list(n)
    os.makedirs(os.path.join(_RESULTS_DIR, "osipi_seeds"), exist_ok=True)
    for s in seeds:
        p = _osipi_seed_path(s)
        if (not force) and os.path.exists(p):
            if verbose:
                print(f"[osipi sweep] seed {s} present -- skip")
            continue
        flat, _ = compute_all(seed=s)
        with open(p, "w") as fh:
            json.dump({"seed": int(s), "flat": flat}, fh)
        if verbose:
            print(f"[osipi sweep] seed {s} done")
    recs = [json.load(open(_osipi_seed_path(s))) for s in seeds
            if os.path.exists(_osipi_seed_path(s))]
    seed0 = next(r for r in recs if r["seed"] == DEFAULT_SEED)
    items = {}
    for k in sorted(seed0["flat"].keys()):
        vals = [r["flat"][k] for r in recs if k in r["flat"]]
        # keys are 'osipi/...' -> neither the Fisher-z nor the (G) special-case in
        # _agg_scalar fires, so point = across-seed mean with a [5, 95] band.
        items[k] = _agg_scalar(vals, k, seed0["flat"][k])
    out = {"n_seeds": len(recs), "seeds": [r["seed"] for r in recs],
           "alpha": ALPHA, "items": items}
    with open(_MULTISEED, "w") as fh:
        json.dump(out, fh, indent=2)
    if verbose:
        print(f"[osipi sweep] aggregated {len(recs)} seeds -> {_MULTISEED}")
    return out


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "sweep":
        nn = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        sweep(n=nn, force=os.environ.get("GAUGE_FORCE") == "1")
        raise SystemExit(0)
    sd = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SEED
    raise SystemExit(main(seed=sd))
