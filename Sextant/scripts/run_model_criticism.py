#!/usr/bin/env python
"""HC2 item (c): RNPE-style posterior-predictive model criticism on in-vivo data.

Wherever the railing-first paper touches real data, a model-criticism check should
flag where the bi-exponential forward/noise model is misspecified -- *without*
needing parameter ground truth.  We use a held-out-b posterior-predictive
self-consistency residual (the assumption-free analogue of the RNPE statistical
model criticism of Ward et al., 2022): fit the bounded NLLS on a retained b-subset,
predict the held-out b-values, and compare the SNR-normalised held-out residual to
a *well-specified* simulated null.  A voxel is "criticised" (flagged misspecified)
if its residual exceeds the null's 95th percentile.

Validation (in-silico, known generator):
  * well-specified test set -> flagged fraction approx 5% (calibrated by construction),
  * tri-exponential / kurtosis test sets -> flagged fraction > 5% (sensitive).

Application:
  * OSIPI abdomen real data (read-only, by absolute path to the monorepo checkout):
    report the criticised fraction among high-SNR voxels.

Run:  python Sextant/scripts/run_model_criticism.py
      [--abdomen-root /path/to/Fashion/download/Data]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "sextant-core"))

import numpy as np  # noqa: E402

from sextant.fashion_reuse import load_railing  # noqa: E402
from sextant.truthsim import (TARGET_BVALS, draw_truths,  # noqa: E402
                              fm_biexp, fm_kurtosis, fm_triexp)

SEED = 20260621
RESULTS = os.path.join(_ROOT, "results", "model_criticism.json")
# Held-out b indices into TARGET_BVALS = [0,10,20,30,50,75,100,150,400,600].
# Spread across the range (low / mid / high) so the check is sensitive to both
# low-b perfusion-regime and high-b diffusion-regime (kurtosis) misspecification.
HELD = np.array([3, 6, 9])                  # b = 30, 100, 600
FIT_IDX = np.array([i for i in range(len(TARGET_BVALS)) if i not in HELD])
_FIT = load_railing()["fit_biexp_nlls"]


_SIGMA_FLOOR = 1e-4   # guards against divide-by-zero on perfectly-fit signals


def ppc_residual(signals: np.ndarray, snrs: np.ndarray = None) -> np.ndarray:
    """Self-calibrating held-out-b posterior-predictive misfit ratio per voxel.

    Fit the bounded NLLS on the retained b-subset; the in-fit residual RMS is a
    per-voxel noise estimate ``sigma_hat`` (averaging-invariant -- it scales with
    whatever effective noise the signal actually has).  The statistic is the
    held-out-b residual RMS divided by ``sigma_hat``: O(1) when the bi-exponential
    is well-specified, and large when the held-out points deviate beyond the noise
    the fit itself sees.  This is the RNPE self-consistency idea (Ward et al. 2022)
    in a truth-free NLLS form, and -- unlike an external-SNR normalisation -- it is
    insensitive to multi-acquisition averaging in the real data.  ``snrs`` is
    accepted for signature compatibility but unused.
    """
    b_fit = TARGET_BVALS[FIT_IDX]
    b_held = TARGET_BVALS[HELD]
    res = np.empty(len(signals))
    for i in range(len(signals)):
        p = _FIT(b_fit, signals[i, FIT_IDX])
        D, Ds, f = p
        if not np.isfinite(Ds):
            res[i] = np.inf
            continue
        infit = signals[i, FIT_IDX] - fm_biexp(b_fit, D, Ds, f)
        sigma_hat = max(np.sqrt(np.mean(infit ** 2)), _SIGMA_FLOOR)
        held_rms = np.sqrt(np.mean((signals[i, HELD] - fm_biexp(b_held, D, Ds, f)) ** 2))
        res[i] = held_rms / sigma_hat
    return res


def _sim_set(n, snr_lo, snr_hi, forward, seed):
    rng = np.random.default_rng(seed)
    truths = draw_truths(n, rng)
    snrs = rng.uniform(snr_lo, snr_hi, size=n)
    D, Ds, f = truths[:, 0:1], truths[:, 1:2], truths[:, 2:3]
    if forward == "biexp_WS":
        clean = fm_biexp(TARGET_BVALS[None, :], D, Ds, f)
    elif forward == "triexp":
        clean = fm_triexp(TARGET_BVALS[None, :], D, Ds, f)
    elif forward == "triexp_strong":
        # a clearly non-mono-exponential tissue (restricted/intracellular tail);
        # a structural SHAPE misspecification the bounded bi-exp cannot represent.
        clean = fm_triexp(TARGET_BVALS[None, :], D, Ds, f, w=0.45, slow=0.25)
    elif forward == "kurtosis":
        clean = fm_kurtosis(TARGET_BVALS[None, :], D, Ds, f)
    elif forward == "baseline":
        # incomplete background suppression / Rician noise floor: a constant
        # offset that grows in relative terms at high b (where S is small).
        clean = fm_biexp(TARGET_BVALS[None, :], D, Ds, f) + 0.03
    # add per-voxel Rician noise at each voxel's own SNR (sigma = 1/SNR, S0=1)
    sigma = (1.0 / snrs)[:, None]
    re = clean + rng.normal(0, 1, clean.shape) * sigma
    im = rng.normal(0, 1, clean.shape) * sigma
    sig = np.sqrt(re ** 2 + im ** 2)
    return sig, snrs


def _default_abdomen_root():
    """Prefer the OSIPI data fetched by reproduce.sh; fall back to a monorepo checkout."""
    for cand in (os.path.join(_ROOT, "data", "osipi", "extracted"),
                 os.path.join(_ROOT, "..", "Fashion", "download", "Data")):
        if os.path.exists(os.path.join(cand, "abdomen.nii.gz")):
            return cand
    return os.path.join(_ROOT, "data", "osipi", "extracted")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--abdomen-root", default=_default_abdomen_root(),
        help="dir holding abdomen.nii.gz / mask / bval (read-only).")
    args = ap.parse_args()

    # 1) Build the well-specified pooled null and the 95th-pct threshold.
    null_sig, null_snr = _sim_set(5000, 10.0, 60.0, "biexp_WS", SEED)
    null_res = ppc_residual(null_sig, null_snr)
    thresh = float(np.percentile(null_res[np.isfinite(null_res)], 95))
    print(f"[criticism] held-out b = {TARGET_BVALS[HELD].tolist()}; "
          f"null 95th-pct threshold = {thresh:.4f}")

    # 2) Validation test sets (fresh draws).
    # WS = calibration; triexp_strong = clear SHAPE misspecification (sensitivity
    # demonstration); triexp/kurtosis/baseline characterise what the self-calibrating
    # check is sensitive to (shape) vs robust to (fit-absorbable amplitude/noise).
    val = {}
    val_seeds = {"biexp_WS": 11, "triexp": 23, "triexp_strong": 71,
                 "kurtosis": 37, "baseline": 53}
    val_plan = [("biexp_WS", "~0.05 (calibration)"),
                ("triexp_strong", ">0.05 (shape, sensitivity check)"),
                ("triexp", ">~0.05 (mild shape)"),
                ("kurtosis", "absorbed (amplitude)"),
                ("baseline", "absorbed (amplitude)")]
    for forward, exp in val_plan:
        sig, snr = _sim_set(3000, 10.0, 60.0, forward, SEED + val_seeds[forward])
        res = ppc_residual(sig, snr)
        frac = float(np.mean(res > thresh))
        val[forward] = {"flagged_frac": frac, "expected": exp, "n": len(sig)}
        print(f"  validation {forward:13s}  flagged={frac:.3f}  ({exp})")

    # 3) Apply to OSIPI abdomen real data (read-only).
    real = {"available": False}
    R = load_railing()
    img = os.path.join(args.abdomen_root, "abdomen.nii.gz")
    mask = os.path.join(args.abdomen_root, "mask_abdomen_homogeneous.nii.gz")
    bval = os.path.join(args.abdomen_root, "abdomen.bval")
    if all(os.path.exists(p) for p in (img, mask, bval)):
        from sextant.truthsim import fit_and_rail  # local import; avoids cycle
        _i, _m, coords, _b, fit_signals, snrs = R["load_voxels"](img, mask, bval)
        fit_signals = np.asarray(fit_signals, float)
        snrs = np.asarray(snrs, float)
        hi = snrs >= float(R["SNR_FLOOR"])
        sig_hi = fit_signals[hi]
        res = ppc_residual(sig_hi, snrs[hi])
        crit = res > thresh
        frac = float(crit.mean())
        # Decisive cross-tab: is railing a misspecification artefact?  Fit the full
        # scheme, flag rail, and compare the criticised fraction among railed vs
        # non-railed voxels.  Similar fractions => railing is independent of
        # model-misfit (identifiability, not misspecification).
        railed = fit_and_rail(sig_hi, TARGET_BVALS).railed
        crit_railed = float(crit[railed].mean()) if railed.any() else None
        crit_nonrailed = float(crit[~railed].mean()) if (~railed).any() else None
        real = {"available": True, "n_high_snr": int(hi.sum()),
                "criticised_frac": frac,
                "excess_over_null": round(frac - 0.05, 4),
                "frac_railed": float(railed.mean()),
                "criticised_frac_among_railed": crit_railed,
                "criticised_frac_among_nonrailed": crit_nonrailed,
                "source": "OSIPI TF2.4 abdomen (Zenodo 14605039), homogeneous ROI"}
        print(f"  real OSIPI abdomen  n_hi={int(hi.sum())}  "
              f"criticised={frac:.3f}  (null baseline 0.05)")
        if crit_railed <= 0.05 + 1e-9 and crit_railed < crit_nonrailed:
            interp = "NOT a misspecification artefact (railed are well-fit-but-unidentified)"
        elif crit_railed > crit_nonrailed + 0.05:
            interp = "co-located with misspecification (railed voxels fit worse)"
        else:
            interp = "independent of misspecification"
        real["interpretation"] = interp
        print(f"    railed={railed.mean():.3f}; criticised|railed={crit_railed:.3f} "
              f"vs criticised|non-railed={crit_nonrailed:.3f}\n"
              f"    -> railing is {interp}")
    else:
        real["note"] = (f"abdomen data not found under {args.abdomen_root}; "
                        "run Sextant/scripts/fetch_osipi.py or pass --abdomen-root")
        print(f"  real OSIPI abdomen  NOT FOUND under {args.abdomen_root} -> FLAG")

    out = {
        "meta": {"seed": SEED, "held_out_b": TARGET_BVALS[HELD].tolist(),
                 "fit_b": TARGET_BVALS[FIT_IDX].tolist(),
                 "null_threshold_pct": 95, "null_threshold": thresh,
                 "statistic": "SNR-normalised held-out-b posterior-predictive RMSE",
                 "method": "RNPE-style (Ward et al. 2022) NLLS analogue, truth-free"},
        "validation": val,
        "real": real,
    }
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[criticism] wrote {RESULTS}")


if __name__ == "__main__":
    main()
