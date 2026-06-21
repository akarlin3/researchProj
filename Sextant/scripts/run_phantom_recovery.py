#!/usr/bin/env python
"""HC2 item (a): known-truth parameter-recovery test for boundary-railing.

The boundary-railing claim needs no ground truth.  Where ground truth *is*
available -- a digital reference object with known (D, D*, f) -- we can go
further and ask whether railing is a faithful flag for genuine D*
unrecoverability.  This is the strongest closure the railing-first paper can make
on the in-vivo parameter-truth gap short of a real phantom scan, and it is
explicit that the substrate is a digital phantom (synthetic signal), not real
tissue.

Three known-truth analyses (all fit with Fashion's exact bounded NLLS):

  B1  Brain digital phantom (Ryghog 2014 / Federau 2012 reference values, committed
      at phantoms/brain/ground_truth/diffusive_groundtruth.json): per-tissue railing
      rate and D* recovery error at the phantom's own 17-point b-scheme.

  B2  f-controlled identifiability sweep: fix (D, D*), vary the perfusion fraction
      f; show the railed fraction and the relative D* recovery error rise together
      as f falls -- railing tracks unrecoverability, validated against known truth.

  B3  Railing-as-flag validity: over the trained-NPE prior, compare D* recovery
      error for railed vs non-railed voxels and score "railed" as a classifier of
      large recovery error (AUC, Spearman) -- a quantitative known-truth check that
      railing flags the voxels whose D* the fit cannot recover.

Run:  python Sextant/scripts/run_phantom_recovery.py
"""
from __future__ import annotations

import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "sextant-core"))

import numpy as np  # noqa: E402

from sextant.truthsim import (TARGET_BVALS, add_noise, draw_truths,  # noqa: E402
                              fit_and_rail, fm_biexp)

SEED = 20260621
RESULTS = os.path.join(_ROOT, "results", "phantom_recovery.json")
_PHANTOM_GT = os.path.join(
    _ROOT, "..", "Fashion", "phantoms", "brain", "ground_truth",
    "diffusive_groundtruth.json")
_PHANTOM_BVAL = os.path.join(
    _ROOT, "..", "Fashion", "phantoms", "brain", "ground_truth", "diffusive.bval")


def _render_truth(D, Ds, f, bvals, snr, m, rng):
    clean = fm_biexp(bvals[None, :], D, Ds, f)
    clean = np.repeat(clean, m, axis=0)
    return add_noise(clean, snr, rng, model="rician")


def auc_binary(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUC of `scores` predicting boolean `labels` (Mann-Whitney; ties at 0.5)."""
    pos = scores[labels]
    neg = scores[~labels]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), float)
    ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    csum = np.cumsum(counts)
    avg = {i: (csum[i] - (counts[i] - 1) / 2.0) for i in range(len(counts))}
    ranks = np.array([avg[i] for i in inv])
    r_pos = ranks[labels].sum()
    return float((r_pos - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg)))


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else float("nan")


# --------------------------------------------------------------------------- #
def b1_brain_phantom(rng):
    gt = json.load(open(_PHANTOM_GT))
    bvals = np.loadtxt(_PHANTOM_BVAL).astype(float)
    tissues = gt["tissues"]
    M = 500
    out = []
    for i, tis in enumerate(tissues):
        D, f, Ds = gt["D"][i], gt["f"][i], gt["Dstar"][i]
        # CSF has f=0, D*=0 (no perfusion); clamp D* into the fit box for rendering
        Ds_render = max(Ds, 3e-3)
        sig = _render_truth(D, Ds_render, f, bvals, 40.0, M,
                            np.random.default_rng(SEED + i))
        fit = fit_and_rail(sig, bvals)
        # recovery error only meaningful where a perfusion compartment exists
        rec_err = (np.nan if f == 0 else
                   float(np.median(np.abs(fit.dstar - Ds) / Ds)))
        out.append({
            "tissue": tis, "D_true": D, "Dstar_true": Ds, "f_true": f,
            "snr": 40.0, "n": M,
            "frac_railed": float(fit.railed.mean()),
            "frac_upper": float(fit.railed_hi.mean()),
            "frac_lower": float(fit.railed_lo.mean()),
            "median_rel_dstar_error": rec_err,
        })
        print(f"  B1 {tis:3s}  D*={Ds*1e3:.0f}e-3 f={f:.3f}  railed={out[-1]['frac_railed']:.3f} "
              f"(up={out[-1]['frac_upper']:.3f})  relErr(D*)="
              f"{'n/a' if f==0 else f'{rec_err:.2f}'}")
    return {"bvals": bvals.tolist(), "source": "Ryghog2014/Federau2012 (committed)",
            "tissues": out}


def b2_f_sweep(rng):
    D, Ds = 1.0e-3, 20e-3            # typical abdominal perfusion truth
    f_grid = [0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.45]
    M = 600
    rows = []
    for j, f in enumerate(f_grid):
        sig = _render_truth(D, Ds, f, TARGET_BVALS, 20.0, M,
                            np.random.default_rng(SEED + 100 + j))
        fit = fit_and_rail(sig, TARGET_BVALS)
        rel = np.abs(fit.dstar - Ds) / Ds
        rows.append({"f": f, "n": M, "snr": 20.0,
                     "frac_railed": float(fit.railed.mean()),
                     "median_rel_dstar_error": float(np.median(rel))})
        print(f"  B2 f={f:.2f}  railed={rows[-1]['frac_railed']:.3f}  "
              f"relErr(D*)={rows[-1]['median_rel_dstar_error']:.2f}")
    return {"D_true": D, "Dstar_true": Ds, "snr": 20.0, "rows": rows}


def b3_flag_validity(rng):
    out = {}
    for snr in (20.0, 40.0):
        N = 3000
        truths = draw_truths(N, np.random.default_rng(SEED + 200 + int(snr)))
        D, Ds, f = truths[:, 0], truths[:, 1], truths[:, 2]
        clean = fm_biexp(TARGET_BVALS[None, :], D[:, None], Ds[:, None], f[:, None])
        sig = add_noise(clean, snr, np.random.default_rng(SEED + 300 + int(snr)))
        fit = fit_and_rail(sig, TARGET_BVALS)
        err = np.abs(fit.dstar - Ds)              # absolute D* recovery error
        rel = err / Ds
        railed = fit.railed
        large = err > np.median(err)              # "large recovery error" label
        out[str(int(snr))] = {
            "n": N,
            "median_abs_err_railed": float(np.median(err[railed])) if railed.any() else None,
            "median_abs_err_nonrailed": float(np.median(err[~railed])) if (~railed).any() else None,
            "median_rel_err_railed": float(np.median(rel[railed])) if railed.any() else None,
            "median_rel_err_nonrailed": float(np.median(rel[~railed])) if (~railed).any() else None,
            "frac_railed": float(railed.mean()),
            "auc_railed_predicts_large_error": auc_binary(railed.astype(float), large),
            "spearman_railed_vs_abs_err": spearman(railed.astype(float), err),
        }
        o = out[str(int(snr))]
        print(f"  B3 SNR{int(snr)}  err(railed)={o['median_abs_err_railed']:.4f} vs "
              f"err(non)={o['median_abs_err_nonrailed']:.4f}  "
              f"AUC={o['auc_railed_predicts_large_error']:.3f}")
    return out


def main():
    rng = np.random.default_rng(SEED)
    print("[phantom] B1 brain digital phantom (known truth):")
    b1 = b1_brain_phantom(rng)
    print("[phantom] B2 f-controlled identifiability sweep (known truth):")
    b2 = b2_f_sweep(rng)
    print("[phantom] B3 railing-as-flag validity over the prior (known truth):")
    b3 = b3_flag_validity(rng)

    out = {
        "meta": {"seed": SEED,
                 "fit": "Fashion fit_biexp_nlls (read-only reuse)",
                 "rail_thresholds": {"lower": 0.0033, "upper": 0.1485},
                 "substrate": "DIGITAL reference phantom + truth-controlled sim "
                              "(synthetic signal; NOT real tissue)"},
        "B1_brain_phantom": b1,
        "B2_f_sweep": b2,
        "B3_flag_validity": b3,
    }
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[phantom] wrote {RESULTS}")


if __name__ == "__main__":
    main()
