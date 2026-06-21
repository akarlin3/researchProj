#!/usr/bin/env python
"""HC2 item (b): in-silico misspecification isolation for boundary-railing.

Question answered: is the NLLS D* boundary-railing we report on real abdominal
data an intrinsic pathology of the bounded estimator under weak D*
identifiability, or could it be an artefact of *simulator-reality mismatch*
(the forward / noise model not matching tissue)?

Design (all truth-controlled; ground truth known at every voxel):
  * Draw N truths over the trained NPE BoxUniform prior.
  * Render at SNR in {10,20,40,100} under one well-specified forward model
    (bi-exponential + Rician -- exactly the model the NLLS assumes) and three
    deliberate misspecifications (tri-exponential tissue tail; diffusional
    kurtosis; non-Rician 4-coil chi noise).
  * Fit every voxel with Fashion's exact bounded NLLS and flag rail at Fashion's
    exact thresholds (reused read-only via sextant.truthsim / fashion_reuse).
  * Report the railed fraction per (forward model x SNR), bootstrap CI, and the
    stratification by true-f and true-D* tercile.

Reading:
  * Well-specified railing > 0 with ZERO misspecification isolates the railing
    as estimator/identifiability pathology, not forward-model misfit.
  * If misspecification does not materially inflate the railed fraction, the
    real-data railing cannot be attributed to simulator-reality mismatch.
  * Honesty gate: whatever the numbers are, they are written from this run.

Run:  python Sextant/scripts/run_misspecification_isolation.py
"""
from __future__ import annotations

import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "sextant-core"))

import numpy as np  # noqa: E402

from sextant.truthsim import (TARGET_BVALS, draw_truths, fit_and_rail,  # noqa: E402
                              render)

SEED = 20260621
N = 2000
SNRS = [10.0, 20.0, 40.0, 100.0]
FORWARDS = ["biexp_WS", "triexp", "kurtosis", "noise_chi"]
N_BOOT = 2000
RESULTS = os.path.join(_ROOT, "results", "misspecification_isolation.json")


def boot_ci(mask: np.ndarray, n_boot: int, rng: np.random.Generator):
    n = len(mask)
    idx = rng.integers(0, n, size=(n_boot, n))
    fracs = mask[idx].mean(axis=1)
    return float(np.percentile(fracs, 2.5)), float(np.percentile(fracs, 97.5))


def tercile_rates(railed: np.ndarray, value: np.ndarray):
    """Railed fraction within low/mid/high terciles of `value` (true f or D*)."""
    q1, q2 = np.percentile(value, [100 / 3, 200 / 3])
    lo = value <= q1
    hi = value > q2
    mid = ~lo & ~hi
    return {
        "low": float(railed[lo].mean()) if lo.any() else None,
        "mid": float(railed[mid].mean()) if mid.any() else None,
        "high": float(railed[hi].mean()) if hi.any() else None,
    }


def main():
    rng = np.random.default_rng(SEED)
    truths = draw_truths(N, rng)            # shared truths across all cells
    f_true, ds_true = truths[:, 2], truths[:, 1]

    # The curated real abdominal high-SNR ROI is concentrated in the hard
    # low-perfusion / high-pseudodiffusion corner; flag those truths so the
    # in-silico railing magnitude can be read in the same regime as the real data.
    corner = (f_true < 0.10) & (ds_true > 75e-3)

    cells = []
    print(f"[misspec] N={N} truths, b={TARGET_BVALS.tolist()}, seed={SEED}")
    print(f"[misspec] hard-corner (f<0.10 & D*>75e-3) holds {int(corner.sum())} truths")
    for forward in FORWARDS:
        for snr in SNRS:
            # independent noise stream per cell (deterministic, seeded)
            cell_rng = np.random.default_rng(
                SEED + 1000 * FORWARDS.index(forward) + int(snr))
            sig = render(truths, TARGET_BVALS, snr, cell_rng, forward=forward)
            fit = fit_and_rail(sig, TARGET_BVALS)
            ci = boot_ci(fit.railed, N_BOOT, np.random.default_rng(SEED + 7))
            cell = {
                "forward": forward, "snr": snr,
                "frac_railed": float(fit.railed.mean()),
                "frac_lower": float(fit.railed_lo.mean()),
                "frac_upper": float(fit.railed_hi.mean()),
                "ci": [ci[0], ci[1]],
                "frac_railed_hard_corner": (
                    float(fit.railed[corner].mean()) if corner.any() else None),
                "by_true_f_tercile": tercile_rates(fit.railed, f_true),
                "by_true_dstar_tercile": tercile_rates(fit.railed, ds_true),
            }
            cells.append(cell)
            corner_str = (f"{cell['frac_railed_hard_corner']:.3f}"
                          if cell['frac_railed_hard_corner'] is not None else "n/a")
            print(f"  {forward:10s} SNR{int(snr):>3}  railed={cell['frac_railed']:.3f} "
                  f"[{ci[0]:.3f},{ci[1]:.3f}]  lo/hi={cell['frac_lower']:.3f}/"
                  f"{cell['frac_upper']:.3f}  corner={corner_str}  "
                  f"f(lo/hi)={cell['by_true_f_tercile']['low']:.2f}/"
                  f"{cell['by_true_f_tercile']['high']:.2f}")

    # Misspecification deltas: misspec railed minus well-specified, per SNR.
    ws = {c["snr"]: c["frac_railed"] for c in cells if c["forward"] == "biexp_WS"}
    deltas = {}
    for forward in FORWARDS:
        if forward == "biexp_WS":
            continue
        deltas[forward] = {
            str(int(c["snr"])): round(c["frac_railed"] - ws[c["snr"]], 4)
            for c in cells if c["forward"] == forward}
    max_abs_delta = max(abs(v) for d in deltas.values() for v in d.values())

    out = {
        "meta": {"seed": SEED, "n_truths": N, "n_boot": N_BOOT,
                 "bvals": TARGET_BVALS.tolist(), "snrs": SNRS,
                 "forwards": FORWARDS,
                 "prior": "trained NPE BoxUniform (D U, D* logU, f U)",
                 "fit": "Fashion fit_biexp_nlls (read-only reuse)",
                 "rail_thresholds": {"lower": 0.0033, "upper": 0.1485}},
        "cells": cells,
        "misspec_delta_vs_WS": deltas,
        "max_abs_delta_vs_WS": round(max_abs_delta, 4),
    }
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n[misspec] well-specified railing (zero misspecification) by SNR: "
          f"{ {int(k): round(v,3) for k,v in ws.items()} }")
    print(f"[misspec] max |railed(misspec) - railed(WS)| across all cells = {max_abs_delta:.3f}")
    print(f"[misspec] wrote {RESULTS}")


if __name__ == "__main__":
    main()
