"""CP3 follow-up: does the post-conformal calibration divergence persist for an
*efficient* estimator (the MAF posterior), not just the segmented reference?

The CP2 gate measured the divergence on Caliper's deliberately over-confident,
homoscedastic segmented estimator; the manuscript flagged the magnitude as
estimator-specific and named this as the natural next experiment. Here we re-run the
same pre-registered gate (matched scan-time, matched CRLB(D*), split-conformal, the
same Delta_sharp / Delta_cond metrics and paired double-bootstrap;
:func:`vernier.gate.run_gate_general`) on the torch MAF posterior, and on the
reference estimator under *identical* train/cal/test splits, so the only thing that
changes is the estimator.

Deterministic (seeded). Writes results/maf_gate.{txt,json}.

    KMP_DUPLICATE_LIB_OK=TRUE python experiments/maf_gate.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")  # torch/numpy OpenMP coexistence
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from vernier import _paths
from vernier.feasibility import THRESH_COND, THRESH_SHARP, select_matched_crlb
from vernier.gate import run_gate_general
from vernier.schemes import CANDIDATE_POOL

_paths.add_caliper()
from caliper.estimator_reference import ReferenceIVIMEstimator  # noqa: E402
from caliper.forward import sample_params  # noqa: E402

SNR = 33.0
N = 6000
SEED = 0
N_BOOT = 1000
RESULTS = Path(__file__).resolve().parent.parent / "results"


def fmt(res):
    L = [f"=== {res['label']} (SNR={res['snr']:g}, n={res['n']}, "
         f"train/cal/test={res['splits']['train']}/{res['splits']['cal']}/{res['splits']['test']}, "
         f"bootstrap={res['n_boot']}) ==="]
    hdr = f"{'scheme':>16} {'CRLB_D*':>8} {'rawCov':>7} {'cov':>6} {'width':>9} {'ECE':>6} {'hiCov':>7}"
    L += [hdr, "-" * len(hdr)]
    for s in res["schemes"]:
        L.append(f"{s['name']:>16} {s['crlb_dstar']:>8.2f} {s['cov_dstar_raw']:>7.3f} "
                 f"{s['cov_dstar']:>6.3f} {s['width_dstar']:>9.3f} {s['ece_dstar']:>6.3f} "
                 f"{s['cond_cov_high_dstar']:>7.3f}")
    cs, cc = res["delta_sharp_ci"], res["delta_cond_ci"]
    L += ["",
          f"  Delta_sharp = {res['delta_sharp']:.4f}  CI [{cs[0]:.4f}, {cs[1]:.4f}]  (thresh {THRESH_SHARP})",
          f"  Delta_cond  = {res['delta_cond']:.4f}  CI [{cc[0]:.4f}, {cc[1]:.4f}]  (thresh {THRESH_COND})",
          f"  VERDICT: {res['verdict']}"]
    return "\n".join(L)


def main():
    from caliper.estimator_maf import MAFPosterior

    params = sample_params(N, np.random.default_rng(SEED))
    matched = select_matched_crlb(CANDIDATE_POOL, params, SNR, tol=0.10, min_keep=3)
    print("matched-CRLB(D*) schemes:", [s.name for s in matched], "\n")

    ref = run_gate_general(lambda s: ReferenceIVIMEstimator(bvalues=s.b), matched,
                           label="REFERENCE (segmented) -- cross-check under MAF splits",
                           n=N, snr=SNR, seed=SEED, n_boot=N_BOOT)
    print(fmt(ref), "\n")
    maf = run_gate_general(lambda s: MAFPosterior(n_bvalues=s.n_b, seed=SEED), matched,
                           label="MAF posterior (efficient estimator)",
                           n=N, snr=SNR, seed=SEED, n_boot=N_BOOT)
    print(fmt(maf), "\n")

    persists = maf["verdict"] == "PASS"
    print(f">>> Reference verdict: {ref['verdict']}  |  MAF verdict: {maf['verdict']}")
    print(f">>> Divergence persists for the efficient MAF estimator: {persists}")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "maf_gate.txt").write_text(fmt(ref) + "\n\n" + fmt(maf) +
                                          f"\n\nPersists for MAF: {persists}\n")
    (RESULTS / "maf_gate.json").write_text(json.dumps(
        {"reference": ref, "maf": maf, "persists_for_maf": persists}, indent=2) + "\n")
    print(f"\n>>> wrote {RESULTS/'maf_gate.txt'} and maf_gate.json")


if __name__ == "__main__":
    main()
