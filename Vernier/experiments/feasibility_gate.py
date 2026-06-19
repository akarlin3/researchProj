"""CP2 feasibility gate -- the hard halt that decides whether Vernier exists.

Runs the pre-registered Experiment A: among matched-scan-time candidate schemes,
select a CRLB(D*)-matched subset (within +/-10%), then ask whether they diverge in
post-conformal D* calibration (sharpness / high-D* conditional coverage) by more
than the pre-set thresholds, with a bootstrap CI excluding zero.

Deterministic (seeded). Publication-independent (Caliper only). Prints the report
and writes it to ``Vernier/results/feasibility_gate.txt`` plus a JSON summary.

    python experiments/feasibility_gate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Make the in-tree ``vernier`` package importable when run as a script.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from vernier import _paths
from vernier import crlb
from vernier.feasibility import format_gate, run_gate, select_matched_crlb
from vernier.schemes import CANDIDATE_POOL

_paths.add_caliper()
from caliper.forward import sample_params, synthetic_cohort  # noqa: E402

PRIMARY_SNR = 33.0
ROBUSTNESS_SNR = (25.0, 50.0)
N = 8000
SEED = 0
N_BOOT = 2000
RESULTS = Path(__file__).resolve().parent.parent / "results"


def _verify_pairing(snr: float, seed: int) -> bool:
    """Confirm matched-n_b schemes share identical params AND noise (paired)."""
    from vernier.schemes import MATCHED_SCANTIME
    a = synthetic_cohort(n=200, bvalues=MATCHED_SCANTIME[0].b, snr=snr, noise="rician", seed=seed)
    b = synthetic_cohort(n=200, bvalues=MATCHED_SCANTIME[1].b, snr=snr, noise="rician", seed=seed)
    return bool(np.allclose(a.params, b.params))


def main() -> None:
    from vernier.schemes import CANDIDATE_POOL as POOL

    rng = np.random.default_rng(SEED)
    params = sample_params(N, rng)

    paired = _verify_pairing(PRIMARY_SNR, SEED)
    print(f"[paired-design check] matched-n_b schemes share true params: {paired}")
    if not paired:
        raise SystemExit("paired-design assumption violated -- aborting")

    # Pre-registered selection: CRLB(D*) within +/-10% of the candidate-pool median.
    cr = {s.name: float(crlb.expected_crlb(s, params, PRIMARY_SNR)[2]) for s in POOL}
    print("\nCandidate-pool CRLB(D*) @ SNR={:g}:".format(PRIMARY_SNR))
    for k, v in cr.items():
        print(f"  {k:>20} {v:7.3f}")
    matched = select_matched_crlb(POOL, params, PRIMARY_SNR, tol=0.10, min_keep=3)
    print("\nMatched-CRLB(D*) schemes selected (within +/-10% of median):")
    for s in matched:
        print(f"  {s.name:>20}  CRLB(D*)={cr[s.name]:.3f}  n_b={s.n_b}  "
              f"low-b<=50:{int(np.sum((s.b<=50)&(s.b>0)))}  high-b>=200:{int(np.sum(s.b>=200))}")

    reports = []
    summary = {"primary_snr": PRIMARY_SNR, "n": N, "seed": SEED, "n_boot": N_BOOT,
               "matched_schemes": [s.name for s in matched],
               "crlb_dstar": {s.name: cr[s.name] for s in matched}, "runs": {}}

    for snr in (PRIMARY_SNR, *ROBUSTNESS_SNR):
        gr = run_gate(matched, n=N, snr=snr, seed=SEED, n_boot=N_BOOT)
        rep = format_gate(gr)
        print("\n" + rep)
        reports.append(rep)
        summary["runs"][f"snr_{snr:g}"] = {
            "verdict": gr.verdict,
            "delta_sharp": gr.delta_sharp, "delta_sharp_ci": gr.delta_sharp_ci,
            "delta_cond": gr.delta_cond, "delta_cond_ci": gr.delta_cond_ci,
            "schemes": [
                {"name": s.name, "crlb_dstar": s.crlb_dstar,
                 "cov_dstar_raw": s.cov_dstar_raw, "cov_dstar": s.cov_dstar,
                 "width_dstar": s.width_dstar, "ece_dstar": s.ece_dstar,
                 "cond_cov_high_dstar": s.cond_cov_high_dstar}
                for s in gr.schemes
            ],
        }

    primary = summary["runs"][f"snr_{PRIMARY_SNR:g}"]["verdict"]
    summary["primary_verdict"] = primary

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "feasibility_gate.txt").write_text(
        f"[paired-design check] matched-n_b schemes share true params: {paired}\n\n"
        + "\n\n".join(reports) + "\n"
    )
    (RESULTS / "feasibility_gate.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n>>> PRIMARY VERDICT (SNR={PRIMARY_SNR:g}): {primary}")
    print(f">>> wrote {RESULTS/'feasibility_gate.txt'} and feasibility_gate.json")


if __name__ == "__main__":
    main()
