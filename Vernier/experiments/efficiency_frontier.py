"""CP3 Experiment B -- the decision-value-per-scan-minute efficiency frontier.

PROVISIONAL (Minos lens). Runs *after* the feasibility gate has PASSED. Sweeps
schemes of different scan-time (7/11/15/22 b-values) and scores each scheme's
conformal-corrected D* posterior on a fixed treat/spare/escalate decision (Minos
utility), reporting decision-value-per-scan-minute -- the protocol-recommendation
deliverable. The decision config (thresholds, cost asymmetry) is fixed once from
the prior so the comparison is about the *protocol*, not the decision.

Deterministic (seeded). Writes results/efficiency_frontier.{txt,json}.

    python experiments/efficiency_frontier.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from vernier import _paths
from vernier.decision import decision_value, format_frontier, make_decision_config
from vernier.schemes import EFFICIENCY_FRONTIER

_paths.add_caliper()
from caliper.forward import sample_params  # noqa: E402

SNR = 33.0
N = 8000
SEED = 0
RESULTS = Path(__file__).resolve().parent.parent / "results"


def main() -> None:
    # Fixed decision config from the prior D* terciles (same for every scheme).
    params = sample_params(N, np.random.default_rng(SEED))
    cfg = make_decision_config(params[:, 2])

    results = [decision_value(s, n=N, snr=SNR, seed=SEED, cfg=cfg) for s in EFFICIENCY_FRONTIER]
    report = format_frontier(results, cfg)
    print(report)

    summary = {
        "snr": SNR, "n": N, "seed": SEED, "provisional": True, "lens": "minos",
        "decision_config": {"t1": cfg.t1, "t2": cfg.t2, "k_under": cfg.k_under, "k_over": cfg.k_over},
        "no_scan_prior_baseline_utility": results[0].baseline_utility,
        "finding": ("decision utility improves monotonically and saturates with scan-time; "
                    "no protocol beats the no-scan prior for D* (identifiability wall)"),
        "schemes": [
            {"name": r.name, "n_b": r.n_b, "scan_minutes": r.scan_minutes,
             "width_dstar": r.width_dstar, "cond_cov_high_dstar": r.cond_cov_high_dstar,
             "mean_utility": r.mean_utility, "baseline_utility": r.baseline_utility}
            for r in results
        ],
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "efficiency_frontier.txt").write_text(report + "\n")
    (RESULTS / "efficiency_frontier.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n>>> wrote {RESULTS/'efficiency_frontier.txt'} and efficiency_frontier.json")


if __name__ == "__main__":
    main()
