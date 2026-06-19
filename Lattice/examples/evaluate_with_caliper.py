"""Score a UQ method on Lattice with Caliper -- the one-way consumption demo.

Lattice does not depend on Caliper; this *example* optionally imports it to show
the canonical coverage / ECE / sharpness scorecard. If Caliper is not importable
it degrades gracefully (the DRO and adapter still work without it).

Run: python examples/evaluate_with_caliper.py
"""

import sys
from pathlib import Path

import numpy as np

import lattice
from examples.evaluate_demo import ReferenceNLLSQuantiles  # reuse the estimator


def _try_import_caliper():
    try:
        from caliper import metrics  # type: ignore
        return metrics
    except ImportError:
        # Caliper is a sibling project, not a dependency. Try the monorepo path.
        here = Path(__file__).resolve()
        for parent in here.parents:
            cand = parent / "Caliper"
            if (cand / "caliper" / "metrics.py").exists():
                sys.path.insert(0, str(cand))
                try:
                    from caliper import metrics  # type: ignore
                    return metrics
                except ImportError:
                    return None
        return None


def main() -> None:
    metrics = _try_import_caliper()
    cohort = lattice.make_cohort("biexp", n=300, snr=40, seed=lattice.DEFAULT_SEED)
    est = ReferenceNLLSQuantiles(cohort.bvalues, sd=[0.11e-3, 6.0e-3, 0.022])
    q = est.predict_quantiles(cohort.signals, lattice.DEFAULT_QUANTILE_LEVELS)
    payload = lattice.to_scorer_inputs(cohort, q)

    if metrics is None:
        print("Caliper not available -- skipping the canonical scorecard.")
        print("The adapter produced scorer-ready inputs:")
        print(f"  y_true {payload['y_true'].shape}, q_pred {payload['q_pred'].shape}, "
              f"levels={list(payload['q_levels'])}, params={payload['param_names']}")
        print("Install/locate Caliper to score coverage / ECE / sharpness.")
        return

    scores = metrics.score_quantiles(**payload, alpha=0.10)
    print(metrics.format_scorecard(scores, title="Lattice biexp cohort scored by Caliper"))


if __name__ == "__main__":
    main()
