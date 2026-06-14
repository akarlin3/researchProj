"""CP0 GATE 0 -- forward-model sanity gate.

Recovering parameters from CLEAN (noise-free) synthetic signals must return the
known truth. If it does not, the generator/estimator is broken and CP1 must not
proceed. Prints per-parameter relative-error statistics over a parameter grid
and a PASS/FAIL verdict.

Run:  python scripts/sanity_forward.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from gauge.forward import ivim_signal, DEFAULT_B_VALUES
from gauge.estimators import fit_nlls
from gauge.cohort import D_RANGE, DSTAR_RANGE, F_RANGE

TOL = 1e-2  # max allowed relative error on clean data


def main():
    rng = np.random.default_rng(20260613)
    n = 400
    D = rng.uniform(*D_RANGE, n)
    Dstar = rng.uniform(*DSTAR_RANGE, n)
    f = rng.uniform(*F_RANGE, n)
    truth = np.stack([D, Dstar, f], axis=1)

    est = np.empty_like(truth)
    for i in range(n):
        s = ivim_signal(DEFAULT_B_VALUES, D[i], Dstar[i], f[i], S0=1.0)  # CLEAN
        e = fit_nlls(s, DEFAULT_B_VALUES)
        est[i] = (e["D"], e["Dstar"], e["f"])

    rel = np.abs(est - truth) / truth

    print("=" * 70)
    print("CP0 / GATE 0 -- forward-model sanity (clean, noise-free recovery)")
    print("=" * 70)
    print(f"b-values ({DEFAULT_B_VALUES.size}): "
          f"{DEFAULT_B_VALUES.astype(int).tolist()}")
    print(f"grid samples: {n}   seed: 20260613")
    print(f"ranges: D={D_RANGE} Dstar={DSTAR_RANGE} f={F_RANGE} (mm^2/s, -, -)")
    print("-" * 70)
    print(f"{'param':>6} | {'median rel err':>15} | {'p95 rel err':>12} | "
          f"{'max rel err':>12}")
    print("-" * 70)
    names = ("D", "D*", "f")
    overall_max = 0.0
    for j, nm in enumerate(names):
        med, p95, mx = (np.median(rel[:, j]), np.percentile(rel[:, j], 95),
                        rel[:, j].max())
        overall_max = max(overall_max, mx)
        print(f"{nm:>6} | {med:>15.2e} | {p95:>12.2e} | {mx:>12.2e}")
    print("-" * 70)
    verdict = "PASS" if overall_max < TOL else "FAIL"
    print(f"overall max relative error: {overall_max:.2e}   "
          f"tolerance: {TOL:.0e}   -> GATE 0 forward-model: {verdict}")
    print("=" * 70)
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
