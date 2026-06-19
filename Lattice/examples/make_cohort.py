"""Generate one cohort per family and print a summary.

Run: python examples/make_cohort.py
"""

import numpy as np

import lattice


def main() -> None:
    print(f"Lattice {lattice.__version__}  |  b-values: {len(lattice.DEFAULT_BVALUES)} points")
    print(f"seed={lattice.DEFAULT_SEED}  param order={lattice.PARAM_NAMES}")
    print("-" * 72)
    for family in lattice.FAMILIES:
        c = lattice.make_cohort(family, n=1000, snr=40, seed=lattice.DEFAULT_SEED)
        D, Dstar, f = c.params.mean(axis=0)
        extra = dict(zip(c.extra_names, c.extra[0])) if c.extra_names else {}
        print(
            f"  {family:22s} n={len(c):4d}  signals={c.signals.shape}  "
            f"mean(D,D*,f)=({D:.2e},{Dstar:.2e},{f:.2f})  extra={extra}"
        )
    print("\nAll five families generated from one shared (D, D*, f) draw per seed.")


if __name__ == "__main__":
    main()
