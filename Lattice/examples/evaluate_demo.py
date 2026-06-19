"""Score a simple UQ method on Lattice using only numpy/scipy (Caliper-free).

Defines a reference estimator that fits IVIM by NLLS and emits Gaussian quantile
intervals with fixed per-parameter SDs -- deliberately mis-scaled, so the raw
coverage departs from nominal (the calibration question Lattice exists to probe).
Scoring here uses Lattice's tiny dependency-free helpers; the canonical scorer
(coverage / ECE / sharpness) is Caliper -- see evaluate_with_caliper.py.

Run: python examples/evaluate_demo.py   (needs the 'selfcheck' extra: scipy)
"""

import numpy as np

import lattice
from lattice.evaluate import central_interval, interval_coverage, mean_sharpness


class ReferenceNLLSQuantiles:
    """Fit (D, D*, f) by NLLS; report Gaussian quantiles at fixed SDs."""

    def __init__(self, bvalues, sd):
        self.bvalues = np.asarray(bvalues, float)
        self.sd = np.asarray(sd, float)  # per-parameter (D, D*, f)

    def predict_quantiles(self, signals, q_levels):
        from scipy.stats import norm
        from lattice.selfcheck import fit_biexp_nlls
        point = np.array([fit_biexp_nlls(self.bvalues, s) for s in signals])  # (n,3)
        z = norm.ppf(np.asarray(q_levels))                                    # (L,)
        return point[:, :, None] + self.sd[None, :, None] * z[None, None, :]  # (n,3,L)


def main() -> None:
    cohort = lattice.make_cohort("biexp", n=300, snr=40, seed=lattice.DEFAULT_SEED)
    est = ReferenceNLLSQuantiles(cohort.bvalues, sd=[0.11e-3, 6.0e-3, 0.022])
    q = est.predict_quantiles(cohort.signals, lattice.DEFAULT_QUANTILE_LEVELS)

    lo, hi = central_interval(q, lattice.DEFAULT_QUANTILE_LEVELS, alpha=0.10)
    cov = interval_coverage(cohort.params, lo, hi)
    sharp = mean_sharpness(lo, hi)

    print("Reference NLLS + fixed-SD Gaussian intervals on Lattice (nominal 0.90)")
    print("-" * 60)
    for name, cvi, shi in zip(cohort.param_names, cov, sharp):
        print(f"  {name:6s} coverage={cvi:.3f}  (gap {cvi-0.90:+.3f})  width={shi:.3e}")
    print("\nRaw coverage departs from nominal -> the calibration problem Lattice probes.")
    print("Score this with the canonical ruler via examples/evaluate_with_caliper.py.")


if __name__ == "__main__":
    main()
