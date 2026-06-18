"""Caliper conformal demo -- numpy only, one command, fixed seeds.

    synthetic IVIM  ->  over-confident reference estimator  ->  CQR + Mondrian
                    ->  calibration ruler + D* tercile table

Run:  python examples/conformal_demo.py      (no torch required)

Every number printed here is produced by this run; the README quotes this
output verbatim. The story in three movements:

1. The raw reference estimator is over-confident (reported quantiles too narrow).
2. Marginal CQR restores *pooled* coverage to nominal -- but not *conditional*
   coverage: stratified by true D* tercile, the well-identified low-D* tercile
   ends up over-covered while the poorly-identified high-D* tercile stays
   under-covered. One global correction cannot serve a steeply heteroscedastic
   problem.
3. Mondrian (group-conditional) CQR restores per-tercile coverage -- but only by
   inflating the high-D* interval width. That width ratio is the price of
   conditional validity at the identifiability wall.
"""
from __future__ import annotations

import numpy as np

from caliper import conformal as C
from caliper import metrics as M
from caliper.estimator_reference import ReferenceIVIMEstimator
from caliper.forward import PARAM_NAMES, synthetic_cohort

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
ALPHA = 0.10           # 90% central intervals
SNR = 40.0
DSTAR = PARAM_NAMES.index("Dstar")
STRATUM_NAMES = {0: "low-D*", 1: "mid-D*", 2: "high-D*"}


def main() -> None:
    print("Caliper conformal demo -- synthetic IVIM -> reference -> CQR/Mondrian "
          "-> ruler\n")

    # 1. synthetic, PHI-free cohorts (in-repo generator, fixed seeds)
    cal = synthetic_cohort(n=4000, snr=SNR, seed=1)
    test = synthetic_cohort(n=9000, snr=SNR, seed=2)
    print(f"cohorts: cal(SNR{SNR:.0f}, n={len(cal)}), test(SNR{SNR:.0f}, "
          f"n={len(test)}); b-values={test.bvalues.size}; nominal coverage "
          f"{1 - ALPHA:.3f}\n")

    est = ReferenceIVIMEstimator()
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)

    # 2. raw scorecard (conditional coverage conditioned on the true params)
    raw = M.score_quantiles(test.params, q_test, LEVELS, alpha=ALPHA,
                            param_names=PARAM_NAMES, conditioning=test.params)
    print(M.format_scorecard(raw, title="RAW reference estimator (over-confident)"))
    print()

    # 3. marginal CQR
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, cal.params)
    q_cqr = cq.apply(q_test)
    cor = M.score_quantiles(test.params, q_cqr, LEVELS, alpha=ALPHA,
                            param_names=PARAM_NAMES, conditioning=test.params)
    print(M.format_scorecard(cor, title="MARGINAL CQR"))

    # 3b. CP1 headline: marginal |coverage - nominal| pre/post CQR
    print("\n=== marginal coverage, pre/post CQR (nominal 0.900) ===")
    print(f"{'param':>7} {'raw cov':>8} {'raw|gap|':>9} {'cqr cov':>8} {'cqr|gap|':>9}")
    for r, c in zip(raw, cor):
        print(f"{r.name:>7} {r.coverage:>8.3f} {abs(r.coverage_gap):>9.3f} "
              f"{c.coverage:>8.3f} {abs(c.coverage_gap):>9.3f}")

    # 4. the conditional-coverage result: D* terciles, three methods
    strata = M.tercile_groups(test.params[:, DSTAR])          # true D* terciles
    groups_cal = M.tercile_groups(cal.params[:, DSTAR])
    mq = C.MondrianConformalQuantile(LEVELS).calibrate(q_cal, cal.params, groups_cal)
    q_mond = mq.apply(q_test, strata)                          # oracle-stratified

    def dstar_interval(q):
        return M.central_interval(q[:, DSTAR, :], LEVELS, ALPHA)

    methods = {"raw": q_test, "marginal-CQR": q_cqr, "Mondrian-CQR": q_mond}
    per_method = {}
    for name, q in methods.items():
        lo, hi = dstar_interval(q)
        per_method[name] = C.conditional_coverage_by_strata(
            test.params[:, DSTAR], lo, hi, strata)

    print("\n" + C.format_strata_table(
        per_method, stratum_names=STRATUM_NAMES,
        title="D* coverage & mean interval width by tercile", nominal=1 - ALPHA))

    # 5. honest conclusion, in the toolkit's own terms
    pooled = per_method["marginal-CQR"]
    mond = per_method["Mondrian-CQR"]
    width_ratio = mond[2].mean_width / mond[0].mean_width
    print("\n--- reading the result honestly ---")
    print(f"* Marginal CQR restores POOLED D* coverage to "
          f"{M.empirical_coverage(test.params[:, DSTAR], *dstar_interval(q_cqr)):.3f} "
          f"(nominal {1 - ALPHA:.3f}),")
    print(f"  but CONDITIONAL coverage is not delivered: low-D* over-covers "
          f"({pooled[0].coverage:.3f}) while high-D* under-covers "
          f"({pooled[2].coverage:.3f}).")
    print(f"* Mondrian CQR equalizes per-tercile coverage "
          f"({mond[0].coverage:.3f}/{mond[1].coverage:.3f}/{mond[2].coverage:.3f}) "
          f"only by inflating width:")
    print(f"  high-D* intervals are {width_ratio:.2f}x the low-D* width.")
    print("* Conformal guarantees marginal coverage unconditionally; conditional "
          "coverage")
    print("  costs sharpness, and at high D* (the identifiability wall) the trade "
          "is steep.")


if __name__ == "__main__":
    main()
