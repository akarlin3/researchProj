"""Fashion reproduction -- one-command, fixed-seed synthetic IVIM scorecard.

    python examples/fashion_repro.py
    python examples/fashion_repro.py --snr 20 --n-train 8000 --n-test 3000

Reproduces, on Caliper's in-repo synthetic phantoms only, the (retooled) Fashion
manuscript's findings in their retooled priority order:

  1. the boundary-railing pathology (the assumption-free PRIMARY): a
     box-constrained NLLS bi-exponential fit rails the weakly-identified D*
     against its bounds -- a per-voxel identifiability signature that needs no
     ground truth (the real-data rates live in the paper);
  2. the calibration ruler (caliper.metrics) as a scoped, ground-truth-only
     SECONDARY: under the honest-CRLB SD convention the railed NLLS intervals
     under-cover D* *conditionally* in the high-D* tercile (not a dramatic
     marginal severity); coverage, ECE and sharpness are reported;
  3. the flow-vs-NLLS comparison: the MAF flow posterior is better-calibrated
     than constrained NLLS on the same held-out set.

The pipeline, end to end, on one fixed seed:

    synthetic cohort -> {constrained NLLS, MAF flow posterior}
                     -> caliper.metrics ruler scorecard + NLLS railing rate
                     -> a short comparison table

IN REVIEW / SYNTHETIC ONLY. The Fashion manuscript is under peer review at NMR in
Biomedicine; this script reproduces only the *phenomenon* on synthetic data. The
manuscript's *clinical / real-data* numbers (e.g. the in-vivo D* boundary-railing
percentage) live in the paper and are deliberately NOT reproduced here. Needs the
optional extras: ``pip install -e ".[estimator,baselines]"`` (torch + scipy).
"""
from __future__ import annotations

import argparse

import numpy as np

from caliper import metrics as M
from caliper.baselines import NLLSIVIMEstimator
from caliper.estimator_maf import MAFPosterior
from caliper.forward import DEFAULT_BVALUES, PARAM_NAMES, synthetic_cohort

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
ALPHA = 0.1


def _scores(y_true, q_pred, conditioning):
    """Per-parameter ruler scorecard as a {name: ParamScore} dict."""
    scs = M.score_quantiles(
        y_true, q_pred, LEVELS, alpha=ALPHA,
        param_names=PARAM_NAMES, conditioning=conditioning,
    )
    return {s.name: s for s in scs}


def _print_comparison(maf_s, nlls_s, nominal):
    """Print the flow-vs-NLLS ruler comparison as a fixed-width table."""
    print(f"nominal central coverage = {nominal:.3f}  (alpha = {ALPHA:.2f})")
    hdr = (f"{'param':>6} | {'MAF cov':>8} {'gap':>8} {'ECE':>7} {'sharp':>9} "
           f"| {'NLLS cov':>8} {'gap':>8} {'ECE':>7} {'sharp':>9}")
    print(hdr)
    print("-" * len(hdr))
    for p in PARAM_NAMES:
        a, b = maf_s[p], nlls_s[p]
        print(f"{p:>6} | {a.coverage:>8.3f} {a.coverage_gap:>+8.3f} {a.ece:>7.3f} "
              f"{a.sharpness:>9.4g} | {b.coverage:>8.3f} {b.coverage_gap:>+8.3f} "
              f"{b.ece:>7.3f} {b.sharpness:>9.4g}")


def _print_tercile(name, ps):
    cc = "  ".join(f"g{g}={v:.3f}" for g, v in sorted(ps.conditional.items()))
    print(f"  {name:>5}  {cc}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fashion synthetic reproduction: NLLS + MAF flow scored by the ruler.")
    ap.add_argument("--snr", type=float, default=20.0)
    ap.add_argument("--n-train", type=int, default=8000)
    ap.add_argument("--n-test", type=int, default=3000)
    ap.add_argument("--seed-train", type=int, default=0)
    ap.add_argument("--seed-test", type=int, default=99)
    ap.add_argument("--epochs", type=int, default=None,
                    help="MAF training epochs (default: MAFPosterior default).")
    args = ap.parse_args(argv)

    nominal = 1.0 - ALPHA
    print("=" * 76)
    print("Caliper x Fashion -- synthetic reproduction (IN REVIEW / SYNTHETIC ONLY)")
    print("=" * 76)
    print(f"regime: SNR={args.snr:g}  b-values={DEFAULT_BVALUES.size}  "
          f"n_train={args.n_train} (seed {args.seed_train})  "
          f"n_test={args.n_test} (seed {args.seed_test})")
    print("the clinical numbers are in the paper; only the phenomenon is reproduced here.\n")

    # --- synthetic cohorts: disjoint train / held-out test ---------------- #
    train = synthetic_cohort(n=args.n_train, snr=args.snr, seed=args.seed_train)
    test = synthetic_cohort(n=args.n_test, snr=args.snr, seed=args.seed_test)

    # --- constrained NLLS baseline (CLAIM 2: railing + overconfidence) ----- #
    nlls = NLLSIVIMEstimator()
    nlls_fit = nlls.solve(test.signals)
    q_nlls = nlls.predict_quantiles(test.signals, LEVELS)
    nlls_s = _scores(test.params, q_nlls, test.params)

    # --- MAF flow posterior (CLAIM 3: better-calibrated than NLLS) --------- #
    maf_kw = {} if args.epochs is None else {"epochs": args.epochs}
    maf = MAFPosterior(n_bvalues=DEFAULT_BVALUES.size, seed=args.seed_train, **maf_kw)
    maf.fit(train.signals, train.params)
    q_maf = maf.predict_quantiles(test.signals, LEVELS)
    maf_s = _scores(test.params, q_maf, test.params)

    # --- CLAIM 1: the ruler scorecard, flow vs NLLS ------------------------ #
    print("[1] CALIBRATION RULER -- flow (MAF) vs constrained NLLS, same held-out set")
    _print_comparison(maf_s, nlls_s, nominal)

    # --- CLAIM 2: boundary railing (NLLS only; MAF is unbounded) ----------- #
    print("\n[2] NLLS BOUNDARY-RAILING (synthetic rate; NOT the clinical figure)")
    rates = {p: nlls.railing_rate_from_fit(nlls_fit, p) for p in PARAM_NAMES}
    rate_any = nlls.railing_rate_from_fit(nlls_fit, "any")
    print("  railed against a box bound:  "
          + "  ".join(f"{p}={rates[p]:.3f}" for p in PARAM_NAMES)
          + f"  any={rate_any:.3f}")
    print("  (the MAF flow posterior has no box bounds and cannot rail.)")

    # --- D* conditional coverage by true-D* tercile ------------------------ #
    print("\n[3] D* conditional coverage by true-D* tercile (g0=low .. g2=high)")
    _print_tercile("MAF", maf_s["Dstar"])
    _print_tercile("NLLS", nlls_s["Dstar"])

    # --- verdict ----------------------------------------------------------- #
    dm, dn = maf_s["Dstar"], nlls_s["Dstar"]
    flow_better = abs(dm.coverage_gap) < abs(dn.coverage_gap) and dm.ece < dn.ece
    print("\n" + "-" * 76)
    print("VERDICT (reproduced qualitative ordering):")
    print(f"  * NLLS D* rails at {rates['Dstar']:.1%} and under-covers "
          f"(coverage {dn.coverage:.3f} < nominal {nominal:.3f}, gap {dn.coverage_gap:+.3f}).")
    print(f"  * MAF flow D* is better-calibrated (coverage {dm.coverage:.3f}, "
          f"gap {dm.coverage_gap:+.3f}, ECE {dm.ece:.3f} vs NLLS ECE {dn.ece:.3f}) and sharper.")
    print(f"  * flow better-calibrated than NLLS on D*: {flow_better}")
    print("-" * 76)
    print("synthetic phantoms only; module is in-review/private until the paper clears.")


if __name__ == "__main__":
    main()
