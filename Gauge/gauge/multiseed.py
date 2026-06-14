"""Gauge-CI -- multi-seed confidence bands for the conditional / finite-sample
quantities (Checkpoints C-E).

Re-runs the deterministic pipeline under multiple RNG seeds and aggregates the
*empirical* (E) quantities into mean + SD + empirical [5, 95] percentile bands.
Seed ``20260613`` is **seed index 0** and its point estimates are preserved for
the data-resampling-only (non-NN) quantities; for the torch-NN-derived quantities
the point estimate is the **across-seed mean** (the B.0-probe decision: the legacy
seed-0 init is a favourable outlier, so the init epistemic spread must be sampled,
not pinned).

Design
------
* ``collect_seed(seed)`` runs the stages for one seed (cohort + NLLS/HGB/MCMC +
  per-seed-init NN bases via ``build_predictions``; the conditional-attack method
  bank; robustness CP0/CP1; the marginal Table-1 sweep) and returns a flat dict of
  scalar (E)/(G) values (+ grids + per-seed cell counts).
* Each seed's raw quantities persist to ``results/seeds/<seed>.json`` and are
  skipped if already present -- the resumable / parallel-burst hook.
* ``aggregate(...)`` reads every seed json and writes ``results/multiseed.json``
  (tidy) + ``results/multiseed_report.txt`` (human). The marginal (G) items live
  in a clearly separated MC-precision block. ``r`` is aggregated in Fisher-z.

Run:  python -m gauge.multiseed                 # 16 seeds, resume-aware
      python -m gauge.multiseed --n 8 --force
      python -m gauge.multiseed --aggregate-only
"""
import argparse
import json
import os
import time

import numpy as np

from gauge.cohort import DEFAULT_SEED, DSTAR_RANGE
from gauge import conditional_attack as ca
from gauge.conditional import mondrian_split
from gauge.conformal import interval_width, empirical_coverage
from gauge import robustness as rob
from gauge import benchmark as bench
from gauge.baselines import PARAM_NAMES

_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
_SEED_DIR = os.path.join(_RESULTS_DIR, "seeds")
A = ca.REP_ALPHA                     # 0.10 -- the level all (E) numbers use
DSTAR = ca.DSTAR

# Distinct seeds appended after seed-0. Fixed list (no Date/random in scripts).
# Distinct, unrelated integers (avoid aliasing the internal +7/+11/+321 offsets).
EXTRA_SEEDS = [20260614, 20270101, 19990412, 20260720, 271828, 314159,
               40520000, 42424242, 13371337, 90210, 20261612,
               77777, 20240229, 20251231, 555555, 60221408,
               888, 19720101, 11235813, 31415926]


# --------------------------------------------------------------------------- #
# Which leaf quantities are torch-NN-derived (point estimate = across-seed mean).
# Everything else is data-resampling-only (point estimate = seed-0, byte-identical
# to the committed reports). A leaf is NN-derived iff its method/arm uses the MDN /
# PNN / DeepEnsemble base.
# --------------------------------------------------------------------------- #
_NN_METHODS = {"conformalized-MDN", "MDN+LCP/features", "MDN+CondConf/Gibbs",
               "raw-MDN"}
_NN_MARG_ARMS = {"raw:PNN-Gaussian", "raw:MDN-DeepEnsemble",
                 "raw:DeepEnsemble-Point", "conformalized:PNN-Gaussian",
                 "conformalized:MDN-DeepEnsemble",
                 "conformalized:DeepEnsemble-Point"}


def _is_nn(key):
    if key == "attack/width_crlb_r":
        return True                                  # computed on MDN+LCP widths
    if key.startswith("attack/method/"):
        # method names themselves contain '/', so match the full method path prefix
        # rather than positional splitting.
        return any(f"attack/method/{nm}/" in key for nm in _NN_METHODS)
    if key.startswith("marginalG/"):
        return key.split("/")[1] in _NN_MARG_ARMS
    return False


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


# --------------------------------------------------------------------------- #
# One seed -> flat dict of (E)/(G) scalars (+ grids, counts).
# --------------------------------------------------------------------------- #
def collect_seed(seed, verbose=True):
    t0 = time.time()
    R, b, a, feat_cal, feat_train, feat_test, train_extra, names = \
        ca._load_inputs(seed=seed)
    assert abs(a - A) < 1e-12
    dtrue = R["test_true"][:, DSTAR]
    regime, reg_edges = ca._regime_from_true(dtrue)
    test_snr = R["test_snr"]
    snr_levels = sorted(set(int(s) for s in R["meta"]["snr_grid"]))

    methods = ca.build_methods(R, feat_cal, feat_train, feat_test, train_extra, a)
    methods.pop("_lcp_bandwidth", None)
    cal_true_d = R["cal_true"][:, DSTAR]
    methods["split (Mondrian/SNR)"] = mondrian_split(
        R["nlls_cal"][:, DSTAR], cal_true_d, R["cal_snr"],
        R["nlls_test"][:, DSTAR], R["test_snr"], a)
    methods["raw-MDN"] = (R[f"MDN-DeepEnsemble_test_lo_{a}"][:, DSTAR],
                          R[f"MDN-DeepEnsemble_test_hi_{a}"][:, DSTAR])

    cc = {name: ca.conditional_coverage(lo, hi, dtrue, regime, test_snr, snr_levels)
          for name, (lo, hi) in methods.items()}

    flat = {}
    grids = {}
    for name, c in cc.items():
        flat[f"attack/method/{name}/hi_marg"] = float(c["hi_marg"])
        flat[f"attack/method/{name}/hi_worst"] = float(c["hi_worst"])
        hi_row = np.asarray(c["grid"])[ca.N_REGIME - 1]          # hi-D* tercile row
        for k, s in enumerate(snr_levels):
            flat[f"attack/method/{name}/hi_cell/SNR{s}"] = float(hi_row[k])
        grids[name] = np.asarray(c["grid"]).tolist()

    # routing error (non-NN)
    rt = ca.routing_analysis(R)
    flat["attack/routing/misroute_pct"] = float(100.0 * (1 - rt["hi_sensitivity"]))
    flat["attack/routing/hi_sensitivity"] = float(rt["hi_sensitivity"])

    # CRLB(D*)/tercile-width per tercile -- terciles recomputed from THIS seed
    crlb_vox = ca.crlb_per_voxel(R, b)
    finite = np.isfinite(crlb_vox)
    bin_w = np.diff(np.concatenate([[DSTAR_RANGE[0]], reg_edges, [DSTAR_RANGE[1]]]))
    tnames = ["lo", "mid", "hi"]
    for r in range(ca.N_REGIME):
        m = (regime == r) & finite
        med = float(np.median(crlb_vox[m]))
        flat[f"attack/crlb_over_width/{tnames[r]}"] = float(med / bin_w[r])
        flat[f"attack/crlb_median/{tnames[r]}"] = med
    dstars, crlb_abs, _ = ca.crlb_sweep(b, snr_levels)
    lo_abs = float(np.mean([crlb_abs[s][0] for s in snr_levels]))
    hi_abs = float(np.mean([crlb_abs[s][-1] for s in snr_levels]))
    flat["attack/crlb_abs_growth"] = hi_abs / lo_abs

    # width-CRLB log-log r on the MDN+LCP widths (NN-derived)
    lo_b, hi_b = methods["MDN+LCP/features"]
    w = interval_width(lo_b, hi_b)
    mfit = finite & np.isfinite(w)
    flat["attack/width_crlb_r"] = float(np.corrcoef(
        np.log(crlb_vox[mfit] + 1e-12), np.log(w[mfit] + 1e-12))[0, 1])

    # smoking gun: split(Mondrian/D-hat*) by OBSERVED stratum vs by TRUE tercile
    e_hat = rt["edges_hat"]
    obs_str = np.digitize(R["nlls_test"][:, DSTAR], e_hat)
    lo_md, hi_md = methods["split (Mondrian/D-hat*)"]
    for r in range(ca.N_REGIME):
        flat[f"attack/smoking_by_observed/{tnames[r]}"] = float(empirical_coverage(
            lo_md[obs_str == r], hi_md[obs_str == r], dtrue[obs_str == r]))
        flat[f"attack/smoking_by_true/{tnames[r]}"] = float(empirical_coverage(
            lo_md[regime == r], hi_md[regime == r], dtrue[regime == r]))

    # ---- robustness CP0 (per-break) + CP1 (acquisition) + latent (all non-NN) --
    P = rob.compute_all(seed=seed, verbose=False)
    for name, r in P["cp1"].items():
        flat[f"robustness/acq/{name}/hi_marg"] = float(r["hi_marg"])
        flat[f"robustness/acq/{name}/hi_worst"] = float(r["hi_worst"])
        flat[f"robustness/acq/{name}/crlb_over_width"] = float(r["crlb_over_width"])
        flat[f"robustness/acq/{name}/marg_dstar"] = float(r["marg_dstar"])
    for name, r in P["cp0"].items():
        for j, p in enumerate(("D", "Dstar", "f")):
            flat[f"robustness/break/{name}/naive/{p}"] = float(r["cov_naive"][j])
            flat[f"robustness/break/{name}/weighted/{p}"] = float(r["cov_weighted"][j])
    flat["robustness/latent/marginal"] = float(P["latent"]["marginal"])
    flat["robustness/latent/hi_marg"] = float(P["latent"]["hi_marg"])
    flat["robustness/latent/hi_worst"] = float(P["latent"]["hi_worst"])

    # ---- marginal Table-1 coverage @a=0.10 (G items, MC precision) -------------
    arms, M = bench.evaluate(R)
    for (arm, j, aa), met in M.items():
        if abs(aa - A) < 1e-12:
            flat[f"marginalG/{arm}/{PARAM_NAMES[j]}"] = float(met["coverage"])

    counts = np.asarray(cc["conformalized-MDN"]["hi_counts"])
    meta = {
        "seed": int(seed), "alpha": A,
        "reg_edges": [float(x) for x in reg_edges],
        "snr_levels": snr_levels,
        "n": {"test": int(R["meta"]["sizes"]["test"]),
              "hi_tercile": int(counts.sum()),
              "hi_worst_cell": int(counts.min()),
              "hi_cells": [int(x) for x in counts],
              "robustness_test": int(rob.N_TEST)},
        "wall_s": round(time.time() - t0, 1),
    }
    if verbose:
        print(f"[multiseed] seed {seed} collected ({meta['wall_s']:.0f}s)")
    return {"meta": meta, "flat": flat, "grids": grids}


# --------------------------------------------------------------------------- #
# Run / resume
# --------------------------------------------------------------------------- #
def seed_list(n):
    return [DEFAULT_SEED] + EXTRA_SEEDS[:max(0, n - 1)]


def _seed_path(seed):
    return os.path.join(_SEED_DIR, f"{int(seed)}.json")


def run(n=16, force=False, verbose=True):
    os.makedirs(_SEED_DIR, exist_ok=True)
    seeds = seed_list(n)
    for s in seeds:
        path = _seed_path(s)
        if (not force) and os.path.exists(path):
            if verbose:
                print(f"[multiseed] seed {s} present -- skip")
            continue
        rec = collect_seed(s, verbose=verbose)
        with open(path, "w") as fh:
            json.dump(_jsonable(rec), fh)
    return seeds


# --------------------------------------------------------------------------- #
# Aggregation
# --------------------------------------------------------------------------- #
def _load_records(seeds):
    recs = []
    for s in seeds:
        p = _seed_path(s)
        if os.path.exists(p):
            with open(p) as fh:
                recs.append(json.load(fh))
    return recs


def _agg_scalar(vals, key, seed0_val):
    vals = np.asarray(vals, float)
    nn = _is_nn(key)
    if key == "attack/width_crlb_r":                 # Fisher-z aggregation
        z = np.arctanh(np.clip(vals, -0.999999, 0.999999))
        mean = float(np.tanh(z.mean()))
        lo, hi = (float(np.tanh(q)) for q in np.percentile(z, [5, 95]))
        sd = float(np.tanh(z.std(ddof=1) + z.mean()) - mean)  # indicative
    else:
        mean = float(vals.mean())
        lo, hi = (float(q) for q in np.percentile(vals, [5, 95]))
        sd = float(vals.std(ddof=1)) if vals.size > 1 else 0.0
    # Checkpoint-D rule change: the SHIP point for EVERY banded (E) item is the
    # across-seed mean (a frozen seed-0 point can land outside its own [5,95] band
    # and so cannot ship). seed-0 is retained as the point ONLY for the (G)
    # marginal block (marginalG/*); the determinism gate checks seed-0 separately.
    is_g = key.startswith("marginalG/")
    point = float(seed0_val) if is_g else mean
    return {"point": point,
            "point_kind": "seed-0 (G)" if is_g else "across-seed mean (E)",
            "mean": mean, "sd": sd, "lo5": lo, "hi95": hi,
            "vmin": float(vals.min()), "vmax": float(vals.max()),
            "n_seeds": int(vals.size), "nn_derived": bool(nn),
            "seed0": float(seed0_val)}


def aggregate(seeds, verbose=True):
    recs = _load_records(seeds)
    if not recs:
        raise RuntimeError("no seed records found in results/seeds/")
    seed0 = next((r for r in recs if r["meta"]["seed"] == DEFAULT_SEED), recs[0])
    keys = sorted(seed0["flat"].keys())
    agg = {}
    for k in keys:
        vals = [r["flat"][k] for r in recs if k in r["flat"]]
        agg[k] = _agg_scalar(vals, k, seed0["flat"][k])

    # per-cell hi-D* heatmap aggregation (for the D-checkpoint 0.90-crossing flags)
    out = {
        "n_seeds": len(recs),
        "seeds": [r["meta"]["seed"] for r in recs],
        "alpha": A,
        "n_per_seed": seed0["meta"]["n"],
        "items": agg,
    }
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    with open(os.path.join(_RESULTS_DIR, "multiseed.json"), "w") as fh:
        json.dump(_jsonable(out), fh, indent=2)
    _write_report(out, recs)
    if verbose:
        print(f"[multiseed] aggregated {len(recs)} seeds -> results/multiseed.json"
              f" + results/multiseed_report.txt")
    return out


def _fmt(it):
    star = "*" if it["nn_derived"] else " "
    return (f"{it['point']:.3f} [{it['lo5']:.3f}, {it['hi95']:.3f}] "
            f"(mean {it['mean']:.3f} sd {it['sd']:.3f}){star}")


def _write_report(out, recs):
    L = []
    def w(*a):
        L.append(" ".join(str(x) for x in a))
    A_ = out["alpha"]
    w("=" * 100)
    w("GAUGE-CI -- MULTI-SEED CONFIDENCE BANDS (Checkpoints C-E)")
    w("=" * 100)
    w(f"n_seeds = {out['n_seeds']}   seeds[0]={DEFAULT_SEED} (index 0)   alpha={A_}")
    w(f"per-seed n: test={out['n_per_seed']['test']}, hi-D* tercile="
      f"{out['n_per_seed']['hi_tercile']}, hi-D* worst-cell="
      f"{out['n_per_seed']['hi_worst_cell']}, robustness test="
      f"{out['n_per_seed']['robustness_test']}")
    w("point = seed-0 for data-resampling (non-NN) items; across-seed MEAN for "
      "torch-NN-derived items (marked *).")
    w("bands = empirical [5, 95] percentile across seeds; r aggregated in Fisher-z.")
    it = out["items"]

    def block(title, keys):
        w(""); w("-" * 100); w(title); w("-" * 100)
        for label, key in keys:
            if key in it:
                w(f"  {label:<46} {_fmt(it[key])}")

    # ---- (E) Table 2: hi-D* per method --------------------------------------
    t2 = []
    for m in ["raw-MDN", "CQR (plain)", "split (Mondrian/SNR)",
              "CQR (Mondrian/SNR)", "conformalized-MDN"]:
        t2.append((f"{m} hi-D* marg", f"attack/method/{m}/hi_marg"))
        t2.append((f"{m} hi-D* worst-SNR", f"attack/method/{m}/hi_worst"))
    block("(E) TABLE 2 -- high-D* tercile coverage per method", t2)

    # ---- (E) Fig 4: 11 label-free methods -----------------------------------
    f4 = []
    for m in ["CQR (plain)", "conformalized-MDN", "split (Mondrian/D-hat*)",
              "CQR (Mondrian/D-hat*)", "split (LCP/features)", "CQR (LCP/features)",
              "split (CondConf/Gibbs)", "CQR (CondConf/Gibbs)",
              "richer-CQR (signal+proxies)", "MDN+LCP/features",
              "MDN+CondConf/Gibbs"]:
        f4.append((f"{m} hi-D* marg", f"attack/method/{m}/hi_marg"))
        f4.append((f"{m} hi-D* worst", f"attack/method/{m}/hi_worst"))
    block("(E) FIG 4 -- 11 label-free hi-D* methods (marg + worst-SNR)", f4)

    # ---- (E) identifiability -------------------------------------------------
    block("(E) IDENTIFIABILITY -- routing, CRLB ratios, width-CRLB r", [
        ("routing misroute %", "attack/routing/misroute_pct"),
        ("CRLB/tercile-width  lo", "attack/crlb_over_width/lo"),
        ("CRLB/tercile-width  mid", "attack/crlb_over_width/mid"),
        ("CRLB/tercile-width  hi  (>1 = wall)", "attack/crlb_over_width/hi"),
        ("abs CRLB growth (x)", "attack/crlb_abs_growth"),
        ("width-CRLB log-log r", "attack/width_crlb_r"),
        ("smoking gun: hi by OBSERVED stratum", "attack/smoking_by_observed/hi"),
        ("smoking gun: hi by TRUE tercile", "attack/smoking_by_true/hi"),
    ])

    # ---- (E) Table 4: acquisition -------------------------------------------
    t4 = []
    for s in ["clinical (11 b)", "CRLB-optimal (11 b)", "dense (22 b)"]:
        t4.append((f"{s} hi-D* marg", f"robustness/acq/{s}/hi_marg"))
        t4.append((f"{s} hi-D* worst", f"robustness/acq/{s}/hi_worst"))
        t4.append((f"{s} CRLB/width (>1=wall)", f"robustness/acq/{s}/crlb_over_width"))
    block("(E) TABLE 4 -- acquisition: hi-D* coverage + resolution ratio by scheme", t4)

    # ---- (E) Table 3: robustness breaks -------------------------------------
    t3 = []
    for sc in ["SNR shift (low)", "prior shift (harder tissue)", "tri-exp misspec"]:
        for p in ("D", "Dstar", "f"):
            t3.append((f"{sc} naive {p}", f"robustness/break/{sc}/naive/{p}"))
            t3.append((f"{sc} weighted {p}", f"robustness/break/{sc}/weighted/{p}"))
    t3 += [("latent hi-D* marg", "robustness/latent/hi_marg"),
           ("latent hi-D* worst", "robustness/latent/hi_worst")]
    block("(E) TABLE 3 -- exchangeability breaks (naive vs weighted) + latent", t3)

    # ---- (G) MC-precision block ---------------------------------------------
    w(""); w("=" * 100)
    w("(G) MARGINAL COVERAGE -- finite-sample GUARANTEED; band = MC precision, "
      "NOT uncertainty about the guarantee")
    w("=" * 100)
    g = sorted(k for k in it if k.startswith("marginalG/"))
    for k in g:
        w(f"  {k.split('/',1)[1]:<46} {_fmt(it[k])}")
    w("=" * 100)
    with open(os.path.join(_RESULTS_DIR, "multiseed_report.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args()
    if args.aggregate_only:
        aggregate(seed_list(args.n))
        return 0
    seeds = run(n=args.n, force=args.force)
    aggregate(seeds)
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    raise SystemExit(main())
