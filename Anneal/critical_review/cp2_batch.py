"""CP2 (C2) full ensemble re-run + classification + equilibration-filtered Weibull refit.

For beta in {0.110, 0.130} at N=256 (the cells C2 targets), re-run ALL 300 logged seeds with
field dumping (tau reproduces the ensemble bit-for-bit), classify each run's pre-death object,
and refit the censored Weibull on the equilibration-filtered set.

Per-run classification (c_struct = 0.08 spatial-contrast gate; head count = dominant spatial
Fourier wavenumber of rho_k, drift-invariant):
  - formed                : max rho_std over [t_skip, tau] > c_struct  (a chimera ever existed)
  - structured_life_frac  : fraction of [t_skip, tau] with rho_std > c_struct
  - last-10% window stats : struct_frac, single_head_frac (m*==1 among structured snaps),
                            dominant_head
  - category:
      never_formed       : not formed
      degenerate_channel : formed but structured_life_frac < 0.5  (chimera dissolved long
                           before the rho_std<0.04 criterion fired -> NON-canonical death channel)
      multi_head         : structured to death, dominant pre-death head count >= 2
      canonical_single   : structured to death, dominant pre-death head count == 1

Filters for the refit:
  unfiltered     : all runs (must reproduce the published k for the cell)
  chimera_death  : drop never_formed + degenerate_channel (keep single+multi-head chimera deaths)
  single_arc     : keep only canonical_single

Decision (spec): high canonical/chimera-death fraction AND filtered k>1 survives -> validated;
else STOP & report contamination.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RING = os.path.join(ROOT, "anneal-hazard")
sys.path.insert(0, RING)
sys.path.insert(0, HERE)

from cp2_fields import (run_field, rho_profile, detect_death_ring,  # noqa: E402
                        EPS_STD, DT_HOLD, T_MAX)
from src.survival import fit_all  # noqa: E402

OUT = os.path.join(HERE, "cp2_validation")
os.makedirs(OUT, exist_ok=True)
ENS = os.path.join(RING, "results", "ensemble_N256.csv")

C_STRUCT = 0.08
T_SKIP = 50.0
LAST_FRAC = 0.10
MAXM = 12
CELLS = [(0.110, 256), (0.130, 256)]


def head_count(rho):
    p = np.abs(np.fft.rfft(rho - rho.mean())) ** 2
    if len(p) <= 1:
        return 0
    return 1 + int(np.argmax(p[1:MAXM + 1]))


def classify(beta, N, seed, tau_ens):
    t, rmean, rstd, th, P = run_field(beta, N, seed)
    tau, ev = detect_death_ring(t, rstd, EPS_STD, DT_HOLD, T_MAX)
    di = int(np.searchsorted(t, tau, side="right")) - 1
    di = min(max(di, 0), len(t) - 1)
    live = (t >= T_SKIP) & (t <= t[di])
    if not live.any():
        live = t <= t[di]
    rstd_live = rstd[live]
    max_rstd = float(rstd_live.max()) if rstd_live.size else 0.0
    formed = max_rstd > C_STRUCT
    structured_life_frac = float((rstd_live > C_STRUCT).mean()) if rstd_live.size else 0.0

    win = np.where((t >= tau - LAST_FRAC * tau) & (t <= t[di]))[0]
    if win.size == 0:
        win = np.array([di])
    struct_mask = rstd[win] > C_STRUCT
    heads = [head_count(rho_profile(th[i], P)) for i, s in zip(win, struct_mask) if s]
    struct_frac10 = float(struct_mask.mean())
    single_head_frac = float(np.mean([h == 1 for h in heads])) if heads else 0.0
    dominant_head = Counter(heads).most_common(1)[0][0] if heads else 0

    if not formed:
        cat = "never_formed"
    elif structured_life_frac < 0.5:
        cat = "degenerate_channel"
    elif dominant_head >= 2:
        cat = "multi_head"
    else:
        cat = "canonical_single"

    rec = {
        "beta": beta, "N": N, "seed": seed, "tau": float(tau), "tau_ens": float(tau_ens),
        "tau_match": bool(abs(tau - tau_ens) < 1e-6), "event": int(ev),
        "formed": int(formed), "max_rho_std": max_rstd,
        "structured_life_frac": structured_life_frac,
        "struct_frac_last10": struct_frac10, "single_head_frac_last10": single_head_frac,
        "dominant_head": int(dominant_head), "n_struct_snaps_last10": len(heads),
        "category": cat,
    }
    return rec, rstd.astype(np.float32), rmean.astype(np.float32)


def _worker(args):
    return classify(*args)


def load_cell_seeds(beta, N):
    out = []
    with open(ENS) as f:
        for r in csv.DictReader(f):
            if abs(float(r["beta"]) - beta) < 1e-9 and int(r["N"]) == N:
                out.append((int(r["run_index"]), int(r["seed"]), float(r["tau"])))
    return sorted(out)


def fit_subset(records, keep):
    sub = [r for r in records if keep(r)]
    if len(sub) < 10:
        return {"n": len(sub), "k": float("nan"), "k_lo": float("nan"),
                "k_hi": float("nan"), "ci_excl_1": False, "lrt_p": float("nan")}
    tau = np.array([r["tau"] for r in sub]); ev = np.array([r["event"] for r in sub])
    fs = fit_all(tau, ev)
    klo, khi = fs.weibull.k_ci
    return {"n": len(sub), "n_events": int(fs.n_events), "k": float(fs.weibull.k),
            "k_lo": float(klo), "k_hi": float(khi),
            "ci_excl_1": bool(np.isfinite(klo) and klo > 1.0), "lrt_p": float(fs.lrt_p)}


def main():
    n_workers = max(1, (os.cpu_count() or 2) - 1)
    all_blob = {"c_struct": C_STRUCT, "last_frac": LAST_FRAC, "cells": {}}
    for beta, N in CELLS:
        seeds = load_cell_seeds(beta, N)
        print(f"\n{'='*78}\nCELL beta={beta} N={N}: {len(seeds)} runs ({n_workers} workers)\n{'='*78}",
              flush=True)
        args = [(beta, N, s, te) for _, s, te in seeds]
        records = []
        traces = {"seed": [], "rho_std": [], "rho_mean": []}
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            for rec, rs, rm in ex.map(_worker, args, chunksize=4):
                records.append(rec)
                traces["seed"].append(rec["seed"])
                traces["rho_std"].append(rs); traces["rho_mean"].append(rm)
        # provenance: all tau reproduce?
        nmis = sum(1 for r in records if not r["tau_match"])
        cats = Counter(r["category"] for r in records)
        ntot = len(records)
        print(f"  tau mismatches: {nmis}/{ntot}")
        print(f"  categories: " + "  ".join(f"{k}={v} ({100*v/ntot:.1f}%)"
                                              for k, v in cats.most_common()))
        canonical_frac = cats.get("canonical_single", 0) / ntot
        chimera_death_frac = (cats.get("canonical_single", 0) + cats.get("multi_head", 0)) / ntot
        degenerate_frac = (cats.get("degenerate_channel", 0) + cats.get("never_formed", 0)) / ntot
        print(f"  canonical(single-arc) frac = {canonical_frac:.3f} | "
              f"chimera-death(1+multi) frac = {chimera_death_frac:.3f} | "
              f"degenerate/never frac = {degenerate_frac:.3f}")

        f_un = fit_subset(records, lambda r: True)
        f_ch = fit_subset(records, lambda r: r["category"] in ("canonical_single", "multi_head"))
        f_sa = fit_subset(records, lambda r: r["category"] == "canonical_single")
        for lab, fr in [("unfiltered ", f_un), ("chimera-death", f_ch), ("single-arc ", f_sa)]:
            print(f"  k[{lab}] n={fr['n']:>3} k={fr['k']:.3f} "
                  f"CI=[{fr['k_lo']:.3f},{fr['k_hi']:.3f}] CI>1={fr['ci_excl_1']} "
                  f"LRTp={fr['lrt_p']:.1e}")

        # save per-run records + traces
        with open(os.path.join(OUT, f"cp2_runs_b{beta}_N{N}.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(records[0].keys()))
            w.writeheader(); w.writerows(records)
        np.savez_compressed(
            os.path.join(OUT, f"cp2_traces_b{beta}_N{N}.npz"),
            seed=np.array(traces["seed"]),
            rho_std=np.array(traces["rho_std"], dtype=object),
            rho_mean=np.array(traces["rho_mean"], dtype=object),
            decimate=10, dt=0.05)

        all_blob["cells"][f"b{beta}_N{N}"] = {
            "n": ntot, "tau_mismatches": nmis, "categories": dict(cats),
            "canonical_single_frac": canonical_frac,
            "chimera_death_frac": chimera_death_frac,
            "degenerate_frac": degenerate_frac,
            "fit_unfiltered": f_un, "fit_chimera_death": f_ch, "fit_single_arc": f_sa,
        }
    with open(os.path.join(OUT, "cp2_summary.json"), "w") as f:
        json.dump(all_blob, f, indent=2)
    print(f"\n[saved] {os.path.join(OUT, 'cp2_summary.json')}")


if __name__ == "__main__":
    main()
