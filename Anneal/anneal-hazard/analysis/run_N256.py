"""CHECKPOINT A — extend the ring hazard run to N=256 (analysis add-on).

Reuses the EXACT CP3/CP4 pipeline and the frozen β-run protocol:
  - model + integrator: src.ring_fast.integrate_ring_fast (numba), same as CP3
  - IC sampler: src.ring_model.make_ring_ic (same ic config)
  - death: rho_std < eps_std(0.04) sustained dt_hold(50), censored at T_max(12000)
  - fit: src.survival.fit_all -> Weibull with profile-likelihood 95% CI on k (as CP4)
  - M = 300, all 5 β ∈ {0.110,0.115,0.120,0.125,0.130}, N = 256 only

Seeding: continues CP3's GLOBALLY-UNIQUE scheme seed = seed_base + global_index.
CP3 consumed indices 0..4499 (5β × 3N × 300 = seeds 20260609..20265108). The N=256
block continues AFTER that: seed = seed_base + 4500 + (beta_idx*M + run_index)
= 20265109..20266608. Disjoint by construction; asserted against ensemble.csv at launch.

Outputs (new files; the frozen CP4 artifacts are never touched):
  results/ensemble_N256.csv   one row per run (same schema as ensemble.csv)
  results/cp_fits_N256.json   per-β Weibull fit (cp4_fits.json-compatible keys)

Usage:
  python analysis/run_N256.py --pilot 10     # time/size a small batch first
  python analysis/run_N256.py --full         # the real M=300 run
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.config_io import load_config
from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast
from src.ring_detector import detect_death_ring, precollapse_stats_ring, dwell_stat_ring
from src import survival as S

RESULTS = os.path.join(ROOT, "results")
CFG = load_config(os.path.join(ROOT, "config.yaml"))

R = CFG["r"]
DT = CFG["dt"]
T_MAX = CFG["T_max"]
DT_HOLD = CFG["dt_hold"]
DECIMATE = CFG["output"]["decimate"]
EPS_STD = CFG["eps_std"]
STOP_EPS = 0.03                 # identical to CP3 (smallest re-sweep eps)
DWELL_BAND_HI = 0.10            # identical to CP3
IC = CFG["ic"]
SEED_BASE = CFG["seed_base"]
BETA_SWEEP = [float(b) for b in CFG["beta_sweep"]]
M_FULL = CFG["M_full"]
N_TARGET = 256

# CP3 consumed exactly this many globally-unique indices (the ORIGINAL 3-N design):
SEED_OFFSET = len(BETA_SWEEP) * len(CFG["N_candidates"]) * M_FULL   # = 4500


def run_one(task):
    """IDENTICAL to src.cp3_ensemble.run_one (copied verbatim so the protocol matches)."""
    beta, N, run_index, seed = task
    P = max(1, int(round(R * N)))
    rng = np.random.default_rng(seed)              # FRESH, from this run's seed only
    th0 = make_ring_ic(N, IC, rng)
    res = integrate_ring_fast(th0, N, P, beta, DT, T_MAX, DECIMATE, seed,
                              stop_eps=STOP_EPS, dt_hold=DT_HOLD)
    tau, ev = detect_death_ring(res.t, res.rho_std, EPS_STD, DT_HOLD, T_MAX)
    st = precollapse_stats_ring(res.t, res.rho_std, res.rho_mean, tau, ev, 50.0, T_MAX)
    dwell = dwell_stat_ring(res.t, res.rho_std, tau, ev, band_hi=DWELL_BAND_HI)
    return {
        "condition": f"b{beta:.3f}_N{N}", "beta": beta, "N": N, "P": P,
        "run_index": run_index, "seed": seed, "tau": tau, "event": ev,
        "dwell_stat": dwell, "rho_std_plateau": st["rho_std_mean"],
        "collapse_rho_mean": float(res.rho_mean[-1]),
    }


def build_tasks(M):
    """N=256 tasks with seeds continuing CP3's global counter (collision-free)."""
    tasks = []
    for beta_idx, beta in enumerate(BETA_SWEEP):
        for run_index in range(M):
            seed = SEED_BASE + SEED_OFFSET + beta_idx * M + run_index
            tasks.append((beta, N_TARGET, run_index, seed))
    return tasks


def _existing_seeds():
    path = os.path.join(RESULTS, "ensemble.csv")
    seeds = set()
    if os.path.exists(path):
        with open(path) as f:
            for row in csv.DictReader(f):
                seeds.add(int(row["seed"]))
    return seeds


def _fit_beta(rows_beta):
    tau = np.array([r["tau"] for r in rows_beta], float)
    ev = np.array([r["event"] for r in rows_beta], int)
    fit = S.fit_all(tau, ev)
    klo, khi = fit.weibull.k_ci
    return {
        "n": fit.n, "n_died": fit.n_events, "n_cens": fit.n - fit.n_events,
        "median_tau": float(np.median(tau)),
        "max_tau": float(tau.max()), "frac_censored": float((ev == 0).mean()),
        "weibull_k": fit.weibull.k, "weibull_k_lo": klo, "weibull_k_hi": khi,
        "weibull_lambda": fit.weibull.lam, "exp_lambda": fit.exp_lambda,
        "lrt_stat": fit.lrt_stat, "lrt_p": fit.lrt_p,
        "ci_excludes_1": bool(klo == klo and klo > 1.0),   # klo==klo guards NaN
    }


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--pilot", type=int, metavar="K",
                   help="run K runs per β (timing/size probe), write nothing permanent")
    g.add_argument("--full", action="store_true", help="run the real M=300 ensemble")
    args = ap.parse_args()

    M = args.pilot if args.pilot else M_FULL
    mode = f"PILOT (M={M}/β)" if args.pilot else f"FULL (M={M}/β)"
    tasks = build_tasks(M)

    # --- safety: distinct + disjoint from the frozen CP4 seeds -------------------
    seeds = [t[3] for t in tasks]
    assert len(set(seeds)) == len(seeds), "N=256 SEEDS NOT DISTINCT — aborting"
    clash = set(seeds) & _existing_seeds()
    assert not clash, f"N=256 seeds collide with ensemble.csv: {sorted(clash)[:5]}…"

    n_workers = max(1, (os.cpu_count() or 2) - 1)
    print(f"[run_N256] {mode}: {len(tasks)} runs, N={N_TARGET}, P={max(1,round(R*N_TARGET))}, "
          f"{n_workers} workers, T_max={T_MAX:.0f}, eps_std={EPS_STD}, dt_hold={DT_HOLD:.0f}",
          flush=True)
    print(f"[run_N256] seeds {seeds[0]}..{seeds[-1]} (disjoint from CP4 "
          f"{SEED_BASE}..{SEED_BASE+SEED_OFFSET-1})", flush=True)

    t0 = time.time()
    rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for row in ex.map(run_one, tasks, chunksize=2):
            rows.append(row)
            done += 1
            if done % max(1, len(tasks) // 10) == 0:
                el = time.time() - t0
                print(f"  {done}/{len(tasks)}  ({el:.0f}s, {el/done*1000:.0f} ms/run)",
                      flush=True)
    wall = time.time() - t0

    # --- per-β fits + table ------------------------------------------------------
    print(f"\n[run_N256] {mode} DONE in {wall:.0f}s ({wall/len(tasks)*1000:.0f} ms/run avg)\n",
          flush=True)
    header = (f"{'beta':>6} {'N':>4} | {'n_died':>6} {'n_cens':>6} {'cens%':>6} | "
              f"{'med τ':>7} {'max τ':>7} | {'k̂ (95% CI profile-lik)':>26}")
    print(header)
    print("-" * len(header))
    fits = {}
    for beta in BETA_SWEEP:
        rb = [r for r in rows if abs(r["beta"] - beta) < 1e-9]
        f = _fit_beta(rb)
        fits[f"b{beta:.3f}_N{N_TARGET}"] = {"beta": beta, "N": N_TARGET, **f}
        ci = f"{f['weibull_k']:.3f} [{f['weibull_k_lo']:.3f}, {f['weibull_k_hi']:.3f}]"
        print(f"{beta:>6.3f} {N_TARGET:>4} | {f['n_died']:>6} {f['n_cens']:>6} "
              f"{f['frac_censored']*100:>5.1f}% | {f['median_tau']:>7.0f} {f['max_tau']:>7.0f} | "
              f"{ci:>26}")

    worst = max(fits.values(), key=lambda v: v["frac_censored"])
    print(f"\n[run_N256] worst censoring: β={worst['beta']:.3f} at "
          f"{worst['frac_censored']*100:.1f}%  (CP4 worst cell was β0.110/N32 ≈ 9.3%)")

    if args.full:
        # write ensemble rows (same schema as ensemble.csv)
        fields = ["condition", "beta", "N", "P", "run_index", "seed", "tau", "event",
                  "dwell_stat", "rho_std_plateau", "collapse_rho_mean"]
        cpath = os.path.join(RESULTS, "ensemble_N256.csv")
        with open(cpath, "w", newline="") as fcsv:
            w = csv.DictWriter(fcsv, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        jpath = os.path.join(RESULTS, "cp_fits_N256.json")
        meta = {"_wall_seconds": round(wall, 1), "_ms_per_run": round(wall/len(tasks)*1000, 1),
                "_n_runs": len(tasks), "_seed_range": [seeds[0], seeds[-1]],
                "_protocol": {"eps_std": EPS_STD, "dt_hold": DT_HOLD, "T_max": T_MAX,
                              "M": M, "r": R, "dt": DT}}
        with open(jpath, "w") as fjson:
            json.dump({"_meta": meta, **fits}, fjson, indent=2)
        print(f"\n[run_N256] wrote {cpath}\n[run_N256] wrote {jpath}")
    else:
        print("\n[run_N256] pilot mode: nothing written. Re-run with --full for M=300.")


if __name__ == "__main__":
    main()
