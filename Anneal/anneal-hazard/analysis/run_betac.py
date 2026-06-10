"""Independent measurement of the ring chimera existence boundary beta_c — ladder driver.

The manuscript cites beta_c ~ 0.13-0.14 from the literature; this campaign measures it
independently with the SAME validated pipeline as CP3 (numba integrate_ring_fast, the
frozen IC recipe, and the frozen death criterion rho_std < 0.04 sustained 50 t.u.,
T_max = 12000), on a beta ladder crossing the boundary from below.

PRE-STATED ENDPOINT (frozen before the ladder ran; design probes only, see
paper/revision-data-gated/results_betac.json):
  Per (beta, N) cell, P_persist = fraction of runs that
    (a) ESTABLISH a chimera: max rho_std over the live window [t_skip=50, t_di] > 0.08,
        where t_di is the last decimated sample <= tau (exact critical_review/cp2_batch.py
        "formed" convention, C_STRUCT=0.08, T_SKIP=50; fallback window [0, tau] if empty), AND
    (b) survive past the fixed horizon T* = 400 t.u. (the pooled median lifetime of the
        last pre-registered known-good cell beta=0.130 at N in {64,128,256}: medians
        422/404/345; the post-boundary IC-transient decay floor seen in design probes at
        beta >= 0.145 has median ~ 230 << T*).
  beta_c(N) = midpoint of a 4-parameter logistic fit to P_persist(beta) (binomial MLE,
  free floor/ceiling), bootstrap CI over runs. Fitting lives in analyze_betac.py; this
  driver only simulates and writes the per-run CSV.

Ladder: beta in {0.1250, 0.1275, ..., 0.1550} (step 0.0025), N in {64, 128, 256},
M = 100 runs per cell -> 3900 runs.
Seeds: GLOBALLY UNIQUE, seed = SEED_BASE + gidx with SEED_BASE = 20270000, disjoint from
all prior campaigns (ensemble*.csv span 20260609..20268108; verified at launch).

Usage:
  python3 analysis/run_betac.py --verify     # RNG/reproducibility gate (run first)
  python3 analysis/run_betac.py              # full ladder -> results/betac_ladder.csv
  python3 analysis/run_betac.py --extend b1,b2,...   # append extra beta rungs (same N, M;
                                             seeds continue from max seed in the CSV + 1)
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from src.config_io import load_config            # noqa: E402
from src.ring_model import make_ring_ic          # noqa: E402
from src.ring_fast import integrate_ring_fast    # noqa: E402
from src.ring_detector import detect_death_ring  # noqa: E402

RESULTS = os.path.join(ROOT, "results")
CSV_PATH = os.path.join(RESULTS, "betac_ladder.csv")

CFG = load_config(os.path.join(ROOT, "config.yaml"))
R, DT, T_MAX, DT_HOLD = CFG["r"], CFG["dt"], CFG["T_max"], CFG["dt_hold"]
DECIMATE, EPS_STD = CFG["output"]["decimate"], CFG["eps_std"]
IC = CFG["ic"]
STOP_EPS = 0.03            # same early-stop as cp3_ensemble / run_extra_N

SEED_BASE = 20270000       # disjoint from 20260609..20268108 (all ensemble*.csv)
BETA_LADDER = [round(0.1250 + 0.0025 * i, 4) for i in range(13)]   # 0.1250 .. 0.1550
N_LADDER = [64, 128, 256]
M_CELL = 100

C_STRUCT = 0.08            # establishment gate (critical_review/cp2_batch.py convention)
T_SKIP = 50.0

FIELDS = ["condition", "beta", "N", "P", "run_index", "seed",
          "tau", "event", "established", "max_rho_std_early"]


def run_one(task):
    """One ladder run: frozen CP3 integration + death detection + cp2 establishment gate."""
    beta, N, run_index, seed = task
    P = max(1, int(round(R * N)))
    rng = np.random.default_rng(seed)              # FRESH, from this run's seed only
    th0 = make_ring_ic(N, IC, rng)
    res = integrate_ring_fast(th0, N, P, beta, DT, T_MAX, DECIMATE, seed,
                              stop_eps=STOP_EPS, dt_hold=DT_HOLD)
    tau, ev = detect_death_ring(res.t, res.rho_std, EPS_STD, DT_HOLD, T_MAX)
    # establishment: EXACT cp2_batch.py 'formed' windowing (live window [T_SKIP, t_di])
    t, rstd = res.t, res.rho_std
    di = int(np.searchsorted(t, tau, side="right")) - 1
    di = min(max(di, 0), len(t) - 1)
    live = (t >= T_SKIP) & (t <= t[di])
    if not live.any():
        live = t <= t[di]
    max_rstd = float(rstd[live].max()) if live.any() else 0.0
    established = int(max_rstd > C_STRUCT)
    return {"condition": f"b{beta:.4f}_N{N}", "beta": beta, "N": N, "P": P,
            "run_index": run_index, "seed": seed, "tau": float(tau), "event": int(ev),
            "established": established, "max_rho_std_early": max_rstd}


def build_tasks(betas, seed_start):
    tasks, gidx = [], 0
    for beta in betas:
        for N in N_LADDER:
            for run_index in range(M_CELL):
                tasks.append((float(beta), int(N), int(run_index), int(seed_start + gidx)))
                gidx += 1
    return tasks


def all_prior_seeds():
    seeds = set()
    for path in glob.glob(os.path.join(RESULTS, "ensemble*.csv")):
        with open(path) as f:
            for row in csv.DictReader(f):
                seeds.add(int(row["seed"]))
    return seeds


def verify():
    tasks = build_tasks(BETA_LADDER, SEED_BASE)
    seeds = [t[3] for t in tasks]
    assert len(set(seeds)) == len(seeds), "SEEDS NOT DISTINCT"
    assert seeds == list(range(SEED_BASE, SEED_BASE + len(seeds))), "seed layout unexpected"
    prior = all_prior_seeds()
    assert not (set(seeds) & prior), "SEEDS COLLIDE WITH PRIOR CAMPAIGNS"
    print(f"[verify] {len(seeds)} tasks, {len(set(seeds))} distinct seeds "
          f"({seeds[0]}..{seeds[-1]}), disjoint from {len(prior)} prior seeds  OK")
    # reproducibility: same seed -> identical tau & establishment
    t0 = (0.14, 64, 0, 998001)
    r1 = run_one(t0); r2 = run_one(t0)
    assert r1["tau"] == r2["tau"] and r1["max_rho_std_early"] == r2["max_rho_std_early"], \
        "NON-REPRODUCIBLE for identical seed"
    print(f"[verify] same seed -> identical tau ({r1['tau']:.2f}) and gate value "
          f"({r1['max_rho_std_early']:.4f})  OK")
    # independence: different seeds -> different IC
    ic_a = make_ring_ic(64, IC, np.random.default_rng(998001))
    ic_b = make_ring_ic(64, IC, np.random.default_rng(998002))
    assert not np.array_equal(ic_a, ic_b), "DIFFERENT SEEDS GAVE IDENTICAL IC"
    b = run_one((0.14, 64, 1, 998002))
    print(f"[verify] different seeds -> different IC "
          f"(max|dtheta0|={np.max(np.abs(ic_a - ic_b)):.3f}), tau {r1['tau']:.1f} vs "
          f"{b['tau']:.1f}  OK")
    # parallel == serial for a small sample
    sample = tasks[:6] + tasks[1950:1956]
    with ProcessPoolExecutor(max_workers=4) as ex:
        par = list(ex.map(run_one, sample))
    for task, row in zip(sample, par):
        ser = run_one(task)
        assert ser["tau"] == row["tau"], f"PARALLEL != SERIAL for seed {task[3]}"
    print(f"[verify] {len(sample)} parallel runs reproduce serially (no fork RNG bleed)  OK")
    print("[verify] ALL RNG CHECKS PASSED")


def run_ladder(betas, seed_start, append):
    tasks = build_tasks(betas, seed_start)
    seeds = [t[3] for t in tasks]
    assert len(set(seeds)) == len(seeds), "SEEDS NOT DISTINCT — aborting"
    prior = all_prior_seeds()
    assert not (set(seeds) & prior), "SEEDS COLLIDE WITH PRIOR CAMPAIGNS — aborting"
    if append and os.path.exists(CSV_PATH):
        with open(CSV_PATH) as f:
            existing_csv = {int(r["seed"]) for r in csv.DictReader(f)}
        assert not (set(seeds) & existing_csv), "SEEDS COLLIDE WITH EXISTING LADDER CSV"

    nw = max(1, (os.cpu_count() or 2) - 1)
    print(f"[run_betac] {len(tasks)} runs: beta in {betas}, N in {N_LADDER}, M={M_CELL}; "
          f"{nw} workers, T_max={T_MAX:.0f}; seeds {seeds[0]}..{seeds[-1]}", flush=True)

    t0 = time.time()
    rows, done = [], 0
    with ProcessPoolExecutor(max_workers=nw) as ex:
        for row in ex.map(run_one, tasks, chunksize=4):
            rows.append(row)
            done += 1
            if done % max(1, len(tasks) // 12) == 0:
                el = time.time() - t0
                print(f"  {done}/{len(tasks)}  ({el:.0f}s, {el / done * 1000:.0f} ms/run)",
                      flush=True)
    wall = time.time() - t0
    print(f"[run_betac] DONE in {wall:.0f}s ({wall / len(tasks) * 1000:.0f} ms/run)", flush=True)

    mode = "a" if (append and os.path.exists(CSV_PATH)) else "w"
    with open(CSV_PATH, mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if mode == "w":
            w.writeheader()
        w.writerows(rows)
    print(f"[run_betac] CSV -> {CSV_PATH} ({mode}, {len(rows)} rows)")

    # quick per-cell table
    print(f"\n{'beta':>7} {'N':>4} | {'P_est':>6} {'P(tau>400)':>10} {'med tau':>8}")
    print("-" * 46)
    for beta in betas:
        for N in N_LADDER:
            sub = [r for r in rows if abs(r["beta"] - beta) < 1e-9 and r["N"] == N]
            est = np.mean([r["established"] for r in sub])
            pp = np.mean([(r["established"] == 1) and (r["tau"] > 400.0) for r in sub])
            med = np.median([r["tau"] for r in sub])
            print(f"{beta:>7.4f} {N:>4} | {est:>6.2f} {pp:>10.2f} {med:>8.0f}")
    return wall


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--extend", type=str, default=None,
                    help="comma-separated extra beta rungs to append (seeds continue)")
    args = ap.parse_args()
    if args.verify:
        verify()
        return
    if args.extend:
        betas = [round(float(b), 4) for b in args.extend.split(",")]
        with open(CSV_PATH) as f:
            mx = max(int(r["seed"]) for r in csv.DictReader(f))
        run_ladder(betas, mx + 1, append=True)
        return
    run_ladder(BETA_LADDER, SEED_BASE, append=False)


if __name__ == "__main__":
    main()
