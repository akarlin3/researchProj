"""CHECKPOINT 3 — full ensemble (beta_sweep x N_candidates, M_full runs each).

RNG: each run gets a GLOBALLY-UNIQUE seed (seed_base + global_index). The worker creates a
FRESH np.random.default_rng(seed) *inside* the run from that seed only (no global/inherited
state), so multiprocessing cannot duplicate or correlate RNG streams. `--verify` asserts this
before the real run (distinct seeds; same seed -> identical tau; different seeds -> different IC).

Outputs:
  results/ensemble.csv                     one row per run {condition,beta,N,run_index,seed,
                                           tau,event,dwell_stat,rho_std_plateau,collapse_rho_mean}
  results/traces/cond_b{beta}_N{N}.npz     decimated rho_std & rho_mean traces (float32) for the
                                           eps re-sweep (t reconstructed as arange*decimate*dt)
"""
from __future__ import annotations

import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from src.config_io import load_config
from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast
from src.ring_detector import detect_death_ring, precollapse_stats_ring, dwell_stat_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
TRACES = os.path.join(RESULTS, "traces")

CFG = load_config(os.path.join(ROOT, "config.yaml"))
R = CFG["r"]
DT = CFG["dt"]
T_MAX = CFG["T_max"]
DT_HOLD = CFG["dt_hold"]
DECIMATE = CFG["output"]["decimate"]
EPS_STD = CFG["eps_std"]
STOP_EPS = 0.03            # = smallest CP4 eps-sweep value, so traces re-sweep over {0.03..0.06}
DWELL_BAND_HI = 0.10       # terminal-descent band (eyeballed at CP3 pre-flight)
IC = CFG["ic"]


def build_tasks():
    """Deterministic task list with globally-unique seeds."""
    tasks = []
    gidx = 0
    for beta in CFG["beta_sweep"]:
        for N in CFG["N_candidates"]:
            for run_index in range(CFG["M_full"]):
                seed = CFG["seed_base"] + gidx
                tasks.append((float(beta), int(N), int(run_index), int(seed)))
                gidx += 1
    return tasks


def run_one(task):
    beta, N, run_index, seed = task
    P = max(1, int(round(R * N)))
    rng = np.random.default_rng(seed)              # FRESH, from this run's seed only
    th0 = make_ring_ic(N, IC, rng)
    res = integrate_ring_fast(th0, N, P, beta, DT, T_MAX, DECIMATE, seed,
                              stop_eps=STOP_EPS, dt_hold=DT_HOLD)
    tau, ev = detect_death_ring(res.t, res.rho_std, EPS_STD, DT_HOLD, T_MAX)
    st = precollapse_stats_ring(res.t, res.rho_std, res.rho_mean, tau, ev, 50.0, T_MAX)
    dwell = dwell_stat_ring(res.t, res.rho_std, tau, ev, band_hi=DWELL_BAND_HI)
    row = {
        "condition": f"b{beta:.3f}_N{N}", "beta": beta, "N": N, "P": P,
        "run_index": run_index, "seed": seed, "tau": tau, "event": ev,
        "dwell_stat": dwell, "rho_std_plateau": st["rho_std_mean"],
        "collapse_rho_mean": float(res.rho_mean[-1]),
    }
    return row, res.rho_std.astype(np.float32), res.rho_mean.astype(np.float32)


def verify():
    tasks = build_tasks()
    seeds = [t[3] for t in tasks]
    assert len(set(seeds)) == len(seeds), "SEEDS NOT DISTINCT"
    assert seeds == list(range(CFG["seed_base"], CFG["seed_base"] + len(seeds))), "seed layout unexpected"
    print(f"[verify] {len(seeds)} tasks, {len(set(seeds))} distinct seeds "
          f"({seeds[0]}..{seeds[-1]})  OK")
    # reproducibility: same seed -> identical tau
    t0 = (0.12, 64, 0, 999001)
    r1 = run_one(t0)[0]; r2 = run_one(t0)[0]
    assert r1["tau"] == r2["tau"], "NON-REPRODUCIBLE for identical seed"
    print(f"[verify] same seed -> identical tau ({r1['tau']:.2f})  OK")
    # independence: different seeds -> different IC and (generically) different tau
    a = run_one((0.12, 64, 0, 999001))[0]
    b = run_one((0.12, 64, 1, 999002))[0]
    ic_a = make_ring_ic(64, IC, np.random.default_rng(999001))
    ic_b = make_ring_ic(64, IC, np.random.default_rng(999002))
    assert not np.array_equal(ic_a, ic_b), "DIFFERENT SEEDS GAVE IDENTICAL IC"
    print(f"[verify] different seeds -> different IC (max|Δθ0|={np.max(np.abs(ic_a-ic_b)):.3f}), "
          f"tau {a['tau']:.1f} vs {b['tau']:.1f}  OK")
    # run a tiny parallel batch and re-check those seeds reproduce serially
    sample = tasks[:8] + tasks[1500:1508]
    with ProcessPoolExecutor(max_workers=4) as ex:
        par = list(ex.map(run_one, sample))
    for task, (row, _, _) in zip(sample, par):
        ser = run_one(task)[0]
        assert ser["tau"] == row["tau"], f"PARALLEL != SERIAL for seed {task[3]}"
    print(f"[verify] {len(sample)} parallel runs reproduce serially (no fork RNG bleed)  OK")
    print("[verify] ALL RNG CHECKS PASSED")


def main():
    if "--verify" in sys.argv:
        verify()
        return

    os.makedirs(TRACES, exist_ok=True)
    tasks = build_tasks()
    # safety: distinct seeds (the one failure that forces a full re-run)
    seeds = [t[3] for t in tasks]
    assert len(set(seeds)) == len(seeds), "SEEDS NOT DISTINCT — aborting"

    n_workers = max(1, (os.cpu_count() or 2) - 1)
    print(f"CP3 ensemble: {len(tasks)} runs over {len(CFG['beta_sweep'])}x{len(CFG['N_candidates'])} "
          f"conditions, M={CFG['M_full']}, {n_workers} workers, T_max={T_MAX}", flush=True)

    rows = []
    traces = {}  # condition -> dict(rho_std=[], rho_mean=[], seed=[], run_index=[])
    t_start = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for row, rs, rm in ex.map(run_one, tasks, chunksize=4):
            rows.append(row)
            c = row["condition"]
            d = traces.setdefault(c, {"rho_std": [], "rho_mean": [], "seed": [], "run_index": []})
            d["rho_std"].append(rs); d["rho_mean"].append(rm)
            d["seed"].append(row["seed"]); d["run_index"].append(row["run_index"])
            done += 1
            if done % 300 == 0:
                el = time.time() - t_start
                print(f"  {done}/{len(tasks)}  ({el:.0f}s, {el/done*1000:.0f} ms/run)", flush=True)

    # write CSV
    fields = ["condition", "beta", "N", "P", "run_index", "seed", "tau", "event",
              "dwell_stat", "rho_std_plateau", "collapse_rho_mean"]
    csv_path = os.path.join(RESULTS, "ensemble.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # write per-condition trace npz (object arrays of variable-length float32 traces)
    for c, d in traces.items():
        np.savez_compressed(
            os.path.join(TRACES, f"cond_{c}.npz"),
            rho_std=np.array(d["rho_std"], dtype=object),
            rho_mean=np.array(d["rho_mean"], dtype=object),
            seed=np.array(d["seed"]), run_index=np.array(d["run_index"]),
            decimate=DECIMATE, dt=DT,
        )

    el = time.time() - t_start
    print(f"\nDONE {len(rows)} runs in {el:.0f}s. CSV -> {csv_path}\n", flush=True)

    # (beta,N) -> {n_died, n_censored, median tau} table
    import numpy as _np
    print(f"{'beta':>6} {'N':>4} | {'n_died':>7} {'n_cens':>7} | {'median tau':>10} | {'rho_std':>7}")
    print("-" * 60)
    for beta in CFG["beta_sweep"]:
        for N in CFG["N_candidates"]:
            sub = [r for r in rows if abs(r["beta"] - beta) < 1e-9 and r["N"] == N]
            taus = _np.array([r["tau"] for r in sub]); evs = _np.array([r["event"] for r in sub])
            nd = int(evs.sum()); nc = len(sub) - nd
            print(f"{beta:>6.3f} {N:>4} | {nd:>7} {nc:>7} | {_np.median(taus):>10.0f} | "
                  f"{_np.mean([r['rho_std_plateau'] for r in sub]):>7.3f}")


if __name__ == "__main__":
    main()
