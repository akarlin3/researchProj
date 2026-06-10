"""Regenerate decimated rho_std/rho_mean traces for an ensemble cell from its logged seeds
(deterministic; reproduces tau bit-for-bit). Used to fill the N=192/256 cells that were not
trace-saved, so the CP3 alternative-criterion re-detection can run on the SAME trajectories.

Usage: python3 regen_traces.py <beta> <N> <ensemble_csv>
Saves critical_review/cp3_criterion/cond_b{beta}_N{N}.npz {rho_std,rho_mean,seed,decimate,dt}.
"""
from __future__ import annotations

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RING = os.path.join(ROOT, "anneal-hazard")
sys.path.insert(0, RING)

from src.config_io import load_config  # noqa: E402
from src.ring_model import make_ring_ic  # noqa: E402
from src.ring_fast import integrate_ring_fast  # noqa: E402
from src.ring_detector import detect_death_ring  # noqa: E402

CFG = load_config(os.path.join(RING, "config.yaml"))
R = CFG["r"]; DT = CFG["dt"]; T_MAX = CFG["T_max"]; DT_HOLD = CFG["dt_hold"]
DECIMATE = CFG["output"]["decimate"]; EPS_STD = CFG["eps_std"]
STOP_EPS = 0.03; IC = CFG["ic"]
OUT = os.path.join(HERE, "cp3_criterion")
os.makedirs(OUT, exist_ok=True)


def run_one(args):
    beta, N, seed, tau_ens = args
    P = max(1, int(round(R * N)))
    th0 = make_ring_ic(N, IC, np.random.default_rng(seed))
    res = integrate_ring_fast(th0, N, P, beta, DT, T_MAX, DECIMATE, seed,
                              stop_eps=STOP_EPS, dt_hold=DT_HOLD)
    tau, ev = detect_death_ring(res.t, res.rho_std, EPS_STD, DT_HOLD, T_MAX)
    return seed, res.rho_std.astype(np.float32), res.rho_mean.astype(np.float32), \
        float(tau), int(ev), abs(tau - tau_ens) < 1e-6


def main(beta, N, ens_csv):
    seeds = []
    with open(ens_csv) as f:
        for r in csv.DictReader(f):
            if abs(float(r["beta"]) - beta) < 1e-9 and int(r["N"]) == N:
                seeds.append((beta, N, int(r["seed"]), float(r["tau"])))
    seeds.sort(key=lambda x: x[2])
    print(f"regen beta={beta} N={N}: {len(seeds)} runs")
    out = {"seed": [], "rho_std": [], "rho_mean": []}
    nmis = 0
    nw = max(1, (os.cpu_count() or 2) - 1)
    with ProcessPoolExecutor(max_workers=nw) as ex:
        for seed, rs, rm, tau, ev, match in ex.map(run_one, seeds, chunksize=4):
            out["seed"].append(seed); out["rho_std"].append(rs); out["rho_mean"].append(rm)
            nmis += (0 if match else 1)
    print(f"  tau mismatches: {nmis}/{len(seeds)}")
    p = os.path.join(OUT, f"cond_b{beta:.3f}_N{N}.npz")
    np.savez_compressed(p, seed=np.array(out["seed"]),
                        rho_std=np.array(out["rho_std"], dtype=object),
                        rho_mean=np.array(out["rho_mean"], dtype=object),
                        decimate=DECIMATE, dt=DT)
    print(f"  [saved] {p}")


if __name__ == "__main__":
    main(float(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
