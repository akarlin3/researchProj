"""Diagnostic: wide (beta, r) grid at fixed small N to locate a pocket with BOTH a clear
chimera (rho_std healthy) AND substantial deaths within T_max. Prints frac|rho_std grid.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast as integrate_ring
from src.ring_detector import detect_death_ring, precollapse_stats_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    nseed = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    T_max = float(sys.argv[3]) if len(sys.argv) > 3 else 6000.0
    betas = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [0.12, 0.14, 0.16, 0.18]
    rs = [float(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [0.10, 0.15, 0.22, 0.30, 0.35]
    eps = 0.04; dt_hold = 50.0; dt = 0.05
    seeds = [20260609 + k for k in range(nseed)]
    ic = {"incoherent_frac": 0.5, "coherent_scale": 0.05}

    print(f"(beta,r) grid  N={N} nseed={nseed} T_max={T_max} eps={eps}", flush=True)
    print("cell = death_frac | rho_std | median(died)", flush=True)
    hdr = "beta\\r | " + " ".join(f"{r:>22.2f}" for r in rs)
    print(hdr); print("-" * len(hdr), flush=True)
    out = {}
    for beta in betas:
        cells = []
        for r in rs:
            P = max(1, int(round(r * N)))
            taus, evs, plats = [], [], []
            for seed in seeds:
                rng = np.random.default_rng(seed)
                th0 = make_ring_ic(N, ic, rng)
                res = integrate_ring(th0, N, P, beta, dt, T_max, decimate=10, seed=seed,
                                     stop_eps=eps, dt_hold=dt_hold)
                tau, ev = detect_death_ring(res.t, res.rho_std, eps, dt_hold, T_max)
                st = precollapse_stats_ring(res.t, res.rho_std, res.rho_mean, tau, ev, 50.0, T_max)
                taus.append(tau); evs.append(int(ev)); plats.append(st["rho_std_mean"])
            taus = np.array(taus); evs = np.array(evs)
            frac = float(evs.mean()); plat = float(np.nanmean(plats))
            medd = float(np.median(taus[evs == 1])) if evs.sum() else float("nan")
            out[f"{beta},{r}"] = {"frac": frac, "rho_std": plat, "median_died": medd,
                                  "taus": taus.tolist(), "events": evs.tolist()}
            cells.append(f"{frac:.2f}|{plat:.3f}|{medd if not np.isnan(medd) else 0:5.0f}")
        print(f"{beta:5.3f} | " + " ".join(f"{c:>22}" for c in cells), flush=True)
    with open(os.path.join(RESULTS, f"ring_grid_N{N}.json"), "w") as f:
        json.dump({"N": N, "T_max": T_max, "nseed": nseed, "cells": out}, f, indent=2)
    print(f"\njson -> ring_grid_N{N}.json", flush=True)


if __name__ == "__main__":
    main()
