"""Characterise chimera lifetime vs N in the ring model (the headline relationship).
Runs nseed seeds per N, records (tau, event), and reports death fraction, median, and
the KM-style median. Saves per-seed taus to JSON and a lifetime-vs-N figure with an
exponential ECDF overlay at one N (memoryless check). Picks the experiment's N/T_max.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast as integrate_ring
from src.ring_detector import detect_death_ring, precollapse_stats_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")


def main():
    beta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.12
    r = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
    Ns = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [32, 40, 48, 56, 64]
    nseed = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    T_max = float(sys.argv[5]) if len(sys.argv) > 5 else 5000.0
    eps = 0.04
    dt_hold = 50.0
    dt = 0.05
    seeds = [20260609 + k for k in range(nseed)]
    ic = {"incoherent_frac": 0.5, "coherent_scale": 0.05}

    print(f"lifetime vs N  beta={beta} r={r} eps={eps} dt_hold={dt_hold} "
          f"T_max={T_max} dt={dt} nseed={nseed}", flush=True)
    print(f"{'N':>4} {'P':>4} | {'died':>6} {'frac':>5} | {'med(all)':>8} {'med(died)':>9} "
          f"| {'mean(died)':>10} | {'rho_std':>7}", flush=True)
    print("-" * 78, flush=True)

    data = {}
    for N in Ns:
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
        taus = np.array(taus, float); evs = np.array(evs, int)
        ndied = int(evs.sum()); frac = ndied / nseed
        med_all = float(np.median(taus))
        med_died = float(np.median(taus[evs == 1])) if ndied else float("nan")
        mean_died = float(np.mean(taus[evs == 1])) if ndied else float("nan")
        plat = float(np.nanmean(plats))
        data[N] = {"P": P, "taus": taus.tolist(), "events": evs.tolist(),
                   "n_died": ndied, "frac": frac, "median_all": med_all,
                   "median_died": med_died, "mean_died": mean_died, "rho_std": plat}
        print(f"{N:4d} {P:4d} | {ndied:2d}/{nseed:<3d} {frac:5.2f} | {med_all:8.0f} {med_died:9.0f} "
              f"| {mean_died:10.0f} | {plat:7.3f}", flush=True)

    meta = {"beta": beta, "r": r, "eps": eps, "dt_hold": dt_hold, "T_max": T_max,
            "dt": dt, "nseed": nseed, "seeds": seeds, "by_N": data}
    jpath = os.path.join(RESULTS, f"ring_lifetime_vs_N_b{beta}_r{r}.json")
    with open(jpath, "w") as f:
        json.dump(meta, f, indent=2)

    # figure: median lifetime vs N (log-y) + ECDF vs exponential at the largest fully-dying N
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    Ns_arr = np.array(Ns)
    med_died_arr = np.array([data[N]["median_died"] for N in Ns])
    ax[0].semilogy(Ns_arr, med_died_arr, "o-", color="C0")
    ax[0].set_xlabel("N"); ax[0].set_ylabel("median lifetime (deaths)")
    ax[0].set_title(f"lifetime vs N  (β={beta}, r={r})")
    ax[0].grid(True, which="both", alpha=0.3)
    for N in Ns:
        ax[0].annotate(f"{data[N]['frac']*100:.0f}% died", (N, data[N]["median_died"]),
                       fontsize=7, textcoords="offset points", xytext=(0, 6), ha="center")

    # pick the N with the most deaths for the ECDF/exponential check
    bestN = max(Ns, key=lambda n: data[n]["n_died"])
    dd = np.array(data[bestN]["taus"])[np.array(data[bestN]["events"]) == 1]
    if len(dd) >= 3:
        dd_sorted = np.sort(dd)
        ecdf = np.arange(1, len(dd_sorted) + 1) / len(dd_sorted)
        ax[1].step(dd_sorted, ecdf, where="post", color="C3", label=f"empirical (N={bestN}, n={len(dd)})")
        lam = 1.0 / dd.mean()
        tt = np.linspace(0, dd_sorted.max() * 1.05, 200)
        ax[1].plot(tt, 1 - np.exp(-lam * tt), "k--", lw=1, label="exponential (same mean)")
        ax[1].set_xlabel("lifetime τ"); ax[1].set_ylabel("CDF")
        ax[1].set_title(f"lifetime CDF vs exponential (N={bestN})")
        ax[1].legend(fontsize=8)
    fig.tight_layout()
    figpath = os.path.join(RESULTS, f"ring_lifetime_vs_N_b{beta}_r{r}.png")
    fig.savefig(figpath, dpi=120)
    plt.close(fig)
    print(f"\njson -> {jpath}\nfig  -> {figpath}", flush=True)


if __name__ == "__main__":
    main()
