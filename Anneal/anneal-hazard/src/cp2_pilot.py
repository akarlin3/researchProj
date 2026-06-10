"""CHECKPOINT 2 — pilot (M=30) at a representative condition.

Reports: died vs censored, median tau, lifetime histogram, and an eps-sensitivity table
(re-running the detector over eps_std in {0.03,0.04,0.05,0.06} on the STORED rho_std traces,
no re-simulation). Integrates with stop_eps = smallest sweep eps so all eps are re-detectable.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.config_io import load_config
from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast
from src.ring_detector import detect_death_ring, precollapse_stats_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")

EPS_SWEEP = [0.03, 0.04, 0.05, 0.06]


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    beta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.12
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    M = int(sys.argv[3]) if len(sys.argv) > 3 else cfg["M_pilot"]

    r = cfg["r"]; P = max(1, int(round(r * N)))
    dt = cfg["dt"]; T_max = cfg["T_max"]; dt_hold = cfg["dt_hold"]
    decimate = cfg["output"]["decimate"]
    eps_primary = cfg["eps_std"]
    stop_eps = min(EPS_SWEEP)  # strictest -> latest death -> trace covers all sweep eps
    seeds = [cfg["seed_base"] + k for k in range(M)]

    traces = []  # (t, rho_std) per run for eps re-sweep
    taus, events, plats = [], [], []
    for k, seed in enumerate(seeds):
        rng = np.random.default_rng(seed)
        th0 = make_ring_ic(N, cfg["ic"], rng)
        res = integrate_ring_fast(th0, N, P, beta, dt, T_max, decimate, seed,
                                  stop_eps=stop_eps, dt_hold=dt_hold)
        traces.append((res.t, res.rho_std))
        tau, ev = detect_death_ring(res.t, res.rho_std, eps_primary, dt_hold, T_max)
        st = precollapse_stats_ring(res.t, res.rho_std, res.rho_mean, tau, ev, 50.0, T_max)
        taus.append(tau); events.append(ev); plats.append(st["rho_std_mean"])

    taus = np.array(taus); events = np.array(events)
    ndied = int(events.sum()); ncens = M - ndied
    med_all = float(np.median(taus))
    med_died = float(np.median(taus[events == 1])) if ndied else float("nan")

    # eps-sensitivity: re-detect on stored traces
    eps_table = []
    for e in EPS_SWEEP:
        tt, ee = [], []
        for (t, rs) in traces:
            tau, ev = detect_death_ring(t, rs, e, dt_hold, T_max)
            tt.append(tau); ee.append(ev)
        tt = np.array(tt); ee = np.array(ee)
        eps_table.append({"eps": e, "n_died": int(ee.sum()),
                          "median_died": float(np.median(tt[ee == 1])) if ee.sum() else float("nan"),
                          "median_all": float(np.median(tt))})

    # histogram
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    dd = taus[events == 1]
    ax[0].hist(dd, bins=20, color="C0", alpha=0.85, edgecolor="k", lw=0.4)
    ax[0].axvline(med_died, color="C3", ls="--", label=f"median={med_died:.0f}")
    ax[0].set_xlabel("lifetime τ (deaths)"); ax[0].set_ylabel("count")
    ax[0].set_title(f"CP2 pilot  β={beta}, N={N}, r={r}  (M={M})\n"
                    f"{ndied} died / {ncens} censored")
    ax[0].legend(fontsize=8)

    es = np.array([row["eps"] for row in eps_table])
    em = np.array([row["median_died"] for row in eps_table])
    en = np.array([row["n_died"] for row in eps_table])
    ax[1].plot(es, em, "o-", color="C2")
    for e, m_, n_ in zip(es, em, en):
        ax[1].annotate(f"{n_}/{M}", (e, m_), fontsize=8, textcoords="offset points", xytext=(0, 6))
    ax[1].set_xlabel("ε_std (detector threshold)"); ax[1].set_ylabel("median lifetime (deaths)")
    ax[1].set_title("ε-sensitivity"); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    figpath = os.path.join(RESULTS, f"cp2_pilot_b{beta}_N{N}.png")
    fig.savefig(figpath, dpi=130); plt.close(fig)

    out = {"beta": beta, "N": N, "P": P, "M": M, "r": r, "T_max": T_max,
           "eps_primary": eps_primary, "dt_hold": dt_hold,
           "n_died": ndied, "n_censored": ncens, "median_all": med_all,
           "median_died": med_died, "chimera_rho_std": float(np.nanmean(plats)),
           "eps_sensitivity": eps_table, "figure": figpath}
    with open(os.path.join(RESULTS, f"cp2_pilot_b{beta}_N{N}.json"), "w") as f:
        json.dump(out, f, indent=2)

    print("=" * 64)
    print(f"CHECKPOINT 2 pilot  —  β={beta}  N={N}  r={r}  M={M}  T_max={T_max}")
    print("=" * 64)
    print(f"  died {ndied}/{M}  censored {ncens}/{M}  ({100*ncens/M:.0f}% censored)")
    print(f"  median τ (all) = {med_all:.0f}   median τ (deaths) = {med_died:.0f}")
    print(f"  chimera rho_std plateau = {np.nanmean(plats):.3f}  (sync floor ~0.018, ε={eps_primary})")
    print(f"\n  ε-sensitivity (re-detected on stored traces):")
    print(f"  {'ε_std':>6} {'n_died':>7} {'median(died)':>13} {'median(all)':>12}")
    for row in eps_table:
        print(f"  {row['eps']:>6.2f} {row['n_died']:>5}/{M:<2} "
              f"{row['median_died']:>13.0f} {row['median_all']:>12.0f}")
    print(f"\n  figure -> {figpath}")


if __name__ == "__main__":
    main()
