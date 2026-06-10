"""CP2 calibration: re-run a few N=256 seeds with field dumping, verify tau reproduces
the ensemble bit-for-bit, and visualise the spatial structure (space-time rho_k heatmap,
pre-death profiles, rho_k histogram) to choose the coherent/incoherent threshold."""
from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cp2_fields import run_field, rho_profile, detect_death_ring, EPS_STD, DT_HOLD, T_MAX

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "cp2_validation")
os.makedirs(OUT, exist_ok=True)


def ens_tau(beta, N):
    d = {}
    with open(os.path.join(ROOT, "anneal-hazard", "results", "ensemble_N256.csv")) as f:
        for r in csv.DictReader(f):
            if abs(float(r["beta"]) - beta) < 1e-9 and int(r["N"]) == N:
                d[int(r["seed"])] = (int(r["run_index"]), float(r["tau"]), int(r["event"]))
    return d


def calibrate(beta, N, n_show=3):
    tab = ens_tau(beta, N)
    seeds = sorted(tab, key=lambda s: tab[s][0])[:n_show]
    print(f"\n=== beta={beta} N={N} ===")
    fig, axes = plt.subplots(n_show, 3, figsize=(13, 3.2 * n_show))
    for row, seed in enumerate(seeds):
        ri, tau_ens, ev = tab[seed]
        t, rmean, rstd, th, P = run_field(beta, N, seed)
        tau_re, ev_re = detect_death_ring(t, rstd, EPS_STD, DT_HOLD, T_MAX)
        match = abs(tau_re - tau_ens) < 1e-6
        print(f"  run {ri} seed {seed}: tau_ens={tau_ens:.1f} tau_reint={tau_re:.1f} "
              f"event={ev_re} match={match} P={P} n_dec={len(t)}")
        di = int(np.searchsorted(t, tau_re, side="right")) - 1
        # full rho_k field
        field = np.array([rho_profile(th[i], P) for i in range(len(t))])  # (n_dec, N)
        # last-10% window
        t0 = tau_re * 0.9
        win = (t >= t0) & (t <= t[di])

        ax = axes[row, 0]
        im = ax.imshow(field.T, aspect="auto", origin="lower", cmap="magma",
                       extent=[t[0], t[di], 0, N], vmin=0, vmax=1)
        ax.axvline(tau_re, color="cyan", lw=1)
        ax.axvline(t0, color="white", ls=":", lw=1)
        ax.set_title(f"rho_k(t,k)  seed{seed} tau={tau_re:.0f}", fontsize=8)
        ax.set_xlabel("t"); ax.set_ylabel("k")
        plt.colorbar(im, ax=ax, fraction=0.04)

        ax = axes[row, 1]
        idxs = np.where(win)[0]
        for i in idxs[:: max(1, len(idxs) // 6)]:
            ax.plot(field[i], color="gray", alpha=0.5, lw=0.8)
        ax.plot(field[idxs].mean(0), color="crimson", lw=2, label="mean last-10%")
        ax.set_title("pre-death rho_k profiles", fontsize=8)
        ax.set_xlabel("k"); ax.set_ylim(0, 1.02); ax.legend(fontsize=7)

        ax = axes[row, 2]
        alive = (t >= 50) & (t < t0)
        ax.hist(field[alive].ravel(), bins=40, range=(0, 1), density=True, alpha=0.5,
                label="alive", color="steelblue")
        if win.any():
            ax.hist(field[idxs].ravel(), bins=40, range=(0, 1), density=True, alpha=0.5,
                    label="last-10%", color="crimson")
        ax.set_title("rho_k histogram", fontsize=8)
        ax.set_xlabel("rho_k"); ax.legend(fontsize=7)
    fig.tight_layout()
    p = os.path.join(OUT, f"calib_b{beta}_N{N}.png")
    fig.savefig(p, dpi=130)
    print(f"  [fig] {p}")


if __name__ == "__main__":
    calibrate(0.130, 256, n_show=3)
    calibrate(0.110, 256, n_show=3)
