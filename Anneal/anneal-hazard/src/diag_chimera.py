"""Diagnostic: visualise clean chimeras at small beta and measure their lifetimes.
Plots r1,r2 over the full T_max for a list of (beta, A) at N, prints tau per run.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config_io import load_config, make_params
from .model import make_initial_conditions, integrate
from .detector import detect_death

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    N = 256
    cases = [(0.02, 0.20), (0.05, 0.20), (0.07, 0.20),
             (0.08, 0.20), (0.06, 0.30), (0.07, 0.30)]
    seed = cfg["seed_base"]
    Tmax = cfg["T_max"]

    fig, axes = plt.subplots(len(cases), 1, figsize=(11, 2.0 * len(cases)), sharex=True)
    for ax, (beta, A) in zip(axes, cases):
        p = make_params(A, beta)
        rng = np.random.default_rng(seed)
        th0 = make_initial_conditions(N, cfg["ic"], rng)
        res = integrate(theta0=th0, p=p, N=N, dt=cfg["dt"], T_max=Tmax,
                        decimate=cfg["output"]["decimate"], seed=seed)
        tau, ev = detect_death(res.t, res.r1, res.r2, cfg["eps"], cfg["dt_hold"], Tmax)
        m = (res.t >= 100) & (res.t < (tau if ev else Tmax))
        rinc = np.minimum(res.r1, res.r2)[m]
        ax.plot(res.t, res.r1, lw=0.5, color="C0")
        ax.plot(res.t, res.r2, lw=0.5, color="C3")
        if ev:
            ax.axvline(tau, color="k", ls="--", lw=0.8)
        ax.set_ylim(0, 1.03)
        ax.set_ylabel(f"β={beta}\nA={A}", fontsize=8)
        lbl = f"r_incoh={rinc.mean():.3f}±{rinc.std():.3f}   τ={tau:.0f}{'' if ev else ' (CENSORED)'}"
        ax.text(0.99, 0.08, lbl, transform=ax.transAxes, ha="right", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))
        print(f"beta={beta:.2f} A={A:.2f} N={N}: r_incoh={rinc.mean():.3f}±{rinc.std():.3f} "
              f"tau={tau:.1f} event={ev}")
    axes[0].set_title(f"Clean-chimera search  N={N}, seed={seed}  (r1=C0, r2=C3)")
    axes[-1].set_xlabel("t")
    fig.tight_layout()
    out = os.path.join(RESULTS, "diag_chimera.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    main()
