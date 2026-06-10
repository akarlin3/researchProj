"""CP2 calibration pass 2: high-res single-run heatmaps + numeric arc-count statistics
across thresholds and windows, to choose a robust classifier."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cp2_fields import (run_field, rho_profile, detect_death_ring, count_arcs_circular,
                        EPS_STD, DT_HOLD, T_MAX)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cp2_validation")


def adaptive_thresh(rho, lo=15, hi=85):
    a, b = np.percentile(rho, [lo, hi])
    return 0.5 * (a + b)


def arc_stats(field, t, P, t0, t1, thresh_mode, min_arc_frac=0.06):
    N = field.shape[1]
    min_w = max(2, int(round(min_arc_frac * N)))
    win = np.where((t >= t0) & (t <= t1))[0]
    ncoh, ninc, cfrac, single = [], [], [], []
    for i in win:
        rho = field[i]
        thr = adaptive_thresh(rho) if thresh_mode == "adaptive" else float(thresh_mode)
        coh = rho > thr
        a = count_arcs_circular(coh, min_w)
        b = count_arcs_circular(~coh, min_w)
        ncoh.append(a); ninc.append(b); cfrac.append(coh.mean())
        single.append(a == 1 and b == 1)
    return {
        "n": len(win),
        "median_ncoh": float(np.median(ncoh)) if win.size else float("nan"),
        "median_ninc": float(np.median(ninc)) if win.size else float("nan"),
        "frac_single_arc": float(np.mean(single)) if win.size else float("nan"),
        "mean_coh_frac": float(np.mean(cfrac)) if win.size else float("nan"),
    }


def analyse(beta, N, seed, tau_ens):
    t, rmean, rstd, th, P = run_field(beta, N, seed)
    tau, ev = detect_death_ring(t, rstd, EPS_STD, DT_HOLD, T_MAX)
    di = int(np.searchsorted(t, tau, side="right")) - 1
    field = np.array([rho_profile(th[i], P) for i in range(len(t))])
    print(f"\nseed {seed} beta {beta}: tau={tau:.1f} (ens {tau_ens:.1f}) P={P} n_dec={len(t)} "
          f"rho_std[alive max]={rstd[(t>=50)&(t<0.9*tau)].max():.3f}")
    windows = {"mature[0.5,0.9]tau": (0.5 * tau, 0.9 * tau),
               "last10%[0.9,1.0]tau": (0.9 * tau, t[di]),
               "midlife[0.3,0.6]tau": (0.3 * tau, 0.6 * tau)}
    for wlabel, (a, b) in windows.items():
        for tm in ("adaptive", 0.6, 0.7):
            s = arc_stats(field, t, P, a, b, tm)
            print(f"   {wlabel:20s} thr={str(tm):8s} n={s['n']:4d} "
                  f"med_ncoh={s['median_ncoh']:.1f} med_ninc={s['median_ninc']:.1f} "
                  f"single-arc%={100*s['frac_single_arc']:5.1f} cohfrac={s['mean_coh_frac']:.2f}")
    # high-res heatmap
    fig, ax = plt.subplots(figsize=(11, 4))
    im = ax.imshow(field.T, aspect="auto", origin="lower", cmap="magma",
                   extent=[t[0], t[di], 0, N], vmin=0, vmax=1)
    ax.axvline(0.9 * tau, color="cyan", ls=":", lw=1)
    ax.set_title(f"rho_k(t,k)  beta={beta} N={N} seed={seed} tau={tau:.0f}")
    ax.set_xlabel("t"); ax.set_ylabel("ring index k")
    plt.colorbar(im, ax=ax, fraction=0.03, label="rho_k")
    fig.tight_layout()
    p = os.path.join(OUT, f"hires_b{beta}_s{seed}.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    print(f"   [fig] {p}")


if __name__ == "__main__":
    # clean-looking and suspicious runs from pass 1
    analyse(0.130, 256, 20266309, 1405.5)   # row1 clean
    analyse(0.130, 256, 20266310, 792.0)    # row2 fine-striped (suspicious)
    analyse(0.110, 256, 20265109, 434.0)    # short
    analyse(0.110, 256, 20265110, 1439.5)   # long
