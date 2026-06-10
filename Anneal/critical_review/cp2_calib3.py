"""CP2 calibration pass 3 — ground-truth the classifier.

Two robust, threshold-light discriminators:
  (1) spatial CONTRAST gate: rho_std(t) (exact spatial std of rho_k). A coherent/incoherent
      coexistence (a live chimera) has rho_std ~ 0.15-0.28; a dissolved near-uniform 'limbo'
      state has rho_std ~ 0.04-0.08; death is rho_std<0.04. -> 'structured' if rho_std>c_struct.
  (2) HEAD COUNT via the spatial Fourier spectrum of rho_k (rotation/drift invariant): the
      dominant non-zero wavenumber m* of (rho_k - mean). A single coherent arc -> m*=1; a
      q-headed (multi-arc) chimera -> m*=q. Threshold-free.
We cross-check m* against a contrast-gated contiguous-arc count.

Plots rho_std(t) + individual snapshots (with m*) for 4 known seeds to fix c_struct & confirm
the head count behaves (310 should read degenerate pre-death; 110/309/109 structured)."""
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


def head_count(rho, max_m=12):
    """Dominant non-zero spatial wavenumber of rho_k (drift-invariant)."""
    f = np.fft.rfft(rho - rho.mean())
    p = np.abs(f) ** 2
    if len(p) <= 1:
        return 0
    m = 1 + int(np.argmax(p[1:max_m + 1]))
    return m


def arc_count_adaptive(rho, min_arc_frac=0.08):
    N = len(rho)
    a, b = np.percentile(rho, [15, 85])
    thr = 0.5 * (a + b)
    coh = rho > thr
    minw = max(2, int(round(min_arc_frac * N)))
    return count_arcs_circular(coh, minw), count_arcs_circular(~coh, minw)


SEEDS = [(0.130, 20266309, 1405.5), (0.130, 20266310, 792.0),
         (0.110, 20265109, 434.0), (0.110, 20265110, 1439.5)]

fig, axes = plt.subplots(len(SEEDS), 4, figsize=(15, 3.0 * len(SEEDS)))
for row, (beta, seed, tau_ens) in enumerate(SEEDS):
    t, rmean, rstd, th, P = run_field(beta, 256, seed)
    tau, ev = detect_death_ring(t, rstd, EPS_STD, DT_HOLD, T_MAX)
    di = int(np.searchsorted(t, tau, side="right")) - 1
    field = np.array([rho_profile(th[i], P) for i in range(len(t))])

    # rho_std trace
    ax = axes[row, 0]
    ax.plot(t[:di + 1], rstd[:di + 1], lw=0.8, color="navy")
    for c in (0.04, 0.08, 0.10):
        ax.axhline(c, ls=":", lw=0.8, color="gray")
    ax.axvline(0.9 * tau, color="orange", lw=1)
    ax.set_title(f"b{beta} s{seed} tau={tau:.0f}\nrho_std(t)", fontsize=8)
    ax.set_xlabel("t"); ax.set_ylim(0, 0.32)

    # snapshots: one mid-life, one early-last10%, one just-before-death
    picks = {"midlife": int(0.5 * di), "last10-start": int(0.9 * di), "pre-death": di - 2}
    for col, (lab, i) in enumerate(picks.items(), start=1):
        i = max(0, min(i, di))
        rho = field[i]
        m = head_count(rho)
        nc, ni = arc_count_adaptive(rho)
        ax = axes[row, col]
        ax.plot(rho, color="crimson", lw=1)
        ax.set_ylim(0, 1.02)
        ax.set_title(f"{lab} t={t[i]:.0f}\nrho_std={rstd[i]:.3f} m*={m} arcs={nc}/{ni}",
                     fontsize=8)
        ax.set_xlabel("k")
fig.tight_layout()
p = os.path.join(OUT, "calib3_snapshots.png")
fig.savefig(p, dpi=130)
print(f"[fig] {p}")

# numeric: structured-fraction of last-10% and head-count distribution there
print(f"\n{'seed':>10} {'beta':>5} {'tau':>7} | structured-frac & head-count over last-10% "
      f"(gate rho_std>0.08)")
for beta, seed, tau_ens in SEEDS:
    t, rmean, rstd, th, P = run_field(beta, 256, seed)
    tau, ev = detect_death_ring(t, rstd, EPS_STD, DT_HOLD, T_MAX)
    di = int(np.searchsorted(t, tau, side="right")) - 1
    win = np.where((t >= 0.9 * tau) & (t <= t[di]))[0]
    field = np.array([rho_profile(th[i], P) for i in win])
    struct = rstd[win] > 0.08
    heads = [head_count(field[j]) for j in range(len(win)) if struct[j]]
    sf = float(struct.mean())
    h1 = float(np.mean([h == 1 for h in heads])) if heads else 0.0
    from collections import Counter
    hc = Counter(heads)
    print(f"{seed:>10} {beta:>5} {tau:>7.0f} | struct_frac={sf:.2f}  "
          f"single-head|struct={h1:.2f}  head_dist={dict(sorted(hc.items()))}")
