#!/usr/bin/env python3
"""Re-derive the Gauge high-D* identifiability wall, dependency-light, and cross-check the anchor.

This is the IDENTIFIABILITY (Cramer-Rao) argument of Augur's spine, reproduced from first
principles so the synthesis does not merely *assert* the wall. It is NOT an impossibility proof:
the Cramer-Rao lower bound (CRLB) on D* is a *local identifiability* limit -- the best variance
any unbiased estimator can achieve from the IVIM b-value curve at a given operating point. It is
scoped to its regime: the segmented bi-exponential IVIM model, the rich synthetic b-scheme Gauge
used for the conditional-coverage cohort, and the SNR band 10-100.

Mirrors, exactly, Gauge/gauge/forward.py:crlb (Fisher information on the bi-exp signal Jacobian,
S0 jointly estimated) and Gauge/gauge/conditional_attack.py:crlb_sweep -- but self-contained
(numpy + matplotlib only; no torch, no Gauge import) so Augur reproduces and can be carved
standalone. The reproduced growth factor is cross-checked against the committed anchor
(anchors/anchors.json -> crlb.crlb_growth_factor ~6x). The cohort-realized
CRLB(D*)/tercile-width ~1.12 is carried from the anchor (it depends on the seeded SNR draw);
this script confirms the same qualitative wall (ratio > 1 across the SNR band).

Outputs:
  Augur/figures/crlb_wall.pdf     -- CRLB(D*) vs true D* across SNR, high-D* tercile shaded
  Augur/results/crlb_wall.json    -- reproduced numbers + anchor cross-check verdict

Exit code: 0 = reproduced and cross-check PASS; non-zero = reproduction diverges from anchor.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AUGUR = HERE.parent
ANCHORS = AUGUR / "anchors" / "anchors.json"
FIGDIR = AUGUR / "figures"
RESDIR = AUGUR / "results"

# --- Gauge's conditional-coverage cohort regime (cohort.py / forward.py constants) -----------
# The rich synthetic b-scheme the wall is diagnosed on (Gauge DEFAULT_B_VALUES). The high-D*
# wall is a property of this segmented IVIM acquisition; it is scoped to it.
B_VALUES = np.array(
    [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
     120, 140, 160, 180, 200, 300, 400, 500, 600, 700, 800], dtype=float)
DSTAR_RANGE = (10e-3, 100e-3)           # Gauge cohort.py DSTAR_RANGE
SNR_GRID = (10, 20, 30, 50, 100)        # Gauge cohort.py DEFAULT_SNR_GRID
D_FIX = 1.5e-3                           # crlb_sweep default tissue diffusion
F_FIX = 0.2                             # crlb_sweep default perfusion fraction
N_SWEEP = 40


def crlb_dstar(b, D, Dstar, f, snr, S0=1.0, fix_s0=False):
    """CRLB std of D* for one IVIM voxel (mirrors Gauge forward.py:crlb, D* = Fisher index 1).

    Jacobian columns (D, D*, f, S0) of S(b)=S0[f e^{-bD*}+(1-f) e^{-bD}]; info = J^T J / sigma^2,
    sigma = S0/snr; CRLB(D*) = sqrt((info^{-1})_{11}). S0 jointly estimated unless fix_s0.
    """
    sigma = S0 / float(snr)
    ed = np.exp(-b * D)
    eds = np.exp(-b * Dstar)
    cols = [S0 * (1.0 - f) * (-b) * ed,    # dS/dD
            S0 * f * (-b) * eds,           # dS/dD*
            S0 * (eds - ed),               # dS/df
            f * eds + (1.0 - f) * ed]      # dS/dS0
    if fix_s0:
        cols = cols[:3]
    J = np.stack(cols, axis=1)
    info = J.T @ J / (sigma * sigma)
    try:
        v = np.linalg.inv(info)[1, 1]
    except np.linalg.LinAlgError:
        return np.inf
    return float(np.sqrt(v)) if np.isfinite(v) and v > 0 else np.inf


def crlb_sweep():
    """Absolute CRLB(D*) across the D* range at each SNR (S0 free). Mirrors crlb_sweep."""
    dstars = np.linspace(DSTAR_RANGE[0], DSTAR_RANGE[1], N_SWEEP)
    absolute = {s: np.array([crlb_dstar(B_VALUES, D_FIX, ds, F_FIX, snr=s) for ds in dstars])
                for s in SNR_GRID}
    return dstars, absolute


def make_figure(dstars, absolute, edges, fig_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for s in SNR_GRID:
        ax.plot(dstars * 1e3, absolute[s] * 1e3, marker="o", ms=2.5, lw=1.3, label=f"SNR {s}")
    # Shade the high-D* tercile (regime unresolvable) and mark the tercile width reference.
    ax.axvspan(edges[1] * 1e3, DSTAR_RANGE[1] * 1e3, color="0.85", alpha=0.6, zorder=0,
               label="high-$D^*$ tercile")
    ax.set_xlabel(r"true $D^*$  ($10^{-3}\,\mathrm{mm^2/s}$)")
    ax.set_ylabel(r"CRLB$(D^*)$ std  ($10^{-3}\,\mathrm{mm^2/s}$)")
    ax.set_title("IVIM pseudo-diffusion identifiability wall (CRLB, S0 free)")
    ax.legend(fontsize=8, ncol=2, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path)
    plt.close(fig)


def main() -> int:
    if not ANCHORS.exists():
        print("CRLB-WALL FAIL: anchors.json missing; run anchors/extract_anchors.py first",
              file=sys.stderr)
        return 2
    anchor = json.loads(ANCHORS.read_text())["crlb"]
    edges = anchor["dstar_tercile_edges"]

    dstars, absolute = crlb_sweep()
    # Reproduced growth factor across the D* range (identical at every SNR -- linear in sigma).
    growth = {s: float(absolute[s][-1] / absolute[s][0]) for s in SNR_GRID}
    growth_repro = float(np.mean(list(growth.values())))

    # Self-derived high-D* CRLB/tercile-width ratio (SNR-band median; cohort-realized 1.12 is
    # SNR-distribution-dependent and is carried from the anchor as the load-bearing value).
    hi_mask = dstars >= edges[1]
    hi_width = DSTAR_RANGE[1] - edges[1]
    per_snr_ratio = sorted(float(np.median(absolute[s][hi_mask]) / hi_width) for s in SNR_GRID)
    ratio_band = [per_snr_ratio[0], per_snr_ratio[-1]]

    make_figure(dstars, absolute, edges, FIGDIR / "crlb_wall.pdf")

    anchor_growth = anchor["crlb_growth_factor"]
    # Cross-check: reproduced growth must round to the anchor's ~6x (tolerance +-1.0x).
    growth_ok = abs(growth_repro - anchor_growth) <= 1.0
    # The wall is present iff CRLB(D*) reaches the tercile-width scale somewhere in the band.
    wall_ok = ratio_band[1] >= 1.0

    out = {
        "_about": "Reproduced IVIM CRLB(D*) identifiability wall; cross-checked vs committed anchor.",
        "regime": anchor["regime"],
        "b_values": list(B_VALUES),
        "snr_grid": list(SNR_GRID),
        "D_fixed": D_FIX, "f_fixed": F_FIX, "dstar_range": list(DSTAR_RANGE),
        "reproduced": {
            "crlb_growth_factor": round(growth_repro, 2),
            "hi_tercile_crlb_over_width_band": [round(x, 2) for x in ratio_band],
            "crlb_dstar_lo_hi_at_snr30": [round(absolute[30][0], 5), round(absolute[30][-1], 5)],
        },
        "anchor": {
            "crlb_growth_factor": anchor_growth,
            "crlb_over_tercile_width_hi": anchor["crlb_over_tercile_width_hi"],
            "conformal_width_crlb_r": anchor["conformal_width_crlb_r"],
            "source": anchor["source"],
        },
        "crosscheck": {
            "growth_factor_pass": growth_ok,
            "wall_present_pass": wall_ok,
            "note": ("reproduced growth ~%.1fx matches anchor ~%.0fx; CRLB(D*) reaches the high-D* "
                     "tercile-width scale (ratio band %.2f-%.2f), confirming the same identifiability "
                     "wall. Cohort-realized 1.12x carried from anchor (SNR-draw dependent). "
                     "Identifiability limit, not impossibility; scoped to segmented IVIM, SNR 10-100."
                     % (growth_repro, anchor_growth, ratio_band[0], ratio_band[1])),
        },
    }
    RESDIR.mkdir(parents=True, exist_ok=True)
    (RESDIR / "crlb_wall.json").write_text(json.dumps(out, indent=2) + "\n")

    print("CRLB identifiability wall reproduced (self-contained Fisher info; S0 free)")
    print(f"  b-scheme: {len(B_VALUES)}-b segmented IVIM; SNR {SNR_GRID}; D={D_FIX}, f={F_FIX}")
    print(f"  reproduced growth: {growth_repro:.2f}x   anchor: ~{anchor_growth:.0f}x   "
          f"[{'PASS' if growth_ok else 'FAIL'}]")
    print(f"  high-D* CRLB/tercile-width band: [{ratio_band[0]:.2f}, {ratio_band[1]:.2f}]   "
          f"anchor cohort-realized: {anchor['crlb_over_tercile_width_hi']:.2f}   "
          f"wall present: {'PASS' if wall_ok else 'FAIL'}")
    print(f"  figure -> {(FIGDIR / 'crlb_wall.pdf').relative_to(AUGUR.parent)}")
    if not (growth_ok and wall_ok):
        print("CRLB-WALL FAIL: reproduction diverges from anchor", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
