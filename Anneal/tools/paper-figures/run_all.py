"""
Driver — regenerate Figures 1 & 7 and the report, deterministically, from
figures.config.json.

Steps:
  CP0  cp0_select.py            audit + pin-checked Fig.1 run selection
  F1   fig1_trace.mjs (node)    regenerate the example trace from the seed
       fig1.py                  plot Figure 1 + caption
  F7   fig7_curves.py           SN/Hopf series+numeric, homoclinic bisection (slow)
       fig7.py                  plot Figure 7 + caption
  APP  figA.py, figB.py         appendix figures (beta=0.10 slice; noise test)
  F10  fig10_ring.py            ring summary (anneal-hazard/results/)
  F12  fig12_crosssystem.py     cross-system comparison (committed records)
  REV  fig_corner_map.py, fig_mech_probe.py, fig_betac.py,
       fig_depth_aging.py, fig_xlabel.py
                                revision appendix figures (committed data)
  RPT  make_report.py           FIGURES_REPORT.md

Usage:
  python3 tools/paper-figures/run_all.py [--skip-curves]
    --skip-curves   reuse cached paper_figures/fig7_curves.json (skip the ~4-5
                    min homoclinic trace) when only plots/report changed.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "paper_figures")

# Paper figure <- source result figure (stem, no extension). Figures 1, 7, A and
# B are generated directly into paper_figures/ by their own scripts; the
# single-panel figures below are produced by the per-tool analysis pipelines and
# assembled here (previously a manual copy step). The two-panel composites
# (fig3, fig5, fig8) are built by compose.py — see compose.PAIR_MAP. All source
# titles are CP-codename-free (publication style).
FIG_MAP = {
    "fig2": "absorption_results/tau_old_vs_new",
    "fig4": "absorption_results/geometric_p",
    "fig6": "transient_results/transient_decider",
    "fig9": "reduced_results/cp4_scatter",
}


def run(cmd, **kw):
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=ROOT, **kw)


def assemble_figs():
    """Copy the per-tool source figures into paper_figures/figN.{pdf,png}.
    Scripted (not manual) so stale copies cannot drift from their sources."""
    print("\n# assembling paper figures from source pipelines")
    for fig, stem in FIG_MAP.items():
        for ext in ("pdf", "png"):
            src = os.path.join(ROOT, f"{stem}.{ext}")
            dst = os.path.join(OUT, f"{fig}.{ext}")
            if os.path.exists(src):
                shutil.copyfile(src, dst)
                print(f"  {fig}.{ext} <- {stem}.{ext}")
            else:
                print(f"  [warn] missing source {stem}.{ext} (run its analysis pipeline first)")


def main():
    skip_curves = "--skip-curves" in sys.argv
    os.makedirs(OUT, exist_ok=True)

    run([sys.executable, os.path.join(HERE, "cp0_select.py")])

    # Figure 1: deterministic trace (node) → plot
    with open(os.path.join(OUT, "fig1_trace.json"), "w") as f:
        run(["node", os.path.join(HERE, "fig1_trace.mjs")], stdout=f)
    run([sys.executable, os.path.join(HERE, "fig1.py")])

    # Figure 7: curves (slow) → plot
    cache = os.path.join(OUT, "fig7_curves.json")
    if skip_curves and os.path.exists(cache):
        print("\n[skip-curves] reusing cached fig7_curves.json")
    else:
        run([sys.executable, os.path.join(HERE, "fig7_curves.py")])
    run([sys.executable, os.path.join(HERE, "fig7.py")])

    # Appendix figures: beta=0.10 slice (A) and finite-size noise test (B).
    run([sys.executable, os.path.join(HERE, "figA.py")])
    run([sys.executable, os.path.join(HERE, "figB.py")])

    # Ring summary figure (reads anneal-hazard/results/).
    run([sys.executable, os.path.join(HERE, "fig10_ring.py")])

    # Cross-system comparison (Fig. 12): mean-field plateau vs ring scaling and
    # the Weibull shapes of both campaigns, from committed run records and fits.
    run([sys.executable, os.path.join(OUT, "fig12_crosssystem.py")])

    # Revision appendix figures: corner-generality map, mechanism probe, and
    # ring existence-boundary measurement. Each plots from committed data; the
    # heavy regeneration commands are documented in the scripts' docstrings
    # (tools/reduced-ode/corner_map_data.py, tools/noise-test/mech_probe.py +
    # tools/manifold-probe/ws_decomposition.mjs, anneal-hazard/analysis/run_betac.py).
    run([sys.executable, os.path.join(OUT, "fig_corner_map.py")])
    run([sys.executable, os.path.join(OUT, "fig_mech_probe.py")])
    run([sys.executable, os.path.join(OUT, "fig_betac.py")])

    # Contribution-raising items: depth-onset reanalysis (committed corner data
    # via tools/reduced-ode/depth_aging_fit.py) and the ring graze/absorption
    # null (committed re-run records via tools/xlabel/run_xlabel.py).
    run([sys.executable, os.path.join(OUT, "fig_depth_aging.py")])
    run([sys.executable, os.path.join(OUT, "fig_xlabel.py")])

    # Assemble single-panel Figures 2, 4, 6, 9 from their source pipelines.
    assemble_figs()

    # Compose the two-panel Figures 3, 5, 8 from their committed source panels.
    sys.path.insert(0, HERE)
    import compose

    compose.compose_all()

    run([sys.executable, os.path.join(HERE, "make_report.py")])
    print("\nDone. Artifacts in paper_figures/.")


if __name__ == "__main__":
    main()
