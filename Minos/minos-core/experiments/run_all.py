"""Minos-Core experiment driver.

Seeded, config-driven sweep over the calibration knob ``tau`` and the shift knob
``delta``. Reproduces all four checkpoint gates from a clean seed, prints every
number that RESULTS.md cites, and writes the four publication figures as vector PDFs.

Run from the project root:  ``python experiments/run_all.py``
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from minos.config import MinosConfig  # noqa: E402
from minos.diagnostics import central_interval_coverage, ece  # noqa: E402
from minos.gate import (  # noqa: E402
    detection_auc,
    expected_utility_gated,
    gate_signal,
    gate_threshold,
    votg,
)
from minos.generative import make_population, realise  # noqa: E402
from minos.seeding import make_rng  # noqa: E402
from minos.voi import (  # noqa: E402
    evpi,
    expected_utility,
    posterior_eu_curve,
    value_of_error_bar,
    voc,
)

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(os.path.dirname(HERE), "figures")

RUN_CFG = MinosConfig(n_voxels=1_000_000)
TAUS = np.round(np.arange(0.5, 2.0001, 0.05), 3)
DELTAS = np.round(np.arange(0.0, 2.5001, 0.25), 3)
LEVELS = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95])
DELTA_TEST = 1.5

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def hr(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    cfg = RUN_CFG
    base = make_population(cfg, make_rng(cfg.seed))
    out = {}

    hr("CONFIG")
    print(f"seed={cfg.seed}  n_voxels={cfg.n_voxels}")
    print(f"t1={cfg.t1} t2={cfg.t2} k_under={cfg.k_under} k_over={cfg.k_over} s={cfg.s}")
    print(f"prior weights={cfg.mix_weights} means={cfg.mix_means} stds={cfg.mix_stds}")
    print(f"shift alpha={cfg.alpha} beta={cfg.beta}  gate q={cfg.q_gate} "
          f"g*={gate_threshold(cfg):.4f}  delta_test={DELTA_TEST}")

    # ---- GATE 1: EU ordering + degenerate EVPI ----------------------------------
    hr("GATE 1 — decision core")
    eu_point = expected_utility("point", base, cfg)
    eu_post = expected_utility("posterior", base, cfg, tau=1.0)
    eu_oracle = expected_utility("oracle", base, cfg)
    evpi1 = evpi(base, cfg, tau=1.0)
    print(f"EU_point    = {eu_point:+.6f}")
    print(f"EU_posterior= {eu_post:+.6f}")
    print(f"EU_oracle   = {eu_oracle:+.6f}")
    print(f"EVPI-analog = {evpi1:+.6f}")
    assert eu_oracle >= eu_post >= eu_point, "EU ordering violated"
    cfg_deg = cfg.replace(n_voxels=200_000, s=1e-4)
    base_deg = make_population(cfg_deg, make_rng(cfg_deg.seed))
    evpi_deg = evpi(base_deg, cfg_deg, tau=1.0)
    print(f"EVPI-analog (degenerate s=1e-4) = {evpi_deg:.3e}")
    assert evpi_deg < 1e-2, "EVPI did not vanish at point-mass limit"
    print("GATE 1 PASS: EU_oracle >= EU_posterior >= EU_point ; EVPI -> 0 at point mass")
    out.update(eu_point=eu_point, eu_post=eu_post, eu_oracle=eu_oracle,
               evpi=evpi1, evpi_deg=evpi_deg)

    # ---- GATE 2: VoC + value-of-error-bar + calibration -------------------------
    hr("GATE 2 — Value of Calibration")
    voc_curve = np.array([voc(base, cfg, t) for t in TAUS])
    post_eu_curve = posterior_eu_curve(base, cfg, TAUS)
    veb = value_of_error_bar(base, cfg)
    i_one = int(np.argmin(np.abs(TAUS - 1.0)))
    print(f"VoC(tau=1)        = {voc_curve[i_one]:+.3e}")
    print(f"VoC(tau=0.5)      = {voc(base, cfg, 0.5):+.6f}")
    print(f"VoC(tau=2.0)      = {voc(base, cfg, 2.0):+.6f}")
    print(f"argmin VoC at tau = {TAUS[int(np.argmin(voc_curve))]:.2f}")
    print(f"value of error bar (EU_post(1) - EU_point) = {veb:+.6f}")
    assert abs(voc_curve[i_one]) < 1e-6
    assert int(np.argmin(voc_curve)) == i_one, "VoC not minimised at tau=1"
    assert voc_curve.min() >= -5e-5
    assert veb > 0
    cov = {lvl: central_interval_coverage(base, cfg, lvl, tau=1.0) for lvl in (0.5, 0.8, 0.9)}
    print("central-interval coverage @ tau=1 :", {k: round(v, 4) for k, v in cov.items()})
    ece_cal = ece(base, cfg, LEVELS, tau=1.0)
    ece_over = ece(base, cfg, LEVELS, tau=0.6)
    ece_under = ece(base, cfg, LEVELS, tau=1.5)
    print(f"ECE tau=1.0 -> {ece_cal:.5f} | tau=0.6 -> {ece_over:.5f} | tau=1.5 -> {ece_under:.5f}")
    print("GATE 2 PASS: VoC(1)=0 minimum, VoC>0 away from 1, value-of-error-bar>0")
    out.update(voc_05=float(voc(base, cfg, 0.5)), voc_20=float(voc(base, cfg, 2.0)),
               veb=veb, ece_cal=ece_cal, ece_over=ece_over, ece_under=ece_under, cov=cov)

    # ---- GATE 3: shift + trust-gate + VoTG --------------------------------------
    hr("GATE 3 — trust-gate")
    votg_curve = np.array([votg(base, cfg, delta=d) for d in DELTAS])
    auc_curve = np.array([detection_auc(base, cfg, d) for d in DELTAS])
    votg0 = votg(base, cfg, delta=0.0)
    votg_t = votg(base, cfg, delta=DELTA_TEST)
    auc0 = detection_auc(base, cfg, 0.0)
    auc_t = detection_auc(base, cfg, DELTA_TEST)
    eu_post_shift = expected_utility("posterior", base, cfg, delta=DELTA_TEST, shift=True)
    eu_gated_shift = expected_utility_gated(base, cfg, delta=DELTA_TEST, shift=True)
    print(f"VoTG(delta=0)          = {votg0:+.6f}")
    print(f"VoTG(delta={DELTA_TEST})        = {votg_t:+.6f}")
    print(f"detection AUC(delta=0) = {auc0:.4f}")
    print(f"detection AUC(delta={DELTA_TEST}) = {auc_t:.4f}")
    print(f"posterior regret @shift= {-eu_post_shift:+.6f}")
    print(f"gated regret     @shift= {-eu_gated_shift:+.6f}")
    assert abs(votg0) < 0.02
    assert votg_t > 0.05
    assert (-eu_gated_shift) < (-eu_post_shift), "gate did not reduce regret under shift"
    assert auc_t > 0.7
    assert abs(auc0 - 0.5) < 0.01
    print("GATE 3 PASS: VoTG(0)~0, VoTG(delta)>0, gated regret<posterior regret, AUC>0.7")
    out.update(votg0=votg0, votg_t=votg_t, auc0=auc0, auc_t=auc_t,
               regret_post_shift=-eu_post_shift, regret_gated_shift=-eu_gated_shift)

    # ---- Figures ----------------------------------------------------------------
    hr("FIGURES")
    _fig_a(cfg, base, post_eu_curve)
    _fig_b(voc_curve, evpi1)
    _fig_c(cfg, base, votg_curve, auc_curve)
    _fig_d(cfg, base)
    for name in ("fig_a_regret_vs_tau", "fig_b_voc_evpi", "fig_c_gate_roc_votg",
                 "fig_d_utility_bars"):
        print("wrote", os.path.join("figures", name + ".pdf"))

    hr("ALL GATES PASS")
    return out


def _fig_a(cfg, base, post_eu_curve):
    point_regret = -expected_utility("point", base, cfg)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(TAUS, -post_eu_curve, lw=2, color="#1f4e79", label="posterior policy")
    ax.axhline(point_regret, ls="--", color="#a0a0a0", label="point policy (ignores error bar)")
    ax.axvline(1.0, ls=":", color="#cc5500", alpha=0.8, label="calibrated (tau=1)")
    ax.set_xlabel("calibration quality  tau  (1 = calibrated)")
    ax.set_ylabel("decision regret  =  -E[U]")
    ax.set_title("(a) Decision regret vs calibration quality")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_a_regret_vs_tau.pdf"))
    plt.close(fig)


def _fig_b(voc_curve, evpi1):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(TAUS, voc_curve, lw=2, color="#1f4e79", label="VoC(tau)")
    ax.axhline(0.0, color="#404040", lw=0.8)
    ax.axvline(1.0, ls=":", color="#cc5500", alpha=0.8, label="calibrated (tau=1)")
    ax.axhline(evpi1, ls="--", color="#8000a0",
               label=f"EVPI-analog @tau=1 = {evpi1:.3f}")
    ax.set_xlabel("calibration quality  tau")
    ax.set_ylabel("value of calibration  (utility)")
    ax.set_title("(b) Value of Calibration and the EVPI-analog")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_b_voc_evpi.pdf"))
    plt.close(fig)


def _fig_c(cfg, base, votg_curve, auc_curve):
    n = base.theta.shape[0]
    mask = np.arange(n) < (n // 2)
    _, w = realise(base, cfg, delta=DELTA_TEST, shift=mask)
    g = gate_signal(w, cfg)
    thr = np.quantile(g, np.linspace(0, 1, 200))
    tpr = np.array([(g[mask] > t).mean() for t in thr])
    fpr = np.array([(g[~mask] > t).mean() for t in thr])

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(10, 4))
    axl.plot(fpr, tpr, lw=2, color="#1f4e79")
    axl.plot([0, 1], [0, 1], ls="--", color="#a0a0a0")
    axl.scatter([(g[~mask] > gate_threshold(cfg)).mean()],
                [(g[mask] > gate_threshold(cfg)).mean()],
                color="#cc5500", zorder=5, label=f"operating point (q={cfg.q_gate})")
    axl.set_xlabel("false-positive rate")
    axl.set_ylabel("true-positive rate")
    axl.set_title(f"(c1) Trust-gate ROC at delta={DELTA_TEST}\nAUC={detection_auc(base, cfg, DELTA_TEST):.3f}")
    axl.legend(frameon=False, fontsize=9, loc="lower right")

    axr.axhline(0.0, color="#404040", lw=0.8)
    axr.plot(DELTAS, votg_curve, lw=2, marker="o", ms=3, color="#1f4e79")
    axr.set_xlabel("distribution shift  delta")
    axr.set_ylabel("VoTG(delta)  (utility recovered)")
    axr.set_title("(c2) Value of the Trust-Gate vs shift")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_c_gate_roc_votg.pdf"))
    plt.close(fig)


def _fig_d(cfg, base):
    conds = [("in-distribution", 0.0, False), (f"shift delta={DELTA_TEST}", DELTA_TEST, True)]
    policies = ["point", "posterior", "gated", "oracle"]
    colors = ["#a0a0a0", "#1f4e79", "#2e8b57", "#cc5500"]
    vals = np.zeros((len(conds), len(policies)))
    for i, (_, d, sh) in enumerate(conds):
        vals[i, 0] = expected_utility("point", base, cfg, delta=d, shift=sh)
        vals[i, 1] = expected_utility("posterior", base, cfg, delta=d, shift=sh)
        vals[i, 2] = expected_utility_gated(base, cfg, delta=d, shift=sh)
        vals[i, 3] = expected_utility("oracle", base, cfg)

    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(conds))
    width = 0.2
    for j, pol in enumerate(policies):
        ax.bar(x + (j - 1.5) * width, vals[:, j], width, label=pol, color=colors[j])
    ax.set_xticks(x)
    ax.set_xticklabels([c[0] for c in conds])
    ax.set_ylabel("expected utility  E[U]")
    ax.set_title("(d) Expected utility by policy")
    ax.legend(frameon=False, fontsize=9, ncol=4, loc="lower center")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_d_utility_bars.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    summary = main()
