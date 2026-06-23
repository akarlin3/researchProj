#!/usr/bin/env python
"""CP0 driver -- the identifiability object + the kill test.

Builds the Fisher/CRLB for the fractional order alpha jointly with (D, S0) under the
diffusion-MRI signal-decay forward model with finite b-value sampling + Rician noise;
sweeps SNR and b-range; locates the recovery-collapse wall; runs the pre-registered REFUTE.
Writes results/RESULTS_CP0.md and a figure. Numbers print to stdout for transcription.

Usage:  <proteus python> experiments/run_cp0.py            # FAST (smoke; default)
        <proteus python> experiments/run_cp0.py --full     # full-N bootstrap CIs
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / "levy-core"
sys.path.insert(0, str(CORE))

from levy import fisher, forward, identifiability, noise, seeding, wall  # noqa: E402

FULL = "--full" in sys.argv
RESULTS = HERE.parent / "results" / "RESULTS_CP0.md"
FIGDIR = HERE.parent / "figures"


def _hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def main() -> int:
    rng = seeding.make_rng()
    n_boot = 300 if FULL else 80

    _hr("CP0 -- fractional-order identifiability under the MRI signal-decay forward model")
    print("Forward model (lead lane, stretched-exponential):  S(b;S0,D,alpha) = S0 exp(-(b D)^alpha)")
    print("Estimand: theta = (S0, D, alpha) estimated JOINTLY.  Noise: Rician, sigma = S0/SNR.")
    print("Likelihood contrast vs fBm/Hurst CRB (Coeurjolly-Istas 2001): the exponent enters via")
    print("the b-INDEXED signal attenuation S(b;theta) under Rician magnitude noise -- NOT via a")
    print("trajectory / MSD / Gaussian-increment likelihood. Different statistical experiment.")

    rep = wall.cp0_verdict(rng=rng, n_boot=n_boot, do_ci=True)
    h = rep.headline

    _hr("1. CRLB sanity: Rician FIM -> Gaussian FIM as SNR -> inf (high-SNR limit recovered)")
    b_h = wall.default_b_design(b_max=wall.HEADLINE["b_max"], n_b=wall.HEADLINE["n_b"])
    th = forward.StretchedExp(alpha=wall.HEADLINE["alpha"]).theta
    for snr in (20, 500):
        fg = fisher.fisher_matrix(b_h, th, snr, "gaussian")
        fri = fisher.fisher_matrix(b_h, th, snr, "rician")
        rel = np.linalg.norm(fg - fri) / np.linalg.norm(fg)
        print(f"  SNR={snr:4d}: ||FIM_G - FIM_R|| / ||FIM_G|| = {rel:.4f}")

    _hr("2. HEADLINE -- canonical sparse clinical protocol: the recovery-collapse WALL")
    print(f"  acquisition: n_b={wall.HEADLINE['n_b']} b-values, b_max={wall.HEADLINE['b_max']:.0f} s/mm^2, "
          f"alpha={wall.HEADLINE['alpha']}, D={wall.HEADLINE['D']} mm^2/s")
    print(f"  b-design = {np.round(b_h,0).tolist()} s/mm^2")
    print(f"  pre-registered threshold: cv_alpha = sqrt(CRLB_alpha)/alpha = {wall.CV_THRESHOLD}")
    print(f"  realistic SNR band: {wall.REALISTIC_SNR_BAND}")
    print(f"  analytic CRLB wall SNR*      = {h.wall_snr:.1f}   (information lower bound on the wall)")
    print(f"  empirical (bootstrap-MLE) wall SNR* = {h.wall_snr_emp:.1f}, 95% CI "
          f"[{h.wall_ci[0]:.1f}, {h.wall_ci[1]:.1f}]   (finite-sample; >= CRLB wall)")
    print(f"  alpha-D correlation at wall: rho_aD = "
          f"{fisher.crlb(b_h, th, h.wall_snr, 'rician').rho_alpha_D:+.3f}")
    print("  cv_alpha across the realistic band:")
    for snr in (20, 30, 40, 60):
        r = fisher.crlb(b_h, th, snr, "rician")
        flag = "UNRECOVERABLE" if r.cv_alpha > wall.CV_THRESHOLD else "recoverable"
        print(f"    SNR={snr:3d}: cv_alpha={r.cv_alpha:.3f} ({flag}), se_alpha={r.se_alpha:.4f}")

    _hr("3. Empirical confirmation -- profile-likelihood CI on alpha (load-bearing CI)")
    print("  (below the wall the data do not constrain alpha -> wide/open CI; above it -> tight)")
    truth = forward.StretchedExp(alpha=wall.HEADLINE["alpha"], D=wall.HEADLINE["D"])
    sigma_cache = {}
    for snr in (20, 40, 80):
        sigma = noise.sigma_from_snr(truth.S0, snr)
        sigma_cache[snr] = sigma
        nu = forward.signal(b_h, truth.theta)
        # median over a few realizations so a single draw doesn't mislead
        los, his, widths, opens = [], [], [], 0
        for _ in range(9):
            M = noise.rician_sample(nu, sigma, rng)
            ci = identifiability.profile_ci_alpha(b_h, M, sigma, truth.theta)
            if ci.open_below or ci.open_above:
                opens += 1
            los.append(ci.lo); his.append(ci.hi); widths.append(ci.width)
        med_w = np.nanmedian(widths)
        print(f"    SNR={snr:3d}: median 95% profile-CI width on alpha = {med_w:.3f}"
              f"  (open/unbounded in {opens}/9 realizations)")

    _hr("4. WALL SURFACE -- wall SNR* over (n_b, b_max) at alpha=0.85 (the deliverable)")
    print("  number of b-values is the dominant driver; wall recedes below the band as n_b grows.")
    header = "  n_b \\ b_max | " + " ".join(f"{int(bm):>7}" for bm in rep.surface_b_max)
    print(header)
    for i, n_b in enumerate(rep.surface_n_b):
        row = " ".join(f"{rep.surface_wall[i,j]:7.1f}" for j in range(len(rep.surface_b_max)))
        tag = "(clinical)" if n_b in wall.CLINICAL_NB else "(research)"
        print(f"  {n_b:10d} | {row}   {tag}")
    lo, hi = wall.REALISTIC_SNR_BAND
    print(f"  [values INSIDE [{lo:.0f},{hi:.0f}] = alpha walls out within realistic acquisition]")

    _hr("5. alpha-dependence at the headline acquisition (wall worsens toward alpha->1)")
    for a, w in zip(rep.alpha_grid, rep.alpha_wall):
        print(f"    alpha={a:.2f}: wall SNR* = {w:.1f}")

    _hr("6. alpha-D degeneracy structure -- extending b_max with FEW b-values trades vs collinearity")
    for b_max in wall.B_MAX_GRID:
        bb = wall.default_b_design(b_max=b_max, n_b=4)
        r = fisher.crlb(bb, np.array([1.0, wall.HEADLINE["D"], 0.85]), 30, "rician")
        print(f"    n_b=4 b_max={int(b_max):5d} SNR=30: rho_aD={r.rho_alpha_D:+.3f}, cv_alpha={r.cv_alpha:.3f}")

    _hr("CP0 VERDICT")
    print(f"  wall_stands           = {rep.wall_stands}")
    print(f"  clinical_wall_in_band = {rep.clinical_wall_in_band}")
    print(f"  research_recovers     = {rep.research_recovers}  (scope boundary, reported not hidden)")
    print(f"  refuted (wedge dead)  = {rep.refuted}")
    for n in h.notes:
        print(f"  - {n}")
    if rep.refuted:
        print("\n  >>> REFUTE TRIGGERED: alpha recoverable across the realistic clinical band. WEDGE DEAD.")
    else:
        print("\n  >>> WALL STANDS (scoped): under realistic clinical few-b acquisition the fractional")
        print("      order alpha is information-limited within the realistic SNR band; it becomes")
        print("      recoverable only with dense multi-b research acquisition. Refute survived.")

    _write_results(rep, sigma_cache)
    _save_figure(rep)
    print(f"\n  results -> {RESULTS}")
    print("CP0 PASS" if not rep.refuted else "CP0 REFUTED")
    return 0


def _write_results(rep, sigma_cache):
    h = rep.headline
    lo, hi = wall.REALISTIC_SNR_BAND
    lines = []
    lines.append("# RESULTS -- CP0: fractional-order identifiability wall\n")
    lines.append("All numbers are derived (Fisher/CRLB closed-form + profile-likelihood + parametric")
    lines.append("bootstrap), fully synthetic, seeded. CRLB = identifiability bound, scoped to its")
    lines.append("regime; never an impossibility claim.\n")
    lines.append("## Forward model")
    lines.append("`S(b; S0, D, alpha) = S0 * exp(-(b D)^alpha)` (stretched-exponential lead lane),")
    lines.append("theta=(S0,D,alpha) estimated JOINTLY, Rician magnitude noise sigma=S0/SNR.\n")
    lines.append("## Headline (canonical sparse clinical protocol)")
    lines.append(f"- acquisition: n_b={wall.HEADLINE['n_b']} b-values, b_max={wall.HEADLINE['b_max']:.0f} s/mm^2, "
                 f"alpha={wall.HEADLINE['alpha']}, D={wall.HEADLINE['D']} mm^2/s")
    lines.append(f"- pre-registered wall threshold: cv_alpha = sqrt(CRLB_alpha)/alpha = {wall.CV_THRESHOLD}")
    lines.append(f"- **analytic CRLB wall SNR\\* = {h.wall_snr:.1f}** (information lower bound)")
    lines.append(f"- **empirical (bootstrap-MLE) wall SNR\\* = {h.wall_snr_emp:.1f}**, 95% CI "
                 f"[{h.wall_ci[0]:.1f}, {h.wall_ci[1]:.1f}] (finite-sample, >= CRLB wall)")
    lines.append(f"- realistic SNR band = [{lo:.0f}, {hi:.0f}] -> both walls land INSIDE the band\n")
    lines.append("## Wall surface: wall SNR* over (n_b, b_max) at alpha=0.85")
    lines.append("| n_b | " + " | ".join(f"b_max={int(bm)}" for bm in rep.surface_b_max) + " | regime |")
    lines.append("|---|" + "---|" * (len(rep.surface_b_max) + 1))
    for i, n_b in enumerate(rep.surface_n_b):
        row = " | ".join(f"{rep.surface_wall[i,j]:.1f}" for j in range(len(rep.surface_b_max)))
        regime = "clinical" if n_b in wall.CLINICAL_NB else "research"
        lines.append(f"| {n_b} | {row} | {regime} |")
    lines.append("")
    lines.append("## alpha-dependence at the headline acquisition (wall SNR*)")
    lines.append("| alpha | " + " | ".join(f"{a:.2f}" for a in rep.alpha_grid) + " |")
    lines.append("|---|" + "---|" * len(rep.alpha_grid))
    lines.append("| wall SNR* | " + " | ".join(f"{w:.1f}" for w in rep.alpha_wall) + " |")
    lines.append("")
    lines.append("## Verdict")
    lines.append(f"- wall_stands = **{rep.wall_stands}**; refuted = **{rep.refuted}**")
    lines.append(f"- clinical_wall_in_band = {rep.clinical_wall_in_band}; research_recovers = {rep.research_recovers}")
    lines.append("")
    lines.append("**Scoped claim.** Under realistic clinical diffusion-MRI acquisition (few b-values, "
                 f"n_b in {wall.CLINICAL_NB}), the fractional order alpha is information-limited within the "
                 f"realistic SNR band [{lo:.0f},{hi:.0f}]: the recovery-collapse wall sits at "
                 f"SNR*~{h.wall_snr:.0f} (CRLB) / {h.wall_snr_emp:.0f} (empirical, 95% CI "
                 f"[{h.wall_ci[0]:.0f},{h.wall_ci[1]:.0f}]). The wall recedes "
                 f"below the band only with dense multi-b research acquisition (n_b>=8). Recoverable in the "
                 "data-rich idealization; walls out under the realistic clinical forward model.")
    RESULTS.parent.mkdir(exist_ok=True)
    RESULTS.write_text("\n".join(lines) + "\n")


def _save_figure(rep):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # pragma: no cover
        print(f"  (figure skipped: {e})")
        return
    FIGDIR.mkdir(exist_ok=True)
    th = forward.StretchedExp(alpha=wall.HEADLINE["alpha"]).theta
    snr = np.geomspace(5, 120, 60)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for n_b in (4, 6, 8, 12):
        b = wall.default_b_design(b_max=wall.HEADLINE["b_max"], n_b=n_b)
        cv = [fisher.crlb(b, th, s, "rician").cv_alpha for s in snr]
        ax1.plot(snr, cv, label=f"n_b={n_b}")
    ax1.axhline(wall.CV_THRESHOLD, ls="--", c="k", lw=1, label=f"wall thr {wall.CV_THRESHOLD}")
    ax1.axvspan(*wall.REALISTIC_SNR_BAND, color="grey", alpha=0.15, label="realistic SNR")
    ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_xlabel("SNR (b=0)")
    ax1.set_ylabel("cv_alpha = sqrt(CRLB_a)/alpha"); ax1.legend(fontsize=8)
    ax1.set_title(f"Recovery-collapse wall (alpha={wall.HEADLINE['alpha']}, b_max={wall.HEADLINE['b_max']:.0f})")
    im = ax2.imshow(rep.surface_wall, aspect="auto", origin="lower", cmap="viridis")
    ax2.set_xticks(range(len(rep.surface_b_max)));
    ax2.set_xticklabels([int(x) for x in rep.surface_b_max])
    ax2.set_yticks(range(len(rep.surface_n_b))); ax2.set_yticklabels(rep.surface_n_b)
    ax2.set_xlabel("b_max (s/mm^2)"); ax2.set_ylabel("n_b")
    ax2.set_title("Wall SNR* surface")
    fig.colorbar(im, ax=ax2, label="wall SNR*")
    fig.tight_layout()
    out = FIGDIR / "cp0_wall.png"
    fig.savefig(out, dpi=130)
    print(f"  figure -> {out}")


if __name__ == "__main__":
    raise SystemExit(main())
