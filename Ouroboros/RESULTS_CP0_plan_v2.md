# Checkpoint 0 — Recon + Plan (v2 recalibrated CNSNS assessment)

Run-then-write. Numbers come only from RESULTS files. This plan enumerates the frozen
rules and edit sites; **no manuscript edits happen in this checkpoint.**

## 1. Frontmatter / TODO state (verbatim, `manuscript/manuscript.tex`)

- L11  `\newcommand{\zenododoi}{<PENDING>}`  ← keep intact, do NOT fabricate.
- L21  `\author{Avery Karlin\corref{cor1}}`
- L22  `\ead{averykarlin3@gmail.com}`
- L23  `\cortext[cor1]{Corresponding author.}`
- L24  `\address{Independent Researcher}`  ← replace with dual affiliation.
- L25-26 ORCID TODO  ← keep intact, do NOT fabricate.
- L27-28 affiliation `% DECISION` TODO  ← remove (resolved).
- L361,363 Zenodo `<PENDING>` TODO + sentence  ← keep intact.

## 2. Tikhonov λ schedule (the oracle assist)

`ouroboros_mitigation.py:174` `lambda_map = {100:0.01, 80:0.1, 60:1.0, 40:10.0, 30:30.0,
20:100.0, 10:500.0}` and `ouroboros_cp3_fairgrid.py:46`
`tikhonov_lambda(snr) = 10**(3.2 - 0.05*snr)`. Picking λ **requires knowing the SNR** —
this is the oracle assist flagged by the assessment. The fair-grid table
(`RESULTS_mitigation_fairgrid.md`) used this SNR-tuned λ; under it Tikhonov beats weak-form
at α=0.9 ([<10,10] dB vs weak [15,20]).

## 3. A(α) and target-derivative code

- A(α): `ouroboros_noise_analysis.py:compute_analytic_factors`, A(α)=h^{-2α}‖w(α)‖₂².
  Existing factors in `RESULTS_noise_amplification.md` §2 (A(0.5)≈127.1, A(0.9)≈7189.5).
- Reconcile: `RESULTS_Aalpha_threshold_reconcile.md` already shows A(α) rises 8.76 dB/step
  but pointwise brackets rise only 5 dB/step → **−3.76 dB/step shortfall**, and explicitly
  flags "the exact signal-amplitude factor is not computed here ... flagged as a gap." CP2
  fills exactly this gap.
- Target derivative: `gl_derivative_time(u_clean, dt, α)` at the true order, on the clean
  trajectory `solve_fractional_system(α, Nt, dt, Nx, y0)` with the canonical y0
  (Gaussian p0, c0, n0; Nx=50, Nt=500, T=5, dt≈0.01002, k_start=20).

## 4. FROZEN no-SNR-knowledge Tikhonov rule (CP1)

**PRIMARY = GCV-selected λ (noisy-data-only).** For each realization, choose a single
global λ by minimizing the Generalized Cross-Validation functional of the second-difference
Tikhonov smoother S(λ)=(I+λ D2ᵀD2)⁻¹ computed on the **noisy trajectory only** — no SNR, no
clean data, no true order. Computed in the eigenbasis of L=D2ᵀD2 (one eigendecomposition,
O(n) per λ): GCV(λ) ∝ RSS(λ)/tr(I−S(λ))², with RSS(λ)=Σᵢ rᵢ(λ)²‖Cᵢ‖², rᵢ=λμᵢ/(1+λμᵢ),
Cᵢ the eigen-coefficients pooled over all 150 field columns. λ grid = `logspace(-1, 3.5, 60)`.

**SECONDARY robustness = fixed λ = 28.2** = 10^(3.2−0.05·35), the SNR-tuned schedule's value
at the grid midpoint SNR=35 dB — the single best noise-agnostic guess if forced to commit
to one λ.

Justification: GCV is the textbook noisy-data-only criterion for the smoothing parameter;
it is the fair deployable analog of the SNR-tuned schedule and gives Tikhonov its best
honest shot. The fixed-midpoint λ is the alternative "commit to one λ" baseline. **This rule
is frozen now; its honest output is reported in CP1 even if Tikhonov no longer beats
weak-form at α=0.9 — or if it still does (honesty guard 1).**

Everything else identical to the fair grid: grid {10,15,…,60} dB, 500 realizations,
seed `int(snr)+42`, clean-derivative true-order oracle scoring, smooth-then-pointwise-select.
The run re-emits weak-form and SNR-tuned-Tikhonov cells in the same harness as self-checks
(must reproduce weak [15,20] and SNR-tuned Tikhonov [<10,10] at α=0.9).

## 5. CP2 plan

Compute RMS of the true target fractional derivative D^α_t x on the clean trajectory, per α
on the {0.3,…,1.0} grid (same as A(α)), over the SINDy evaluation range [k_start:]. Convert
to power dB (10·log10 RMS²) and report the per-0.2-step rise. Compare to the 3.76 dB/step
shortfall: if the signal-power rise ≈ 3.76 dB/step it quantitatively accounts for the
A(α)-vs-threshold gap (keep the offset claim); if much smaller, soften the claim (honesty
guard 2). Emit a new subsection into `RESULTS_noise_amplification.md`.

## 6. Manuscript edit sites (CP3/CP4, numbers only from CP1/CP2)

- **Frontmatter** L21-28: dual affiliation block; remove DECISION TODO; keep ORCID/Zenodo.
- **Abstract** L31: trim for concision (ceiling-vs-deployable, benchmark, Lyapunov caution).
- **§3.3 A(α) caveat** L199: update/soften the "8.75-vs-5 dB/step" offset sentence per CP2;
  cite new RESULTS subsection.
- **§3.4 mitigation** L234, L239, Table L221: update Tikhonov text per CP1 fair-λ result.
- **§4 Discussion** L327 (novelty/contribution), L329 (remedy framing): significance
  reframe to identifiability-limits/cautionary map + explicit practitioner takeaway;
  fair-λ Tikhonov update; add generality limitation (two low-D systems, dimension-agnostic
  A(α), higher-D/real data = future work).
- **§4 biology** L346-354 table + surrounding: lightly de-emphasize tumor framing; keep the
  claim-constraint table.
- **Conclusion** L357: align with fair-λ Tikhonov + significance reframe.
- **Highlights** L45-46: align with fair-λ result + reframe.

## 7. Compute estimate

- CP1: 4 methods (weak, tikhonov-snrtuned, tikhonov-gcv, tikhonov-fixed) × 3 α × 11 SNR =
  132 cells × 500 trials, Pool(8). GCV adds cheap eigenbasis ops. ~5–15 min.
- CP2: 8 α × one clean solve + one GL derivative. < 1 min.

Proceeding to CP1.
