# Gnomon CP3 verdict — reproduce-or-refute of Fashion's ruler

**Verdict: PARTIAL — the diagnostic and the mechanism reproduce; two headline
*marginal* coverage magnitudes do not, and the reason is a substantive,
under-documented methodological choice (not fabrication).**

Run: `KMP_DUPLICATE_LIB_OK=TRUE PYTHONPATH=Gnomon python -m gnomon.reproduce`
→ machine record: [`results/reproduction.json`](results/reproduction.json) (seed 20260620).
This is a **CP3 hard halt**: it changes the retool, so Gnomon does **not** auto-proceed
to packaging the numbers as "reproduced."

## Scorecard (frozen targets vs the clean rebuild)

| Target | Claimed | Rebuilt (95% boot CI) | Tol | Result |
|--------|---------|-----------------------|-----|--------|
| **T1** NLLS D\* railing (real OSIPI abdomen) | **54.7%** | **54.2%** [52.0, 56.4] (rail_tol 1e‑3; 56.2% at 1e‑2; 58.7% SNR>25) | ±5 pp | **REPRODUCES** ✅ |
| **T3c** D\* coverage, MCMC **quantile** | **0.94** | **0.90** [0.89, 0.91] | ±0.05 | **REPRODUCES** ✅ |
| **T3c′** D, f coverage, MCMC quantile | **0.94** | D 0.95, f 0.96 | ±0.05 | **REPRODUCES** ✅ |
| **T4** MAF flow vs railed NLLS (ECE/sharp/cov) | directional | ECE 0.069<0.121, sharp 0.112<0.181, cov 0.979>0.763 (all CIs exclude 0) | sign+CI | **REPRODUCES** ✅ |
| **T3a** D\* coverage, Laplace **SD** | **0.30** | **0.80** [0.78, 0.82] | ±0.10 | **DIVERGES** ✗ |
| **T3b** D\* coverage, MCMC **SD** | **0.67** | **0.90** [0.89, 0.92] | ±0.10 | **DIVERGES** ✗ |

4 of 6 reproduce, including the single most load-bearing one (the 54.7% railing rate,
on the **real open data**, sha256-verified). The divergence is confined to the two
*severe Gaussian-overconfidence* marginal numbers.

## What reproduces — and why that matters

- **The railing diagnostic is real.** An independent NLLS, on the open OSIPI abdomen
  ROI, rails D\* in **54.2%** of voxels — Fashion's 54.7% sits inside the bootstrap CI.
  This is the paper's per-voxel motivating fact, and it survives a clean rebuild.
- **The qualitative mechanism is real.** Symmetric/Gaussian intervals under-cover D\*,
  the **shape-correct quantile interval recovers** it (T3c), and the **amortized flow
  beats the railed NLLS** on ECE, sharpness, and coverage simultaneously (T4).
- **The failure is located where Fashion/Gauge say it is.** Per true-D\*-tercile, the
  Gaussian under-coverage concentrates in **high D\***: Laplace 0.91 / 0.86 / **0.63**;
  MCMC-SD 0.95 / 0.95 / **0.81**. The "identifiability wall" reproduces.

## The divergence — which numbers, by how much, likely cause

Fashion's **marginal** Gaussian coverages (Laplace 0.30, MCMC-SD 0.67) are *much* more
severe than the clean rebuild's (0.80, 0.90). Two compounding, **under-documented**
choices explain the gap — both run, not asserted
([`scripts/divergence_diagnostic.py`](scripts/divergence_diagnostic.py)):

1. **Cohort regime (≈ half the gap).** Fashion's pooled 0.30/0.67 imply a cohort
   concentrated in the hard **high-D\*/low-perfusion** corner. Gnomon draws truths
   across the Lattice physiological prior, so pooling *dilutes* the wall: pooled
   Laplace 0.80 vs high-D\*-tercile **0.63**. Fashion's "3 pancreas truths" are never
   listed in its prose — the exact regime is a **completeness gap**.
2. **Railed-voxel uncertainty convention (the rest of the gap).** A railed/unidentified
   D\* has *no local information*; the **honest** asymptotic covariance therefore makes
   its interval **wide** (→ over-cover). Fashion's baseline is "overconfident by
   design" — it must instead assign those voxels **near-zero SD**. Swapping Gnomon to
   that floored convention drops Laplace D\* coverage to pooled **0.68** / high-tercile
   **0.41** — and the floored *pooled* number (0.68) lands essentially on Fashion's
   MCMC-SD claim (0.67). This is precisely the **CRLB/uncertainty assumption** the
   reviewers flagged as undocumented.

**Reading:** the phenomenon is genuine, but the dramatic headline *marginal* numbers
are **not regime- or convention-robust**. They are an artifact of (1) an unstated hard
cohort and (2) an unstated overconfident treatment of unidentified voxels. Presented
as marginal coverages without those caveats, they are not independently reproducible —
which is consistent with the "internal inconsistency / incompleteness" the manuscript
was returned for. The honest, reproducible version of the same claim is **conditional**
(per-D\*-tercile) coverage with the uncertainty convention stated.

## Implication for the retool (why this is a hard halt)

This changes how the numbers should be presented; it does **not** sink the paper:

- **Keep, strengthened:** the railing result (T1 — real-data, exact, the strongest
  asset), the quantile-interval fix (T3c), and the flow-vs-railed-NLLS comparison (T4).
- **Re-frame, do not re-quote:** replace the headline **marginal** D\* coverages
  (0.30/0.67) with **per-D\*-tercile (conditional)** coverages, which reproduce and are
  regime-robust, OR pin the exact truths + uncertainty convention so the marginal
  numbers become reproducible.
- **Document the flagged assumption:** state explicitly that the "overconfident"
  baseline floors the SD of railed/unidentified D\* (otherwise the honest CRLB
  over-covers there) — this is the Cramér–Rao/uncertainty item the reviewers wanted.

Gnomon is the clean technical core for that retool: an independent implementation that
shares no code with Fashion, plus the complete methods write-up
([`docs/METHODS.md`](docs/METHODS.md)). It does **not** become a standalone paper.
