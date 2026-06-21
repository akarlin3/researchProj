# Reproduction targets (narrative) — frozen before running

The authoritative, machine-readable version is [`gnomon/manifest.py`](gnomon/manifest.py)
(`manifest.validate()` enforces internal consistency). This file is the human-readable
companion. **All targets are read from Fashion's prose (`*.md`), never its source.**
Tolerances are **frozen now**, before any rebuild runs (guardrail 4); CP3 compares
against them and may not loosen them afterward.

## T1 — NLLS D\* boundary-railing rate (open real data)

- **Claimed:** **54.7%** of 1618 high-SNR ROI voxels have a boundary-railed NLLS D\*
  estimate. — *Fashion/REVIEWER_RESPONSE_R2.md:39,50*
- **Substrate:** OSIPI TF2.4 abdomen scan, Zenodo **14605039** (open, CC),
  download-on-demand. Clinical-sparse 8-b.
- **Tolerance:** point estimate within **±5 pp** of 54.7% (i.e. [49.7%, 59.7%]),
  **or** 54.7% inside the voxel bootstrap 95% CI. Reported at `rail_tol` 1e-3 and 1e-2.
- **Note:** the "high-SNR ROI" / 1618-voxel selection is under-specified in Fashion's
  prose; Gnomon documents its own threshold and reports railing-rate sensitivity to it.
- **Context (not an independent target):** Caliper's *synthetic* NLLS D\* railing at
  SNR 20 is ~9.5% — the synthetic-vs-real gap is itself part of the methods story.

## T3 — the headline coverage table (synthetic 9-cell set)

Design: 3 ground-truth conditions × SNR {10, 20, 40} × 200 noise realizations,
nominal 0.95, parameter D\*.

| Sub-target | Interval kind | Claimed coverage | Tolerance | Source |
|---|---|---|---|---|
| **T3a** | Laplace Gaussian **SD** | **0.30** (severely overconfident) | ±0.10 | README.md:64 |
| **T3b** | MCMC Gaussian **SD** | **0.67** (overconfident) | ±0.10 | README.md:65 |
| **T3c** | MCMC 2.5/97.5 **quantile** | **0.94** (≈ nominal — *the headline*) | ±0.05 | README.md:66 |
| **T3c′** | MCMC quantile, **D and f** | **≈0.94** each | ±0.05 | README.md:69 |

The mechanism Gnomon must reproduce: the right *shape* (quantiles), not a wider SD,
recovers D\* coverage. D and f are already near-nominal — only D\* needed the fix.

## T4 — MAF flow vs boundary-railed NLLS (ECE & sharpness *behavior*)

- Fashion's prose **defines** ECE and sharpness but states **no numeric values**
  (README.md:35-37). The target is therefore the **direction/behavior**, not a point.
- **Gate:** on a held-out synthetic set, the MAF flow has **lower D\* ECE** and is
  **sharper** than the railed NLLS baseline, at **coverage ≥ the baseline**; each
  gap's bootstrap CI excludes 0.
- **Context (not independent):** Caliper documents NLLS D\* coverage 0.786 vs flow
  0.875; ECE 0.075 vs 0.016; sharpness 151.7 vs 48.96.

## Verdict logic (CP3, hard halt)

- **REPRODUCES** — every gated target within tolerance (with bootstrap CIs):
  Fashion's numbers are real; its rejection was *presentational*. → CP4.
- **DOES NOT REPRODUCE** — any gated target misses: *substantive* inconsistency.
  Emit a precise divergence report (which numbers, by how much, likely cause). **STOP.**
