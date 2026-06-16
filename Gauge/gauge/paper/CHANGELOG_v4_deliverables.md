# Manuscript integration changelog — three deliverables (v4 sign-off)

Integrates the three deferred experiment deliverables into `gauge_v3_revised.tex`.
**Every numeral is pulled verbatim from a gated, seeded results file and survives
`consistency.py` (GATE 3: PASS, 34/34 verbatim + 36/36 16-seed bands).** No
hand-typed numbers. The central characterize-don't-solve claim is unchanged and was
strengthened only within its existing scope (label-free, bounded — no universal
quantifier).

## Falsification gate — NOT tripped
- **03-2 Bayesian-shrinkage:** `VERDICT: A — DISSOCIATION / WALL HOLDS` (confirming).
  Deciding numbers (single-seed report): shrinkage hi-D* RMSE 24.15e-3 vs MCMC
  20.54e-3; shrinkage hi-D* conformalized worst-SNR cell 0.762 vs nominal 0.90.
  Sanity gates pass (no width collapse 1.14×; routing not fixed) → not an
  illegitimate collapse.
- **02-3 Arm 1 alt-forward-model:** `BRANCH A — WALL IS GENERAL / CIRCULARITY BROKEN`
  (confirming). Continuity gate exact (max |ΔS| = 0.00e+00).
- **04-3 Arm 2:** envelope with zero-deviation recovery + fires-before-failure.

All Branch A. No abstract rescope required.

## Provenance / method decisions
- **PR #19 (shrinkage dissociation) merged to `main` first** (per sign-off), so
  `gauge/dissociation.py` + results are canonical before integration.
- **16-seed dissociation bands were already produced & committed** by PR #19, folded
  into the main benchmark sweep: `results/multiseed.json` carries n_seeds=16 with all
  `dissoc/*` keys and `consistency.py` already asserts them. The bands are
  deterministic, gated, and gate-passing; re-running the sweep would only reproduce
  identical-but-NN-nondeterministically-drifting numbers (a seed-0 spot-check
  reproduced VERDICT A: 24.17e-3 / 0.756, within the banded NN tolerance — which is
  exactly why NN-derived numbers ship as bands, not byte-traced). We therefore cite
  the committed gated bands rather than recompute.
- **Shipped numbers are the across-seed MEANS** from `multiseed.json` /
  `altmodel_multiseed.json` (the (E) convention), NOT the single-seed report values.
  e.g. Arm-1 surrogate-A recal hi-D*eff ships as 0.796 [0.748, 0.854] (16-seed mean),
  not the report's seed-0 0.850.

## Text edits A–J (`gauge_v4_easy_fixes.md`) — idempotency: ALL ALREADY APPLIED
Merged via PR #18 (`9909cce`) before this run; verified present in the canonical
`.tex`. No re-application. (Spec file itself is not in-repo; #18's merged A–J is the
authoritative applied set.)

## Claims touched (old → new, source)

| # | Location | Old → New | Source |
|---|----------|-----------|--------|
| 1 | Abstract — Results | (b-value/monitor sentence) → +clause: biased shrinkage buys point precision not hi-D* coverage; gap replicates under non-bi-exp generator; off-model envelope charted | 03-2 A + 02-3 Arm1 A + 04-3 Arm2 |
| 2 | Abstract — Conclusion | "unrecoverable by the label-free methods we tested" → "… now including a deliberately biased shrinkage estimator, and replicating under a non-bi-exponential generator …" | 03-2 A + 02-3 Arm1 A |
| 3 | Intro — contribution (i) | "a property of the **bi-exponential model class** and the data" → "a property of the **perfusion-estimation problem** and the data … persists for a biased shrinkage estimator and replicates under a non-bi-exp generator … off-model envelope" | 03-2 A + 02-3 Arm1 A + 04-3 Arm2 |
| 4 | §4.8 OSIPI — scope sentence | "intrinsic to the **bi-exponential model** and the data" → "intrinsic to the **perfusion-estimation problem** … replicates even when the ground truth is non-bi-exponential (velocity-dispersion)" | 02-3 Arm1 A |
| 5 | **NEW §4.9** `sec:dissociation-shrinkage` + **Table `tab:dissociation`** | (none) → point-vs-coverage dissociation: shrinkage lowers lo-tercile RMSE 15.1 [14.7,15.6] but not hi-D* (24.6 [24.1,25.3] vs 20.5), hi-D* marg 0.827 [0.801,0.845] / worst-SNR 0.721 [0.663,0.765]; no collapse (1.15×); misroutes 34% | 03-2 A; bands `dissoc/*` in `multiseed.json` |
| 6 | **NEW §4.10** `sec:altmodel` | (none) → non-bi-exp replication: recal hi-D*eff 0.796 [0.748,0.854] (A) / 0.800 [0.756,0.845] (B); marg restored 0.898; log-normal 0.797; continuity CV=0 exact + 0.798; monitor at chance (AUC 0.501) = hidden channel | 02-3 Arm1 A; `altmodel_multiseed.json` |
| 7 | **NEW §4.11** `sec:envelope` | (none) → off-model envelope: hi-deviation D* cov 0.769 (tri-exp) / 0.770 (stretched), f 0.791 (dispersion); zero-dev recovery 0.903; monitor AUC weak/family-dependent (tri-exp 0.532) | 04-3 Arm2; `altmodel_multiseed.json` |
| 8 | **NEW Fig `fig:altmodel`** | (none) → wired existing gated `figures/altmodel.pdf` (A: Arm-1 tercile coverage; B: Arm-2 degradation; C: Arm-2 monitor AUC) | altmodel experiment |
| 9 | Discussion | (observable-vs-latent ¶) → +sentence: three independent stresses corroborate information-not-artifact split | 03-2 + 02-3 Arm1 + 04-3 Arm2 |
| 10 | Limitations | "tri-exponential … schematic probes" → "off-model behavior now **characterized** (non-bi-exp replication + validity envelope); schematic deviations not a tissue model; in-vivo unverifiability **unchanged**" | 02-3 Arm1 + 04-3 Arm2 |
| 11 | `consistency.py` | +3 verbatim CHECKS (dissociation Branch-A verdict + 2 sanity gates). `dissoc/*` + `altmodel/*` band assertions already present (PRs #19/#20–22). | gating |

## Build / gate
- `gauge/paper/build.sh` → GATE 3 **PASS**; `tectonic` compiles clean; 0 undefined
  refs/citations, 0 `??`. `gauge_v3_revised_R2.pdf` rebuilt (20 → 23 pp, +3).

## Flags for human review
- **Central-claim wording changed** (bounded re-generalization to "perfusion-estimation
  problem"): `Gauge_MRM_cover_letter` may need a matching pass — **not edited** here.
- **`gauge_v3_revised_FIXED.tex` / `_R2_FIXED.pdf`** are a divergent parallel lineage
  (149-line diff from the canonical source, not built by `build.sh`) — **left
  untouched**; update separately if it is a live submission artifact.
- The single tri-exp robustness probe (§4.4) is retained; §4.11 generalizes it rather
  than replacing the text — confirm that is the intended relationship.
