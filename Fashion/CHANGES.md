# CHANGES.md — friction-remediation run (`fix/friction-remediation`)

Addresses the Universal Friction Engine findings (HC1–HC6, CS1–CS4, G1–G4) on the
Project Fashion IVIM uncertainty-quantification manuscript. Honesty-first: results
that weakened/tempered/corrected a claim are in **ADVERSE_RESULTS.md**; open items
in **FLAGS.md**. Built on `fix/assessor-remediation` (manuscript + frozen numbers).

One commit per checkpoint. New numbers are additive; the numbers-freeze gate
(`paper/check_numbers.sh`) passes (no baseline count decreased).

## CP0 — enumerate & diagnose
`paper/FRICTION_DIAGNOSIS.md`: compute/data inventory (sbi installed in an isolated
venv → NPE retrains runnable; trained `.pt` absent → retrain; in-vivo raw not local)
and a per-item runnable-now / -with-compute / needs-data / manuscript-only table.

## CP1 — simulation-budget sweep (HC4/CS3)  [executed]
`npe/run_cp1_budget_sweep.py` retrains the NPE at budgets {50k,100k,250k,500k,1M},
recipe otherwise fixed (set/NSF/log-D*/seed 0/box prior/clinical-sparse), and
recomputes per-SNR median D* claimed/achieved SD-to-CRLB ratios via the *identical*
efficiency audit. Outputs `npe/cp1_budget_sweep.csv`, supplementary figure
`figS6_budget_sweep`. Summariser validated against the frozen baseline (reproduces
Table 2 D* 0.084/0.160/0.376/0.671 and below-floor 0.86→0.52).
**Verdict: CONFIRMATORY — budget-invariant.** D* claimed ratio and below-floor
fraction (~0.69 overall) are flat across the full 20× range (50k→1M); the
overconfidence is not simulation starvation. Added as a fourth ablation
(manuscript Limitations, Supplementary Methods S6, Figure S6).

## CP2 — bias-aware floor (HC1/CS1)  [done]
`npe/run_cp2_bias_aware_floor.py` + `npe/cp2_bias_aware_floor.csv`. Concedes
(Hero & Fessler 1993) that the unbiased CRLB does not bound a biased estimator;
re-audits the claimed D* SD against a van-Trees Bayesian floor under the actual
log-uniform prior. Overconfidence **survives**: 78%(SNR10)→92%(SNR100) of flagged
points stay below the bias-aware floor (58.8% of all grid pts vs 69.5% unbiased);
prior explains 8–22%. Manuscript: Hero-Fessler cited (ref 28), unbiased-CRLB ratios
demoted to a diagnostic, CRLB-independent evidence foregrounded. (CRLB recompute
matches the frozen `crlb_sd` column to 1e-14.)

## CP3 — NLLS baseline survivorship (HC5)  [done]
`npe/run_cp3_railing_audit.py`. Abdominal S4 D* spread now reported BOTH ways —
railed-included (all 1618): SD 0.27 / IQR 0.21; railed-excluded (733, OLD): SD 0.41
/ IQR 0.38. Brain N=500 NLLS coverage documented as already railed-included (failed
fits counted not-covered). Supplement S4 caption + main-text clause updated; OLD
numbers retained and labelled.

## CP4 — OOD gate decorrelation (HC6)  [executed]
`npe/run_cp45_decorrelate_control.py`. Re-scores the self-consistency gate in
simulation against independent (parameter-recovery) targets. **ADVERSE:** the
AUC ~0.99 reproduces only against the shared-fit held-out residual; against true
D* parameter-recovery error the AUC collapses toward chance. Gate reframed as a
noise/abstention triage signal. [final numbers pending 500k → ADVERSE_RESULTS.md]

## CP5 — in-silico control + RNPE (HC2/CS2)  [executed]
`npe/run_cp45_decorrelate_control.py`. In-silico control (correct forward model,
truth known): in-distribution parameter coverage ≈ nominal even stratified →
**signal-domain coverage ≠ parameter calibration** (TEMPERED framing). RNPE-style
model criticism with an SNR-matched well-specified null applied to the in-vivo
self-consistency residuals. Real in-vivo parameter-truth gap FLAGGED as
needs-external-data. [final numbers pending 500k]

## CP6 — novelty positioning vs the graveyard (G1–G4/CS4)  [done]
Manuscript Discussion gains a positioning paragraph conceding documented
amortized-posterior overconfidence (Hermans 2022; Cannon 2022; Ward 2022; Frazier
2020) and ensemble/recalibration degradation under shift (Ovadia 2019), then
sharpens the contribution (bias-aware information-floor diagnostic in IVIM;
aggregate-vs-pointwise in qMRI). CS4 subtle(sim)-vs-gross(in-vivo) scoping added;
recalibration safeguard noted to inherit Ovadia's caveat. Six verified references
(23–28) added to `manuscript.tex` and `refs.bib`; the original 21 + dataset intact.

## Build & gates
`tectonic manuscript.tex` and `supplement.tex` compile clean; new refs resolve.
`paper/check_numbers.sh`: PASS (no baseline number decreased; all changes additive).
Bibliography: 21 original + 1 dataset (osipidata, prior run) + 6 verified graveyard
refs = 28 entries; all 6 cited inline (23–28); no fabricated identifiers.

## Final disposition — NOT submission-ready

This run hardened the science and the honesty of the manuscript; it did **not**
make it submission-ready. Classifying the structural items by what the experiments
actually did:

**Genuinely downgraded by experiment**
- **OOD self-consistency gate (CP4/HC6).** The headline AUC 0.99 is largely a
  shared-fit/shared-noise artefact: against an independent parameter-recovery
  target the AUC collapses to ~0.60. The gate is downgraded from "detects the
  unobservable miscalibration" to a noise/abstention triage heuristic. This is the
  one claim materially weakened by experiment.

**Reframed / tempered (claim survived the test, framing corrected)**
- **Bias-aware floor (CP2/HC1).** Conceptual concession (Hero & Fessler): the
  unbiased CRLB is not a floor for a biased estimator, so the unbiased-CRLB ratios
  are demoted to a diagnostic. But the effect **survived** the van-Trees bias-aware
  floor (78–92% of flagged points), so the result is reframed, not downgraded.
- **In-silico control (CP5/HC2).** In-distribution parameter calibration holds
  (coverage ≈ nominal, even stratified), so "pointwise miscalibrated" is tempered to
  "information-floor-inefficient and transfer-fragile," and the in-vivo limb is
  scoped as a held-out-signal/OOD check.
- **NLLS railing (CP3/HC5).** Both railed-included and railed-excluded spreads now
  reported; the brain-coverage survivorship premise was corrected (already
  railed-included).
- **Novelty positioning (CP6).** Reframed to concede prior work and sharpen the
  delta; no result changed.

**Strengthened by experiment**
- **Budget sweep (CP1/HC4).** Budget-invariance across 50k–1M confirms the
  "intrinsic" claim as a fourth ablation.

**Open — cannot be closed by CC (see FLAGS.md)**
- Real in-vivo parameter ground truth (HC2/CS2): needs phantom or
  repeated-acquisition data. The in-vivo limb remains a held-out-signal/OOD check.
- In-vivo 3-way b-split gate decorrelation (HC6): downloadable raw data + `.pt`
  required; only the simulation decorrelation was run.
- Per-protocol recalibration experiment: proposed, not implemented; inherits
  Ovadia's shift-degradation caveat.

Figures 1–4 and S1–S5 were not regenerated (committed CSV numbers unchanged; new
analyses are additive); figS6 is newly generated. A camera-ready pass should
regenerate figures from current models and resolve the FLAGS.md typeset items.
