# Checkpoint 0 — Friction-remediation enumerate & diagnose

Continuation of the Project Fashion manuscript remediation, addressing the
**Universal Friction Engine** findings (HC1–HC6, CS1–CS4, G1–G4). These target
the *scientific claims* (CRLB below-floor argument, NLLS survivorship, OOD-gate
correlation, simulation budget, bias-aware floor, novelty positioning), not the
formatting/provenance items closed in `DIAGNOSIS.md` (assessor remediation).

Branch `fix/friction-remediation`, based on `fix/assessor-remediation` (carries
`paper/manuscript.tex`, `paper/supplement.tex`, `paper/refs.bib` (21 entries +
dataset), figures, and `NUMBERS_FROZEN.txt`).

## Compute / data inventory (verified this run)

| Resource | Status |
|---|---|
| CPU (Apple M4, 10-core) | available |
| GPU | none |
| `sbi` 0.26.1 + torch 2.12.1 + numpy 2.3.5 (py3.11) | **installed** in isolated venv `$CLAUDE_JOB_DIR/venv-npe` (inherits proteus torch; needs `KMP_DUPLICATE_LIB_OK=TRUE`). NPE **retrains are runnable**. |
| Trained NPE posteriors (`npe/*.pt`) | **absent** (gitignored). Committed CSVs were produced from models that no longer exist ⇒ any per-model recompute must **retrain**. |
| Committed analysis CSVs (efficiency_map, cp3_crlb_compare, g_ood_gating_voxels, f2_realdata, figS4) | present; per-voxel / per-grid-point. |
| In-vivo raw data (`download/Data/brain.*`, `*abdomen*`) | **not present locally** (fetched via `zenodo_get` record 14605039). In-vivo *recompute* needs download + a trained `.pt`. |
| Canonical NPE recipe | `train_npe.py --mode set --budget 500000 --epochs 200 --log-dstar --seed 0` (NSF, boxuniform, clinical_sparse 8-b = {0,50,100,200,400,600,800,1000}). Verified: summariser reproduces committed Table 2 (D* claimed 0.084/0.160/0.376/0.671; below-floor 0.86/0.76/0.63/0.52). |

## Runnability table (friction items)

| Item | Friction claim under attack | Runnability | Plan |
|---|---|---|---|
| **HC4 / CS3** (CP1) | "below-floor D* overconfidence is intrinsic" untested vs simulation starvation | **runnable-with-compute** | Retrain NPE at budgets {50k,100k,250k,500k,1M}; recompute per-SNR median D* claimed/achieved ratios. `run_cp1_budget_sweep.py`. |
| **HC1 / CS1** (CP2) | unbiased CRLB does not bound a biased estimator ⇒ "below floor ⇒ overconfident" is a category error (Hero–Fessler) | **runnable-now** (analysis) | van-Trees Bayesian / prior-regularized floor under the actual log-uniform prior; surviving below-floor fraction. `run_cp2_bias_aware_floor.py`. |
| **HC5** (CP3) | in-vivo NLLS spread/coverage silently excludes the 54.7% boundary-railed voxels | **runnable-now** (abdominal, from committed per-voxel CSV); brain coverage verified from code | Report railed-included vs excluded SD/IQR side by side; document brain inclusion policy. `run_cp3_railing_audit.py`. |
| **HC6** (CP4) | gate AUC 0.99 inflated by shared-fit correlation | **runnable-with-compute** (sim, truth-based); 3-way in-vivo split = needs-data | Re-score gate vs independent targets (parameter-recovery error, CI miscoverage) in simulation. `run_cp45_decorrelate_control.py`. |
| **HC2 / CS2** (CP5) | signal-domain coverage ≠ parameter calibration; in vivo confounds width vs misspecification | **runnable-with-compute** (in-silico control, RNPE on sim null + committed in-vivo scores); real-parameter-truth in vivo = **needs-external-data** | In-silico parameter coverage (no misfit) + RNPE-style model criticism. FLAG phantom/repeat-acquisition gap. |
| **G1–G4 / CS4** (CP6) | novelty vs documented amortized-overconfidence / recalibration-degradation literature | **manuscript-only** | Concede Hermans 2022 / Cannon–Ward 2022 / Frazier 2020 / Ovadia 2019; sharpen contribution; subtle-vs-gross scoping; add 6 verified refs. |

## Hard-constraint posture (unchanged from assessor run)
- Numbers frozen to `NUMBERS_FROZEN.txt`; `paper/check_numbers.sh` gate (no baseline count may decrease).
- Legitimate recomputes (CP1 budget, CP2 bias-aware, CP3 railed-included, CP4 decorrelated AUC) carry OLD + NEW, both labelled.
- Adverse / claim-weakening outcomes recorded verbatim in `ADVERSE_RESULTS.md`; do not spin.
- One commit per checkpoint; no force-push.

## Out of scope (cannot be closed by CC)
- **HC2/CS2 real in-vivo parameter ground truth** (phantom or repeated-acquisition): needs external data CC cannot manufacture. Mitigated by the in-silico control; FLAGGED.
