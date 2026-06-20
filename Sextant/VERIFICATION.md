# VERIFICATION — pre-coding gate (CP0) and registered thresholds

This file records the human-reviewable facts established **before** any analysis
code ran, mirroring `Minos/VERIFICATION.md`.

## CP0.1 — Data source is open and clean-IP

* **Primary cohort:** OSIPI TF2.4 open data, Zenodo record
  [14605039](https://zenodo.org/records/14605039), **CC-BY-4.0**. Archive
  `OSIPI_TF24_data_phantoms.zip`, MD5 `e7b3fe1d811a7a45c5aaf6c604c82793`
  (verified on download). The human-abdominal scan (`Data/abdomen.nii.gz`,
  144×144×21×104) is Philips 3T DWI, sourced from the IVIMNET repository
  (`Data/abdomen_readme.txt`). 12 unique b-values {0,10,20,30,40,50,75,100,150,
  250,400,600} with replicates; analysis uses Fashion's 10-pt evaluation subset.
* **Clean-IP confirmation:** `pancData3` is **absent from the working tree and
  from all git history** (no path, no commit message). The only `MSK` strings in
  the monorepo are (a) upstream open-source algorithm code (the MRI-QAMPER IVIM
  MATLAB package shipped inside the public OSIPI collection — an algorithm, not
  patient data) and (b) exclusion notes / synthetic "pancreatic-anchor" priors in
  Minos. No `.nii/.nii.gz/.dcm` is committed anywhere. Sextant is born inside the
  monorepo and never adds pancData3 — clean by construction.
* **Independent replication (pending sign-off):** TCGA-LIHC liver DWI (TCIA, DOI
  10.7937/K9/TCIA.2016.IMMQW8UQ, **CC BY 3.0**), a DICOM-verified 4-b-value liver
  series (0/50/500/800) — human-abdominal, non-pancreatic, non-MSK, downloadable
  via the NBIA REST API. Requires a license/posture sign-off before download.

## CP0.2 — Reuse is read-only

The boundary-railing computation is **not reimplemented**. Sextant loads
Fashion's exact `fit_biexp_nlls`, `load_voxels`, `TARGET_BVALS`, `SNR_FLOOR`,
`DSTAR_LOWER_RAIL`, `DSTAR_UPPER_RAIL` from `Fashion/npe/run_s4_figure.py` (and
the wide-bounds variant from `run_crlb_sampling_bound.py`) by AST extraction — the
source file imports torch + the npe package, so it cannot be imported directly,
but the dependency-light railing definitions are exec'd verbatim. The ruler is
loaded read-only from `Fashion/uq/calib.py`. `tests/test_fashion_reuse.py` pins
the constants so drift fails loudly.

## CP0.3 — Replication thresholds (registered BEFORE running)

The railing metric is the fraction of high-SNR (SNR ≥ 8) voxels whose NLLS D\*
estimate is within ε of a fit bound, using Fashion's exact bounds/thresholds;
bootstrap = voxel-level resample, **B = 5000**, **seed 20260613**, 95% percentile
CI.

* **REPLICATES (general claim):** point estimate ≥ **0.30** *and* bootstrap 95% CI
  lower bound ≥ **0.20**. "Strong" if point ∈ [0.40, 0.70] (brackets the original
  54.7%) or CI overlaps 0.547.
* **DOES NOT REPLICATE → scope down:** CI upper bound < 0.20 → the claim becomes
  "railing occurs under conditions X (SNR / b-scheme / organ)", reported honestly
  as a scoping result, not a failure.

These thresholds were fixed here before `scripts/run_railing.py` was executed.

## CP0.4 — Outcome (filled after running; numbers from `results/railing_results.json`)

| Cohort | n (high-SNR) | railed | 95% CI | verdict |
|---|---:|---:|---|---|
| `abdomen_homogeneous` | 1 618 | 0.5470 | [0.5222, 0.5711] | REPLICATES-STRONG |
| `abdomen_full` | 19 652 | 0.4781 | [0.4710, 0.4851] | REPLICATES-STRONG |

Both clear the registered thresholds. The independent TCGA-LIHC replication is the
CP3 hard-halt decision (license sign-off).
