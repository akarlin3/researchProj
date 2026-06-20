# RESULTS — CP3: the human-abdominal replication gate

Pre-registered thresholds (see `VERIFICATION.md`): REPLICATES if point ≥ 0.30 and
95% bootstrap CI lower bound ≥ 0.20; SCOPE-DOWN if CI upper < 0.20.

## CP3a — internal generalisation (OSIPI full-abdomen ROI) — DONE

The original analysis used a curated *homogeneous* ROI. CP3a re-runs railing on the
**entire** abdominal ROI (`mask_full_temp`, 28 053 voxels → 19 652 high-SNR), a
broader and harder test on the same open data.

| metric | value |
|---|---|
| high-SNR voxels | 19 652 |
| railed (tight) | **0.4781 (47.81%)** |
| 95% bootstrap CI | **[0.4710, 0.4851]** |
| lower / upper rail | 0.046 / 0.432 |
| railed (wide bounds) | 0.2743 |
| **verdict** | **REPLICATES-STRONG** |

Railing generalises strongly beyond the curated ROI: ~48% of high-SNR fits across
the whole abdomen rail to a bound, with a tight CI, dominated by the upper
(high-D\*) wall, and surviving generous wide bounds (27%). This clears the
registered thresholds.

## CP3b — independent cohort (TCGA-LIHC liver DWI) — PENDING SIGN-OFF

A truly independent human-abdominal cohort was located and DICOM-verified:

* **TCGA-LIHC** (liver hepatocellular carcinoma), TCIA, DOI
  10.7937/K9/TCIA.2016.IMMQW8UQ, **CC BY 3.0**.
* Verified at the DICOM level: a real 4-b-value liver DWI series
  **b = 0 / 50 / 500 / 800 s/mm²** (spans the IVIM low+high range), human in-vivo,
  `BodyPartExamined: LIVER`.
* Non-pancreatic, non-MSK. Downloadable via the NBIA REST API (the same mechanism
  Gauge already uses): `getSeries?Collection=TCGA-LIHC&Modality=MR` →
  `getImage?SeriesInstanceUID=…`.

This requires a human license/posture sign-off (CC BY 3.0 + TCIA Data Usage
Policy) before download — the CP3 hard halt. Kidney fallbacks if needed: TCGA-KIRP
(3-b, CC BY 3.0), CPTAC-CCRCC (CC BY 4.0, thinner multi-b).

**Verdict so far:** the primary claim REPLICATES strongly on the open OSIPI
human-abdominal data (original ROI and full-abdomen generalisation). Independent
cross-dataset confirmation (TCGA-LIHC) is staged and awaits sign-off; if declined,
the claim is reported as demonstrated across the full OSIPI human-abdominal
acquisition with independent confirmation noted as future work (honest scoping, not
a failure).
