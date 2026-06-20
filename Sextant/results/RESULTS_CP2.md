# RESULTS — CP2: boundary-railing as the primary analysis

Seed 20260613, 5000-bootstrap. All numbers from `results/railing_results.json`
(regenerate with `python scripts/run_railing.py`). Cohort: OSIPI open
human-abdominal IVIM (CC-BY-4.0), curated homogeneous ROI = the original cohort on
which Fashion first computed the railing figure.

## Primary result (original cohort `abdomen_homogeneous`)

* High-SNR voxels (SNR ≥ 8): **1 618** (exactly matches Fashion's reported count).
* **NLLS D\* railing fraction: 0.5470 (54.70%)** — reproduces Fashion's 54.7%.
* 95% bootstrap CI: **[0.5222, 0.5711]**.
* Rail direction: lower bound 15.6%, **upper bound 39.1%** (high-D\* wall dominates).
* Wide-bounds sensitivity: **0.3214** — still a substantial minority with generous
  bounds, so railing is *not* an artefact of the tight prior box.

## When and why it rails (characterisation)

Railing fraction by SNR stratum (homogeneous ROI):

| SNR band | n | railed | lower | upper |
|---|---:|---:|---:|---:|
| [8, 15) | 633 | 0.542 | 0.013 | 0.529 |
| [15, 30) | 930 | 0.531 | 0.224 | 0.308 |
| [30, 60) | 55 | 0.873 | 0.655 | 0.218 |

Railing **does not vanish with SNR** — it is pervasive across all strata, which is
the signature of *structural* weak identifiability of D\*, not of measurement noise
alone. At low SNR the optimiser is pushed to the **upper** D\* bound; in
high-SNR homogeneous tissue (genuinely near-zero perfusion) D\* collapses to the
**lower** bound. Either way the fit terminates on a boundary rather than at an
interior optimum.

## Why this is the stronger primary claim

The fraction of fits that rail is computed entirely from the optimiser's output
and the data — **no ground truth, no noise-model trust assumption**. It therefore
cannot be "overextended", which directly answers the reviewer critique that demoted
the calibration ruler. The ruler is retained as a scoped secondary diagnostic
(`sextant.ruler`), explicitly limited to data with ground truth.
