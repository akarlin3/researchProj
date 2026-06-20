# ASSUMPTIONS — read-only reuse wiring & dependency posture

Mirrors `Minos/future/ASSUMPTIONS.md`: the single place that records exactly what
Sextant borrows from sibling subrepos and how.

## Read-only reuse map

| Sextant uses | from | file | how |
|---|---|---|---|
| `fit_biexp_nlls`, `load_voxels`, `TARGET_BVALS`, `SNR_FLOOR`, `DSTAR_LOWER_RAIL`, `DSTAR_UPPER_RAIL` | Fashion | `Fashion/npe/run_s4_figure.py` | AST extraction (`sextant.fashion_reuse.load_railing`) — exec'd verbatim with `numpy/scipy/nibabel` only; torch/npe never imported |
| `fit_biexp_wide`, `WIDE_LOW/HIGH`, `_RAIL_TOL`, `_rail_fraction` | Fashion | `Fashion/npe/run_crlb_sampling_bound.py` | AST extraction (`load_wide`) |
| `coverage`, `ece`, `sharpness_rel`, `LEVELS` (secondary ruler) | Fashion | `Fashion/uq/calib.py` | `importlib` load by path (`sextant.ruler`) — numpy/scipy only |
| open-data fetch + provenance manifest pattern | Gauge | `Gauge/scripts/fetch_osipi.py` | mirrored in `scripts/fetch_osipi.py` |
| seed convention `20260613` | Gauge / Minos | — | `sextant.seeding.GLOBAL_SEED` |

Sextant **never writes** into the Fashion or Gauge trees. If Fashion's source
drifts, extraction fails loudly (missing name) or `tests/test_fashion_reuse.py`
catches it (pinned constants).

## Runtime dependencies (synthetic + analysis)

`numpy ≥ 1.23`, `scipy ≥ 1.9`, `nibabel ≥ 5.0` (NIfTI I/O). No torch required —
the reused railing definitions are torch-free. Tested on the `proteus` conda env
(Python 3.11, numpy 2.3.5, scipy 1.17.1, nibabel 5.4.2).

## Clean-IP posture

**Resolved at CP0; confirmed clean. No `pancData3`, no MSK clinical data, in tree
or history.** All imaging is download-on-demand under `data/` (git-ignored); only
provenance manifests and seeded result JSON are committed. See `VERIFICATION.md`.

## Data assumptions

* The OSIPI abdomen acquisition is N = 1 subject; population inference is **not**
  claimed from it. The railing fraction is a within-subject voxel statistic with a
  voxel-bootstrap CI; the cross-cohort comparison (homogeneous vs full vs, pending,
  TCGA-LIHC) is what supports generality.
* SNR is estimated from b-value replicates (Fashion's `load_voxels`); voxels with
  no usable replicates default to SNR = 30 (Fashion convention, reused verbatim).
