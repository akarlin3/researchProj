# CP3 validation -- repeatability scale-calibration on ACRIN-6698 (PROVISIONAL)

n=76 test-retest tumors, level=0.1 (90% intervals), seed=20260613, conformal backend: caliper.conformal.conformal_offset.
Analytic repeat-coverage target for a correctly measurement-scaled interval: **0.755** (NOT 0.90).

## Headline (single-voxel SNR = 20, eff SNR = SNR x sqrt(n_vox))

| param | coverage [BCa 95%] | scale ratio R=σ_impl/σ_rep | z-disp | σ_implied | σ_repeat |
|---|---|---|---|---|---|
| D | 0.263 [0.158,0.355] | 0.247 | 3.95 | 4.92e-05 | 1.99e-04 |
| Dstar | 0.797 [0.661,0.864] | 19.544 | 0.05 | 8.66e-02 | 4.43e-03 |
| f | 0.421 [0.303,0.526] | 0.469 | 2.64 | 2.77e-02 | 5.91e-02 |

## Sensitivity: coverage(D) vs single-voxel SNR (honesty curve)

| single-voxel SNR | 10 | 20 | 30 | 50 | 100 |
|---|---|---|---|---|---|
| coverage(D) | 0.434 | 0.263 | 0.158 | 0.105 | 0.053 |

Crossing (single-voxel SNR where coverage(D)==0.755): outside grid.

Gauge Sec 4.2.2 baseline (RANK): D Spearman r=+0.60 [0.42,0.72]; D* null. Echo measures SCALE (coverage), an independent axis.

## GATE (parameter D): **LETHE**
- coverage_D = 0.263, CI [0.158,0.355]; under=True, saturated=False, discriminative=True
- R(D) = 0.247  (R<1: interval narrower than repeatability; R>1: wider)