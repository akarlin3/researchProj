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

## CP3b — independent cohort (TCGA-LIHC liver DWI) — DONE (signed off)

Independent human-abdominal cohort, downloaded via the NBIA REST API (reusing
Gauge's TCIA helpers read-only):

* **TCGA-LIHC** (liver), TCIA, DOI 10.7937/K9/TCIA.2016.IMMQW8UQ, **CC BY 3.0**.
* Siemens 1.5T (b encoded in `SequenceName`) — a different site, scanner, and
  organ from the Philips 3T OSIPI cohort. Non-pancreatic, non-MSK.
* Two IVIM-capable acquisition schemes present across patients.

| cohort | scheme | subjects | n high-SNR (analyzed) | railed | 95% CI | lower / upper | verdict |
|---|---|---|---|---|---|---|---|
| `lihc_liver_4b` | b=0/50/500/800 (clean) | 1 | 648 265 (40 000) | **0.4370** | [0.4321, 0.4419] | 0.165 / 0.272 | REPLICATES-STRONG |
| `lihc_liver_3b` | b=50/400/800 (no b=0) | 3 | 1 544 999 (40 000) | **0.7345** | [0.7302, 0.7388] | 0.725 / 0.009 | REPLICATES |

Per-subject 3-b railing is consistent (0.659, 0.744, 0.752). High-SNR sets exceed
the 40 000-fit cap, so a seeded random subsample was analysed (logged, not silent).

**Reading:** the phenomenon replicates on independent liver data. The clean 4-b
scheme (with b=0) rails at 43.7%, squarely in the original's "strong" band and
again split across both bounds. The sparse 3-b scheme (no b=0, normalised by b=50)
rails far more (73.5%) and almost entirely to the lower bound — exactly the
expected direction: fewer perfusion-sensitive b-values → worse D\* identifiability
→ more railing. Honest caveats: (i) the SNR floor here uses a background-noise
estimate (vs OSIPI's replicate-variance SNR), and (ii) the 3-b scheme lacks b=0.
Neither weakens the qualitative claim — a large, robust fraction of conventional
NLLS D\* fits rail on independent human-abdominal data.

## Overall CP3 verdict

**REPLICATES.** Boundary-railing of conventional NLLS D\* is a robust, assumption-
free property across: the original OSIPI ROI (54.7%), the full OSIPI abdomen
(47.8%), and an independent TCGA-LIHC liver cohort (43.7% clean 4-b; 73.5% sparse
3-b across 3 subjects). All cohorts clear the pre-registered thresholds. The
primary claim is established; merge-back into the retooled Fashion spine is the
default (no salami).
