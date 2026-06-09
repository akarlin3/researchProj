# Tier 2 — in-vivo breast DWI (ACRIN-6698)

**Claim (bounded):** the detector flags the genuine synthetic → in-vivo distribution
shift, and the per-unit OOD flag tracks loss of ADC test-retest trustworthiness.
**D\*/f are not validated in vivo here** — they stay synthetic until the liver data
(Tier 3). See `VERIFICATION.md` for the pre-coding manifest check.

## Design

| | Arm 1 — matched-scheme (primary) | Arm 2 — imputation → dense (secondary) |
|---|---|---|
| Estimator | none (estimator-free) | dense 10-b µGUIDE (reused from Tier 1) |
| Feature space | ACRIN 4-b b0-normalized signal (3 informative decay ratios) | length-10 embedding + signal residual |
| Detectors | Family 2 only (kNN density-ratio / conformal) | Family 1 (Mahalanobis on embedding) **and** Family 2 (residual) |
| Why this split | Family 1 needs the length-10 embedder, which does **not** exist at the bare 4-b scheme — we state this asymmetry rather than forcing a degenerate 4-b embedder | the two-family comparison lives here |

**Imputation (Arm 2):** per unit, a quick segmented IVIM fit from the 4 real points
is forward-projected onto the 6 missing dense-grid b-values; the 4 real values are
spliced back in. **Labelled circularity:** the imputed points are generated *by the
IVIM model*, so a low Arm-2 score there can be the imputer fitting the model rather
than true in-distribution-ness. Hence Arm 2 is secondary, and the robustness readout
is ρ(Arm-1, Arm-2): agreement ⇒ missing-input handling isn't driving the result;
divergence ⇒ imputation is injecting IVIM-model artifacts.

**Trust reference (headline):** ADC test-retest repeatability, computed
mono-exponentially straight from the 4-b data (a clean external reference — *not*
from µGUIDE). Default unit: whole-tumor ROI mean ADC (~per patient), matching the
trial's coffee-break design and provided segmentations. Coupling = Spearman ρ between
the OOD score (test scan) and the per-unit ADC repeatability (|ΔADC| or, primary,
within-subject CoV).

**Controls:** AUROC on two proxy label axes (synthetic-ID vs real-in-vivo; QA-pass vs
QA-fail) since there is no in-vivo ground truth; degeneracy-FPR (low-f, low-ADC) to
confirm the detector isn't just flagging degenerate units; and an ADC-level partial
correlation to confirm the coupling isn't merely both quantities tracking ADC level.

## Running

```bash
export PYTHONPATH=$(pwd)

# Synthetic validation harness (default; no external data) — validates the machinery
# and the SIGN/logic of the coupling, NOT the in-vivo result.
python3 -m sibyl.experiments.tier2 --results-dir results/tier2 --dense-epochs 80

# Arm-1 only (fast, no µGUIDE):
python3 -m sibyl.experiments.tier2 --no-arm2

# Real ACRIN units (after building a UnitTable from the download, see below):
python3 -m sibyl.experiments.tier2 --units path/to/acrin_units.npz

# Tests
python3 -m pytest tests/ -q
```

## Synthetic validation harness — canonical result (800 units, 80-epoch dense)

The harness gives each in-vivo stand-in unit a latent "acquisition quality" that
drives **both** its SNR and its corruption probability, so OOD-ness and ADC
unrepeatability emerge from the *same* per-unit cause (mirroring the real structure:
bad acquisitions are both off-distribution and unrepeatable). Nothing is hand-set.

| Quantity | Value |
|---|---|
| **Arm 1** score: ID vs in-vivo | 0.13 vs 1.58 |
| **Arm 1** AUROC sim-to-real / QA | **0.84** / **0.90** |
| **Coupling** ρ(score, ADC CoV) | **0.46** (p=2e-42) |
| Coupling ρ(score, \|ΔADC\|) | 0.31 |
| Partial ρ(score, CoV \| ADC level) | **0.44** (coupling survives controlling ADC) |
| Degeneracy-FPR, ID-only (low-f / low-ADC; baseline 0.05) | 0.06 / 0.08 |
| **Arm 2** AUROC sim-to-real (F1 / F2) | 0.74 / 0.77 (**below** Arm 1 → imputation dilutes detection) |
| **Robustness** ρ(Arm1, Arm2 F1 / F2) | 0.64 / 0.83 (agreement) |
| Family agreement ρ(F1, F2) | 0.62 |

Note the harness cohort wCV (~19 %) is deliberately wide — it spans very low SNR to
exercise the full trust spectrum — and is **not** comparable to the real ACRIN
ROI-mean wCV (~4.8–5.4 %). This artifact validates logic and sign only.

## Real-data path

1. Download the ACRIN-6698 **test/retest subset** (and the whole-tumor segmentations
   + ancillary QA spreadsheet) via NBIA Data Retriever (`.tcia` manifest). Fully
   public, CC BY 4.0.
2. Build a `UnitTable`:
   - If you convert the DWI to NIfTI (4D + `.bval` + ROI mask), use
     `sibyl.data.acrin_ingest.build_units_from_nifti(cases)` — fully runnable.
   - To go straight from DICOM, implement the series-discovery sketch in
     `build_units_from_acrin_dicom` on top of the tested pure-numeric core
     (`reconcile_b0`, `roi_mean_per_b`, `assemble_4b_vector`).
   - Join QA flags with `load_qa_flags(xlsx, col_map=...)` — set `col_map` to the real
     header strings from the ancillary XLSX.
3. Run `run_tier2(units_path=...)`. The synthetic-ID reference, both arms, the trust
   coupling and the controls all apply unchanged; the `is_synth_id` axis becomes the
   sim-to-real label and `qa_pass` the analyzability positive control.
```
