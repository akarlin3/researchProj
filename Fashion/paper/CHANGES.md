# CHANGES — assessor-remediation run

Branch `fix/assessor-remediation`. One commit per checkpoint. Every edit verified against the
number-freeze gate (`paper/check_numbers.sh`) — no reported result value changed.

All manuscript/supplement line references are to the files as of the final commit.

| CP | Issue (assessor) | Fix | File / location |
|----|------------------|-----|-----------------|
| 0 | No editable source existed | Reconstructed canonical LaTeX (manuscript + supplement + refs.bib) verbatim from the latest merged MRM PDF; froze all numbers | `paper/manuscript.tex`, `paper/supplement.tex`, `paper/refs.bib`, `Fashion/NUMBERS_FROZEN.txt` |
| 1 | Data Availability missing (Format FAIL, Medium) | Added a Data Availability Statement pointing to the Zenodo reproducibility archive (`doi:10.5281/zenodo.20649669`); lists code, calibration grid, efficiency map; notes OSIPI offer | `manuscript.tex` — Data availability section (Declarations block) |
| 2 | In-vivo datasets unnamed/uncited (Provenance FAIL, Medium) | Named + cited both in-vivo datasets to their single recoverable source, the OSIPI TF2.4 IVIM-MRI Code Collection data (Gurney-Champion et al., Zenodo 2025, `doi:10.5281/zenodo.14605039`); new reference `[22]`; markers at Methods Part 2, Results (brain + abdomen), Limitations, Fig 4 caption; accession named in Supp S3 and S4 | `manuscript.tex` (ref [22]; 5 in-text markers), `supplement.tex` (S3, S4) |
| 3 | N=500 vs ~2,000 unreconciled (Norm gap / R1 Minor) | Verified in code (`npe/run_f_realdata.py`, `npe/run_g_ood_gating.py`): same brain image, same gray-matter mask, same held-out-b partition, same RNG seed (42), differ only in `--n-voxels`. The 500 is **not** a strict subset of the 2,000 for plausible voxel counts, so stated the accurate relationship: overlapping random draws from the same gray-matter ROI, sized differently because a per-voxel ROC (Supp S3) needs more voxels than a coverage curve (Fig 4B) | `manuscript.tex` Methods Part 2; `supplement.tex` S3 |
| 4 | R2 contribution/positioning (prose only) | Reframed to the transferable diagnostic (information-floor / CRLB auditing as a general check for learned qMRI posteriors; aggregate-passes/pointwise-fails as the portable lesson; NPE as worked example). Added a positioning paragraph with one-sentence deltas vs Casali 2026 [19], µGUIDE/Jallais & Palombo 2024 [20], Manzano-Patrón 2025 [21]. Stated explicit scope early (end of Introduction). Zero new numbers, only existing citations | `manuscript.tex` — Introduction |
| 5 | Prose density / long sentences (R3 Minor; MRM return cause) | Split four of the longest compound sentences (Introduction + Discussion) at em-dash/semicolon joins; meaning and numbers preserved exactly | `manuscript.tex` — Introduction, Discussion (amortized-posteriors, limitations) |
| 6 | Venue retarget → MRI (Elsevier) | Title page venue swap (MRM → Magnetic Resonance Imaging (Elsevier)); Elsevier declarations block (CRediT, Declaration of competing interest, Funding, Data availability, Declaration of generative AI — reusing the existing AI-disclosure text); reference style → Elsevier bracketed numeric ([N] in-text + numbered list); fresh cover letter to EIC John C. Gore | `manuscript.tex` (title, declarations, citations); `paper/cover_letter.tex` (new) |

## What was deliberately *not* done

- **No fabrication.** Every identifier (repro DOI, dataset accession) was recovered from the
  repository; **no `[[TODO]]` placeholder was needed** in the manuscript.
- **No numerical change.** All coverage values, CRLB ratios, percentages, Ns, AUCs, and the
  21 reference entries are byte-for-byte as transcribed from the source PDF.
- **Reference style:** in-text citations are Elsevier bracketed numeric and the list is
  numbered, but the 21 entries' internal text was kept verbatim rather than re-flowed into
  the `elsarticle-num` micro-format (initials-first), to honor the "content unchanged"
  constraint. `elsarticle.cls` was verified to build under `tectonic`; a production
  Elsevier-class typesetting can be generated from the committed `refs.bib` without touching
  entry content. The structured abstract was preserved (which `elsarticle`'s abstract
  environment does not natively support), another reason the working `article` build was kept.

## Build

`cd Fashion/paper && make` → `manuscript.pdf`, `supplement.pdf`, `cover_letter.pdf`
(tectonic; clean, no errors, all cross-references resolve).

---

# In-vivo validation (follow-up work, beyond the assessor remediation)

This addresses the highest-severity open item (R2's significance axis: in-vivo evidence).
Unlike CP0–CP6 (prose/formatting, zero new numbers), this is **new analysis and therefore
produces new numbers** — all computed by running code on real open data; nothing fabricated.

| Step | What was done | Artifacts |
|------|---------------|-----------|
| Env + data | Built the NPE venv (`~/ProjectFashion/.venv-npe`: torch 2.12, sbi 0.26.1, nibabel, dcm2niix); downloaded the OSIPI brain+abdomen data (Zenodo 14605039) | (gitignored data/venv) |
| Retrain | Retrained the setB NPE from scratch (set/nsf/log-D*/clinical-sparse/500k budget, seed 0; 99 epochs, matching the original early-stop) | `npe/loss_setB_retrain.json`; `npe_posterior_setB.pt` (gitignored) |
| **Reproduce (option 3)** | A from-scratch retrain **reproduces the manuscript's headline in-vivo numbers**: brain held-out-b coverage **NPE 0.031 / NLLS 0.904** at nominal 0.95 (paper: 0.03 / 0.90); OOD gate **self-consistency AUC 0.993 vs χ² 0.594**, 49% retained, residual SD 11.0 (paper: 0.99 / 0.59) | regenerated `npe/f2_realdata.csv`, `npe/g_ood_gating*.csv` |
| **Cohort (option 1)** | New **multi-subject in-vivo cohort**: open hypovascular-liver IVIM dataset (Dryad 10.5061/dryad.xwdbrv1cg; GE 1.5T; 77 patients, 59 analyzable). Off-scheme held-out-b coverage, NPE vs per-voxel NLLS. **NPE under-covers NLLS in 100% of patients (59/59); cohort-mean 0.58 (NPE) vs 0.78 (NLLS) at nominal 0.95.** | `npe/run_liver_cohort.py`, `npe/make_liver_cohort_figure.py`, `npe/liver_cohort_coverage.csv`, `npe/liver_cohort_summary.json`, `figures/manuscript/figS6_liver_cohort.{png,pdf,csv}` |
| Write-up | New main-text "Multi-subject cohort confirmation" paragraph; Supplementary §S6 + Figure S6; updated Limitations and Data Availability; new reference [23] | `manuscript.tex`, `supplement.tex` |

**Methodological notes / honest caveats (also in the manuscript):**
- In vivo has no voxel-level ground truth, so the cohort compares the two **estimators**
  (NPE vs NLLS), not either against truth. NLLS is itself below nominal in vivo (model
  misfit + inter-series assembly noise); the load-bearing result is the systematic
  NPE-below-NLLS gap, which is robust to the shared noise model.
- The liver dataset stores each b-value as a separate two-volume series; per-series b=0
  normalization, filename-derived b-values, a geometry-consistency filter (8 patients
  excluded), and a robust median-residual SNR (~13) are documented in `run_liver_cohort.py`.

**Reproduce:** install `npe/requirements.txt` (+ `zenodo_get nibabel dcm2niix`) into a venv;
`python utilities/data_simulation/Download_data.py`; `zenodo_get 4408313` + `dcm2niix` per
patient; `python npe/train_npe.py --mode set --budget 500000 --density-estimator nsf
--log-dstar --b-scheme clinical_sparse --output npe/npe_posterior_setB.pt ...`; then
`python npe/run_f_realdata.py`, `npe/run_g_ood_gating.py`, and `npe/run_liver_cohort.py`.

The number gate still PASSES for all pre-existing results (no reported value decreased); the
cohort numbers are additions.
