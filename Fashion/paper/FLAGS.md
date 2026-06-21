# FLAGS — items needing the author's decision or confirmation

No `[[TODO]]` placeholders were inserted in the manuscript: every required identifier was
recovered from the repository. The items below are decisions and provenance confirmations,
not blockers.

## DECISIONS (author's call)

- **`[[DECISION: disclose prior MRM review? author's call]]`** — The cover letter
  (`paper/cover_letter.tex`) deliberately makes **no** mention of the prior *Magnetic
  Resonance in Medicine* review. Whether to proactively disclose that history to the *MRI*
  editor is the author's decision; a flagged source comment marks where a short factual
  paragraph would go if desired.

## CONFIRMATIONS (recovered, but please verify before submission)

- **In-vivo dataset citation [22].** Both in-vivo datasets (brain N=500 / ~2,000-voxel gate;
  abdomen N=1) trace to a single source: `utilities/data_simulation/Download_data.py` fetches
  **only** Zenodo record **14605039** (`OSIPI_TF24_data_phantoms.zip`) and unzips it into
  `download/Data/`, from which the analysis code reads both `brain.*` and `abdomen.*`. Cited as
  *Gurney-Champion O, et al. Data for the OSIPI TF2.4 IVIM-MRI Code Collection. Zenodo; 2025.
  doi:10.5281/zenodo.14605039* (authors/year from the Zenodo record metadata).
  - *Confirm:* (a) the author list/title match the form you want (this is the version DOI;
    a concept DOI may be preferred); (b) if the OSIPI bundle re-distributes the brain/abdomen
    acquisitions from a primary acquisition source, you may wish to additionally cite that
    primary source. The archive is named `..._phantoms.zip` but the manuscript treats brain
    and abdomen as in-vivo acquisitions; the citation is to the archive that the code actually
    consumes.

- **Reproducibility DOI in the Data Availability Statement.** `doi:10.5281/zenodo.20649669`
  was taken from the `Fashion/README.md` DOI badge. Confirm this is the intended public
  archive version for the code, calibration grid, and efficiency map.

## R2 significance axis — substantially addressed (follow-up in-vivo work)

The highest-severity reviewer concern (in-vivo evidence needs to go beyond a single subject)
has now been **substantially addressed** by adding an open **multi-subject liver cohort**
(N = 59 patients; Dryad 10.5061/dryad.xwdbrv1cg), on which the amortized NPE under-covers
held-out b relative to per-voxel NLLS in 100% of patients (main-text "Multi-subject cohort
confirmation"; Supplementary §S6 / Figure S6). A from-scratch retrain also reproduces the
manuscript's single-subject in-vivo numbers. See CHANGES.md ("In-vivo validation").

**What still remains open (honest residual):**
- **No in-vivo reference standard.** In vivo lacks voxel-level ground truth, so the cohort
  compares NPE vs NLLS, not against truth; a cohort with a reference standard (or a controlled
  phantom cohort) is still the ideal closing evidence.
- **NLLS is below nominal in vivo** (~0.78 at 0.95) due to biexponential model misfit and the
  liver dataset's per-series-b0 assembly noise; the result is framed as the NPE-vs-NLLS gap,
  not as NLLS being perfectly calibrated.
- The cohort is liver (GE 1.5T); a same-organ multi-subject brain cohort matching the trained
  context was not openly available at the time of writing (see the dataset search notes).
- Trained checkpoints (`*.pt`) and raw data are gitignored; results require the documented
  retrain + download to reproduce bit-for-bit (the science reproduces; exact values vary with
  training stochasticity).
