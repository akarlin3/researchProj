# Checkpoint 0 — Enumerate & diagnose

Project Fashion manuscript: *"Calibration and Efficiency of Uncertainty Estimates in
Intravoxel Incoherent Motion Imaging: Quantile Intervals, Cross-Paradigm Comparison, and a
Cramér–Rao Audit of Amortized Posteriors"* (Avery Karlin).

## Repo inventory (manuscript-relevant)

| Artifact | Location | Notes |
|---|---|---|
| Manuscript source | **did not exist** | Only flattened PDF exports in `~/Downloads/` (master likely Google Docs). Reconstructed here as `paper/manuscript.tex`. |
| Supplement source | **did not exist** | Reconstructed as `paper/supplement.tex` (figs S1–S5). |
| Bibliography | **did not exist** | Reconstructed as `paper/refs.bib` (21 entries verbatim). |
| Figures | `figures/manuscript/*.pdf` | All 9 present (fig1–4, figS1–S5). |
| Cover letter | **did not exist** | To be authored in CP6. |
| Build script | `paper/Makefile` (new) | tectonic-based. |
| Analysis code | `npe/`, `uq/` | NPE + UQ campaign. |
| In-vivo data loader | `utilities/data_simulation/Download_data.py` | `zenodo_get https://zenodo.org/records/14605039`. |
| Reviewer responses | `REVIEWER_RESPONSE.md`, `REVIEWER_RESPONSE_R2.md` | Reference `Fashion_Manuscript_Merged_MRM.pdf`. |

Source PDFs used for reconstruction (newest):
`~/Downloads/Fashion_Manuscript_Merged_MRM (16).pdf` (13 pp, 2026-06-13),
`~/Downloads/Fashion_Supplementary_Information (5).pdf` (5 pp, 2026-06-12).

## Located identifiers / datasets

- **Open-reproducibility artifact (Zenodo DOI):** `10.5281/zenodo.20649669`
  (URL `https://doi.org/10.5281/zenodo.20649669`).
  Source: `Fashion/README.md` DOI badge. **Recovered.**
- **In-vivo brain dataset** (held-out-b N=500; OOD gate ~2,000 GM voxels):
  fetched in `Download_data.py` via `zenodo_get https://zenodo.org/records/14605039` →
  **Zenodo record 14605039** (OSIPI TF2.4 IVIM reference data; files `brain.nii.gz`,
  `brain_mask_gray_matter.nii.gz`, `brain.bval`). **Record ID recovered**; exact
  authors/title for a polished citation would need a Zenodo metadata fetch — cited by
  record/accession.
- **In-vivo abdominal case** (N=1, Supp Fig S4): described only as "one open IVIM
  acquisition." No source/accession in code or PDF. → `[[TODO: abdominal dataset
  citation]]` + FLAG.

## Diagnosis table

| # | Issue (assessor) | File | Lines (manuscript.tex) | Planned fix | Needs external info? |
|---|---|---|---|---|---|
| 1 | Data Availability Statement missing (Format FAIL, Medium) | manuscript.tex | after Discussion / before References | Add DAS → Zenodo 20649669 (code, calibration grid, efficiency map, OSIPI offer) | **N** (DOI in hand) |
| 2 | In-vivo datasets unnamed/uncited (Provenance FAIL, Medium) | manuscript.tex (Methods, Fig 4 caption); supplement.tex (S3, S4) | Methods Part 2; captions | Name+cite brain (Zenodo 14605039), add bib entry; abdominal → TODO+FLAG | **Partial** (abdominal Y) |
| 3 | N=500 vs ~2,000 unreconciled (Norm gap / R1 Minor) | manuscript.tex (Methods, Robustness); supplement.tex (S3) | held-out-b clause + S3 | State 500 = held-out-b subset of the ~2,000-voxel GM ROI feeding the gate; verify subset in code | **N** (verify in code) |
| 4 | R2 reframe: contribution & positioning (R2) | manuscript.tex (Intro end; Discussion) | end of Intro; new positioning para | Transferable-contribution framing; explicit scope sentence; 1 positioning para w/ one-sentence deltas vs [19][20][21]; zero new numbers | **N** |
| 5 | Prose density / long sentences (R3 Minor; MRM return cause) | manuscript.tex (Intro, Discussion) | longest compound sentences | Split, meaning-preserving | **N** |
| 6 | Venue retarget → MRI (Elsevier) | manuscript.tex title page; new cover_letter.tex | title; declarations; refs | Title swap; fresh cover letter to EIC J.C. Gore; Elsevier declarations block; ref style → Elsevier numeric | **N** |

## Hard-constraint posture

- Numbers frozen to `Fashion/NUMBERS_FROZEN.txt` (multiset over rendered PDFs).
  Gate rule: no baseline number's count may decrease (`paper/check_numbers.sh`).
- 21 bibliography entries transcribed verbatim; only dataset citations to be *added*.
- One commit per checkpoint on branch `fix/assessor-remediation`.

## Out of scope for this run (carried to final summary)

The highest-severity assessor item — **R2's significance axis (needs new in-vivo
validation data)** — is **out of scope**. Checkpoint 4 mitigates the positioning but does
not close it.
