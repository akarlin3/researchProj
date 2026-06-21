# Manuscript source — *Calibration and Efficiency of Uncertainty Estimates in IVIM*

Editable LaTeX source for the IVIM uncertainty-quantification manuscript and its
supplement.

## Provenance

This LaTeX source was **reconstructed from the merged submission PDF**
(`Fashion_Manuscript_Merged_MRM`) and supplement (`Fashion_Supplementary_Information`)
because no editable source had previously been committed to the repository. All prose,
numerical results, tables, and the 21 references were transcribed verbatim from that PDF.
`NUMBERS_FROZEN.txt` (repo root of `Fashion/`) records every reported quantity at the
reconstruction baseline so subsequent prose edits can be checked for numerical drift.

## Files

| File | What it is |
|---|---|
| `manuscript.tex` | Main manuscript |
| `supplement.tex` | Supplementary Information (figures S1–S5) |
| `refs.bib`       | Bibliography database (21 references + dataset entries) |
| `Makefile`       | Build via `tectonic` |

Figures are read from `../figures/manuscript/*.pdf` (committed).

## Build

```sh
make            # builds manuscript.pdf and supplement.pdf
# or
tectonic manuscript.tex
tectonic supplement.tex
```

`tectonic` fetches LaTeX packages on first run (needs network once). The bibliography is
rendered from an inline `thebibliography` block so the build needs no separate BibTeX pass;
`refs.bib` is maintained in parallel as the structured citation record.
