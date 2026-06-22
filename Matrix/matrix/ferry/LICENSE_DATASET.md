# Ferry dataset — license & provenance record (clean-IP gate)

Ferry grounds Matrix's substrate on a **public** RT dataset with a **verified open license**.
**No patient data or image blobs are committed to this repository** — the loader
(`dataset.py`) downloads by script from the TCIA NBIA REST API into a git-ignored cache and
keeps only tiny derived label/dose grids (also git-ignored). This file records exactly what is
used and under what terms.

## Dataset

- **Name:** Pancreatic-CT-CBCT-SEG — *Pancreatic CT, CBCT, and structures with delivered dose
  for locally advanced pancreatic cancer ablative radiotherapy*
- **Host:** The Cancer Imaging Archive (TCIA), accessed via the NBIA REST API
- **DOI:** `10.7937/TCIA.ESHQ-4D90`
- **Landing page:** https://www.cancerimagingarchive.net/collection/pancreatic-ct-cbct-seg/
- **Version used:** **Version 2** (2022-08-23) — the version that adds the **RTDOSE** objects
  (Version 1 had CT + RTSTRUCT only).

## License (verified)

- **Data license:** **Creative Commons Attribution 4.0 International (CC BY 4.0)** — the Data
  License field on the TCIA collection page, corroborated by the *Scientific Data* data
  descriptor: *"licensed under a Creative Commons Attribution 4.0 International License."*
- **Terms:** redistribution and reuse (including **commercial** use) permitted **with
  attribution**; no non-commercial or share-alike restriction; **no registration / data-use
  agreement gate**.
- License text: https://creativecommons.org/licenses/by/4.0/

## What Ferry pulls (and what it does NOT)

- Pulls, for **one patient** (default `Pancreas-CT-CB_001`): the **RTSTRUCT** (clinician
  contours) and **RTDOSE** (delivered 3-D dose grid) series only — a few tens of MB.
- Does **NOT** pull the CT / CBCT image volumes (not needed: contours are rasterised directly
  onto the RTDOSE grid, which shares the structures' Frame of Reference).
- Does **NOT** commit any DICOM, pixel data, or derived grid to git.

## Real vs synthetic (honest ceiling)

- **REAL:** anatomy (target VOI + abdominal OAR contours → twin `labels`) and dose geometry
  (the delivered RTDOSE grid → twin `dose`).
- **SYNTHETIC:** all IVIM/perfusion `(D, D*, f)`, the scan/noise model, the SNR map, and the
  high-D* sub-region — there is no scanner, so there is **no real diffusion/perfusion data**.

## Required attribution / citation

> Hong, J., et al. *Pancreatic-CT-CBCT-SEG* (Version 2) [Data set]. The Cancer Imaging Archive,
> 2022. DOI: 10.7937/TCIA.ESHQ-4D90.
>
> Clark, K., et al. *The Cancer Imaging Archive (TCIA): Maintaining and Operating a Public
> Information Repository.* Journal of Digital Imaging, 26(6), 1045–1057 (2013).

Verify the current license text on the landing page before any redistribution; TCIA records
the license per collection version.
