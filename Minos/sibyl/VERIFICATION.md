# Tier 2 pre-coding verification — ACRIN-6698 public TCIA manifest

The directive required verifying the *actual* public manifest before any coding.
Two independent web-research passes (TCIA collection page + data descriptor +
Partridge 2018 / Newitt 2019 repeatability papers) were run. Verdict: **gate
PASSED**, with three design-shaping nuances folded into the implementation.

## The three required checks

### 1. Are the repeat (test/retest) raw DWI series public — or only ADC maps?
**CONFIRMED: raw multi-b source DWI is public, not ADC-only.** The collection ships
raw acquired DWI (TRACE / MergedMSMB series) *and* derived ADC maps. The test/retest
"coffee-break" repeatability sub-study is a **separate download option in the full
collection**: **71 analyzable same-session pairs** (60 at T0, 11 at T1; 89 consented).
→ No fallback to ADC-maps-only is needed. Implemented: `acrin_ingest` reads the
per-b source signal (TRACE/MSMB) and computes ADC itself.

### 2. Are QA analyzability flags present in the public metadata?
**CONFIRMED: yes**, in the *"Full Collection Ancillary Patient Information"* XLSX
(DWI QC ratings: fat-suppression / artifacts / SNR → poor·moderate·good, plus a
usable-for-ADC pass/fail and an ROI-confidence 1–3). **Not in DICOM.** The one item
not verifiable from public web text is the *exact column header strings* (they live
inside the XLSX). → Implemented: `acrin_ingest.load_qa_flags()` parameterizes the
column names via `col_map` (defaults in `DEFAULT_COL_MAP`), to be confirmed against
the downloaded file.

### 3. How is b0 reconciled across the three 2-b mini-acquisitions?
**SETTLED.** Most sites acquire a single 4-b series (one b0). One site acquired three
separate 2-b series (0/100, 0/600, 0/800); the trial reconciles them into a
**MergedMSMB series = all non-zero b images + a single _averaged_ b0**. So the
canonical 4-b vector is `[mean(b0s), b100, b600, b800]`. → Implemented exactly in
`acrin_ingest.reconcile_b0()` / `assemble_4b_vector()` (averaging when >1 b0 is
present); if the derived TRACE/MSMB series is consumed, the averaging is already done.

## Three findings that shaped the design (all anticipated by the directive's hedges)

- **Voxel-level test/retest is NOT supported.** Test and retest ROIs were drawn
  *independently with no cross-referencing*, with repositioning between scans, and
  are not co-registered. → We **stay ROI-level** (the directive's default unit:
  whole-tumor ROI mean ADC, ~per patient). Voxel-level would require self-registration
  and is intentionally omitted.
- **No published numeric SNR.** "SNR" in ACRIN is only a categorical QA rating. → The
  synthetic ID-reference noise level is a **free parameter** (`ACRIN_REF_SNR`, to be
  matched to the empirical b0 SNR / swept), not hardcoded from literature.
- **ADC convention is fixed:** mono-exponential, log-linear least squares over **all
  four** b-values 0/100/600/800 (Partridge 2018; Newitt 2019). → `adc_monoexp_fit`
  uses all four. Published ROI-mean repeatability benchmark: **wCV ≈ 4.8–5.4 %,
  RC ≈ 0.16×10⁻³ mm²/s, ICC ≈ 0.97**.

## Data access
Fully public (CC BY 4.0) via **NBIA Data Retriever** + `.tcia` manifest. Full
collection ~1.94 TB; the test/retest subset is a much smaller separate download.
385 subjects total (242 primary analysis, 71 test/retest). The imaging is **not**
present in this environment — see `docs/tier2.md` for the download + ingestion path.

## Sources
- TCIA ACRIN-6698 collection: https://www.cancerimagingarchive.net/collection/acrin-6698/
- Data descriptor PDF (Newitt, 2021-05-20): ACRIN-6698-ISPY2-DWI-and-DCE-MRI-Data-Descriptions
- Partridge et al., *Radiology* 2018 (PMC6283325) — ADC definition, primary cohort.
- Newitt et al., *JMRI* 2019;49:1617 (PMC6524146) — test/retest repeatability (wCV/RC/ICC).
