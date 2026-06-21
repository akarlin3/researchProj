# Clean-room / IP statement

Gnomon is an **independent reimplementation** of Fashion's IVIM calibration ruler.
Independence is the entire point: a rebuild that shares no code with Fashion (or
with Caliper, which packages Fashion's method) cannot inherit Fashion's internal
inconsistencies, so a reproduce-or-refute test on it is meaningful.

## What "clean-room" means here

- **The methods are reimplemented from spec, not from source.** The IVIM
  bi-exponential forward model, the box-constrained NLLS railing diagnostic, the
  Laplace and per-voxel MCMC posteriors, the MAF/NPE amortized posterior, and the
  coverage / ECE / sharpness scorer are written from the published equations and
  definitions. The original implementation lives in `gnomon/`.
- **Fashion's *numbers* may be read as reproduction targets — only from its prose.**
  Every target in `gnomon/manifest.py` cites a Fashion `.md` file (README /
  REVIEWER_RESPONSE), never a `.py`. A static test (`test_every_target_cites_
  fashion_prose_not_code`) enforces this.
- **Caliper's ruler module is OFF-LIMITS as an import.** `caliper.metrics`,
  `caliper.baselines`, `caliper.conformal`, and the Fashion-reproduction example are
  Fashion's method as code; importing any of them would defeat the rebuild. This is
  enforced three ways: `gnomon/_paths.py` allows **only** `lattice` and never adds
  Caliper to the path; `_paths.FORBIDDEN = ("caliper",)` with `assert_no_caliper()`;
  and a static AST test (`test_no_gnomon_source_imports_caliper`).
  *Standard metrics may be — and are — re-derived independently from their
  references (Gneiting & Raftery 2007; quantile-calibration literature).*

## Data: synthetic + open only, in tree

- **No clinical / patient / scanner data is committed** — no DICOM, no NIfTI, no
  `pancData3`, no MSK data. (This file *names* those datasets only to state their
  exclusion; a test asserts no data-like files exist in the tree.)
- **Synthetic** cohorts come from the **Lattice** DRO (sibling, MIT, synthetic-only),
  imported **read-only**. Gnomon reimplements its own forward model for the
  estimators/CRLB; it uses Lattice only as the ground-truth *data* substrate.
- **Open real** data — the OSIPI TF2.4 abdomen acquisition (Zenodo **14605039**,
  CC-licensed) — is fetched **on demand** by `gnomon/osipi.py` into a **gitignored**
  `download/` directory with a provenance manifest, and is **never redistributed**
  in-tree (mirroring Lattice's and Lethe's posture). Importing the module touches no
  network.

## One-way dependency

```
Gnomon  ──imports (read-only)──►  Lattice        (synthetic substrate)
Gnomon  ──reads numbers from──►   Fashion *.md*   (reproduction targets, prose only)
Gnomon  ──MUST NOT import──►       Caliper        (Fashion's method as code)
Gnomon  ◄──imported by── (nothing)
```

The `gnomon` package imports only the standard library + numpy at its core
(scipy/torch/requests are optional extras for NLLS-MCMC / MAF / OSIPI fetch). It
imports **nothing** from Caliper or any paper project.
