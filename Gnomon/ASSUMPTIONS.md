# ASSUMPTIONS.md — what is SOLID vs PROVISIONAL in Gnomon

> **Gnomon is a reproduce-or-refute test, not a result yet.** Until CP3 renders the
> verdict, the only solid things are the *scaffold and the contract*: the embedding,
> the clean-room boundary, and the **frozen** targets manifest. Every *number Gnomon
> produces* is PROVISIONAL until the CP3 run, and the verdict itself can go either
> way (REPRODUCES → presentational; DOES NOT REPRODUCE → substantive, hard halt).

Last audited: 2026-06-20. Built from `origin/main` of the `ResearchProj` monorepo.
Machine-readable mirror of the targets: [`gnomon/manifest.py`](gnomon/manifest.py).

## 0. SOLID vs PROVISIONAL

| | Depends on the CP3 run? | Status |
|---|---|---|
| **Embedding** (top-level subrepo, clean history, README registration) | No | **SOLID** |
| **Clean-room boundary** (`_paths.py` lattice-only; Caliper forbidden; static tests) | No | **SOLID** |
| **Targets manifest** (claimed numbers + frozen tolerances + provenance) | No — frozen before running | **SOLID (frozen)** |
| **Completeness checklist** (the Huang-flagged items METHODS.md must cover) | No | **SOLID** |
| **CP2 implementations** (forward/NLLS/Laplace/MCMC/MAF/metrics) | They *produce* the numbers | **SOLID** — 16/16 self-consistency gates pass |
| **Every reproduced number** (railing rate, coverages, ECE/sharpness gaps) | **Yes** | **RENDERED (CP3)** — in [`results/reproduction.json`](results/reproduction.json) |
| **The verdict** (REPRODUCES / DOES NOT REPRODUCE) | **Yes** | **RENDERED (CP3): PARTIAL** — see [`VERDICT.md`](VERDICT.md) |
| **Merge-back as Fashion's clean core** | Only on the REPRODUCES branch | **CONDITIONAL** |

## 1. Fashion — pinned inputs (the targets Gnomon must reproduce)

**Status (PINNED):** Fashion is *in review at MRM (R2)*, **no DOI**. The targets are
read from Fashion's prose only (`README.md`, `REVIEWER_RESPONSE*.md`) and frozen in
`gnomon/manifest.py`. Source `file:line` is recorded per target. If Fashion's
claimed numbers change in revision, the manifest pins are what must be re-frozen.

Known **under-specifications in Fashion's prose** that Gnomon resolves by documented
clean-room choice (and flags as such — these are exactly the completeness gaps):

- **clinical-sparse b-scheme:** Fashion's prose gives the dense 16-b set verbatim but
  never lists the 8-b clinical-sparse values. Gnomon pins an explicit 8-b set
  (`manifest.B_SCHEMES["clinical_sparse"]`) and says so.
- **the 1618-voxel "high-SNR ROI" selection** for the 54.7% number: Gnomon documents
  its own ROI/SNR threshold for the OSIPI abdomen scan; the railing rate's
  sensitivity to that choice is reported (a genuine reproducibility risk).
- **NLLS box bounds / init** and **MCMC sampler settings**: stated explicitly in
  METHODS.md (Fashion surfaced these only late, in reviewer responses).

## 2. Lattice — pinned substrate

**Status (PINNED, read-only):** synthetic cohorts via `lattice.make_cohort` at the
manifest seed. Lattice is MIT, synthetic-only, and on `origin/main`. Gnomon imports
it read-only and reimplements its own forward model for the estimators.

## 3. The reproduction tolerances are frozen

Set **before** running (guardrail 4): T1 ±5 pp; T3a/b ±0.10; T3c ±0.05; T4
directional with bootstrap CI excluding 0. CP3 may not loosen them post-hoc.

## 4. What would invalidate Gnomon's verdict

- A change to Fashion's claimed numbers in revision (re-freeze the manifest).
- A bug in a CP2 implementation that the self-consistency gates miss (mitigated by:
  clean-signal round-trip recovers truth; continuity limits hold; bootstrap CIs).
- Using the wrong substrate/condition (mitigated by pinning b-scheme, SNR, seed).
