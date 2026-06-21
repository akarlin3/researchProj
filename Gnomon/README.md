# Gnomon — a clean-room ruler for *reproduce-or-refute* of Fashion

> **A gnomon is the part of a sundial that casts the shadow — the bare reference
> that tells you the time without trusting the dial's markings.** Gnomon is the
> bare reference for Fashion's calibration ruler: an independent, from-scratch
> rebuild that re-measures Fashion's load-bearing numbers without trusting (or
> importing) Fashion's — or Caliper's — implementation.

*No standalone paper by default — Gnomon is a **verdict** that feeds the Fashion
retool.* Fashion (*"Calibration and Efficiency of Uncertainty Estimates in IVIM …"*,
in review at *MRM*) was returned at review on **methods**: internal inconsistencies,
incompleteness (under-specified dataset IDs, training/fitting detail, the CRLB
assumption), and overextended claims. Gnomon is the clean-slate hedge — a ruler that
**cannot inherit** those inconsistencies because it shares no code with Fashion, and
is documented completely from line one.

## The question Gnomon answers

> *Does a clean, from-scratch rebuild reproduce Fashion's load-bearing numbers?*

- **If it reproduces** (every pinned target within its frozen tolerance): the
  results are real and the rejection's inconsistencies were **presentational**.
  Gnomon becomes the clean reference implementation plus the complete methods
  write-up Fashion lacked.
- **If it does not**: the inconsistencies were **substantive**. Gnomon emits a
  precise divergence report (which numbers, by how much, likely cause) and **stops** —
  that changes the retool and possibly the claims.

This is a **CP3 hard halt either way.**

## The targets (pinned *before* running — see [`gnomon/manifest.py`](gnomon/manifest.py))

All read from Fashion's **prose** (`*.md`), never its source:

| Target | Claimed | Tolerance | Substrate |
|--------|---------|-----------|-----------|
| **T1** NLLS D\* boundary-railing | **54.7%** of 1618 high-SNR ROI voxels | ±5 pp | OSIPI abdomen (Zenodo 14605039, open) |
| **T3a** D\* coverage @0.95, Laplace SD | **0.30** (overconfident) | ±0.10 | synthetic 9-cell |
| **T3b** D\* coverage @0.95, MCMC SD | **0.67** (overconfident) | ±0.10 | synthetic 9-cell |
| **T3c** D\* coverage @0.95, MCMC **quantile** | **0.94** (≈ nominal — headline) | ±0.05 | synthetic 9-cell |
| **T4** MAF flow vs railed NLLS (ECE, sharpness) | *directional* (Fashion states no numbers) | sign + CI≠0 | synthetic held-out |

## Clean-room rules

1. **Independent implementation.** Forward model, NLLS railing, Laplace + MCMC
   posteriors, MAF flow, and the coverage/ECE/sharpness ruler are reimplemented
   from spec. Standard metrics are re-derived from their published definitions.
2. **Caliper's ruler is OFF-LIMITS as an import** — it is Fashion's method as code.
   Enforced at the seam by [`gnomon/_paths.py`](gnomon/_paths.py) (lattice-only;
   `caliper` is in `FORBIDDEN`) and by a static test.
3. **Data substrate:** synthetic from **Lattice** (read-only sibling import) + the
   **open** OSIPI abdomen scan (download-on-demand, gitignored). No proprietary /
   clinical data in tree or history. See [`CLEANROOM.md`](CLEANROOM.md).

## Layout

```
gnomon/
  manifest.py    frozen reproduction targets + tolerances + completeness checklist
  _paths.py      read-only sibling bootstrap — LATTICE ONLY (caliper forbidden)
  forward.py     [CP2] clean-room IVIM bi-exponential forward model + Jacobian/CRLB
  cohort.py      [CP2] Lattice-backed synthetic cohorts (read-only adapter)
  nlls.py        [CP2] box-constrained NLLS + D* boundary-railing diagnostic
  bayes.py       [CP2] Laplace + per-voxel MCMC (SD vs quantile intervals)
  flow.py        [CP2] MAF amortized posterior (NPE)            [optional: flow extra]
  metrics.py     [CP2] coverage / ECE / sharpness, re-derived independently
  bootstrap.py   [CP2] seeded bootstrap CIs for load-bearing numbers
  osipi.py       [CP2] OSIPI abdomen download-on-demand + provenance [data extra]
  reproduce.py   [CP3] run rebuild, compare to manifest, emit the verdict
docs/METHODS.md  the complete methods write-up (closes Fashion's completeness gaps)
tests/           CP1 scaffold gates (import, manifest, clean-room boundary, clean IP)
```

## Status

**CP1 — scaffold + targets manifest (this commit).** Structure mirrors the sibling
subrepos (Minos/Lattice/Datum/Vernier); the manifest is frozen; the clean-room
boundary is enforced and tested. The numerical rebuild (CP2) and the reproduction
verdict (CP3) are next. Everything here is **PROVISIONAL** until CP3 renders the
verdict — see [`ASSUMPTIONS.md`](ASSUMPTIONS.md) and [`VERIFICATION.md`](VERIFICATION.md).

## Setup

```bash
# Run inside the ResearchProj monorepo (Lattice is a sibling, imported read-only).
python -m pytest Gnomon/tests -q          # CP1 gates (pure-python)
# CP2+ extras, when implemented:
#   pip install -e 'Gnomon[flow,data]'    # torch (MAF) + requests (OSIPI fetch)
bash Gnomon/reproduce.sh                   # one-command verdict (CP3)
```

MIT-licensed software (the rebuilt technical core); see [`LICENSE`](LICENSE).
