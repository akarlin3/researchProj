# ASSUMPTIONS.md — the load-bearing assumption of Datum

> **Datum is a benchmark built on a ruler that is still in review.** It packages a
> fixed IVIM uncertainty-calibration task, a curated baseline panel, and reference
> numbers — **under the explicit assumption that Fashion's calibration ruler
> survives to publication as submitted.** Fashion is in review at *NMR in
> Biomedicine* (retooled, boundary-railing-first; resubmitted from MRM). The retool
> demotes the ruler to a **scoped secondary** (ground-truth/synthetic only) reported
> under the **honest CRLB** convention — Datum's substrate is synthetic, so the
> scope holds. Therefore the ruler version is pinned here, and **every reference
> number Datum produces by scoring through that ruler is flagged PROVISIONAL** in
> the table / output where it appears.
>
> If Fashion's ruler shifts in revision (its nominal levels, its coverage/ECE/
> sharpness recipe, or its headline about which interval is "the calibrated one"),
> Datum's reference numbers are **invalidated** and must be re-run. This manifest is
> what makes that re-run cheap: bump the pins in `datum/manifest.py`, run
> `python revalidate.py`, and the gates tell you what still holds. Datum shares
> Minos's applied-half rework gate.

Last audited: 2026-06-19. Audited against the working tree of the `researchProj`
monorepo (GitHub `akarlin3/ResearchProj`), branch built from `main` @ `73d588e`.
Machine-readable mirror of this file: [`datum/manifest.py`](datum/manifest.py).

---

## 0. What is solid regardless, and what is assumption-dependent

| | Depends on Fashion publishing? | Status |
|---|---|---|
| **Task definition** (`datum/task.py` — substrate, levels, metrics, conditioning, baseline set) | No — a fixed spec, runs on synthetic data | **SOLID** |
| **Substrate adapters** (`datum/substrate.py` — Lattice DRO (primary), Gauge bootstrap, OSIPI DRO) | No — read-only reuse of Lattice/Gauge | **SOLID** |
| **Baseline registry** (`datum/baselines.py`) | No — identities/paradigms/sources | **SOLID** |
| **Ruler adapter** (`datum/ruler.py`) + **submission interface** (CP3) | No — the *wiring* is fixed | **SOLID** |
| **Embedding, README, manifest, tests, CI** | No | **SOLID** |
| **Every reference number** — per-baseline coverage / coverage-gap / ECE / sharpness / pinball / interval-score, marginal **and** per-D\* tercile, with bootstrap CIs | **Yes** — produced by scoring through Fashion's ruler | **PROVISIONAL** |

The ruler and the substrate are imported **read-only** (`datum/_paths.py`); Datum
never edits Caliper or Gauge. The dependency is one-way: **nothing imports Datum.**

---

## 1. FASHION — pinned inputs (the calibration ruler Datum is built on)

**Status (PINNED):** *in review at NMR in Biomedicine* (retooled,
boundary-railing-first; resubmitted from MRM) — **NOT finalized**, no manuscript
DOI assigned. Source: `Fashion/paper_retool/` (NMR in Biomedicine, Wiley NJD-v2),
`Gnomon/handoff/CLAIMS_LEDGER.md`.

| key | pinned value | source | role for Datum |
|---|---|---|---|
| `ruler.definition_artifact` | `Fashion/uq/calib.py` | Fashion repo | the **standard** Datum scores against |
| `ruler.symbols` | `coverage`, `ece`, `sharpness_rel` | `Fashion/uq/calib.py` | the metric recipe |
| `ruler.nominal_levels` | `[0.50, 0.68, 0.80, 0.90, 0.95, 0.99]` | `LEVELS` in `calib.py` | frozen confidence levels |
| `ruler.version` | `0.1.0` | `Fashion/pyproject.toml` | package version |
| `ruler.commit` | `f078802` (last touch of `calib.py`) | `git log -1 -- Fashion/uq/calib.py` | code provenance |
| `ruler.code_zenodo` | `10.5281/zenodo.20649669` | `Fashion/README.md` | citable code snapshot |
| `ruler.scope` | scoped secondary — ground-truth/synthetic only | retooled Fashion / Gnomon hand-off | the ruler cannot touch a real scan; Datum's substrate is synthetic, so the scope holds |
| `ruler.convention` | honest CRLB (default; floored rejected) | `Gnomon/docs/METHODS.md` §5b | Datum scores the honest-CRLB intervals Caliper produces |
| `ruler.manuscript_doi` | `None` | — | **`None` ⇒ PROVISIONAL in force** |

**The Fashion assumption Datum relies on (retooled):** that Fashion's calibration
ruler — the coverage/ECE/sharpness triple at the frozen nominal levels — survives
review *as a scoped, ground-truth-only secondary* reported under the honest CRLB.
The retool drops the dramatic *marginal* D\* under-coverage (0.30/0.67) and keeps
the honest, milder, **conditional** finding (high-D\* tercile under-coverage); the
assumption-free primary is boundary railing on real data, which needs no ruler.
Datum operationalizes the scoped ruler as a benchmark standard on synthetic data.
If the ruler's recipe or nominal levels change, **all reference numbers are
invalid** and must be re-run.

**Ruler implementation Datum calls (read-only):** `caliper.metrics.score_quantiles`
(Caliper v0.1.0, MIT) — Caliper packages Fashion's recipe model-agnostically.
This is the **ruler-as-standard, differentiated from Casali**: Casali (2026) is a
model-based UQ *method* documenting D\* overconfidence; Fashion's ruler is the
*standard* that measures any method's honesty, and Datum is the benchmark that
applies it — the revision asset for Fashion's Casali differentiation.

---

## 2. SUBSTRATE — pinned inputs (the data methods are scored on)

Lattice (the intended substrate) is now **built and is the primary substrate** —
task v2 swapped it in for the v1 Gauge bootstrap (the build-order dependency,
resolved). Gauge's cohort is kept runnable as the bootstrap/cross-check; an OSIPI
DRO is the external validation. Lattice's convention matches Gauge's exactly
((D, D\*, f) physical), so the same `datum.convert` mapping applies unchanged.

| key | pinned value | source | role |
|---|---|---|---|
| `substrate.primary.entrypoint` | `lattice.make_cohort` | `Lattice/lattice/cohort.py` | **synthetic** labeled IVIM DRO (primary) |
| `substrate.primary.seed` | `20260619` | `lattice.DEFAULT_SEED` | determinism |
| `substrate.primary.commit` | `cbabffe` | `git log -1 -- Lattice/lattice/cohort.py` | provenance |
| `substrate.primary.version` | `0.1.0` | `Lattice/pyproject.toml` | version pin |
| `substrate.bootstrap` | Gauge cohort, `gauge.cohort.generate_cohort`, seed `20260613`, commit `b4ada17` | `Gauge/gauge/cohort.py` | pre-Lattice bootstrap (still runnable) |
| `substrate.external_validation` | OSIPI TF2.4 DRO, DOI `10.5281/zenodo.14605039` | `Gauge/scripts/fetch_osipi.py`, `Gauge/results/osipi_provenance.json` | synthetic DRO, download-on-demand, git-ignored |

**Guardrail:** synthetic-only. No clinical / in-vivo / MSK data is materialised in
this tree; the OSIPI DRO is synthetic and only its provenance is committed.

---

## 3. The re-validation path (one command)

```
python revalidate.py          # check pins + prove substrate→ruler resolves
python revalidate.py --full   # (CP2+) regenerate every PROVISIONAL reference number
```

When Fashion's ruler locks: set `RULER['manuscript_doi']` in `datum/manifest.py`
(and bump any changed `version`/`commit`), run `python revalidate.py --full`, and
the regenerated reference numbers cease to be provisional. Until then, **no
ruler-dependent number is ever presented as final.**

---

## 4. Citable-release path (documented, NOT executed)

A JOSS / Zenodo release of Datum is gated on Fashion's ruler locking + DOIs (the
same gate as Caliper). Usable internally now under the provisional flags; not
released as a citable artifact until the ruler is final.
