# Random-floor measurement + lift decomposition

**Run date:** 2026-06-14
**Where:** local Mac (M4), off the existing GCS artifacts — **no GCE, no re-search**.
**Reproduce:** `PYTHONPATH=src python -m proteus.floor` (config `floor` block; seed 1729).
**Outputs:** `data/processed/floor_sample.txt`, `data/processed/floor.{json,csv}`.

## TL;DR

We replaced the statistically-thin pilot floor (1 event in 60 = 1.7%, whose own CI
ran ~0–9%) with a properly-sampled uniform draw of **1500 random HQ-clust30 Atlas
proteins**, screened through the **unchanged** S4→S5 pipeline at the enriched run's
widened line (**−1.1587**), and decomposed the enrichment lift.

- **Tightened floor:** above-line rate **0.80%** (12/1500), Wilson 95% CI
  **[0.46%, 1.39%]**. This replaces 1.7%.
- **The enriched lift is real and highly significant:** enriched 7.67% (23/300) vs
  floor 0.80% — Fisher p = 7.3e-11, rate ratio **9.6×** (95% CI 4.8–19.0).
- **DECOMPOSITION VERDICT — the lift is entirely S4 (the fold-class/triad search),
  not S5 (the cleft line):**
  - **S4 / triad rate:** 99.3% (enriched) vs **1.87%** (random) — rate ratio
    **53×** (95% CI 37–77), p ≈ 0. The fold-class search is the whole engine.
  - **S5 / conditional rate (above-line GIVEN a triad — the load-bearing number):**
    7.7% (enriched) vs **42.9%** (random) — random triad-bearers clear the cleft
    line **~5× MORE often** than the fold-class hits do (p = 9e-9; rate ratio
    enriched/floor = 0.18). **S5 adds no positive enrichment beyond "has a triad";
    on real data it selects *against* the fold-class hits.**

## Checkpoint 0 — PRECONDITION REPORT (audit)

| Precondition | Status |
|---|---|
| `screen` / S4 / S5 / `calibrate` intact | ✅ reused untouched; imports + run clean in the `proteus` conda env (numpy/scipy/biotite; fpocket on PATH). |
| Widened line value + source | **−1.1587**, from the enriched sweep `gs://projproteus-fold/atlas-sweep/2026-06-14/funnel.json` (`threshold`). Re-derivable from the controls via `calibrate.recovery_screen` (min over IsPETase/LCC + recovered divergent positives PET46/Cut190/TfCut2 at precision 1.0, divergent-recall 1.0). **Local re-derivation jitters** (−1.149 / −1.156 / −1.161 / −1.168 across runs — fpocket non-determinism, ±~0.01), so we **pin the decision line to −1.1587** for the floor so both sides are judged at the identical line. |
| Atlas HQ clust30 accession universe (same population as the enriched set) | ✅ `highquality_clust30.lookup` — the Foldseek DB's own representative list, **36,986,627 entries**, served Range-enabled at `https://foldcomp.steineggerlab.workers.dev/highquality_clust30.lookup` → `steineggerlab.s3.amazonaws.com/...` (**1,030,274,634 bytes**, `Accept-Ranges: bytes`). This is the SAME population the enriched fold-class search ran against; the floor is drawn uniformly from it with **no fold-class pre-filter**. |
| Foldcomp / fetch path (API was 403 during the GCE run) | ✅ Two paths confirmed. (1) The **foldcomp worker** serves the DB files (`.lookup`/`.index`/data) with HTTP Range. (2) `fetchPredictedStructure/{MGYP}.pdb` is **currently live** (verified 1500/1500 + 39/39 probes = 200) and returns the **same ESMFold V0 models** the GCE foldcomp DB holds — proven structurally identical: re-screening enriched candidates locally reproduces their **catalytic-triad residues exactly** (e.g. 497/489/411, 221/435/334). We fetch by accession from the live API (with retries); the foldcomp worker is the documented fallback. |
| Enriched funnel (comparison set) | ✅ 300 screened → 298 triad+ → 23 above the −1.1587 line (`funnel.json`); conditional 23/298 = 7.72%. |

No precondition missing — proceeded.

## Checkpoint 1 — Fair random sample

- **Uniform draw of `floor.n` = 1500** accessions from the HQ clust30 universe by
  seeded random **byte offsets** across the whole 1.03 GB `.lookup` (window 256 B,
  first complete line per window). Lines are ~27.85 B (near-constant width), so
  uniform-over-bytes ≈ uniform-over-accessions.
- **Seed = `random_seed` = 1729.** Deterministic regardless of fetch parallelism
  (a seeded RNG fixes the offset stream; windows are fetched in parallel then
  walked in stream order). Confirmed reproducible: two independent runs produced
  the identical leading accessions.
- **Guard — same population, no fold-class pre-filter:** these are random Atlas
  proteins, the null. 1500 unique drawn from 1875 offsets.
- Recorded to `data/processed/floor_sample.txt` (accessions + seed + per-hit pLDDT).

## Checkpoint 2 — Fetch + screen through the SAME pipeline

- Fetched all 1500 (100% fetch success; mean pLDDT 0.878, 0–1 scale). Screened
  through the unchanged pipeline at **−1.1587**: S4 geometry → fpocket (only on
  triad-positives) → S5 cleft → composite, scored against the SAME IsPETase/LCC
  anchor as calibration. Parallelized over cores; fpocket runs in an isolated temp
  dir per call (race-free).
- **fpocket is non-deterministic run-to-run** (the widened line wandered −1.149 to
  −1.168 across local runs; individual composites can swing ~±1.0). S4 geometry is
  deterministic (triads reproduce exactly). To make the rate comparison immune to
  this, we **re-screened the enriched 300 in the SAME local pass** (same anchor,
  same line, same fpocket session) as a matched control — see CP3.

### Per-stage funnel (random floor)

| Stage | Count | Rate |
|---|---:|---|
| Fetched & screened | **1500** | 100% fetched (mean pLDDT 0.878) |
| → triad-positive (S4) | **28** | **1.87%** of screened |
| → catalytic pocket (S5, ≤12 Å of Ser OG) | 26 | 92.9% of triad |
| → **above the widened line** (composite ≥ −1.1587) | **12** | **0.80%** of screened; **42.9%** of triad |

Enriched comparator (same line): **300 → 298 triad+ → 23 above-line**
(canonical GCE); local re-screen **300 → 298 triad+ → 27 above-line** (reproduces
canonical within fpocket jitter — 298 triad exact, 27 vs 23 above-line).

## Checkpoint 3 — Statistics + lift decomposition

### 1. Tightened floor (replaces 1.7%)

Above-line rate **12/1500 = 0.80%**, **Wilson 95% CI [0.46%, 1.39%]**. The old
1/60 = 1.7% sat inside this CI but its own CI ran ~0–9% — uninformative. The floor
is now pinned to <1.4% with 95% confidence.

### 2. Corrected significance of the enriched lift

The old claim ("the enriched CI excludes 1.7%") was invalid (1/60's CI ran ~0–9%).
Tested correctly against the properly-sampled floor:

| Comparison | Enriched | Floor | Fisher p | rate ratio (enr/floor), 95% CI |
|---|---|---|---|---|
| above-line / screened | 23/300 = **7.67%** [5.16, 11.24] | 12/1500 = **0.80%** [0.46, 1.39] | **7.3e-11** | **9.58×** [4.82, 19.05] |

Two-proportion z = 7.86, p = 3.8e-15. **The enriched rate significantly exceeds the
properly-sampled floor.**

### 3. DECOMPOSITION — where does the lift come from?

`overall = triad_rate × conditional_rate` (enriched: 0.993 × 0.077 = 0.077 = 23/300).

| Component | Enriched | Floor (random) | Fisher p | rate ratio enr/floor |
|---|---|---|---|---|
| **S4 / triad rate** (triad+ / screened) | 298/300 = **99.3%** [97.6, 99.8] | 28/1500 = **1.87%** [1.29, 2.68] | 6e-305 | **53.2×** [36.9, 76.8] |
| **S5 / conditional rate** (above-line / triad+) | 23/298 = **7.72%** [5.20, 11.31] | 12/28 = **42.86%** [26.51, 60.93] | 3.3e-6 | **0.18×** [0.10, 0.32] |

**Matched local re-screen (jitter-immune, same fpocket session):** enriched
conditional 27/298 = **9.06%** [6.30, 12.86] vs floor 42.86% — p = 1.4e-7, rate
ratio enriched/floor = 0.21 [0.12, 0.37]. The conditional inversion is **not a
cross-machine fpocket artifact**: the floor conditional CI lower bound (26.5%) lies
far above the enriched upper bound (~12.9%) under matched conditions.

> **Triad count caveat (flagged per spec):** the floor yielded **28 triad+ (< 50)**,
> so the conditional estimate's denominator is modest. Raising `floor.n` would
> tighten the floor-conditional CI, **but cannot change the verdict** — the floor
> conditional (43%, CI ≥ 26.5%) already sits entirely above the enriched conditional
> (≤ 12.9%). The direction (S5 does *not* enrich) is secure at n = 1500.

## Checkpoint 4 — Verdict + honest scope

### The decomposition verdict (stated plainly)

**The enrichment-vs-random lift is entirely S4 (the fold-class / triad search). S5
(the cleft line) does not earn its place as a discriminator on real Atlas data.**

- **S4 carries everything.** The fold-class Foldseek search raises the catalytic-triad
  rate from the random **1.87%** to **99.3%** — a **53× enrichment**. That is the
  whole engine of the lift: "the fold-class search finds α/β-hydrolase triads."
- **S5 provides no positive enrichment beyond "has a triad."** Among triad-bearers,
  random Atlas proteins clear the widened cleft line at **42.9%**, the fold-class
  hits at only **7.7–9.1%**. The cleft line passes random triad-bearers **~5× more
  often** (p ≈ 1e-7). On real data the S5 peripherality line is *anti-correlated*
  with the fold-class hits, not discriminating in their favour.

**Why (mechanism, confirmed by the data):** the S5 line is exposure-dominated (0.40
weight on Ser peripherality; druggability/depth down-weighted) — it favours
shallow, surface-exposed, lid-less active sites and penalises deep/buried ones. The
enriched tier is **289/300 (96%) sourced from the 1EA5 anchor — *Torpedo*
acetylcholinesterase, a deep-gorge hydrolase** (11/300 from 1CRL lipase; **zero**
from the IsPETase/LCC/PET46/Cut190/TfCut2 PETase-cutinase anchors). Deep-gorge /
lid-bearing hydrolases bury the catalytic Ser → low peripherality → they **fail** the
exposure-favouring line (only 7.7% pass). Random triad-bearers carry no gorge bias,
so ~43% sit peripheral enough to pass. The cleft line is doing its job (rejecting
buried sites) — but on a fold-class shortlist that is itself dominated by buried-site
AChE-like enzymes, that means S5 mostly *removes* the enriched hits.

### Scope guard (carry it)

This measurement validates **enrichment-vs-random only**. It does **NOT**:

- **Address AChE-branch dominance.** Quantified here and worse than feared: the
  enriched tier is 96% acetylcholinesterase-branch, 0% PETase/cutinase-branch. The
  surfaced leads — and the floor's 12 above-line hits — are **exposed-site serine
  hydrolases of unknown function, not PET candidates**. This is the per-query-tiering
  problem (the top-by-bits tier is swamped by the one anchor with the highest-scoring
  Atlas matches), separate from this measurement.
- **Address PET-specificity.** Neither S4 (triad geometry) nor S5 (a generic
  exposed-cleft score) tests for PET turnover. Precision against the true
  metagenomic negative space — vast non-PET serine-hydrolase fold space the 6-control
  panel never saw — remains **unmeasured**. A specificity filter is a separate piece
  of work.

### Caveats

- **fpocket non-determinism** handled by pinning the line to −1.1587 and re-screening
  the enriched set in the same pass; the headline rates are robust to it.
- **Floor triad+ = 28 (< 50)** — flagged; verdict robust regardless (see CP3).
- The floor's 12 above-line hits are leads at best (exposed-site triad-bearers); the
  single highest-scoring (composite 6.10) out-scores the enriched rank-1 (4.29),
  underscoring that the cleft line is not PET-specific.

## Reproducibility

- **Seed** 1729; **floor.n** 1500; **line** −1.1587 (pinned to the enriched sweep).
- **Universe** `highquality_clust30.lookup` (36,986,627 entries; 1,030,274,634 B).
- **Fetch** `api.esmatlas.com/fetchPredictedStructure/{MGYP}.pdb` (live; structurally
  identical to the GCE foldcomp models — triads reproduce exactly). Foldcomp worker
  Range fallback documented.
- **Anchor** IsPETase+LCC, percentile mode, controls separated (margin ~1.40).
- **Artifacts** `floor_sample.txt`, `floor.{json,csv}`; enriched comparator from
  `gs://projproteus-fold/atlas-sweep/2026-06-14/`.
