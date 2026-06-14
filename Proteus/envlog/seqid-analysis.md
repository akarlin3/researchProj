# Sequence-identity reach analysis — Figure 1 of the methods paper

**Run date:** 2026-06-14
**Where:** local Mac (M4), off the existing GCS sweep artifacts — **no GCE, no re-search, no re-fold**.
**Reproduce:** `conda activate proteus && python scripts/seqid_analysis.py <out_dir> <tmp_dir>`
(deterministic; biotite Smith-Waterman local alignment, BLOSUM62, gap −10/−1; seed-free). Inputs are
the per-query-tiering screened records + cached Atlas PDBs and the control query PDBs (paths set at the
top of the script).
**Central question:** is the above-floor structural signal **near-homolog in sequence**, or do any
sequence-divergent-but-structurally-close hits clear the cleft line? This re-plots the
per-query-tiering PETase-branch gradient against **% sequence identity** instead of Foldseek bits,
and triages the 96 above-line hits into *near-homolog* vs *structure-first*.

## TL;DR

**The above-floor signal is entirely near-homolog. Discovery is cleanly flat — structure-first adds
nothing beyond homology's reach.** Every one of the 96 above-line hits is ≥25% identical to a known
PETase query, with full catalytic-domain coverage (93/96 are ≥30%). The above-line│triad rate rises
**monotonically** with sequence identity (6% → 27% → 70% → 92%), crossing the 42.9% random floor only
at **~40% identity** — squarely in homology territory. Below the **~30% twilight threshold** the
structural cleft signal *collapses to 6.1%* (well under the random floor), so structure does **not**
persist where sequence search fails — they fail together. The genuinely divergent dark tail
(<20% identity) the project targets is **empty**: no screened triad-bearer is that divergent.
The paper's central claim now stands rigorously.

---

## Checkpoint 0 — Precondition report (all inputs recovered; no STOP)

| Precondition | Status |
|---|---|
| Screened **PETASE-branch triad-bearers (294)** + **96 above-line** + anchor map | ✅ `per_query_screen_PETASE.json` (300 screened → 294 triad → 96 above-line), each carrying `best_query` (structural anchor) and `petase_query`. For the PETASE branch the anchor **is** the nearest PETase query by bits (`best_query == petase_query` for all 300). |
| **Hit sequences** (Atlas, keyed by MGYP) | ✅ all 300 fetched ESMFold PDBs cached in `per_query_cache/` (the foldcomp models from the live `fetchPredictedStructure` path). Sequence read from CA ATOM records — no `highquality_clust30.fasta` slice or re-fetch needed. |
| **Query sequences** (PETase/cutinase fold-class queries) | ✅ control PDBs in `structures/`: IsPETase 6EQE (265 aa), LCC_WT 4EB0 (258), PET46 8B4U (269), Cut190 4WFI (258), TfCut2 4CG1 (262). Chain-A protein sequence from CA records. |
| **Alignment tool** | ✅ biotite 1.6.0 (`align_optimal`, BLOSUM62) — the pipeline's own structure library. MMseqs2 is present in conda pkgs but unlinked; biotite Smith-Waterman is sufficient and **independently validated** below against Foldseek's own `fident` (Pearson r = 0.987). parasail / biopython absent. |
| Floor comparator | ✅ `floor.json`: 1500 random screened → 28 triad+ → 12 above-line; conditional **12/28 = 42.86%**. Reused as measured. |

Anchor-class distribution of the 300 (best query by Foldseek bits): **8B4U/PET46 226**, 4CG1/TfCut2 33,
6EQE/IsPETase 27, 4EB0/LCC 13, 4WFI/Cut190 1 — i.e. most PETASE-branch neighbours structurally match
the *divergent archaeal* PET46 best, which is exactly why this branch reaches down toward the twilight
zone at all.

---

## Checkpoint 1 — Sequence identity with coverage (method)

For each screened triad-bearer the hit sequence was locally aligned (Smith-Waterman, BLOSUM62,
gap −10/−1, no terminal penalty) against **each** of the five PETase queries. Per alignment we record:

- **% identity** = identical columns ÷ alignment length (MMseqs/Foldseek `fident` convention);
- **query coverage** = aligned query residues ÷ query length;
- **hit coverage** = aligned hit residues ÷ hit length.

Two identities are reported per hit: to the **structural anchor** (`best_query`) and to the **nearest
PETase query overall** (the max-identity query among those clearing a **50% query-coverage floor**).
The coverage floor is the guard the brief requires — it prevents a short high-identity patch from
masquerading as a homolog. **All 300 hits clear it** (0 low-coverage flags): every alignment spans the
full α/β-hydrolase catalytic domain (median query coverage ≈ 0.93), so these are genuine domain-level
homologs, not local-patch artifacts. Per-hit table: **`data/processed/seqid_per_hit.csv`** (300 rows,
identity + coverage to all five queries, plus Foldseek `fident`).

**Independent cross-check (Foldseek `fident`):** my biotite anchor identities track Foldseek's own
structure-derived identity (`result.m8` col 5) at **Pearson r = 0.987** (mean 36.1% vs 33.8%, median
absolute difference 2.0 pts). The recomputation is sound; the small offset is the expected
sequence-SW-vs-structural-alignment gap.

---

## Checkpoint 2 — The gradient re-plotted against sequence identity (Figure 1 data)

**Above-line │ triad rate, stratified by nearest-PETase % identity** (294 triad-bearers; Wilson 95%
CI; compared to the 42.86% random floor by Fisher exact + rate ratio):

| seq-id bin | n (triad) | above-line | rate | Wilson 95% | vs floor 42.9% |
|---|---:|---:|---:|---|---|
| **<20%** | **0** | 0 | — | — | *empty — no divergent dark tail* |
| 20–30% | 49 | 3 | **6.1%** | [2.1%, 16.5%] | p < 1e-4, RR **0.14×** (far below) |
| 30–40% | 190 | 52 | **27.4%** | [21.5%, 34.1%] | p = 0.119, RR 0.64× (below) |
| 40–60% | 43 | 30 | **69.8%** | [54.9%, 81.4%] | p = 0.029, RR **1.63×** (above) |
| >60% | 12 | 11 | **91.7%** | [64.6%, 98.5%] | p = 0.005, RR **2.14×** (above) |

```
above-line|triad rate vs sequence identity      (floor = 42.9% ┊)
 <20%  (n=0)   ·                                 ┊            empty
20-30% (n=49)  ███ 6.1%                           ┊
30-40% (n=190) ██████████████ 27.4%               ┊
40-60% (n=43)  ███████████████████████████████████┊████ 69.8%
 >60%  (n=12)  ███████████████████████████████████┊████████████ 91.7%
```

- **Monotonic** in sequence identity (Spearman of binned rate is +1).
- **Floor crossover ≈ 40% identity**: the rate sits *below* the random floor through 40% identity
  (27.4% at 30–40%) and only clears it in the 40–60% band (69.8%) — the same place a sequence search
  is already comfortably in the homolog regime.
- **Twilight threshold ≈ 25–30%**: in the 25–30% band (the bin is effectively 25–30%; minimum identity
  across all 294 is **25.1%**) the rate is **6.1%**, *one-seventh* of the random floor. Below the
  twilight zone the structural cleft signal does **not** survive — it disappears alongside sequence
  reliability. There is no "structure rescues where sequence fails" regime in the data.

**Bits ↔ sequence-identity correlation (the confirmation test).** Across the 294 triad-bearers,
Foldseek bits and nearest-PETase identity are strongly positively correlated: **Pearson r = 0.735**
(p = 3.3e-51), **Spearman ρ = 0.732** (p = 1.4e-50). The bits-decay crossover the prior run located at
**~1090 bits** (per-query top-100 cut, where the bits-stratified rate ≈ the floor) maps to a median
nearest-PETase identity of **35.2%** — *above* the 30% twilight line. So the bits gradient and the
identity gradient are the **same gradient**: high bits ≈ high identity, and the above-floor tail lives
at near-homolog identities. **The near-homolog read is confirmed.**

---

## Checkpoint 3 — The decisive triage of the 96 above-line hits

| cell | definition | **count** (primary, biotite SW) | count (robustness, Foldseek `fident`) |
|---|---|---:|---:|
| **Near-homolog** | ≥30% identity to a PETase query (full coverage) | **93** | 77 |
| Twilight | 25–30% identity | **3** | 18 |
| **Structure-first** | <25% identity, structurally close, above-line, high composite | **0** | 1 |

- **Near-homolog (93/96):** identity range 30.5%–84.8% to a PETase query — a plain sequence search
  already reaches every one. Not novel. (Even the top-composite candidate, MGYP000644355661, comp 2.76,
  is 32.6% identical to IsPETase.)
- **Twilight (3/96):** the only hits between 25% and 30% identity — `MGYP001815130755` (25.1%),
  `MGYP001401963747` (28.1%), `MGYP002625649965` (29.8%). **All three anchor to PET46 (8B4U)**, the
  divergent archaeal PETase, at ≥93% query coverage, and all three are the *lowest-composite* above-line
  hits (−0.76, −0.94, −0.83 — barely clearing the −1.1587 line). They are full-domain homologs of the
  most divergent *known* PETase, not novel folds.
- **Structure-first (the one cell that could hold a genuine divergent discovery): EMPTY (0).**
  Under the primary metric no above-line hit is <25% identity to a PETase query. Under Foldseek's
  stricter structural `fident`, exactly **one** hit dips just below the line and must be named:

  | MGYP | anchor | nearest id (biotite) | nearest id (Foldseek) | cov_q | composite | pLDDT | bits |
  |---|---|---:|---:|---:|---:|---:|---:|
  | **MGYP001815130755** | 8B4U / PET46 | 25.1% | **23.4%** | 0.96 | −0.757 | 0.979 | 963 |

  This single borderline case is a **96%-coverage, full-catalytic-domain homolog of PET46** (the
  divergent archaeal positive) sitting at the 23–25% twilight boundary with a composite *barely* above
  the line. It is not a sequence-divergent structure-first discovery — it is a marginal PET46 homolog
  that a sequence search seeded with the archaeal PETase would itself surface. **No hit lies in the
  genuinely divergent dark tail (<20%); that bin is empty.**

---

## Checkpoint 4 — Verdict

**The above-floor structural signal is entirely near-homolog in sequence; the discovery is cleanly
flat and structure-first adds nothing beyond homology's reach.** Re-plotted against sequence identity,
the PETASE-branch gradient is a *homology gradient*: the above-line│triad rate climbs monotonically
with identity and only beats the 42.9% random floor at ≥40% identity, deep in the homolog regime. The
bits-decay the prior run measured is the same curve (bits ↔ identity r = 0.74; the ~1090-bit crossover
sits at ~35% identity). Below the ~25–30% twilight threshold the cleft signal collapses to 6% — *below*
random — so structure does not reach where sequence cannot. The triage is unambiguous: **0** structure-
first finds under a sequence-SW metric, and under the strictest structural-identity metric a **single**
borderline (23.4%) full-coverage homolog of the divergent archaeal positive PET46 — not a novel
divergent fold. The divergent dark tail (<20% identity) the project set out to mine is **empty in the
screened tier**: there is nothing there for a structure-first method to find.

This is the rigorous form of the paper's central claim. The honest-negative methods result holds: a
structure-first screen of metagenomic dark proteins recovers PET-hydrolase-like clefts **only among
sequences a homology search already reaches**, and confers no detectable advantage in the divergent
tail. The thin "crack" is a single twilight-boundary homolog of a known archaeal PETase — too marginal
(barely above the line, ≥23% identical, full domain coverage) to anchor a discovery angle in the paper.

### Scope guard (carried)

**Identity here measures sequence-search *reach*, not PET activity.** A hit being near-homolog says a
sequence search would find it; a hit being structurally close says nothing about PET turnover. The
S4/S5 path tests fold-class + exposed-cleft geometry, not catalysis. No candidate is verified — leads
are prioritized, not validated. This analysis settles *reach* (the paper's claim), not function.

## Caveats

- **Identity metric.** Primary = sequence-only Smith-Waterman (BLOSUM62); robustness = Foldseek
  structural `fident`. They agree at r = 0.987; the structure-first count is 0 (SW) / 1 (fident) — the
  verdict is invariant to the choice (at most one marginal twilight-boundary homolog).
- **This tier is the top-300 PETASE branch** (bits 962–1717), the closest structural neighbours of the
  PETase queries — by construction enriched for homologs. That is the point: it is the *most favourable*
  tier for finding a divergent above-line hit, and even here the <20% dark tail is empty. The full Atlas
  certainly contains <20%-identity α/β-hydrolases, but none are structurally close enough to a PETase
  query to enter this screened tier — which is itself the flat-discovery result.
- **Floor n = 28 (<50)** → wide floor CI [26.5%, 60.9%]; the verdict is robust regardless, since the
  divergent bins sit far below the floor point estimate and the empty <20% bin needs no comparator.
- **fpocket non-determinism** is inherited (composites/line pinned at −1.1587 upstream); identity and
  coverage computed here are fully deterministic.

## Reproducibility

- **Inputs:** `per_query_screen_PETASE.json`, `per_query_cache/*.pdb` (300 Atlas foldcomp models),
  control query PDBs `structures/{6EQE,4EB0,8B4U,4WFI,4CG1}.pdb`, `floor.json`, `result.m8.gz`
  (Foldseek cross-check) — all from `gs://projproteus-fold/atlas-sweep/2026-06-14/`.
- **Aligner:** biotite 1.6.0 `align_optimal` local, BLOSUM62, gap (−10,−1), no terminal penalty.
- **Identity:** identical ÷ alignment length; **coverage:** aligned query residues ÷ query length;
  coverage floor 0.50 for the nearest-PETase call.
- **Artifacts:** this report + `data/processed/seqid_per_hit.csv` (per-hit Figure 1 table), pushed to
  `gs://projproteus-fold/atlas-sweep/2026-06-14/`.
