# Revision notes — Proteus (PLOS Computational Biology)

Two highest-leverage fixes from the publication-odds assessment, executed end to end.
The manuscript source (`proteus.tex`) is **not** in the repo — only the compiled PDF
(`~/Downloads/proteus_manuscript.pdf`). All revised prose is therefore delivered as
**drop-in LaTeX blocks** keyed to section in `revision/manuscript_revisions.tex`
(token template) and `revision/manuscript_revisions.filled.tex` (numbers filled from
the committed analysis JSON by `revision/fill_tokens.py`; syntax-checked with
`tectonic`, builds clean).

**Pinned throughout:** decision line **−1.1587**; size-invariant percentile anchor;
original floor seed **1729**. New random draws use a new explicit seed **20260614**
(block *k* under seed 20260614+*k*). Equivalence margin **δ = 10 pp** (single constant,
`analysis/equivalence_tost.py::DELTA_PP`).

---

## Headline change that propagates everywhere
The properly-powered random floor's exposed-cleft rate is **50.2 % (105/209)**, not
the original **42.9 % (12/28)**. The original was a small-*n* upward fluctuation; the
triad-enrichment (~53-fold) is unchanged (enlarged pool gives a 1.82 % random triad
rate, vs the original 1.87 %). **Every occurrence of "42.9 %"/"~43 %" as the floor
must become "50.2 %"/"~50 %"** — abstract, §2.2, §2.3, Discussion, and the dashed
floor line in Figs 2, 3 (right), and 4. Figs 3-right, 4-left, 4-right are regenerated
below; **Fig 2's floor line still needs updating** (its per-band PET bars are
unchanged; only the dashed line 42.9 → 50.2 and the "crosses near ~40 %" → "~40–45 %"
text move).

---

## Change log (one row per change)

| # | Reviewer concern | What changed | New number(s) + provenance | File / section |
|---|---|---|---|---|
| A1a | R1 — "indistinguishable (Fisher p=0.30)" is a failure-to-reject on n=28, not equivalence | **Enlarged** the random triad-bearing baseline 28 → **209** triad-bearers (11,498 random Atlas proteins through the *unchanged* pipeline at −1.1587) | floor = **105/209 = 50.24 %**, Wilson 95 % **43.5–57.0 %**. seed **20260614** (blocks +k); `analysis/enlarge_floor.py` → `data/processed/floor_enlarged.json` | Results §2.2/§2.3; Methods "Enlarged baseline" |
| A1b | R1 — replace Fisher with a formal equivalence/non-superiority test | **TOST + non-superiority** on d = PET − floor, δ = 10 pp | d = **−17.6 pp** (90 % CI −24.8, −10.3); one-sided 95 % **upper bound on enrichment = −10.3 pp**; non-superiority **p = 2.0×10⁻¹⁰**; two-prop p = 7.2×10⁻⁵, Fisher p = 1.0×10⁻⁴. **Strict two-sided equivalence NOT claimed** — PET is significantly *below* the floor (honest: not merely "not enriched"). `analysis/equivalence_tost.py` → `equivalence.json` | Results §2.3 (BLOCK 2); Methods (BLOCK 7); Fig 4 caption |
| A1c | R1 — regenerate Fig 4 (left) with enlarged-baseline CI | Regenerated Fig 4 left (PETase 32.7 %, AChE 8.8 %, floor 50.2 % ± Wilson CI; non-superiority annotation) | `figures/fig4_left_enlarged.{pdf,png}`; `analysis/make_figures.py` | Fig 4 (left) |
| A1d | consistency — bits-gradient "tail" claim was vs the n=28 floor | Recomputed the bits-stratified tail gradient vs the **enlarged** floor | top-25 = **80 %** vs 50.2 % (Fisher **p=0.0053**, 1.59×, still above); top-50 64 % (p=0.085); top-100 44 %; rank 101–294 = 26.8 % (**p=1.5×10⁻⁶ below floor**). Crosses floor between top-50/top-100 (~1180–1090 bits). `analysis/bits_gradient_floor.py` → `bits_gradient_enlarged.json`; `figures/fig4_right_enlarged.*` | Results §2.3/§2.4; Fig 4 (right) |
| A2a | R1 — divergent-tail depletion may be an ESMFold prediction-error artifact (paper flags but never runs the pLDDT comparison) | **pLDDT confound test** across divergent (<25 %, n=2,100), enlarged baseline (n=209), near-homolog above-line (n=96); global + active-site-local | medians **92.9 / 90.1 / 94.5** (%); divergent **NOT lower** than baseline — slightly higher (MWU **p=1.5×10⁻⁹**, KS 1.4×10⁻¹⁰, Cliff's **δ=+0.25**). `analysis/plddt_confound.py` (pLDDT from cached ESMFold V0 PDBs) → `plddt_confound.json` | Results §2.5 (BLOCK 3); Limitations (BLOCK 6) |
| A2b | R1 — robustness: does depletion persist after matching prediction quality? | **pLDDT-matched re-screen** (subsample divergent triad-bearers to baseline pLDDT histogram, seed 1729) | above-line rate **11.4 % → 15.7 %** (n=89) — far below the ~50 % floor; depletion **persists**. `analysis/plddt_confound.py` | Results §2.5; Limitations |
| A2c | R1 — supplementary pLDDT figure | New supp figure: ECDF + violin of pLDDT by band | `figures/figS_plddt_bands.{pdf,png}`; `analysis/make_figures.py` | new Supp Fig (BLOCK 10) |
| B1 | R2+R3 — reframe as decision-relevant cautionary benchmark, not a methods footnote | Rewrote abstract, intro framing, discussion: "excellent serine-hydrolase **finder**, not a PET-function **discriminator**; homology bound quantified" | prose; BLOCKS 1, 4 | Abstract; Intro; Discussion |
| B2 | R2 — engage adjacent prior art | New Discussion paragraph: convergent 2025 MBE structure-vs-sequence phylogenomics result (same low-pLDDT root cause); reconcile with Foldseek mining + ML PET-hydrolase **successes** as homology-reachable / sequence-filtered | Mutti 2025 (10.1093/molbev/msaf149); Jaito 2026 (10.7717/peerj.20462, 711 lipases); Norton-Baker 2025 (10.1021/acscatal.5c03460). **All Crossref-verified.** | Discussion (BLOCK 5); References (BLOCK 9) |
| B3 | R3 — soften the strongest generalization | Title/abstract/discussion qualified: "**for this operationalization** of PET-hydrolyzing function" / "**with these descriptors**"; headline retained | prose | Title/Abstract/Discussion (BLOCKS 1, 4) |
| B4 | PLOS — missing declarations | Added **Competing interests** + **Funding** statements (kept existing Gen-AI disclosure + Data availability) | BLOCK 8 | Declarations |
| B5 | author instruction — affiliation/byline | Byline set to **Avery Karlin — The Annealing Signet Institute** (Columbia + CU-Boulder bio, ORCID 0000-0003-3848-6782), matching the projGauge paper | BLOCK 0 | Title page |

---

## Scientific bottom line (what moved)
- The central claim is **stronger and honest** after powering it: structure-first
  triage's PET branch is **not enriched above** the random serine-hydrolase floor —
  in fact it sits **significantly below** it (the "indistinguishable, p=0.30" verdict
  was an artifact of the n=28 baseline). Non-superiority is decisive; the largest
  plausible enrichment is negative (−10.3 pp).
- The divergent-tail depletion is **not** an ESMFold artifact: the divergent band is
  **not lower** in pLDDT than the baseline (if anything higher), and the depletion
  survives a pLDDT-matched re-screen. This replaces the manuscript's "merits an
  explicit pLDDT comparison" hand-wave with the actual result.
- The only above-floor signal remains the near-homolog tail (top-25 = 80 %,
  p=0.0053), exactly where sequence search already reaches — the homology-bound thesis
  holds and is quantified.

## Provenance / reproduce
```
PYTHONPATH=src python analysis/enlarge_floor.py        # CP1 data  (network + fpocket; ~45 min)
PYTHONPATH=src python analysis/equivalence_tost.py     # CP1 stats
PYTHONPATH=src python analysis/plddt_confound.py       # CP2       (reads cached ESMFold PDBs)
PYTHONPATH=src python analysis/bits_gradient_floor.py  # tail gradient vs enlarged floor
PYTHONPATH=src python analysis/make_figures.py         # all figures
python revision/fill_tokens.py                         # fill LaTeX numbers
```
Env: conda env `proteus` (numpy/biotite/scipy/matplotlib). `fpocket` on PATH.
Heavy upstream inputs not committed (gitignored): the ESM-Atlas PDB caches
(`sensitivity_cache/`, ~2,100 models) and `branch_partition.csv`; small input
summaries are snapshotted under `data/processed/inputs_snapshot/`.

## Still open / `[NEEDS DATA]`
- **None of the new numbers are `[NEEDS DATA]`** — every statistic traces to a
  committed script run on available data.
- **Author action (not data):** update Fig 2's dashed floor line 42.9 → 50.2 % (bars
  unchanged); paste the drop-in blocks into `proteus.tex`; complete the Gen-AI
  disclosure specifics; replace remaining `[VERIFY]` reference placeholders (the three
  prior-art citations added here are Crossref-verified and ready).
