# CHANGES — NMR in Biomedicine journal format (Wiley NJD-v2) + author block

Branch `worktree-nmr-biomed-format`. **Format only — no new claim, no new number.**
`numbers.tex` is unchanged (the numbers-gate `consistency.py` re-derives it
byte-identically: 174 macros, exit 0). Ports `Fashion/paper_retool/manuscript.tex`
from a neutral `article` class to the official **Wiley NJD-v2** class that NMR in
Biomedicine directs LaTeX authors to.

## What changed
- **Document class.** `article` (ebgaramond) → `\documentclass[AMA,STIX1COL]{WileyNJD-v2}`.
  Dropped the now-class-owned packages (geometry, authblk, hyperref, ebgaramond,
  fontenc/inputenc); kept the project macros (`numbers.tex`, `\Ds`, `\unit`, `\d`).
- **Vendored template files** (the Wiley class is **not on CTAN**, so it is committed
  beside the source for a self-contained `tectonic` build): `WileyNJD-v2.cls` (v0.2),
  `NJDnatbib.sty`, `WileyNJD-AMA.bst`.
- **Author block.** Filled the previously blank `\author`/`\affil`: **Avery Karlin**,
  affiliations (1) Annealing Signet Institute, (2) Dept. of Applied Physics and
  Applied Mathematics, Columbia University, (3) Dept. of Computer Science, University
  of Colorado Boulder; `\corres` + `\email{ak5232@columbia.edu}`; `\presentaddress`
  noting student status at Columbia and CU Boulder; `\authormark{Karlin}`.
- **Front matter → class macros.** Structured abstract block → `\abstract[Summary]{}`
  (272 words, ≤300 limit); added `\keywords{}` (6, journal-required, previously
  missing); added an abbreviations `\footnotetext`; added the `\jnlcitation{}`
  "How to cite" content (auto-emitted by the class).
- **References.** `\bibliographystyle{unsrt}` → AMA numbered/superscript via the
  `AMA` class option + `WileyNJD-AMA.bst` (NMR in Biomedicine house style). All 15
  cited references resolve (no undefined citations).
- **Graphical abstract.** WileyNJD-v2 has no GA macro and the journal takes the GA as
  a separate upload, so the inline GA section was removed from the manuscript and
  re-homed into a standalone `graphical_abstract_card.tex` → `graphical_abstract_card.pdf`
  (image `figures/graphical_abstract.pdf` + the ≤80-word TOC text), preserved intact.
- **Build hygiene.** Added `.gitignore` for LaTeX intermediates. `tectonic` builds
  both `manuscript.pdf` (9 pp) and the GA card cleanly; the WileyNJD-v2 class compiles
  under tectonic's XeTeX (STIX fonts auto-fetched).

---

# CHANGES — needle/novelty (1/3): pre-empt "railing is just known $D^*$ non-identifiability"

Branch `worktree-needle+novelty` (requested: `needle/novelty`). Target paper:
`Fashion/paper_retool/manuscript.tex` (the railing-first IVIM submission for NMR in
Biomedicine). **Framing only — no new claim, no new number, no new result.** Every
number in the manuscript is a `\macro` from the auto-generated `numbers.tex`, which is
unchanged (byte-identical) by this work.

## Goal

Make the novelty **unmissable** so a reviewer cannot default to "bound-railing is just
the already-known weak identifiability of $D^*$" (the R2 → Major risk). The fix is to
state, at every place the contribution is named, that the *phenomenon* is **not** the
claim — the *contribution* is: (i) an assumption-free, optimizer-level statistic that
**neutralizes the overextension objection no calibration ruler can answer**, and (ii)
its **first systematic, pre-registered, cross-cohort** quantification with an
**independent clean-room reproduction**. The known weak-identifiability fact is conceded
outright and the contribution pivots off it.

## Edits (3, prose only)

| # | Location | What it adds | Concedes phenomenon is known? |
|---|----------|--------------|-------------------------------|
| 1 | **Abstract — Purpose** | One sentence: "That $D^*$ is weakly identified is itself long known; what is new is not the phenomenon but an assumption-free, optimizer-level statistic that neutralizes this objection where no calibration ruler can, together with its first systematic, pre-registered, cross-cohort quantification and an independent clean-room reproduction." | Yes ("itself long known") |
| 2 | **Introduction — re-aim paragraph (end of)** | Novelty up front: "We are explicit about what is new. That $D^*$ is weakly identified — and that a bounded fit can consequently pin it to a bound — is long established [koh2011, lemke2011]; the contribution is not that observation but an assumption-free, optimizer-level statistic that answers the overextension objection no calibration ruler can, given here its first systematic, pre-registered, cross-cohort quantification and an independent clean-room reimplementation." | Yes ("long established", cited) |
| 3 | **Discussion — Differentiation** | (a) **New** one-sentence delta vs known $D^*$ instability ("textbook [koh2011, lemke2011]; the contribution here is not that observation but an assumption-free … railing statistic …"); (b) Casali tightened to one sentence ("Casali et al. characterise one supervised network's aleatoric/epistemic uncertainty budget, validated on preclinical mouse brain; this work instead documents an estimator-agnostic, real-data identifiability signature in human-abdominal DWI, their own residual $D^*$ overconfidence being the learned echo of the same wall that rails the conventional fit."). | Yes ("textbook", cited) |

Citations used (`koh2011`, `lemke2011`, `casali2026`) are **reused** keys already present
in the manuscript and `refs.bib` — no new reference was added.

## Final-gate confirmations

- **No new claim / no new number.** The diff adds prose only. No numeric token appears in
  any added line except inside reused citation keys; **no result-number macro** (`\Homo*`,
  `\Full*`, `\Lihc*`, `\cov*`, `\flow*`, `\nlls*`, `\gap*`, `\margMq*`, `\gnomon*`,
  `\osipi*`, `\tcga*`, `\rail*`, `\thr*`, `\fashion*`) was introduced.
- **Novelty sentence present in BOTH abstract and intro.** ✔ (edits 1 and 2).
- **The phenomenon is never asserted as new.** ✔ — every novelty statement explicitly
  concedes weak $D^*$ identifiability is known ("long known" / "long established" /
  "textbook") before attaching "first" only to the *systematic, pre-registered,
  cross-cohort quantification + clean-room reproduction*, never to the phenomenon. This is
  consistent with the Results/Abstract statement that the railing rate "reproduc[es] a
  prior report."
- **Consistency gate (one command).** `python3 consistency.py` → `numbers-gate OK: 126
  macros`, exit 0. `numbers.tex` checksum unchanged (`d4e43c029cc525f5a5597dd478bd1a75`
  before and after; `diff` empty).
- **Build.** `tectonic -X compile manuscript.tex` → exit 0, `manuscript.pdf` written; no
  undefined-citation/reference errors (only pre-existing font fallbacks + minor hbox
  spacing).

## Status

**PROVISIONAL — NOT submission-ready.** Per the prompt, the manuscript is **not** to be
marked submission-ready until the author approves the novelty wording. This is needle-moving
1/3; do not auto-merge.
