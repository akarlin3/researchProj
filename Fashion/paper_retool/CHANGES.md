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
claim — the *contribution* is: (i) an assumption-free, optimiser-level statistic that
**neutralises the overextension objection no calibration ruler can answer**, and (ii)
its **first systematic, pre-registered, cross-cohort** quantification with an
**independent clean-room reproduction**. The known weak-identifiability fact is conceded
outright and the contribution pivots off it.

## Edits (3, prose only)

| # | Location | What it adds | Concedes phenomenon is known? |
|---|----------|--------------|-------------------------------|
| 1 | **Abstract — Purpose** | One sentence: "That $D^*$ is weakly identified is itself long known; what is new is not the phenomenon but an assumption-free, optimiser-level statistic that neutralises this objection where no calibration ruler can, together with its first systematic, pre-registered, cross-cohort quantification and an independent clean-room reproduction." | Yes ("itself long known") |
| 2 | **Introduction — re-aim paragraph (end of)** | Novelty up front: "We are explicit about what is new. That $D^*$ is weakly identified — and that a bounded fit can consequently pin it to a bound — is long established [koh2011, lemke2011]; the contribution is not that observation but an assumption-free, optimiser-level statistic that answers the overextension objection no calibration ruler can, given here its first systematic, pre-registered, cross-cohort quantification and an independent clean-room reimplementation." | Yes ("long established", cited) |
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
marked submission-ready until the author approves the novelty wording. This is needle-mover
1/3; do not auto-merge.
