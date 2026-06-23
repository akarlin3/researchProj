# PROVISIONAL_LEDGER.md — what is held, what moves, and what is load-bearing

> **Bottom line (the GATE E claim):** **No held or provisional number imported from an unpublished
> anchor is load-bearing for Augur's spine.** The spine of §5 (*$D^*$ cross-modally orphaned*) rests
> on **in-repository, reproduced anchors** (the CRLB identifiability wall and the $D^*$ test–retest
> null with its bootstrap CI). Fashion and Minos are **contextual forward-citations**; their final
> numbers can shift in review without changing the spine. This ledger records every imported value,
> whether it moves on acceptance, and the re-check that fires.

Last audited: 2026-06-22 (this finish build). Re-run `bash reproduce.sh` after any anchor update.

---

## 1. Numbers imported from each anchor

| # | Value in Augur | Macro | Source anchor | In manuscript? | Load-bearing for spine? | Moves on anchor acceptance? | Re-check that fires |
|---|---|---|---|---|---|---|---|
| 1 | CRLB(D\*) growth ~**6×**; reproduced **5.7×** | `\numCRLBgrowth`, `\numCRLBgrowthRepro` | **Gauge** (in-repo, reproduced) | §4, §5(1), Fig.1 | **YES** | No — reproduced from first principles this build | Confirm reproduced factor still ≈ anchor on Gauge update |
| 2 | high-D\* CRLB/tercile-width **1.12×** | `\numCRLBtercileRatio` | **Gauge** (in-repo) | §4, §5(1) | YES (cohort-realized; carried) | Possibly (seed/SNR-draw dependent) | Re-extract from `conditional_attack_report.txt`; band check in `crlb_wall.json` |
| 3 | conformal width~CRLB **r=0.77** | `\numWidthCRLBr` | **Gauge** (in-repo) | §4 | Supporting | Unlikely | Re-extract on Gauge update |
| 4 | D\* test–retest **r=−0.17**, p=0.13, n=76 | `\numDstarR`,`\numDstarP`,`\numDstarN` | **Gauge** (in-repo; ACRIN-6698 open) | §5(2), abstract | **YES (the null)** | No — committed seeded result | Re-extract from `invivo_real_provenance.json` |
| 5 | D\* 95% **bootstrap CI [−0.39, +0.05]** | `\numDstarCIlo`,`\numDstarCIhi` | **Gauge** (in-repo, BCa seed 20260613) | §5(2), abstract | **YES (must span zero)** | No | Fisher-z reproduction must still match (`retest_ci.py`) |
| 6 | companion **D r=+0.60**, CI [+0.42,+0.72] | `\numDR`,`\numDCIlo`,`\numDCIhi` | **Gauge** (in-repo) | §5(2) | Contrast | No | Re-extract on Gauge update |
| 7 | Sun 2019 D\*–Ktrans **r=0.389**; f·D\* **0.533** | `\numSunR`,`\numSunComposite` | **External lit** (verified) | §5(3) | Corroborating | No — published | Re-verify quote at submission |
| 8 | Yang 2019 D\* **ICC=0.55**, null | `\numYangICC` | **External lit** (verified) | §5(3) | Corroborating | No — published | Re-verify quote at submission |

**Imported from Minos:** the **theory-half** results only — the gap law `G = (1/6)|z*(λ)|·γ`, the
second-order VoI law `V = ½|EU″|G² = O(γ²)` (Plumbline Prop. 3, the Delphi result), and the
label-free floor `AUC = ½` (Plumbline Thm 2). These are **formulas**, machine-verified and
**publication-independent (SOLID)** — they carry no numeric value into Augur that could move on
Minos acceptance.

**Imported from Fashion:** a **qualitative** trust claim only (skew-aware posterior restores
marginal coverage; residual high-D\* conditional gap). **No number** is imported.

**Imported from Lethe:** the wrong-size limit is stated **qualitatively** ("far too narrow"); the
specific coverage figures (0.263 vs 0.755 target) are **not** carried into Augur's prose, so no
Lethe number is load-bearing here.

---

## 2. Why nothing held is load-bearing

The two release-blocking anchors are **Fashion** and **Minos**:
- **Minos** contributes only SOLID theory formulas → nothing numeric to move.
- **Fashion** contributes only a qualitative §Trust claim → nothing numeric to move.

The load-bearing spine numbers (rows 1–6) are all **Gauge** in-repository results, **reproduced**
this build, and Gauge is *recommended* but **not** a release-blocking anchor. Even so they are
in-repo and reproduced, not "held". Hence: **the spine survives even if Fashion's and Minos's final
numbers shift in review.** This is the explicit GATE E confirmation. (Conversely, if review forced a
Fashion/Minos result that *contradicted* the spine, that would be a scientific revision, not a
moved number — outside this ledger's scope.)

---

## 3. Finalization checklist (forward-cite → published swaps)

When an anchor publishes, do **all** of:

- [ ] **Fashion accepted** → in `paper/refs.bib` change `@unpublished{fashion}` to `@article{...}`
      with the real DOI; set `FASHION.published=true` + DOI in `release_config.json`; update
      `ASSUMPTIONS.md §1`; re-verify the §Trust claim against the published text.
- [ ] **Minos accepted** → swap `@unpublished{minos}` → `@article{...}` with DOI; set
      `MINOS.published=true` + DOI in `release_config.json`; re-confirm the Plumbline theorem/prop
      numbering against the published version.
- [ ] **Gauge accepted** (recommended) → swap `@unpublished{gauge}` → DOI; re-run `reproduce.sh`
      and confirm rows 1–6 unchanged.
- [ ] **Lethe accepted** (recommended) → swap `@unpublished{lethe}` → DOI.
- [ ] **Re-verify CITATIONS.md Tier B** (framing references inherited from Plumbline §7) against
      primary sources.
- [ ] Run `python3 release_gate.py` — it lifts only when **Fashion AND Minos** are both published
      (with DOIs). Then `bash submit.sh` proceeds to the pre-submission checklist.
- [x] **Author list + venue confirmed (GATE G, 2026-06-22): NMR in Biomedicine; Avery Karlin (solo).**
- [ ] Set the **NMR in Biomedicine** journal document class in `augur.tex` (kept generic `article` while HELD).

No PROVISIONAL flag may be cleared except through this checklist.
