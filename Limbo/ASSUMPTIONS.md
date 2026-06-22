# ASSUMPTIONS.md — Limbo scope, distinctness, and status pins

**Limbo is a broad _field review_** of trustworthy uncertainty quantification (UQ) for
quantitative/diffusion body MRI and its decision-use in MR-guided adaptive radiotherapy. Its
contribution is **synthesis, taxonomy, and gap-identification** — not new measurements. Target
venue (author input, CP0): **Physics in Medicine & Biology** (Topical Review; spans both the
quantitative-MRI-methods half and the MR-guided-RT half of the spine).

## Status

| item | value |
|---|---|
| Stage | CP1 complete (taxonomy + verified citation base + gate). CP2/CP3 pending. |
| Standalone value | **Trigger-independent.** Strengthens the *first* PhD application (field command + a citable review); does **not** depend on the author's papers publishing. |
| Submission gating | **NOT publish-gated.** Unlike Augur, Limbo can be submitted on its own scientific merit. |
| Embedding | Embedded subrepo (`Limbo/`) for now; carve out to `projLimbo` later **iff** it grades portfolio-worthy (the Caliper/Lattice precedent, PR #59). |
| Clean IP | Published literature + open material only. **No pancData3 / MSK / patient data** in tree or history. |

## Scope boundary (CP0)

**IN.** (a) Trustworthiness/calibration of *uncertainty estimates* for quantitative & diffusion
body MRI — IVIM/DWI/DKI, DCE/DSC perfusion, relaxometry: CRLB/Fisher, Bayesian posteriors,
bootstrap, conformal prediction, deep-learning UQ (ensembles, MC-dropout, evidential), test–retest
repeatability (QIBA), coverage/calibration metrics. (b) The decision-use of those quantitative
biomarkers + their uncertainty in MR-guided adaptive radiotherapy (MR-Linac / MRgART): response-
adaptive dose, dose painting, gating/triggering, robust re-optimisation under biomarker uncertainty.

**OUT.** Non-quantitative/qualitative MRI; brain-only neuro work; generic computer-vision UQ with no
quantitative-MRI tie; reconstruction UQ that never reaches a parameter map; radiotherapy that is not
MR-guided; and — load-bearing — **the author's own unpublished results as primary content** (those
belong to Augur and the individual papers, not to this field survey).

## Distinctness from Augur (CP0 gate — verdict: SEPARABLE)

Limbo and Augur deliberately share the **trust → value-of-information → action** spine. Everything
else separates cleanly, and the following guardrails keep Limbo from collapsing into Augur:

| axis | **Augur** | **Limbo** |
|---|---|---|
| object surveyed | the author's **own** four papers (Fashion/Minos/Lethe/Gauge) | the **external field** (other groups' literature) |
| document type | perspective / synthesis of one's own arc | field **review** / literature survey |
| spine role | narrative arc; the "D\* cross-modally orphaned" thread | neutral taxonomy any paper maps onto |
| publish-gating | **hard-blocked** until Fashion + Minos publish (`Augur/check_anchors.py`) | **not gated**; standalone now |
| dependency on own arc | total (it *is* the arc) | none (own work = a minority of entries, peer-cited) |
| citation base | 2 external (Tier A) + inherited (Tier B) | a large verified **external** base is the bulk |
| value timing | end-stage (post-publication) | first-application value **now** |

**Collapse-risk guardrails (enforced as editorial rules through CP2/CP3):**
1. The author's own papers appear, if at all, as a *minority* of entries cited on equal footing
   with external work — never as the survey's organising centre. They are tracked separately in
   `CITATIONS.md` ("In-portfolio cross-references"), not asserted as published literature.
2. Limbo does **not** import Augur's "D\* cross-modally orphaned" narrative as its spine; the spine
   is a neutral survey axis.
3. The gap map is the **field's** open problems, with the thesis mentioned only as "where our
   program sits" — honest positioning, not overclaim.

## Buttress

Limbo **absorbs Buttress** (the portfolio-thickener). There is no separate Buttress repo; this
review *is* the thickener — a citable artifact that adds field command to the portfolio.

## Citation discipline (the hard gate)

Verified citations are the dominant risk for any review and the repeat failure mode in this
portfolio (Ouroboros's phantom "Sun et al."; Augur's mis-quoted "r≈0.39"). Therefore:

- every `limbo.bib` entry carries a **resolvable identifier** (DOI / arXiv / stable proceedings URL);
- every entry has a **one-line verified claim** in `CITATIONS.md`, traceable to the source;
- `verify_citations.py` fails the build on any orphan or identifier-less entry (run by `reproduce.sh`
  and the test-suite). CP1 verified all identifiers against primary sources on **2026-06-22**; CP3
  re-runs with `--online` resolvability checks and a final no-drift pass.

Four entries (`koo2016`, `blandaltman1986`, `ling2000`, `vanhoudt2021qib`) carry **thesis-level**
claims (the paper's title/abstract-evident thesis) rather than a quoted number; any verbatim
number from these must be re-pulled from the PDF before use. Flagged in `CITATIONS.md`.
