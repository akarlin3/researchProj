# Fashion retool — boundary-railing-first manuscript (NMR in Biomedicine)

The retooled Fashion manuscript, redrafted with a **boundary-railing-first** spine
by integrating two completed, independently validated inputs — the **Sextant**
boundary-railing primary (PR #29) and the **Gnomon CP4** spine-agnostic hand-off
(PR #31). No science was re-run here; this package *integrates* validated
components and honors the Gnomon claims ledger exactly. Target venue:
**NMR in Biomedicine** (EIC Prof. John R. Griffiths).

> **STATUS: PROVISIONAL — NOT SUBMITTED.** The manuscript is now in the official
> **NMR in Biomedicine (Wiley NJD-v2)** journal format; human review is the next
> step (see *Not done* below). Do not submit as-is.

## Build

```bash
bash build.sh          # numbers-gate (integration) -> tectonic -> manuscript.pdf
python3 consistency.py # numbers-gate only (174 macros, stdlib-only)
tectonic manuscript.tex                # manuscript (Wiley NJD-v2 class)
tectonic graphical_abstract_card.tex   # standalone graphical-abstract card
```

The manuscript compiles with the **WileyNJD-v2** class (Wiley's "New Journal
Design"), the class NMR in Biomedicine directs LaTeX authors to. Because that
class is **not on CTAN**, the three required template files are vendored beside
the source (`WileyNJD-v2.cls`, `NJDnatbib.sty`, `WileyNJD-AMA.bst`) so `tectonic`
builds the paper offline with no extra setup.

## Contents

| File | What it is |
|---|---|
| `manuscript.tex` / `.pdf` | The boundary-railing-first manuscript in the **Wiley NJD-v2** journal class (AMA numbered references). |
| `WileyNJD-v2.cls`, `NJDnatbib.sty`, `WileyNJD-AMA.bst` | Vendored Wiley NMR-in-Biomedicine template files (not on CTAN; tracked so the build is self-contained). |
| `graphical_abstract_card.tex` / `_card.pdf` | Standalone graphical abstract (image + ≤80-word text), uploaded separately — the class has no GA macro. |
| `figures/graphical_abstract.pdf` | The 50×60 mm graphical-abstract image. |
| `numbers.tex` | **Auto-generated** macros; every number traces to a seeded Sextant/Gnomon result JSON. Do not hand-edit. |
| `consistency.py` | **Numbers-gate**: re-derives `numbers.tex` from the validated source JSONs and asserts cross-source consistency; non-zero exit aborts the build. |
| `refs.bib` | Bibliography (reused from Fashion + the TCGA-LIHC dataset citation). |
| `HUANG_COVERAGE_MAP.md` | Each Huang critique → where/how the retool resolves it (CP3). |
| `cover_letter_phenomenon.tex` | Cover letter — leads with the phenomenon. |
| `build.sh` | Gate + compile. |

## Claims ledger → manuscript section map

Every disposition in `Gnomon/handoff/CLAIMS_LEDGER.md` maps to a section; no claim
appears outside the ledger.

| Ledger | Disposition | Manuscript location |
|---|---|---|
| **K1** railing (real OSIPI) | KEEP — primary | §Results/“Primary: boundary-railing across cohorts”; Table (railing) |
| **K2** quantile marginal fix | KEEP — secondary | §Results/secondary ruler (“interval shape…”) |
| **K3** flow beats railed NLLS | KEEP — secondary | §Results/secondary ruler (“amortized posterior…”) |
| **R1** Gaussian under-coverage → conditional | REFRAME | §Results/secondary ruler; Table (conditional coverage) |
| **R2** quantile residual high-$D^*$ wall | REFRAME | §Results/secondary ruler (residual gap 0.81) |
| **R3** below-floor → SD convention | REFRAME | §Methods/CRLB(b); §Results/“The convention reconstructs…” |
| **D1** marginal 0.30 / 0.67 headline | **DROP** | Not claimed; shown only as labelled floored-convention illustration |
| OOD gate / timing / brain held-out-$b$ | OUT OF SCOPE | §Discussion — explicitly **not claimed** |

## Numbers provenance (no fabrication)

The numbers-gate reads only already-seeded, already-CI'd source JSONs and never
recomputes:

- `Sextant/results/railing_results.json` — railing rates + bootstrap CIs (4 cohorts)
- `Sextant/results/{osipi,tcga_lihc}_provenance.json` — dataset IDs / DOIs
- `Gnomon/handoff/conditional_coverage.json` — per-$D^*$-tercile coverage, both SD conventions
- `Gnomon/results/reproduction.json` — keep-set (K1 clean-room rail, K2 quantile, K3 flow)

Cross-source assertions enforced by the gate: every bootstrap CI brackets its
point; Sextant's homogeneous railing matches Fashion's reported 54.7% and the
1618-voxel count; the Sextant (54.7%) and Gnomon clean-room (54.2%) real-data CIs
overlap; Gnomon's reconstruction self-check is green; the honest high-$D^*$ anchor
is ~0.63 (R1); the floored reconstruction is ~0.68 (R3); and the **dropped**
0.30/0.67 marginal headline is guarded against (asserted absent as an honest
coverage macro).

## Clean IP

OSIPI (Zenodo 14605039, CC-BY-4.0) and TCGA-LIHC (TCIA, CC BY 3.0) only —
open, download-on-demand, checksum-verified, never redistributed in-tree.
No `pancData3` / MSK data. Synthetic substrate is the read-only Lattice sibling.

## Not done (deliberately — halt for human review)

- **Not submitted.** No submission action has been taken.
- **Format pack (needle/format) APPLIED — pending author sign-off on the GA.**
  NMR-in-Biomedicine required items now present: a graphical abstract (50x60 mm
  figure `figures/graphical_abstract.pdf` + 78-word text), a single-paragraph
  ~300-word abstract (converted from the structured one), and a generative-AI
  declaration in the manuscript declarations. The GA figure choice/text awaits
  author approval before this is called submission-ready (see root `CHANGES.md`).
- **Venue template APPLIED.** The manuscript is now in the official Wiley NJD-v2
  class with AMA numbered/superscript references (`WileyNJD-AMA.bst`), the `Summary`
  abstract macro (272 words, under the 300 limit), `\keywords`, an abbreviations
  footnote, and the auto-emitted "How to cite" box. Final figure/word-count checks
  against the live author guidelines still warrant an author pass before submission.
- **Cover letter.** The phenomenon-led cover letter is provided; the author
  decides whether/how to disclose prior review history.
- **Author/affiliation block FILLED.** Avery Karlin, with three affiliations —
  Annealing Signet Institute (primary); Department of Applied Physics and Applied
  Mathematics, Columbia University; and Department of Computer Science, University
  of Colorado Boulder (student) — correspondence `ak5232@columbia.edu`. Confirm the
  ORCID/email and primary-affiliation order before submission.
- **OUT-OF-SCOPE items** (OOD gate, timing, brain held-out-$b$) remain unevaluated
  by design; if a future revision wants them, they must be sourced elsewhere.
