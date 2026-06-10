# Nonlinear Dynamics (Springer) reformat — changelog

Venue switch of *Death of a Chimera* from REVTeX (PRE/Chaos) to the **Springer
Nature `sn-jnl` template** for submission to **Nonlinear Dynamics**. This is a
**reformat, not a revision**: no result, number, equation, figure, table, or
reference content changed. The science is untouched.

## New files

| File | Purpose |
|---|---|
| `death_of_a_chimera_nld.tex` | The Springer/ND manuscript (sn-jnl class). |
| `death_of_a_chimera_nld.bib` | BibTeX bibliography (all 19 refs, transcribed verbatim). |
| `sn-jnl.cls` | Official Springer Nature class (v2.x, from the Springer Nature LaTeX template). |
| `sn-mathphys.bst` | Math & Physical Sciences reference style (numbered) — the style used. |
| `sn-basic.bst` | Basic Springer Nature style (numbered) — shipped for one-token swap (see below). |

The REVTeX source (`death_of_a_chimera_aip.tex`) and the plain-`article` source
(`death_of_a_chimera.tex`) are **retained unchanged**.

## Target class — confirmed (not guessed)

- **Class:** `sn-jnl.cls` — the universal Springer Nature journal class (the live
  ND submission page is auth-walled; class + options confirmed from the official
  Springer Nature LaTeX template bundle).
- **Reference style:** **numbered, square brackets, in order of citation** —
  confirmed for *Nonlinear Dynamics* via its citation-style record (`[1]`,
  `[1, 2]`, `[1–4]`, in order of appearance). Implemented with
  `\documentclass[pdflatex,sn-mathphys,Numbered]{sn-jnl}` (Math & Physical
  Sciences numbered style — the disciplinarily-correct fit for a
  Phys. Rev./Chaos/Physica D paper).
- **Alternative:** if the ND submission page specifies the *Basic* style, change
  the option to `sn-basic` (both `.bst` files are shipped). Either way the result
  is numbered `[n]` citations in order — only minor reference-list punctuation
  differs.

## REVTeX → sn-jnl mapping applied

| REVTeX construct | sn-jnl equivalent |
|---|---|
| `\documentclass[aip,preprint]{revtex4-2}` | `\documentclass[pdflatex,sn-mathphys,Numbered]{sn-jnl}` |
| `\author{}` + `\affiliation{}` + `\email{}` + `\thanks{ORCID}` | `\author*[1]{\fnm{Avery} \sur{Karlin}\orcid{0000-0003-3848-6782}}` + `\email{averykarlin3@gmail.com}` + `\affil*[1]{\orgname{Independent Researcher}}` |
| `\begin{abstract}…\end{abstract}` | `\abstract{…}` (text kept byte-for-byte) |
| *(none)* | `\keywords{chimera states, coupled phase oscillators, Kuramoto–Sakaguchi model, finite-size scaling, transient dynamics, Weibull aging hazard}` (new, required by Springer) |
| `\begin{quotation}` lay summary | `\begin{quote}\itshape …` (preserved) |
| `figure`/`table` + booktabs | unchanged (`\includegraphics`, `\graphicspath` preserved) |
| `\appendix` | `\begin{appendices} … \end{appendices}` |
| inline `thebibliography` | `\bibliography{death_of_a_chimera_nld}` (sn-mathphys numbered) |
| `\begin{acknowledgments}` | `\bmhead{Acknowledgments}` (replaced with the Springer AI disclosure) |
| Author Declarations block | `\section*{Declarations}` (Funding / Competing interests / Data and code availability / Author contributions) |

> Note: `sn-jnl`'s `\orcid` macro renders an ORCID-logo image (`Orcidlogo.eps`)
> that is not bundled; it is redefined in the preamble to typeset the iD as linked
> text, so the document compiles self-contained.

## Text edits (the only content-level changes)

1. **AI disclosure** — Acknowledgments replaced with the Springer-compliant
   disclosure (drafting/editing role stated; author sole and accountable; AI not
   an author).
2. **Language tempering** (3 edits, verified as the *only* body differences vs. the
   REVTeX source):
   - section heading "…never lets go" → "…does not absorb"
   - "the only immortal chimeras" → "the only persistent, non-absorbing chimeras"
   - "The stability folklore of the deployment" → "The stability assumption of the deployment"

## Carryover-edit audit (flagged, not invented)

- **"Eq. (3)" → "Eq. (2)" — no such reference exists.** Neither manuscript contains
  any internal "Eq. (3)" (or "Eq. (2)") cross-reference. The only equation refs are
  *external* to Abrams 2008 ("Eq. (12) of Ref. [3]", "Eq. 17/18 of Ref."), which are
  correct. This carryover was already resolved in this repo; **nothing fabricated.**
- **Figure numbering is already monotonic** by first citation
  (Fig. 1→9, then A, B). No reordering needed.
- **Captions are complete in source** — every `\caption{}` is a full, punctuated
  sentence. No truncation at the source level (see build caveat below).

## Verification performed

- Body (Introduction → Conclusion) diffed against the REVTeX source: the **only**
  differences are the three tempering edits above. **Appendices byte-identical.**
- Abstract: **byte-identical** to the current source abstract.
- All **19** `\cite` keys resolve to `.bib` entries; no missing, no orphan entries.
- No REVTeX-only commands remain (`\affiliation`, `\thanks`, `acknowledgments`, etc.).

## Build caveat (could not satisfy locally)

- **No LaTeX toolchain in this environment** (no `pdflatex`/`bibtex`) and **no
  poppler** (could not rasterize the existing PDF). The project is **build-ready**
  but the final **"compile clean" check must be run on Overleaf** (which runs
  BibTeX automatically) or any TeX Live install: place `death_of_a_chimera_nld.tex`,
  `death_of_a_chimera_nld.bib`, `sn-jnl.cls`, and the chosen `.bst` in one folder
  with `paper_figures/`, and compile pdflatex → bibtex → pdflatex × 2.
- On that compile, confirm: monotonic figure numbers, no clipped captions,
  abstract + keywords present, declarations present, all `[n]` references resolved.

## Open item for the author to confirm

- **Reference sub-style** (`sn-mathphys` vs `sn-basic`) against the live ND
  submission page — one-token change; both `.bst` shipped.
- **Graphical abstract:** Springer ND may request one for initial submission. The
  manuscript has no dedicated graphical-abstract file; a key figure (e.g. Fig. 7,
  the (β, A) stability diagram, or Fig. 1) can serve. Flagged, not added.
