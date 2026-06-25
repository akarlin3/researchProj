# Limbo — a field review of trustworthy UQ for quantitative body MRI and its use in adaptive RT

**Status: PROVISIONAL · submission-ready compiled manuscript (CP0–CP3 complete) · not publish-gated.**

> **phiRO retarget (CP4 complete; GATE D confirmed).** The manuscript has been reformatted from the
> IOP `iopjournal` class (the old *Phys. Med. Biol.* target) to **Elsevier `elsarticle`** for
> submission to **Physics and Imaging in Radiation Oncology (phiRO)** — open access, via Editorial
> Manager. The phiRO manuscript is [`limbo_phiro.tex`](limbo_phiro.tex) → `limbo_phiro.pdf`; the
> original [`limbo.tex`](limbo.tex) (IOP) is retained as the **content-identity reference**. The
> reformat changed **zero content**: body prose, all 59 references, and the four verbatim quotes are
> byte-identical (proven by diff; enforced by `tests/test_phiro_format.py`). GATE D (HUMAN)
> confirmed: institutional corresponding email `ak5232@columbia.edu`; title and the three
> affiliations as-is; competing-interest + funding statements as-is; and the pre-existing
> `vanhoudt2021qib` ledger en-dash ("test–retest") corrected to the source-faithful hyphen
> ("test-retest"). Submission is via Editorial Manager; the cover letter is a separate upload.

The manuscript [`limbo.tex`](limbo.tex) → [`limbo.pdf`](limbo.pdf) is typeset for **Physics in
Medicine & Biology** (Topical Review) with IOP's `iopjournal` class; it compiles clean with
`tectonic` (gated on the citation gate), the survey cites all 59 verified entries with zero phantom
prose-cites, and the four thesis-level entries have been re-pulled **verbatim** from source (see
[`CITATIONS.md`](CITATIONS.md), *Verbatim re-pulls*).

Limbo is a **broad field review** — a survey of *the literature's* work on trustworthy uncertainty
quantification (UQ) for **quantitative/diffusion body MRI** (IVIM, DWI/ADC, DKI, DCE perfusion,
relaxometry) and its **decision-use in MR-guided adaptive radiotherapy** (MR-Linac / MRgART). It
organises the field along a **trust → value-of-information → action** axis and maps where the
field's UQ-trust questions remain open.

Its value is **trigger-independent**: it strengthens the *first* PhD application (field command + a
citable paper), and it **absorbs Buttress** (the portfolio-thickener) — there is no separate
Buttress; this review *is* the thickener.

Target venue: **Physics and Imaging in Radiation Oncology (phiRO)** — Elsevier, open access, review
article, via Editorial Manager (retargeted from the CP0 input of *Physics in Medicine & Biology*,
Topical Review; see the phiRO-retarget note above). Submission is via Editorial Manager; the cover
letter is a separate upload; an APC applies (open access).

## What this is — and is not

- It **is** synthesis, taxonomy, and gap-identification across the external literature.
- It makes **no new measurement** and asserts **no experimental result** of its own.
- It is **distinct from Augur**, and that distinctness is the CP0 gate (see below).

## Distinct from Augur (CP0 verdict: separable)

[`Augur/`](../Augur) is a *perspective on the author's own arc* (Fashion / Minos / Lethe / Gauge),
hard-blocked from submission until those papers publish. **Limbo is a field survey of _others'_
work**, not publish-gated, organised around the same trust→action spine used as a neutral
classification axis. The author's own papers appear in Limbo, if at all, only as a *minority* of
entries cited on equal footing with external work — never as the organising centre. Full table and
collapse-risk guardrails in [`ASSUMPTIONS.md`](ASSUMPTIONS.md).

## The hard gate — verified citations, zero phantoms

A review's dominant failure mode is the phantom citation (this portfolio's own history: Ouroboros's
non-existent "Sun et al."; Augur's mis-quoted "r≈0.39"). Limbo makes that mechanical:

- [`limbo.bib`](limbo.bib) — **59 verified entries**, each with a resolvable **DOI / arXiv / stable
  proceedings URL**.
- [`CITATIONS.md`](CITATIONS.md) — a one-line **verified claim** per citekey, traceable to source;
  all identifiers resolved against primary sources on 2026-06-22.
- [`verify_citations.py`](verify_citations.py) — fails the build on any entry without a resolvable
  identifier, or any orphan between bib and ledger. `--online` additionally HEAD-checks resolvability
  (CP3).

```
python3 verify_citations.py           # offline gate: 59 entries, zero unverifiable -> exit 0
python3 verify_citations.py --online  # also confirm each DOI/arXiv resolves (network)
./build.sh                            # citation gate -> compile limbo_phiro.tex (tectonic) -> limbo_phiro.pdf
./build.sh --iop                      # build the retained IOP version (limbo.tex -> limbo.pdf)
./reproduce.sh                        # gate + tests + phiRO compile (green == submission-ready)
```

## Layout

| file | purpose |
|---|---|
| `limbo_phiro.tex` / `limbo_phiro.pdf` | **the phiRO submission manuscript** (Elsevier `elsarticle`; target *Phys. Imaging Radiat. Oncol.*) |
| `limbo.tex` / `limbo.pdf` | the IOP `iopjournal` version (old *Phys. Med. Biol.* target) — retained as the **content-identity reference** |
| `TAXONOMY.md` | the trust → VoI → action survey axis (+ foundations + gap-map seams) |
| `SURVEY.md` | the markdown survey draft the manuscript prose was ported from |
| `limbo.bib` | the verified citation base (59 entries) |
| `CITATIONS.md` | per-citekey verified claim + resolvable identifier (+ the 4 verbatim re-pulls) |
| `verify_citations.py` | the citation gate (offline; `--online` for resolvability; scans both `.tex` manuscripts) |
| `ASSUMPTIONS.md` | scope boundary, distinctness-from-Augur, clean-IP, status pins |
| `build.sh` | gate → compile the manuscript with tectonic (`--iop` to build the IOP version) |
| `reproduce.sh` | one-command re-validation (gate + pytest + phiRO compile) |
| `elsarticle.cls`, `elsarticle-num.bst` | vendored Elsevier class + numbered-Vancouver bib style (LPPL; see `ELSEVIER_CLASS_PROVENANCE.md`) |
| `iopjournal.cls`, `orcid.pdf` | vendored IOP class + asset (LPPL; see `IOP_CLASS_PROVENANCE.md`) |
| `tests/` | gate assertions + phiRO reformat checks (quote-identity, cite-set == IOP, content-identity) |

## Checkpoints

- **CP0 — audit + scope + distinctness (HALT, cleared).** Embedding confirmed (top-level subrepo,
  Augur-minimal pattern); scope boundary set; trust→VoI→action taxonomy fixed; distinctness from
  Augur proven (separable); venue taken as input (PMB).
- **CP1 — framework + verified citation base.** Taxonomy built; 59-entry `.bib`, each with a
  resolvable id + verified claim; gate passes with **zero unverifiable entries**.
- **CP2 — survey + gap map.** Survey drafted by axis with the open-problems map (G1–G4); every claim
  cites a verified entry.
- **CP3 — honest-scope + final citation gate + compiled manuscript.** `--online` re-verification
  (all 59 resolve live), no-drift pass, honest-scope section; **`limbo.tex` typeset for PMB
  (`iopjournal`) and compiled to `limbo.pdf`**; the four thesis-level entries re-pulled verbatim
  from source; Buttress absorbed into the discussion/gap map. Staged for review (no auto-merge).
- **CP4 — phiRO retarget (content-preserving reformat).** IOP `iopjournal` → Elsevier `elsarticle`
  for **Physics and Imaging in Radiation Oncology**; numbered-Vancouver references
  (`elsarticle-num`); IOP back-matter macros split into the phiRO declaration `\section*` blocks
  (CRediT · competing interest · funding · data availability · generative-AI). **Content changed by
  zero** — body prose, all 59 refs, and the 4 quotes byte-identical to `limbo.tex` (proven by diff;
  enforced by `tests/test_phiro_format.py`). GATE D (HUMAN) confirmed: institutional corresponding
  email `ak5232@columbia.edu`; the `vanhoudt2021qib` ledger en-dash corrected to a source-faithful
  hyphen. `limbo_phiro.pdf` builds clean (0 unresolved, 59 refs); gate exit 0 (offline + `--online`);
  `reproduce.sh` green. Submission via Editorial Manager (cover letter = separate upload; APC applies).
  Staged for review (no auto-merge).
