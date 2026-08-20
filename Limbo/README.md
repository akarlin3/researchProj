# Limbo — a field review of trustworthy UQ for quantitative body MRI and its use in adaptive RT

**Status: PROVISIONAL · submission-ready compiled manuscript (v2, RT-first, phiRO) · not publish-gated.**

The manuscript [`limbo.tex`](limbo.tex) → [`limbo.pdf`](limbo.pdf) is typeset for **Physics and
Imaging in Radiation Oncology** (phiRO; Elsevier / ESTRO) as a **Review article**, using Elsevier's
`elsarticle` class; it compiles clean with `tectonic` (gated on the citation gate), the survey cites
all 66 verified entries with zero phantom prose-cites, and every reference carries a resolvable
DOI/arXiv identifier with a one-line verified claim (see [`CITATIONS.md`](CITATIONS.md); the four
foundation entries are re-pulled **verbatim** from source in *Verbatim re-pulls*).

Limbo is a **broad field review** — a survey of *the literature's* work on trustworthy uncertainty
quantification (UQ) for **quantitative/diffusion body MRI** (IVIM, DWI/ADC, DKI, DCE perfusion,
relaxometry) and its **decision-use in MR-guided adaptive radiotherapy** (MR-Linac / MRgART). It is
re-anchored *radiation-oncology-first* — the **MR-Linac dose decision** is the spine, and
quantitative-MRI UQ enters as the input that decision consumes, read along a **trust → value →
action** axis — and it maps where the field's UQ-trust questions remain open for adaptive RT.

Its value is **trigger-independent**: it stands on its own scientific merit and does not depend on
the author's own results publishing.

Target venue: **Physics and Imaging in Radiation Oncology** (phiRO; Elsevier / ESTRO), Review article.

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
./build.sh                            # citation gate -> compile limbo.tex (tectonic) -> limbo.pdf
./reproduce.sh                        # gate + tests + manuscript compile (green == submission-ready)
```

## Layout

| file | purpose |
|---|---|
| `limbo.tex` / `limbo.pdf` | the compiled manuscript (Elsevier `elsarticle`; target *Phys Imaging Radiat Oncol*) |
| `TAXONOMY.md` | the trust → VoI → action survey axis (+ foundations + gap-map seams) |
| `SURVEY.md` | the markdown survey draft the manuscript prose was ported from |
| `limbo.bib` | the verified citation base (66 entries) |
| `CITATIONS.md` | per-citekey verified claim + resolvable identifier (+ the 4 verbatim re-pulls) |
| `verify_citations.py` | the citation gate (offline; `--online` for resolvability; scans `limbo.tex`) |
| `ASSUMPTIONS.md` | scope boundary, distinctness-from-Augur, clean-IP, status pins |
| `build.sh` | gate → compile the manuscript with tectonic |
| `reproduce.sh` | one-command re-validation (gate + pytest + compile) |
| `elsarticle.cls`, `elsarticle-num.bst` | vendored Elsevier class + numbered-Vancouver bib style (LPPL; see `ELSEVIER_CLASS_PROVENANCE.md`) |
| `tests/` | gate assertions (zero-unverifiable, bucket coverage, distinctness documented) |

## Checkpoints

- **CP0 — audit + scope + distinctness (HALT, cleared).** Embedding confirmed (top-level subrepo,
  Augur-minimal pattern); scope boundary set; trust→VoI→action taxonomy fixed; distinctness from
  Augur proven (separable); venue taken as input (PMB).
- **CP1 — framework + verified citation base.** Taxonomy built; 59-entry `.bib`, each with a
  resolvable id + verified claim; gate passes with **zero unverifiable entries**.
- **CP2 — survey + gap map.** Survey drafted by axis with the open-problems map (G1–G4); every claim
  cites a verified entry.
- **CP3 — honest-scope + final citation gate + compiled manuscript.** `--online` re-verification
  (all entries resolve live), no-drift pass; the four foundation entries re-pulled verbatim
  from source. (Original PMB/`iopjournal` build; superseded by v2 below.)
- **v2 — RT-first re-anchor + phiRO.** Re-anchored radiation-oncology-first (the MR-Linac dose
  decision as the spine); §7 rewritten as a plain Discussion/Limitations with all
  portfolio/strategy/positioning prose removed; reference base rebalanced to **66** (added the
  MR-Linac validation / dose-painting-trial / on-device-repeatability literature, demoted the
  over-weight general-ML-UQ and decision-theory tail); journal names abbreviated to ISO/LTWA;
  reformatted IOP→phiRO (`elsarticle`, *Review article*, numbered Vancouver); identity unified to
  `ak5232@columbia.edu`. Gate green offline + online. Staged for the human submission gate (no
  auto-merge).
