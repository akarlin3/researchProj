# Limbo — a field review of trustworthy UQ for quantitative body MRI and its use in adaptive RT

**Status: PROVISIONAL · CP1 complete (taxonomy + verified citation base + hard gate) · not publish-gated.**

Limbo is a **broad field review** — a survey of *the literature's* work on trustworthy uncertainty
quantification (UQ) for **quantitative/diffusion body MRI** (IVIM, DWI/ADC, DKI, DCE perfusion,
relaxometry) and its **decision-use in MR-guided adaptive radiotherapy** (MR-Linac / MRgART). It
organises the field along a **trust → value-of-information → action** axis and maps where the
field's UQ-trust questions remain open.

Its value is **trigger-independent**: it strengthens the *first* PhD application (field command + a
citable paper), and it **absorbs Buttress** (the portfolio-thickener) — there is no separate
Buttress; this review *is* the thickener.

Target venue (CP0 input): **Physics in Medicine & Biology** (Topical Review).

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
python3 verify_citations.py --online  # also confirm each DOI/arXiv resolves (network; CP3)
./reproduce.sh                        # gate + tests
```

## Layout

| file | purpose |
|---|---|
| `TAXONOMY.md` | the trust → VoI → action survey axis (+ foundations + gap-map seams) |
| `limbo.bib` | the verified citation base (59 entries) |
| `CITATIONS.md` | per-citekey verified claim + resolvable identifier (the ledger the gate checks) |
| `verify_citations.py` | the citation gate (offline; `--online` for resolvability) |
| `ASSUMPTIONS.md` | scope boundary, distinctness-from-Augur, clean-IP, status pins |
| `reproduce.sh` | one-command re-validation (gate + pytest) |
| `tests/` | gate assertions (zero-unverifiable, bucket coverage, distinctness documented) |

## Checkpoints

- **CP0 — audit + scope + distinctness (HALT, cleared).** Embedding confirmed (top-level subrepo,
  Augur-minimal pattern); scope boundary set; trust→VoI→action taxonomy fixed; distinctness from
  Augur proven (separable); venue taken as input (PMB).
- **CP1 — framework + verified citation base (this commit).** Taxonomy built; 59-entry `.bib`, each
  with a resolvable id + verified claim; gate passes with **zero unverifiable entries**.
- **CP2 — survey + gap map (pending).** Draft the survey by axis and the open-problems map; every
  claim cites a verified entry.
- **CP3 — honest-scope + final citation gate (pending).** `--online` re-verification, no-drift pass,
  honest-scope section (synthesis + gaps, no new results); update root README; stage PR (no auto-merge).
