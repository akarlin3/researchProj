# Citable-release path (documented, NOT executed)

Datum is usable internally now, under the provisional flags. A **citable** release
is intentionally *not* executed yet — it is gated on the same condition as Caliper:
**Fashion's calibration ruler locking** (MRM acceptance + DOIs). This mirrors the
monorepo discipline (cf. Caliper's deliberately un-gated, non-citable status and
Minos `future/`'s provisional build).

## Why gated

Every Datum reference number is scored on Fashion's in-review ruler. Releasing a
citable benchmark now would publish numbers that can still change in revision. So
the release waits until `datum.manifest.RULER['manuscript_doi']` is assigned.

## The path, when the ruler locks

1. **Lock the pins.** Set `RULER['manuscript_doi']` (and bump any changed
   `version`/`commit`) in `datum/manifest.py`; update `ASSUMPTIONS.md`.
2. **Re-validate.** `python revalidate.py --full` regenerates every reference
   number under the locked ruler; the PROVISIONAL flag clears automatically
   (`is_provisional()` keys off the DOI).
3. **JOSS.** Add a `paper.md` + `paper.bib` (JOSS format) describing the benchmark;
   Datum already ships an OSI-approved licence (MIT), tests, and docs.
4. **Zenodo.** Archive a tagged release for a DOI (as Proteus/Caliper do), and add
   a `CITATION.cff` citing Fashion's ruler (now with its DOI) and the OSIPI DRO
   (DOI `10.5281/zenodo.14605039`).
5. **Cross-link.** Register the Datum DOI badge in the monorepo README, alongside
   the others.

Until then: **do not cite Datum numbers as final.**
