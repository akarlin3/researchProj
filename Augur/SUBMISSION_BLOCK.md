# SUBMISSION_BLOCK.md — Augur is COMPLETE but HELD

**The paper is finished; only submission is held.** Augur is the end-stage synthesis of the
IVIM-UQ program. The manuscript is compiled (`paper/augur.pdf`) and `reproduce.sh` is green against
the in-repository anchors. But submitting it before its load-bearing anchors publish would (a) cite
results that may still change in review, and (b) front-run its own dependencies. So submission — and
only submission — is **held**.

## Reproduction and release are separate

- **Reproduction** (`reproduce.sh`) regenerates every in-repository anchor and is **green** (exit 0),
  independent of submission status. This is *not* gated on publication.
- **Release** (`release_gate.py`, `submit.sh`) is the **HOLD**. It refuses to submit while a
  load-bearing anchor is unpublished. (This block used to live inside `reproduce.sh` via the old
  `check_anchors.py`; it has been **lifted out** into the release path, where it belongs. The old
  `check_anchors.py` is superseded by `release_gate.py`.)

## The rule

Augur **may not be submitted** until both load-bearing anchors are published:

- **Fashion** (the trust / ruler anchor) — in review at *NMR in Biomedicine*.
- **Minos** (the value-of-information / decision anchor) — provisional.

**Lethe** and **Gauge** are strongly recommended to be out as well, since §4 of the manuscript rests
on them — but they are *recommended*, not release-blocking. (Augur's load-bearing spine is the
in-repository, reproduced Gauge anchors; see `PROVISIONAL_LEDGER.md`.)

## How the HOLD is enforced

`release_config.json` pins each anchor's published state. `release_gate.py` exits non-zero (HELD)
while either load-bearing anchor is unpublished, refreshes the `SUBMISSION_HOLD` marker, and names
the unmet conditions. `submit.sh` runs the gate and **halts** with
`HELD — awaiting Fashion + Minos publication`. `tests/test_augur.py` asserts the hold is engaged.
Today all return HELD — the intended state.

## Lifting the HOLD (the only path)

Per `PROVISIONAL_LEDGER.md §3`:

1. When an anchor publishes, set `published=true` **with its real DOI** for Fashion and Minos in
   `release_config.json`, and update `ASSUMPTIONS.md §1`. (`release_gate.py` rejects `published=true`
   with no DOI — no fabricated DOIs.)
2. Swap the `@unpublished{...}` forward-cites in `paper/refs.bib` to the published references.
3. **Re-verify `CITATIONS.md` Tier B** against primary sources.
4. Run `python3 release_gate.py` — the hold lifts only when Fashion **and** Minos are both published.
   Then `bash submit.sh` proceeds to the pre-submission checklist (author list + venue, journal
   class, final compile).

Until every step is done, Augur is a finished-but-held draft, not a submitted paper.
