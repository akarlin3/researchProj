# RELEASE.md — the Matrix submission HOLD

Matrix's manuscript is **submission-ready but HELD**. It is complete, compiles, and
**reproduces green on both substrates** (the synthetic twin and the Ferry real-data
substrate); submission is withheld until the two consumed components it leans on publish.

```
RELEASE  iff  FASHION_PUBLISHED  AND  MINOS_PUBLISHED          (both DOIs in release.json)
otherwise HELD — awaiting Fashion + Minos publication
```

**Forge is not a condition.** It is NOT-BUILT (deferred 2027) and is presented as drop-in
future work; it never holds submission. This is recorded explicitly in `release.json`
under `not_conditions` so it is on the record that Forge does not gate release.

## Why a separate gate (not in `reproduce.sh`)

Reproduction and release are different questions:

- `reproduce.sh` — *"does it still reproduce?"* Validates the **harness** on placeholder
  components (NOT-Fashion / NOT-Minos / NOT-Forge). Publication-agnostic; no Fashion/Minos
  short-circuit. Green today.
- `release_gate.py` — *"may it be submitted yet?"* Encodes the publication HOLD. Held today.

Keeping them decoupled means the science keeps reproducing while submission stays correctly
blocked.

## Usage

```bash
python release_gate.py status   # HELD/RELEASE + unmet conditions; refresh the SUBMISSION_HOLD marker
python release_gate.py submit   # refuse (exit 3) with the held message while HELD; else proceed
python release_gate.py check    # exit 0 iff RELEASE, else 3 (for scripts/CI)
```

## To release (when the papers land)

1. When **Fashion** publishes: in `release.json`, set `FASHION_PUBLISHED.satisfied = true`
   and fill its `doi`. Update `ASSUMPTIONS.md` §Fashion and the `\bibitem{fashion}` in
   `paper/matrix.tex` (drop the `PROVISIONAL-fashion` marker).
2. When **Minos** publishes: the same for `MINOS_PUBLISHED` / `\bibitem{minos}`.
3. Run `python release_gate.py status`. With both satisfied → RELEASE; the `SUBMISSION_HOLD`
   marker is removed. Re-run `bash reproduce.sh` once with the real adapters wired in
   (per `PROMOTION.md`) to confirm every gate stays green, then submit.

The central honest negative (F1: action-suppression ≠ outcome-protection) does **not** depend
on Fashion/Minos being correct, so it survives the drop-in (see the stub ledger in
`PROMOTION.md` / `ASSUMPTIONS.md`).
