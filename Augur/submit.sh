#!/usr/bin/env bash
# Augur SUBMIT path -- gated by the release gate, separate from reproduction.
#
# This is the only script that "submits". It refuses while the hold is active and prints exactly
# which release condition(s) are unmet. It does NOT run the reproduction pipeline (that is
# reproduce.sh's job and is already green); release and reproduction are deliberately separate.
#
#   bash submit.sh             # checks the release gate; HALTS if held
#
# Lifting the hold: set published=true WITH a real DOI for Fashion and Minos in
# release_config.json (per PROVISIONAL_LEDGER.md), then re-run.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
if [ ! -x "$PROT" ]; then PROT="$(command -v python3)"; fi

echo "### Augur submit path"
echo

"$PROT" "$HERE/release_gate.py"
GATE_RC=$?
echo

if [ "$GATE_RC" -ne 0 ]; then
  echo "=============================================================="
  echo " HELD -- awaiting Fashion + Minos publication."
  echo " Submission refused. See SUBMISSION_HOLD and PROVISIONAL_LEDGER.md."
  echo "=============================================================="
  exit 1
fi

# Release gate is CLEAR (both load-bearing anchors published). Final pre-submission checklist.
echo "=============================================================="
echo " RELEASE GATE CLEAR -- proceeding to pre-submission checklist."
echo "=============================================================="
echo "  [ ] Re-verify CITATIONS.md Tier B against primary sources."
echo "  [ ] Swap forward-cites (Fashion, Minos) to published DOIs in paper/refs.bib."
echo "  [ ] Confirm author list + venue (GATE G) and set the journal document class."
echo "  [ ] Run: bash reproduce.sh && bash paper/build.sh   (final clean compile)."
echo
echo ">>> Augur is cleared for submission once the checklist above is complete."
exit 0
