#!/usr/bin/env bash
# One-command re-validation of the Matrix synthetic-twin closed loop.
#
# Runs, in order: the pytest suite -> CP1 (twin + harness) -> CP2 (ruler + gates) ->
# CP3 (Forge-shaped dose stage) -> CP4 (closed loop + honest scope). Each stage is a gate;
# the script reports each and exits non-zero if any fails.
#
# This is also the re-validation contract from ASSUMPTIONS.md: when Fashion / Minos / Forge
# land (or revise), point the relevant interface at the real adapter (one line in
# matrix/loop.py:Interfaces), bump the pins in ASSUMPTIONS.md, and run this once to learn
# what still holds.
#
# Usage:  bash Matrix/reproduce.sh
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"

if [ ! -x "$PROT" ]; then
  echo "FATAL: python interpreter not found at PROT=$PROT (set PROT=... to override)"; exit 2
fi

rc=0
run_stage() {  # name  cmd...
  local name="$1"; shift
  echo; echo "######## $name ########"
  if "$@"; then echo ">>> $name: PASS"; else echo ">>> $name: FAIL"; rc=1; fi
}

export PYTHONPATH="$HERE:${PYTHONPATH:-}"

run_stage "pytest suite"                 "$PROT" -m pytest "$HERE/tests" -q
run_stage "CP1 twin + harness"           "$PROT" "$HERE/verify_cp1.py"
run_stage "CP2 ruler + trust/action gates" "$PROT" "$HERE/verify_cp2.py"
run_stage "CP3 Forge-shaped dose stage"  "$PROT" "$HERE/verify_cp3.py"
run_stage "CP4 closed loop + honest scope" "$PROT" "$HERE/verify_cp4.py"

echo
if [ "$rc" -eq 0 ]; then
  echo "================ reproduce.sh: all gates GREEN ================"
else
  echo "================ reproduce.sh: a gate FAILED (rc=$rc) ================"
fi
exit "$rc"
