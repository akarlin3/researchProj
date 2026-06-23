#!/usr/bin/env bash
# One-command re-validation of the Levy build.
#
# Runs, in order: CP0 (the fractional-order identifiability object + the kill test).
# Each stage is a gate; the script stops at the first failure. Stages not yet built are
# reported as PENDING (not failures) -- CP1+ (joint (alpha,beta) degeneracy; manuscript)
# are Phase-3 / later and appear once built.
#
# Usage:  bash Levy/reproduce.sh            # FAST (smoke; default)
#         FULL=1 bash Levy/reproduce.sh     # full-N bootstrap CIs
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
FULLFLAG=""; [ "${FULL:-0}" = "1" ] && FULLFLAG="--full"

if [ ! -x "$PROT" ]; then
  echo "FATAL: python interpreter not found at PROT=$PROT (set PROT=... to override)"; exit 2
fi

rc=0
run_stage() {  # name  script  [args...]
  local name="$1"; local script="$2"; shift 2
  echo; echo "######## $name ########"
  if [ ! -f "$script" ]; then
    echo ">>> PENDING: $script not built yet"; return 0
  fi
  if "$PROT" "$script" "$@"; then
    echo ">>> $name: PASS"
  else
    echo ">>> $name: FAIL"; rc=1
  fi
}

run_stage "CP0 identifiability object + kill test" "$HERE/verify_cp0.py" $FULLFLAG
run_stage "CP1 joint (alpha,beta) degeneracy"      "$HERE/verify_cp1.py" $FULLFLAG
run_stage "CP2 across-alpha wall robustness"       "$HERE/verify_cp2.py" $FULLFLAG
run_stage "CP-paper manuscript consistency"        "$HERE/paper/consistency.py"

echo
if [ "$rc" -eq 0 ]; then
  echo "================ reproduce.sh: all built stages GREEN ================"
else
  echo "================ reproduce.sh: a stage FAILED (rc=$rc) ================"
fi
exit "$rc"
