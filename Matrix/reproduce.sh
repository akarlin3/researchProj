#!/usr/bin/env bash
# One-command re-validation of the Matrix closed loop, over BOTH substrates.
#
# Runs, in order:
#   pytest suite
#   CP1 (twin + harness) -> CP2 (ruler + gates) -> CP3 (Forge-shaped dose) -> CP4 (closed loop)
#       ... on the SYNTHETIC twin   [offline; regenerates results/RESULTS_CP2.json, RESULTS_CP4.json]
#   loop.py BYTE-IDENTITY gate      [the engine is byte-for-byte unchanged across substrates]
#   Ferry CP1 (drop-in proof)       [offline-capable; byte-unchanged + contract + reproducible]
#   Ferry CP2 (grounded closed loop on REAL anatomy + dose) [regenerates results/RESULTS_FERRY_CP2.json]
#       ... on the FERRY real-data substrate. SKIPPED (not failed) iff the public TCIA
#           dataset is unreachable AND no local cache exists -- a real-data anchor cannot be
#           regenerated without the real data; every other gate still gates the build.
#
# Each stage is a gate; the script reports each and exits non-zero if any *required* gate fails.
#
# SEPARATION OF CONCERNS: this script validates the HARNESS on placeholder components
# (NOT-Fashion / NOT-Minos / NOT-Forge); it is publication-agnostic and contains NO
# Fashion/Minos short-circuit. The submission HOLD (awaiting Fashion + Minos publication) lives
# entirely in the separate release gate -- see release_gate.py / RELEASE.md. Reproduction and
# release are deliberately decoupled: `reproduce.sh` answers "does it still reproduce?", the
# release gate answers "may it be submitted yet?".
#
# Usage:  bash Matrix/reproduce.sh
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"

if [ ! -x "$PROT" ]; then
  echo "FATAL: python interpreter not found at PROT=$PROT (set PROT=... to override)"; exit 2
fi

# The git-blob hash matrix/loop.py shipped with (Matrix PR #56). The synthetic and the Ferry
# grounded runs invoke this exact engine; the substrate swap edits it in no way.
LOOP_PY_BLOB_SHA1="4a34806ac4fa55c0ce5453b9864d37c67abfda92"

rc=0
run_stage() {  # name  cmd...   (required: failure -> rc=1)
  local name="$1"; shift
  echo; echo "######## $name ########"
  if "$@"; then echo ">>> $name: PASS"; else echo ">>> $name: FAIL"; rc=1; fi
}

run_stage_optional() {  # name  cmd...   (exit 2 from cmd -> SKIP, not FAIL)
  local name="$1"; shift
  echo; echo "######## $name ########"
  "$@"; local ec=$?
  if   [ "$ec" -eq 0 ]; then echo ">>> $name: PASS"
  elif [ "$ec" -eq 2 ]; then echo ">>> $name: SKIP (real-data substrate unavailable; anchor not regenerated)"
  else echo ">>> $name: FAIL"; rc=1; fi
}

byte_identity_gate() {  # the engine is byte-for-byte unchanged across substrates
  echo; echo "######## loop.py byte-identity (synthetic == Ferry) ########"
  local got
  got="$(git -C "$HERE" hash-object matrix/loop.py 2>/dev/null \
         || "$PROT" - "$HERE/matrix/loop.py" <<'PY'
import hashlib, sys
d = open(sys.argv[1], "rb").read()
print(hashlib.sha1(b"blob %d\0" % len(d) + d).hexdigest())
PY
        )"
  echo "  git-blob sha1(matrix/loop.py) = $got"
  if [ "$got" = "$LOOP_PY_BLOB_SHA1" ]; then
    echo "  == shipped hash (PR #56): the loop engine is unchanged across both substrates."
    echo ">>> loop.py byte-identity: PASS"
  else
    echo "  != shipped hash $LOOP_PY_BLOB_SHA1"; echo ">>> loop.py byte-identity: FAIL"; rc=1
  fi
}

export PYTHONPATH="$HERE:${PYTHONPATH:-}"

# ---- synthetic twin -------------------------------------------------------------------
run_stage "pytest suite"                    "$PROT" -m pytest "$HERE/tests" -q
run_stage "CP1 twin + harness"              "$PROT" "$HERE/verify_cp1.py"
run_stage "CP2 ruler + trust/action gates"  "$PROT" "$HERE/verify_cp2.py"
run_stage "CP3 Forge-shaped dose stage"     "$PROT" "$HERE/verify_cp3.py"
run_stage "CP4 closed loop + honest scope"  "$PROT" "$HERE/verify_cp4.py"

# ---- the engine is the same bytes on both substrates ----------------------------------
byte_identity_gate

# ---- Ferry real-data substrate --------------------------------------------------------
run_stage          "Ferry CP1 drop-in proof (byte-unchanged + contract)" "$PROT" "$HERE/verify_ferry_cp1.py"
run_stage_optional "Ferry CP2 grounded closed loop (REAL anatomy + dose)" "$PROT" "$HERE/verify_ferry_cp2.py"

echo
if [ "$rc" -eq 0 ]; then
  echo "================ reproduce.sh: all required gates GREEN ================"
else
  echo "================ reproduce.sh: a required gate FAILED (rc=$rc) ================"
fi
exit "$rc"
