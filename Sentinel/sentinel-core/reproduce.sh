#!/usr/bin/env bash
# One-command re-validation of projSentinel.
#   GATE B/C : pytest (enabler runs, Matrix byte-unchanged, baselines sanity-pass)
#   GATE D   : the CP0 separation experiment -> prints the VERDICT (RED).
# Matrix is imported READ-ONLY; the run asserts loop.py byte-identity and aborts if
# Matrix changed. Point $SENTINEL_MATRIX_PATH at the Matrix package if relocated.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
export SENTINEL_MATRIX_PATH="${SENTINEL_MATRIX_PATH:-/Users/averykarlin/researchProj/.claude/worktrees/matrix-subrepo/Matrix}"
export PYTHONPATH="$HERE:${PYTHONPATH:-}"
rc=0

echo "######## Matrix byte-identity anchor ########"
"$PROT" - <<'PY' || rc=1
from sentinel.matrix_bridge import assert_matrix_untouched, LOOP_PY_SHA256
print("loop.py sha256:", assert_matrix_untouched(), "(matches anchor:", LOOP_PY_SHA256[:12]+"...)")
PY

echo; echo "######## GATE B/C: tests ########"
"$PROT" -m pytest -q "$HERE/tests" || rc=1

echo; echo "######## GATE D: CP0 separation verdict ########"
"$PROT" "$HERE/experiments/run_cp0.py" || rc=1

echo; echo "EXIT: $rc"
exit "$rc"
