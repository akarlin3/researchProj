#!/usr/bin/env bash
# Augur one-command REPRODUCTION (green against the in-repository provisional anchors).
#
# Augur makes no new measurement; "reproduce" therefore means: regenerate every in-repository
# anchor deterministically, regenerate the manuscript's numbers.tex from them, and run the test
# suite. This is intentionally SEPARATE from release: the submission HOLD lives in the release
# gate (release_gate.py / submit.sh), NOT here. reproduce.sh is green whether or not the paper is
# clear to submit.
#
# Pipeline:  extract anchors -> CRLB wall (+figure) -> D* retest CI -> D*-Ktrans -> consistency
#            (numbers.tex) -> pytest.   Exit 0 = every anchor regenerated and self-consistent.
#
# Uses the proteus env (system python3 lacks numpy); override with PROT=... if needed.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
if [ ! -x "$PROT" ]; then PROT="$(command -v python3)"; fi

rc=0
step() {  # label  script  [args...]
  local label="$1"; shift
  echo; echo "## ${label}"
  if "$PROT" "$@"; then echo ">>> ${label}: OK"; else echo ">>> ${label}: FAIL"; rc=1; fi
}

echo "### Augur reproduction (green against in-repository anchors)"

step "1. Extract in-repo anchors (Gauge results + CITATIONS -> anchors.json)" \
     "$HERE/anchors/extract_anchors.py"
step "2. CRLB identifiability wall (re-derived; figure + anchor cross-check)" \
     "$HERE/scripts/crlb_wall.py"
step "3. D* test-retest correlation interval (bootstrap CI carried + Fisher-z reproduced)" \
     "$HERE/scripts/retest_ci.py"
step "4. D*-Ktrans cross-modal evidence (external literature, verified)" \
     "$HERE/scripts/dstar_ktrans.py"
step "5. CP4 consistency gate (regenerate numbers.tex; macros + spine invariants)" \
     "$HERE/paper/consistency.py"

echo; echo "## 6. Tests (reproduction artifacts + spine invariants + hold engaged)"
if "$PROT" -m pytest "$HERE/tests" -q; then echo ">>> tests: PASS"; else echo ">>> tests: FAIL"; rc=1; fi

echo
if [ "$rc" -eq 0 ]; then
  echo "================ reproduce.sh: GREEN -- every in-repo anchor regenerated ================"
  echo ">>> Reproduction is independent of release. Submission status is reported separately by"
  echo "    release_gate.py / submit.sh (currently HELD until Fashion + Minos publish)."
else
  echo "================ reproduce.sh: a step FAILED (rc=$rc) ================"
fi
exit "$rc"
