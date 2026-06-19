#!/usr/bin/env bash
# One-command re-validation for Echo. Mirrors Minos/future/reproduce.sh.
#
#   bash reproduce.sh          # default
#   FULL=1 bash reproduce.sh   # full-N harness self-test
#
# CP1 method self-test always runs (SOLID, synthetic). CP2 data check + CP3 validation run
# only if real test-retest data is present (download-on-demand; not committed). CP4
# consistency runs only if the manuscript exists. Stops at the first hard failure.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"
N_SELFTEST=$([ "${FULL:-0}" = "1" ] && echo 50000 || echo 20000)

run_stage() {
  local name="$1"; shift
  echo ">>> ${name}: running"
  if "$@"; then echo ">>> ${name}: PASS"; else
    local rc=$?
    if [ "$rc" = "3" ]; then echo ">>> ${name}: SKIP (data gate not satisfied)"; return 0; fi
    echo ">>> ${name}: FAIL (rc=$rc)"; exit "$rc"
  fi
}

echo "=== Echo re-validation (CP1 -> CP4) ==="
run_stage "CP1 method self-test (SOLID)" "$PY" "$HERE/scripts/run_harness.py" --n "$N_SELFTEST"

echo ">>> CP2 data check"
"$PY" "$HERE/scripts/fetch_invivo.py" --check || true

# CP3 renders a VERDICT (PASS or LETHE) -- both are valid outcomes, so exit 0 = rendered,
# 3 = data gate unsatisfied (-> use Reverb/run_harness), 1 = error. Report the verdict.
echo ">>> CP3 real-data validation (PROVISIONAL): running"
"$PY" "$HERE/scripts/run_validation.py"; cp3=$?
if [ "$cp3" = "3" ]; then
  echo ">>> CP3: SKIP (data gate not satisfied -> Reverb: scripts/run_harness.py)"
elif [ "$cp3" = "0" ]; then
  V=$("$PY" -c "import json,sys;print(json.load(open('$HERE/results/RESULTS_VALIDATION.json'))['gate']['VERDICT'])" 2>/dev/null || echo "?")
  echo ">>> CP3: VERDICT RENDERED = $V"
else
  echo ">>> CP3: ERROR (rc=$cp3)"; exit "$cp3"
fi

if [ -f "$HERE/paper/consistency.py" ] && [ -f "$HERE/paper/lethe.tex" ]; then
  run_stage "CP4 manuscript consistency (Lethe / Echo portion)" "$PY" "$HERE/paper/consistency.py"
else
  echo ">>> CP4 manuscript consistency: SKIP (manuscript absent)"
fi

echo "=== done ==="
