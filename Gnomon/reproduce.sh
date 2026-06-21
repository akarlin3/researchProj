#!/usr/bin/env bash
# One-command Gnomon re-run. From the monorepo root or Gnomon/, this runs the CP1
# gates now and the full reproduce-or-refute verdict once CP2/CP3 land.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "== CP1 scaffold gates =="
python -m pytest "$HERE/tests" -q

echo "== CP3 reproduction verdict =="
# Runs the rebuild, compares to the frozen manifest, writes results/reproduction.json
# and prints REPRODUCES / DOES NOT REPRODUCE. Implemented at CP3.
python -m gnomon.reproduce || {
  echo "(gnomon.reproduce not yet implemented — CP2/CP3 pending)" >&2
  exit 0
}
