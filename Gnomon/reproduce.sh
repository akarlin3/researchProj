#!/usr/bin/env bash
# One-command Gnomon re-run. From the monorepo root or Gnomon/, this runs the CP1
# gates and the full reproduce-or-refute verdict (CP2/CP3, see VERDICT.md).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "== CP1 scaffold gates =="
python -m pytest "$HERE/tests" -q

echo "== CP3 reproduction verdict =="
# Runs the rebuild, compares to the frozen manifest, writes results/reproduction.json
# and prints REPRODUCES / DOES NOT REPRODUCE.
python -m gnomon.reproduce
