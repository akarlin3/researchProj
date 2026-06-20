#!/usr/bin/env bash
# One-command re-run of the Sextant boundary-railing analysis.
#
#   bash Sextant/reproduce.sh
#
# Steps: (1) fetch OSIPI open human-abdominal data (download-on-demand, MD5
# verified, raw arrays git-ignored); (2) run the seeded railing analysis +
# bootstrap CIs on both cohorts; (3) run the test suite; (4) build the paper.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"

echo "== 1. fetch OSIPI open data (CC-BY-4.0; human-abdominal DWI + DRO) =="
"$PROT" "$HERE/scripts/fetch_osipi.py" || { echo "fetch FAILED"; exit 1; }

echo "== 2. boundary-railing analysis (seed 20260613, bootstrap CIs) =="
"$PROT" "$HERE/scripts/run_railing.py" || { echo "analysis FAILED"; exit 1; }

echo "== 3. tests =="
"$PROT" -m pytest -q "$HERE/sextant-core/tests" || { echo "tests FAILED"; exit 1; }

echo "== 4. build paper (consistency gate -> tectonic/pdflatex) =="
if [ -f "$HERE/paper/build.sh" ]; then
  bash "$HERE/paper/build.sh" || { echo "paper build FAILED"; exit 1; }
else
  echo "(paper/build.sh not present yet — CP4 deliverable; skipping)"
fi

echo "== done =="
