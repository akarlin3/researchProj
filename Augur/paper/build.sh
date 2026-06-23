#!/usr/bin/env bash
# Build the Augur manuscript. Regenerates the reproduced, anchor-traced results, runs the CP4
# consistency gate (which rewrites numbers.tex and checks macro coverage + spine invariants),
# then compiles with tectonic. Mirrors Gauge/Minos/Lethe paper/build.sh.
#
# Uses the proteus env (system python3 lacks numpy); override with PROT=... if needed.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUGUR="$(cd "$HERE/.." && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
if [ ! -x "$PROT" ]; then PROT="$(command -v python3)"; fi

echo "== regenerate reproduced anchors (extract -> crlb -> retest -> dstar) =="
"$PROT" "$AUGUR/anchors/extract_anchors.py"   || { echo "extract FAILED"; exit 1; }
"$PROT" "$AUGUR/scripts/crlb_wall.py"         || { echo "crlb_wall FAILED"; exit 1; }
"$PROT" "$AUGUR/scripts/retest_ci.py"         || { echo "retest_ci FAILED"; exit 1; }
"$PROT" "$AUGUR/scripts/dstar_ktrans.py"      || { echo "dstar_ktrans FAILED"; exit 1; }

echo "== CP4 consistency gate (regenerates numbers.tex; checks macros + invariants) =="
"$PROT" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling augur.tex with tectonic =="
  tectonic augur.tex && { echo "built: $HERE/augur.pdf"; exit 0; }
  echo "tectonic failed"; exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling with pdflatex (x2 + bibtex) =="
  pdflatex -interaction=nonstopmode -halt-on-error augur.tex >/dev/null 2>&1
  bibtex augur >/dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode -halt-on-error augur.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error augur.tex && { echo "built augur.pdf"; exit 0; }
  echo "pdflatex failed (see augur.log)"; exit 1
else
  echo "== no LaTeX engine found; consistency gate passed, skipping PDF compile =="
  exit 0
fi
