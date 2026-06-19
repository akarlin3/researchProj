#!/usr/bin/env bash
# Build the Lethe manuscript (Echo portion) with the CP4 consistency gate first.
# Mirrors Minos/future/paper/build.sh.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"

if [ ! -f "$HERE/lethe.tex" ]; then
  echo "lethe.tex not present (see README.md)."
  exit 0
fi

echo ">>> CP4 consistency gate"
"$PY" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  tectonic lethe.tex
else
  pdflatex -interaction=nonstopmode lethe.tex && pdflatex -interaction=nonstopmode lethe.tex
fi
echo ">>> built lethe.pdf"
