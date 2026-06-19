#!/usr/bin/env bash
# Build the Vernier manuscript. Runs the CP4 consistency gate first (regenerates
# numbers.tex from the seeded results and checks vernier.tex traceability), then
# compiles. Prefers `tectonic` (self-contained); falls back to pdflatex x2. If no
# LaTeX engine is present, the consistency gate still runs (it needs no LaTeX).
#
# Usage: bash build.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
[ -x "$PROT" ] || PROT="$(command -v python3 || command -v python)"

echo "== CP4 consistency gate =="
"$PROT" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling with tectonic =="
  tectonic vernier.tex && { echo "built vernier.pdf (tectonic)"; exit 0; }
  echo "tectonic failed"; exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling with pdflatex (x2) =="
  pdflatex -interaction=nonstopmode -halt-on-error vernier.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error vernier.tex && { echo "built vernier.pdf (pdflatex)"; exit 0; }
  echo "pdflatex failed (see vernier.log)"; exit 1
else
  echo "== no LaTeX engine found; consistency gate passed, skipping PDF compile =="
  echo "   (on Overleaf, the ebgaramond+microtype preamble builds as-is)"
  exit 0
fi
