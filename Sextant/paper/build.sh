#!/usr/bin/env bash
# Build the Sextant manuscript: CP4 consistency gate -> tectonic/pdflatex.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"

echo "== CP4 consistency gate =="
"$PROT" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling with tectonic =="
  tectonic sextant.tex && { echo "built sextant.pdf (tectonic)"; exit 0; }
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling with pdflatex (x2) =="
  pdflatex -interaction=nonstopmode -halt-on-error sextant.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error sextant.tex && { echo "built sextant.pdf (pdflatex)"; exit 0; }
else
  echo "== no LaTeX engine found; consistency gate passed =="
  exit 0
fi
