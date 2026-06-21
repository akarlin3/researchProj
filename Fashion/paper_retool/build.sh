#!/usr/bin/env bash
# Build the retooled (boundary-railing-first) Fashion manuscript:
# numbers-gate (integration consistency) -> tectonic/pdflatex.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"

echo "== numbers-gate (integration: Sextant + Gnomon -> numbers.tex) =="
"$PY" "$HERE/consistency.py" || { echo "numbers-gate FAILED"; exit 1; }

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling with tectonic =="
  tectonic manuscript.tex && { echo "built manuscript.pdf (tectonic)"; exit 0; }
  echo "tectonic build FAILED"; exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling with pdflatex + bibtex =="
  pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex >/dev/null 2>&1
  bibtex manuscript >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex && { echo "built manuscript.pdf (pdflatex)"; exit 0; }
  echo "pdflatex build FAILED"; exit 1
else
  echo "== no LaTeX engine found; numbers-gate passed, PDF not built =="
  exit 0
fi
