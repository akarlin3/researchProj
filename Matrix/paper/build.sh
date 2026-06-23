#!/usr/bin/env bash
# Build the Matrix manuscript: run the traceability gate, then compile.
#
#   1. python consistency.py   regenerate numbers.tex from seeded results + verify
#   2. tectonic matrix.tex     compile to matrix.pdf (fallback: pdflatex x2)
#
# The consistency gate MUST pass before a PDF is produced: every \num* macro in
# matrix.tex must be defined and every load-bearing assert must hold.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROT="${PROT:-/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python}"
PY="$PROT"; [ -x "$PY" ] || PY="${PYTHON:-python3}"

echo "== consistency gate (numbers trace to seeded results) =="
"$PY" "$HERE/consistency.py" || { echo "consistency gate FAILED"; exit 1; }

echo "== figures (regenerate from seeded results; non-fatal -- committed PDFs are the fallback) =="
PYTHONPATH="$(dirname "$HERE"):${PYTHONPATH:-}" "$PY" "$HERE/figures/make_figures.py" \
  || echo "  [note] figure regeneration skipped/failed; using committed figures/*.pdf"

cd "$HERE"
if command -v tectonic >/dev/null 2>&1; then
  echo "== compiling matrix.tex with tectonic =="
  tectonic matrix.tex && { echo "built matrix.pdf (tectonic)"; exit 0; }
  exit 1
elif command -v pdflatex >/dev/null 2>&1; then
  echo "== compiling matrix.tex with pdflatex (x2) =="
  pdflatex -interaction=nonstopmode -halt-on-error matrix.tex >/dev/null 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error matrix.tex && { echo "built matrix.pdf (pdflatex)"; exit 0; }
  exit 1
else
  echo "== no LaTeX engine found; consistency gate passed, skipping PDF compile =="
  exit 0
fi
