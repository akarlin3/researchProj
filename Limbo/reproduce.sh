#!/usr/bin/env bash
# reproduce.sh — one-command re-validation of Limbo's contract.
#   1. the citation gate (offline): every entry resolvable + claimed, zero orphans,
#      zero phantom \cite in the manuscript prose
#   2. the test-suite (incl. quote-identity + IOP↔phiRO content-identity)
#   3. the manuscript compile (tectonic limbo_phiro.tex -> limbo_phiro.pdf; phiRO target)
# Pass --online to additionally HEAD-check that every DOI/arXiv resolves (network).
# Limbo is trigger-independent (no upstream block): a green run == submission-ready.
set -euo pipefail
cd "$(dirname "$0")"

ONLINE=""
[[ "${1:-}" == "--online" ]] && ONLINE="--online"

echo "== Limbo citation gate =="
python3 verify_citations.py $ONLINE

echo
echo "== Limbo tests =="
if command -v pytest >/dev/null 2>&1; then
  pytest -q tests/
else
  python3 -m pytest -q tests/
fi

echo
echo "== Limbo manuscript (compile: phiRO/Elsevier target) =="
if command -v tectonic >/dev/null 2>&1; then
  tectonic limbo_phiro.tex
  echo "built limbo_phiro.pdf"
else
  echo "ERROR: tectonic not found; install tectonic to compile limbo_phiro.tex" >&2
  exit 1
fi

echo
echo "Limbo OK — gate clean, tests pass, manuscript compiled (phiRO)."
