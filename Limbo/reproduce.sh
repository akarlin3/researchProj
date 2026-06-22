#!/usr/bin/env bash
# reproduce.sh — one-command re-validation of Limbo's CP1 contract.
#   1. the citation gate (offline): every entry resolvable + claimed, zero orphans
#   2. the test-suite
# Pass --online to additionally HEAD-check that every DOI/arXiv resolves (network; CP3).
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
echo "Limbo CP1 OK."
