#!/usr/bin/env bash
# build.sh — compile the Limbo manuscript, gated on the citation gate.
#   Mirrors the sibling pattern (cf. Gauge build.sh): run the hard gate first,
#   then compile with tectonic. The gate's non-zero exit aborts the build
#   (set -e), so a phantom \cite or unresolved identifier never reaches a PDF.
#   Target: limbo_phiro.tex (Elsevier elsarticle; Physics and Imaging in Radiation
#   Oncology). limbo.tex (IOP/PMB) is retained as the content-identity reference; pass
#   --iop to build it instead.
# Pass --online to also HEAD-check every identifier resolves live.
set -euo pipefail
cd "$(dirname "$0")"

ONLINE=""
SRC="limbo_phiro"
for a in "$@"; do
  [[ "$a" == "--online" ]] && ONLINE="--online"
  [[ "$a" == "--iop" ]] && SRC="limbo"
done

echo "== Limbo citation gate (build precondition) =="
python3 verify_citations.py $ONLINE

echo
echo "== compiling ${SRC}.tex with tectonic =="
if command -v tectonic >/dev/null 2>&1; then
  tectonic "${SRC}.tex"
else
  echo "ERROR: tectonic not found; install tectonic to build ${SRC}.pdf" >&2
  exit 1
fi
echo "built: $(pwd)/${SRC}.pdf"
