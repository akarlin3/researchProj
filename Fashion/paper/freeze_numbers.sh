#!/usr/bin/env bash
# Extract every reported number from the rendered manuscript + supplement PDFs
# into a canonical, sorted multiset. Used by the assessor-remediation final gate to
# prove that prose edits change no numerical result.
#
# Usage:
#   ./freeze_numbers.sh > NUMBERS_FROZEN.txt   # freeze baseline
#   ./freeze_numbers.sh > NUMBERS_NOW.txt      # snapshot after edits
#   diff NUMBERS_FROZEN.txt NUMBERS_NOW.txt    # must be empty
#
# Numbers are extracted from the BUILT PDFs (manuscript.pdf, supplement.pdf), not the
# .tex, so LaTeX-only tokens (package versions, macro args) never pollute the set.
set -euo pipefail
cd "$(dirname "$0")"

extract() {
  # -layout keeps tables aligned; we only care about the numeric tokens.
  pdftotext -layout "$1" - 2>/dev/null \
    | grep -oE '[0-9]+([.,][0-9]+)*([eE][-+]?[0-9]+)?' \
    | sort | uniq -c \
    | sed 's/^ *//'
}

echo "# manuscript.pdf"
extract manuscript.pdf
echo "# supplement.pdf"
extract supplement.pdf
