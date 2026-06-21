#!/usr/bin/env bash
# Final-gate numeric check for the assessor-remediation run.
#
# Re-extracts numbers from the current built PDFs and compares against the frozen
# baseline (NUMBERS_FROZEN.txt). The contract:
#
#   * No baseline number may DECREASE in count  -> a reported result was altered or
#     dropped. This FAILS the gate (exit 1).
#   * Numbers may INCREASE or be ADDED -> these come from content legitimately added by
#     the remediation (the repository DOI, dataset accession, declarations block, the
#     reference-style conversion). They are listed for manual confirmation that each is
#     an identifier, not a research result.
#
# Usage: ./check_numbers.sh   (run after `make` so PDFs are current)
set -euo pipefail
cd "$(dirname "$0")"

FROZEN="../NUMBERS_FROZEN.txt"   # canonical baseline lives at the Fashion/ root
NOW="../NUMBERS_NOW.txt"

./freeze_numbers.sh > "$NOW"

awk '
  # read baseline counts keyed by "<section>\t<number>"
  FNR==NR {
    if ($0 ~ /^#/) { sec=$0; next }
    cnt=$1; sub(/^[0-9]+ /,"",$0); base[sec"\t"$0]=cnt; next
  }
  # read current counts
  {
    if ($0 ~ /^#/) { sec=$0; next }
    cnt=$1; sub(/^[0-9]+ /,"",$0); now[sec"\t"$0]=cnt
  }
  END {
    fail=0
    for (k in base) {
      b=base[k]; n=(k in now)?now[k]:0
      if (n < b) { printf "DECREASED  %-14s  baseline=%d now=%d\n", k, b, n; fail=1 }
    }
    print "----- additions (manual review: must be identifiers, not results) -----"
    for (k in now) {
      b=(k in base)?base[k]:0; n=now[k]
      if (n > b) printf "increased  %-14s  baseline=%d now=%d\n", k, b, n
    }
    if (fail) { print "\nRESULT: FAIL — a reported number changed or was dropped."; exit 1 }
    else      { print "\nRESULT: PASS — no reported number decreased." }
  }
' "$FROZEN" "$NOW"
