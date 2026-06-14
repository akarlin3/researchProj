#!/usr/bin/env bash
# Build the Gauge manuscript and run the GATE 3 consistency check.
# Requires `tectonic` (self-contained LaTeX; fetches packages on first run).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

echo "[1/2] GATE 3 consistency check (numbers trace to gated CP printouts)"
python consistency.py

echo "[2/2] compiling gauge.tex with tectonic"
tectonic gauge.tex
echo "built: $HERE/gauge.pdf"
