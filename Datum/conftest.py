"""Make ``datum`` importable when running ``pytest Datum/tests`` from anywhere.

Datum is normally used via ``pip install -e Datum``; this conftest also lets the
CP1 gates run in-place inside the monorepo without an install step by putting the
Datum root (which contains the ``datum`` package) on ``sys.path``. The read-only
Caliper/Gauge dependencies are resolved separately by ``datum._paths.ensure_deps``.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
