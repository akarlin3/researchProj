"""On-demand fetch of the OSIPI DRO, reusing Gauge's pinned provenance.

This reuses (never reinvents) Gauge's fetcher: it imports the pinned constants
(``ZIP_URL``, ``ZIP_MD5``, ``DRO_MEMBER``, ``ZENODO_DOI``, ``REPO_COMMIT``) and the
``verify_and_summarize`` checker from ``Gauge/scripts/fetch_osipi.py`` -- Gauge is
read-only and is never modified. The 245 MB phantom ZIP is downloaded and
md5-verified, ``DRO.npy`` is extracted into the git-ignored ``Datum/data/osipi/``
cache, and a provenance manifest is written under ``Datum/results/`` (the only
tracked artifact). The DRO itself (raw arrays) is never committed.

    python -m datum.osipi_fetch          # fetch if absent, then verify + provenance
    python -m datum.osipi_fetch --force  # re-download
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

from datum import _paths

_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = _ROOT / "data" / "osipi"
DRO_PATH = DATA_DIR / "DRO.npy"
PROV_PATH = _ROOT / "results" / "osipi_provenance.json"


def _gauge_fetch():
    """Import Gauge's fetch_osipi module (read-only) to reuse its pinned constants."""
    root = _paths.find_monorepo_root()
    if root is None:
        raise RuntimeError("monorepo root not found; cannot reuse Gauge's OSIPI pins")
    scripts = str(root / "Gauge" / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    import fetch_osipi  # noqa: E402  -- Gauge's script, imported read-only
    return fetch_osipi


def _md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def fetch(force: bool = False) -> Path:
    g = _gauge_fetch()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if force or not DRO_PATH.exists():
        zip_path = DATA_DIR / g.ZIP_NAME
        if force or not zip_path.exists():
            print(f"[osipi] downloading {g.ZIP_URL} -> {zip_path}")
            urllib.request.urlretrieve(g.ZIP_URL, zip_path)
        got = _md5(zip_path)
        if got != g.ZIP_MD5:
            raise RuntimeError(f"OSIPI zip md5 mismatch: {got} != {g.ZIP_MD5}")
        print(f"[osipi] zip md5 OK ({got})")
        with zipfile.ZipFile(zip_path) as zf, zf.open(g.DRO_MEMBER) as src, \
                open(DRO_PATH, "wb") as out:
            shutil.copyfileobj(src, out)
        print(f"[osipi] extracted {g.DRO_MEMBER} -> {DRO_PATH}")

    summary = g.verify_and_summarize(str(DRO_PATH))   # reuse Gauge's verifier
    PROV_PATH.parent.mkdir(parents=True, exist_ok=True)
    prov = {
        "artifact": "OSIPI TF2.4 IVIM digital reference object (DRO) -- synthetic",
        "reused_from": "Gauge/scripts/fetch_osipi.py (pinned constants + verifier)",
        "zenodo_doi": g.ZENODO_DOI,
        "zip_url": g.ZIP_URL,
        "zip_md5": g.ZIP_MD5,
        "dro_member": g.DRO_MEMBER,
        "repo_commit": g.REPO_COMMIT,
        "license_data": "CC-BY-4.0",
        "license_code": "Apache-2.0",
        "posture": "download-on-demand; raw DRO NOT committed (data/ git-ignored)",
        "summary": summary,
    }
    PROV_PATH.write_text(json.dumps(prov, indent=2) + "\n")
    print(f"[osipi] provenance -> {PROV_PATH}")
    return DRO_PATH


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Fetch the OSIPI DRO (reuses Gauge pins).")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    args = ap.parse_args(argv)
    fetch(force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
