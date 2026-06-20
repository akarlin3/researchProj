#!/usr/bin/env python
"""Download-on-demand fetcher for the OSIPI TF2.4 open data (human-abdominal DWI).

Mirrors the Gauge fetch/provenance pattern (``Gauge/scripts/fetch_osipi.py``):
  * download the OSIPI Zenodo data archive, verify its MD5;
  * extract ONLY the members Sextant needs (human-abdominal scan + masks + the
    synthetic DRO used by the scoped ruler section);
  * commit a JSON provenance manifest (URLs, DOI, license, checksums, b-scheme,
    voxel counts, access date) — never the raw imaging arrays (``data/`` is
    git-ignored).

Primary cohort for boundary-railing is the OSIPI **human-abdominal** acquisition
(Philips 3T, sourced from the IVIMNET repository), distributed in the same
CC-BY-4.0 Zenodo record. This is the open data on which Fashion's 54.7% railing
figure was originally computed; Sextant re-centres it as the primary claim and
extends it to the full-abdomen ROI.

Run:  python scripts/fetch_osipi.py            (uses a cached zip if present)
      python scripts/fetch_osipi.py --force    (re-download)
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import urllib.request
import zipfile

import numpy as np

# --------------------------------------------------------------------------- #
# Pinned, immutable source: OSIPI code repo tag v0.1.0 + Zenodo data record.
# --------------------------------------------------------------------------- #
REPO = "https://github.com/OSIPI/TF2.4_IVIM-MRI_CodeCollection"
REPO_TAG = "v0.1.0"
ZENODO_RECORD = "14605039"
ZENODO_DOI = "10.5281/zenodo.14605039"
ZIP_URL = (f"https://zenodo.org/records/{ZENODO_RECORD}/files/"
           "OSIPI_TF24_data_phantoms.zip?download=1")
ZIP_NAME = "OSIPI_TF24_data_phantoms.zip"
ZIP_MD5 = "e7b3fe1d811a7a45c5aaf6c604c82793"          # from Zenodo file API

# Human-abdominal members (the boundary-railing cohorts) + the synthetic DRO.
ABDOMEN_MEMBERS = [
    "Data/abdomen.nii.gz",
    "Data/abdomen.bval",
    "Data/abdomen_readme.txt",
    "Data/mask_abdomen_homogeneous.nii.gz",   # curated ROI (original cohort)
    "Data/mask_full_temp.nii.gz",             # full-abdomen ROI (generalisation)
]
DRO_MEMBER = "Utilities/DRO.npy"               # synthetic, for the scoped ruler

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DATA_DIR = os.path.join(_ROOT, "data", "osipi")
_EXTRACT_DIR = os.path.join(_DATA_DIR, "extracted")
_RESULTS_DIR = os.path.join(_ROOT, "results")
_PROVENANCE = os.path.join(_RESULTS_DIR, "osipi_provenance.json")


def _md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def _sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()[:16]


def _find_cached_zip():
    """Return a path to an already-downloaded, MD5-correct zip, if any."""
    candidates = [os.path.join(_DATA_DIR, ZIP_NAME)]
    job_tmp = os.environ.get("CLAUDE_JOB_DIR")
    if job_tmp:
        candidates.append(os.path.join(job_tmp, "tmp", ZIP_NAME))
    for c in candidates:
        if os.path.exists(c) and _md5(c) == ZIP_MD5:
            return c
    return None


def download_zip(force=False):
    os.makedirs(_DATA_DIR, exist_ok=True)
    zip_path = os.path.join(_DATA_DIR, ZIP_NAME)
    if not force:
        cached = _find_cached_zip()
        if cached:
            print(f"[osipi] using cached zip {cached}")
            return cached
    print(f"[osipi] downloading {ZIP_URL} -> {zip_path}")
    urllib.request.urlretrieve(ZIP_URL, zip_path)
    got = _md5(zip_path)
    if got != ZIP_MD5:
        raise RuntimeError(
            f"OSIPI zip md5 mismatch: got {got}, expected {ZIP_MD5}. "
            "Delete data/osipi and re-run.")
    print(f"[osipi] zip md5 OK ({got})")
    return zip_path


def extract(zip_path):
    os.makedirs(_EXTRACT_DIR, exist_ok=True)
    members = ABDOMEN_MEMBERS + [DRO_MEMBER]
    out = {}
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        for m in members:
            if m not in names:
                raise RuntimeError(f"member {m} not in archive")
            dest = os.path.join(_EXTRACT_DIR, os.path.basename(m))
            with zf.open(m) as src, open(dest, "wb") as dst:
                dst.write(src.read())
            out[m] = dest
            print(f"[osipi] extracted {m} -> {dest}")
    return out


def summarize(paths):
    """Sanity-check shapes/units and return a provenance summary (no raw arrays)."""
    import nibabel as nib  # local import; nibabel is a runtime dependency
    abd = paths["Data/abdomen.nii.gz"]
    bval = paths["Data/abdomen.bval"]
    img = nib.load(abd).get_fdata()
    bv = np.loadtxt(bval)
    assert img.ndim == 4, f"abdomen.nii.gz not 4D: {img.shape}"
    assert img.shape[3] == bv.size, f"volumes {img.shape[3]} != bvals {bv.size}"
    uniq, counts = np.unique(bv, return_counts=True)
    masks = {}
    for key in ("Data/mask_abdomen_homogeneous.nii.gz", "Data/mask_full_temp.nii.gz"):
        m = nib.load(paths[key]).get_fdata()
        masks[os.path.basename(key)] = int((m > 0).sum())
    files = {os.path.basename(m): {"sha256_16": _sha256_file(p)} for m, p in paths.items()}
    return {
        "abdomen_shape": list(img.shape),
        "n_diffusion_volumes": int(img.shape[3]),
        "b_values_unique_s_per_mm2": uniq.tolist(),
        "b_value_repeats": {str(int(b)): int(c) for b, c in zip(uniq, counts)},
        "roi_voxel_counts": masks,
        "files": files,
    }


def write_provenance(summary):
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    prov = {
        "dataset": {
            "name": "OSIPI TF2.4 open data -- human-abdominal IVIM/DWI + synthetic DRO",
            "primary_cohort": "in-vivo human ABDOMEN, Philips 3T multi-b DWI "
                              "(sourced from the IVIMNET repository), N=1 subject",
            "kind": "OPEN community-standard data; the human-abdominal scan is the "
                    "cohort on which Fashion's 54.7% NLLS railing was first computed",
            "code_repo": REPO,
            "code_repo_tag": REPO_TAG,
            "code_license": "Apache-2.0",
            "data_record_doi": ZENODO_DOI,
            "data_record_url": f"https://zenodo.org/records/{ZENODO_RECORD}",
            "data_license": "CC-BY-4.0",
            "abdomen_source_note": "Data/abdomen_readme.txt: 'Scanned on a Philips 3T "
                                   "system. Taken from https://github.com/oliverchampion/IVIMNET'",
            "zip_name": ZIP_NAME,
            "zip_md5": ZIP_MD5,
            "members_extracted": ABDOMEN_MEMBERS + [DRO_MEMBER],
            "forward_model": "bi-exponential S/S0 = (1-f)exp(-bD) + f exp(-bD*)",
            "citation_data": (
                "Gurney-Champion O, Rashid I, van der Thiel M, Kuppens D, Voorter P, "
                "van Houdt P, Peterson E, Jalnefjord O. Data to "
                "github.com/OSIPI/TF2.4_IVIM-MRI_CodeCollection. Zenodo; "
                f"doi:{ZENODO_DOI}."),
            "ip_posture": "OPEN/CLEAN. NOT pancData3, NOT MSK clinical data. "
                          "download-on-demand; raw imaging NOT committed (data/ git-ignored).",
            "fetched_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            **summary,
        }
    }
    if os.path.exists(_PROVENANCE):
        try:
            old = json.load(open(_PROVENANCE))
            old.update(prov)
            prov = old
        except Exception:
            pass
    with open(_PROVENANCE, "w") as fh:
        json.dump(prov, fh, indent=1)
    print(f"[osipi] wrote provenance -> {_PROVENANCE}")
    return _PROVENANCE


def fetch(force=False):
    zip_path = download_zip(force=force)
    paths = extract(zip_path)
    summary = summarize(paths)
    write_provenance(summary)
    return _EXTRACT_DIR


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    args = ap.parse_args()
    fetch(force=args.force)
