#!/usr/bin/env python
"""Download-on-demand fetch for the independent TCGA-LIHC liver-DWI replication.

TCGA-LIHC (liver hepatocellular carcinoma), TCIA, DOI
10.7937/K9/TCIA.2016.IMMQW8UQ, **CC BY 3.0** — human in-vivo abdominal (liver)
multi-b DWI. Independent of the OSIPI cohort (different site/scanner/organ),
non-pancreatic, non-MSK. Downloaded via the NBIA REST API, reusing Gauge's TCIA
helpers read-only (``sextant.dicom_io``).

Two IVIM-capable acquisition schemes are present across patients:
  * 4-b ``0/50/500/800`` (``liver_diff ax 0-50-500-800``) — has b=0 (clean).
  * 3-b ``50/400/800``  (``ep2d_diff_b50_400_800``)       — lowest b = 50.

Posture: download-on-demand; assembled arrays live under ``data/`` (git-ignored);
only the provenance manifest is committed.

Run:  python scripts/fetch_tcga_lihc.py
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "sextant-core"))

from sextant.dicom_io import api, assemble_series, get_series_zip   # noqa: E402

COLLECTION = "TCGA-LIHC"
DOI = "10.7937/K9/TCIA.2016.IMMQW8UQ"
LICENSE = "CC BY 3.0"
CITATION = ("Erickson, B. J., Kirk, S., Lee, Y., Bathe, O., Kearns, M., Gerdes, C., "
            "Rieger-Christ, K., & Lemmerman, J. (2016). The Cancer Genome Atlas Liver "
            "Hepatocellular Carcinoma Collection (TCGA-LIHC) [Data set]. "
            "The Cancer Imaging Archive. https://doi.org/10.7937/K9/TCIA.2016.IMMQW8UQ")

# IVIM-capable acquisitions, selected by exact SeriesDescription (avoids the
# derived _ADC maps and post-contrast _C variants).
SCHEMES = {
    "liver_4b_0_50_500_800": ["liver_diff ax 0-50-500-800"],
    "liver_3b_50_400_800": ["ep2d_diff_b50_400_800", "ep2d_diff_b50_400_800 FREE BREATHE"],
}

_DATA_DIR = os.path.join(_ROOT, "data", "tcga_lihc")
_PROVENANCE = os.path.join(_ROOT, "results", "tcga_lihc_provenance.json")


def _sha256_arr(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).astype(np.float32).tobytes()).hexdigest()[:16]


def _select_series():
    series = api("getSeries", Collection=COLLECTION, Modality="MR")
    chosen = {}   # scheme -> {patient: series_dict}
    for scheme, descs in SCHEMES.items():
        chosen[scheme] = {}
        for s in series:
            if s.get("SeriesDescription", "") in descs:
                pid = s["PatientID"]
                # one series per patient per scheme; prefer the larger image count
                prev = chosen[scheme].get(pid)
                if prev is None or int(s.get("ImageCount", 0)) > int(prev.get("ImageCount", 0)):
                    chosen[scheme][pid] = s
    return chosen


def fetch():
    os.makedirs(_DATA_DIR, exist_ok=True)
    chosen = _select_series()
    manifest_series = []
    for scheme, by_pt in chosen.items():
        for pid, s in by_pt.items():
            uid = s["SeriesInstanceUID"]
            print(f"[lihc] {scheme} {pid}: {s.get('SeriesDescription')} ...", flush=True)
            vol, bvals = assemble_series(get_series_zip(uid))
            out = os.path.join(_DATA_DIR, f"{scheme}__{pid}.npz")
            np.savez_compressed(out, signals4d=vol.astype(np.float32), bvals=bvals)
            rec = {
                "scheme": scheme, "patient": pid, "series_uid": uid,
                "series_description": s.get("SeriesDescription"),
                "study_date": s.get("StudyDate", ""),
                "bvals": bvals.tolist(), "shape": list(vol.shape),
                "sha256_16": _sha256_arr(vol),
                "npz": os.path.relpath(out, _ROOT),
            }
            manifest_series.append(rec)
            print(f"       shape={vol.shape} bvals={bvals.tolist()}")
    prov = {
        "dataset": {
            "name": "TCGA-LIHC liver DWI -- independent human-abdominal replication",
            "collection": COLLECTION,
            "source": "TCIA / NBIA",
            "doi": DOI,
            "license": LICENSE,
            "citation": CITATION,
            "modality": "MR DWI (real in-vivo, LIVER); human abdominal",
            "manufacturer_note": "Siemens (b-value encoded in SequenceName *ep_b<NN>t)",
            "ground_truth_ivim": False,
            "ip_posture": "OPEN/CLEAN. NOT pancreatic, NOT MSK. download-on-demand; "
                          "assembled arrays NOT committed (data/ git-ignored).",
            "nbia_api": "https://services.cancerimagingarchive.net/nbia-api/services/v1",
            "fetched_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "n_series": len(manifest_series),
            "series": manifest_series,
        }
    }
    os.makedirs(os.path.dirname(_PROVENANCE), exist_ok=True)
    with open(_PROVENANCE, "w") as fh:
        json.dump(prov, fh, indent=1)
    print(f"[lihc] wrote provenance -> {_PROVENANCE}  ({len(manifest_series)} series)")
    return _DATA_DIR


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    fetch()
