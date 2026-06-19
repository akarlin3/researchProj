"""Provenance manifest writer -- mirrors Gauge's download-on-demand data posture.

POSTURE (identical to Gauge): for any real repeatability dataset, NO pixel data is ever
committed. Only a provenance manifest (dataset name, source, DOI, license, attribution,
per-exam checksums, b-scheme) is written to ``results/``. The raw arrays land under
``data/`` which is git-ignored. Echo redistributes nothing.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

# Pinned public-data descriptor for the primary repeatability source.
ACRIN6698 = {
    "name": "ACRIN-6698 / I-SPY2 Breast DWI (same-day test-retest arm, TrT0/TrT1)",
    "source": "TCIA (The Cancer Imaging Archive), NBIA REST API",
    "collection_page": "https://www.cancerimagingarchive.net/collection/acrin-6698/",
    "doi": "10.7937/tcia.kk02-6d95",
    "license": "CC-BY-4.0",
    "citation": (
        "Newitt, D. C., Partridge, S. C., Zhang, Z., Gibbs, J., Chenevert, T., Rosen, M., "
        "et al. (2021). ACRIN 6698/I-SPY2 Breast DWI [Data set]. The Cancer Imaging Archive. "
        "DOI 10.7937/tcia.kk02-6d95."
    ),
    "modality": "MR DWI (real in-vivo, breast); NO ground-truth IVIM parameters",
    "b_values_s_per_mm2": [0.0, 100.0, 600.0, 800.0],
    "test_retest_arm": "same-day coffee-break repeats (TrT0/TrT1); ~76 tumor pairs",
    "posture": "download-on-demand; pixel data NOT committed (data/ git-ignored); "
               "provenance manifest only",
}


def sha16(arr_bytes: bytes) -> str:
    """16-char truncated SHA256 (matches Gauge's ``sha_signals`` convention)."""
    return hashlib.sha256(arr_bytes).hexdigest()[:16]


def write_manifest(path: str | Path, dataset: dict, exams: list[dict],
                   extra: dict | None = None) -> Path:
    """Write a committed provenance manifest. ``exams`` carry checksums, not pixels."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "echo_provenance_version": 1,
        "dataset": dataset,
        "exams": exams,
        "n_exams": len(exams),
    }
    if extra:
        payload.update(extra)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False))
    return path
