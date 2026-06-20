"""DICOM / TCIA I/O for the independent replication cohort.

Reuses Gauge's NBIA download + DICOM helpers **read-only** (``_api``,
``_get_series_zip``, ``_read_dicoms``, ``_rescaled`` from
``Gauge/scripts/fetch_invivo.py``, loaded by path) and adds:

* a robust b-value extractor (standard DICOM tag 0018,9087, with a Siemens
  ``SequenceName`` ``*ep_b<NN>t`` fallback — TCGA-LIHC encodes b this way);
* a scheme-agnostic series assembler (Gauge's ``assemble_trace`` hard-codes the
  ACRIN-6698 b-scheme; liver DWI uses other schemes).
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np

_BVAL_TAG = (0x0018, 0x9087)     # standard DICOM MR Diffusion b-value
_SEQ_TAG = (0x0018, 0x0024)      # SequenceName (Siemens encodes b here)

_gauge = None


def _monorepo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def gauge_dicom():
    """Load Gauge's fetch_invivo helpers read-only (numpy/urllib/pydicom only)."""
    global _gauge
    if _gauge is None:
        path = _monorepo_root() / "Gauge" / "scripts" / "fetch_invivo.py"
        spec = importlib.util.spec_from_file_location("gauge_fetch_invivo", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _gauge = mod
    return _gauge


def api(endpoint, **params):
    return gauge_dicom()._api(endpoint, **params)


def get_series_zip(series_uid):
    return gauge_dicom()._get_series_zip(series_uid)


def bvalue(ds):
    """Robust b-value: standard tag first, then Siemens SequenceName ``*ep_b<NN>``."""
    if _BVAL_TAG in ds:
        try:
            v = float(ds[_BVAL_TAG].value)
            if v == v:           # not NaN
                return v
        except (TypeError, ValueError):
            pass
    seq = ds.get(_SEQ_TAG, None)
    if seq is not None:
        m = re.search(r"_b(\d+)", str(seq.value))
        if m:
            return float(m.group(1))
    return None


def _zloc(ds):
    if getattr(ds, "SliceLocation", None) is not None:
        return round(float(ds.SliceLocation), 3)
    if "ImagePositionPatient" in ds:
        return round(float(ds.ImagePositionPatient[2]), 3)
    return round(float(getattr(ds, "InstanceNumber", 0)), 3)


def assemble_series(zip_bytes):
    """Assemble a single DWI series ZIP into ((H,W,Z,nb) array, sorted bvals).

    Keeps only slice locations present at *every* b-value (a complete grid).
    """
    g = gauge_dicom()
    by = {}
    rows = cols = None
    for ds in g._read_dicoms(zip_bytes):
        b = bvalue(ds)
        if b is None:
            continue
        rows, cols = int(ds.Rows), int(ds.Columns)
        by[(b, _zloc(ds))] = g._rescaled(ds)
    if not by:
        raise ValueError("no b-value-tagged pixel data in series")
    bvals = sorted({b for (b, _) in by})
    zlocs = sorted({z for (_, z) in by})
    full = [z for z in zlocs if all((b, z) in by for b in bvals)]
    if not full:
        raise ValueError("no slice location complete across all b-values")
    vol = np.full((rows, cols, len(full), len(bvals)), np.nan, float)
    for (b, z), img in by.items():
        if z in full:
            vol[:, :, full.index(z), bvals.index(b)] = img
    return vol, np.asarray(bvals, float)
