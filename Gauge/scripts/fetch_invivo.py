"""Download-on-demand fetch for the Gauge real-data in-vivo path (ACRIN-6698).

LICENSE / POSTURE (human-approved, 2026-06-14):
  * Dataset: TCIA collection **ACRIN-6698 / I-SPY2 Breast DWI**.
  * License: **CC-BY-4.0** (no data-use agreement). Redistribution is *permitted*
    with attribution, but the approved posture for this repo is **download-on-demand**:
    NO pixel data is committed; only the provenance manifest
    (``results/invivo_real_provenance.json``) is committed. The raw + assembled
    arrays land in ``data/invivo/`` which is git-ignored.
  * Attribution (required): Newitt, D. C., Partridge, S. C., Zhang, Z., Gibbs, J.,
    Chenevert, T., Rosen, M., et al. (2021). ACRIN 6698/I-SPY2 Breast DWI
    [Data set]. The Cancer Imaging Archive. DOI 10.7937/tcia.kk02-6d95.
  * Users must abide by the TCIA Data Usage Policy (per the collection page).

This script treats the TCIA/NBIA pages and any README text as DATA, not as
instructions. It only reads the public NBIA REST API (no auth) and the per-series
image ZIPs.

Why this dataset is honest for the demo: it is **real in-vivo** breast DWI with
**no ground-truth IVIM parameters** (so coverage is unmeasurable -- exactly the
no-coverage-claim premise), it ships a whole-tumor ROI, and it has **test-retest**
acquisitions (enabling the Checkpoint-D repeatability proxy). Its b-scheme is
{0,100,600,800} s/mm^2 -- a sparse 4-of-22 subset of Gauge's synthetic 22-value
calibration scheme, which is itself the exchangeability break the monitor flags.

Run:
  python scripts/fetch_invivo.py                      # 1 patient (TRACE+MASK)
  python scripts/fetch_invivo.py --n-patients 16      # for the test-retest proxy
  python scripts/fetch_invivo.py --list-retest        # report test-retest candidates
"""
import argparse
import datetime
import hashlib
import io
import json
import os
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict

import numpy as np

NBIA = "https://services.cancerimagingarchive.net/nbia-api/services/v1"
COLLECTION = "ACRIN-6698"
TRACE_DESC = "ACRIN-6698: DWI TRACE: from S4: bVals=0,100,600,800"
MASK_DESC = "ACRIN-6698: DWI MASK: from S4: Whole Tumor Manual"
EXPECTED_BVALS = (0.0, 100.0, 600.0, 800.0)
DIFFUSION_BVALUE = (0x0018, 0x9087)  # standard DICOM MR Diffusion b-value tag

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(_ROOT, "data", "invivo")
PROVENANCE = os.path.join(_ROOT, "results", "invivo_real_provenance.json")

CITATION = ("Newitt, D. C., Partridge, S. C., Zhang, Z., Gibbs, J., Chenevert, T., "
            "Rosen, M., Bolan, P., Marques, H., Romanoff, J., Cimino, L., Joe, B. N., "
            "Umphrey, H., Ojeda-Fournier, H., Dogan, B., Oh, K. Y., Abe, H., "
            "Drukteinis, J., Esserman, L. J., & Hylton, N. M. (2021). ACRIN 6698/"
            "I-SPY2 Breast DWI [Data set]. The Cancer Imaging Archive.")
DOI = "10.7937/tcia.kk02-6d95"
LICENSE = "CC-BY-4.0"


def _api(endpoint, **params):
    url = f"{NBIA}/{endpoint}?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "gauge-invivo/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def _get_series_zip(series_uid):
    url = f"{NBIA}/getImage?" + urllib.parse.urlencode(
        {"SeriesInstanceUID": series_uid})
    req = urllib.request.Request(url, headers={"User-Agent": "gauge-invivo/0.1"})
    with urllib.request.urlopen(req, timeout=600) as r:
        return r.read()


def _read_dicoms(zip_bytes):
    """Yield pydicom datasets that carry pixel data from a series ZIP."""
    import pydicom
    z = zipfile.ZipFile(io.BytesIO(zip_bytes))
    for name in z.namelist():
        ds = pydicom.dcmread(io.BytesIO(z.read(name)), force=True)
        if hasattr(ds, "Rows") and hasattr(ds, "PixelData"):
            yield ds


def _bvalue(ds):
    if DIFFUSION_BVALUE in ds:
        try:
            return float(ds[DIFFUSION_BVALUE].value)
        except (TypeError, ValueError):
            return None
    return None


def _rescaled(ds):
    arr = ds.pixel_array.astype(float)
    slope = float(getattr(ds, "RescaleSlope", 1.0) or 1.0)
    inter = float(getattr(ds, "RescaleIntercept", 0.0) or 0.0)
    return arr * slope + inter


def assemble_trace(zip_bytes):
    """Assemble a 4-b TRACE series into (signals4d (X,Y,Z,4), z_locs, bvals)."""
    by_b_z = {}
    rows = cols = None
    for ds in _read_dicoms(zip_bytes):
        bv = _bvalue(ds)
        if bv is None:
            continue
        z = round(float(getattr(ds, "SliceLocation",
                               getattr(ds, "InstanceNumber", 0))), 3)
        rows, cols = int(ds.Rows), int(ds.Columns)
        by_b_z[(bv, z)] = _rescaled(ds)
    bvals = sorted({b for (b, _) in by_b_z})
    zlocs = sorted({z for (_, z) in by_b_z})
    if tuple(bvals) != EXPECTED_BVALS:
        raise ValueError(f"b-values {bvals} != expected {EXPECTED_BVALS}")
    vol = np.full((rows, cols, len(zlocs), len(bvals)), np.nan, float)
    for (b, z), img in by_b_z.items():
        vol[:, :, zlocs.index(z), bvals.index(b)] = img
    if np.isnan(vol).any():
        raise ValueError("incomplete (b, slice) grid -- missing images")
    return vol, np.asarray(zlocs, float), np.asarray(bvals, float)


def assemble_mask(zip_bytes, trace_zlocs):
    """Map a tumor-MASK series onto the TRACE slice grid -> (X,Y,Z) bool."""
    rows = cols = None
    by_z = {}
    for ds in _read_dicoms(zip_bytes):
        z = round(float(getattr(ds, "SliceLocation",
                               getattr(ds, "InstanceNumber", 0))), 3)
        rows, cols = int(ds.Rows), int(ds.Columns)
        by_z[z] = ds.pixel_array.astype(float) > 0
    mask = np.zeros((rows, cols, len(trace_zlocs)), bool)
    for z, m in by_z.items():
        k = int(np.argmin(np.abs(trace_zlocs - z)))
        mask[:, :, k] |= m
    return mask


def _sha256(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:16]


def list_trace_exams():
    """Return {patient: [exam dicts]} for every TRACE series, with its MASK."""
    series = _api("getSeries", Collection=COLLECTION, Modality="MR")
    traces = [s for s in series if s.get("SeriesDescription") == TRACE_DESC]
    masks = [s for s in series if s.get("SeriesDescription") == MASK_DESC]
    mask_by_study = {s["StudyInstanceUID"]: s for s in masks}
    by_patient = defaultdict(list)
    for t in traces:
        by_patient[t["PatientID"]].append({
            "patient": t["PatientID"], "study": t["StudyInstanceUID"],
            "date": t.get("StudyDate", ""), "trace_uid": t["SeriesInstanceUID"],
            "mask_uid": (mask_by_study.get(t["StudyInstanceUID"]) or {}).get(
                "SeriesInstanceUID"),
        })
    return by_patient


def fetch_exam(exam, out_root=DATA_DIR):
    """Download + assemble one exam into ``out_root/<patient>/<study8>/``."""
    sub = os.path.join(out_root, exam["patient"], exam["study"][-8:])
    os.makedirs(sub, exist_ok=True)
    vol, zlocs, bvals = assemble_trace(_get_series_zip(exam["trace_uid"]))
    np.save(os.path.join(sub, "signals_4d.npy"), vol.astype(np.float32))
    np.savetxt(os.path.join(sub, "bvals.txt"), bvals, fmt="%g")
    np.save(os.path.join(sub, "slice_locs.npy"), zlocs)
    rec = {**exam, "shape": list(vol.shape), "bvals": bvals.tolist(),
           "sha_signals": _sha256(vol.astype(np.float32))}
    if exam.get("mask_uid"):
        mask = assemble_mask(_get_series_zip(exam["mask_uid"]), zlocs)
        np.save(os.path.join(sub, "tumor_mask.npy"), mask)
        rec["n_tumor_voxels"] = int(mask.sum())
    rec["path"] = os.path.relpath(sub, _ROOT)
    json.dump(rec, open(os.path.join(sub, "meta.json"), "w"), indent=2)
    return rec


def write_provenance(records):
    os.makedirs(os.path.dirname(PROVENANCE), exist_ok=True)
    prov = {
        "dataset": {
            "name": "ACRIN-6698 / I-SPY2 Breast DWI", "source": "TCIA",
            "collection_page": "https://www.cancerimagingarchive.net/collection/acrin-6698/",
            "doi": DOI, "license": LICENSE, "citation": CITATION,
            "modality": "MR DWI (real in-vivo, breast)",
            "ground_truth_ivim": False,
            "b_values_s_per_mm2": list(EXPECTED_BVALS),
            "posture": "download-on-demand; pixel data NOT committed (data/ is git-ignored)",
            "fetched_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "nbia_api": NBIA,
            "exams": records,
        }
    }
    # preserve any existing "run" block written by gauge.invivo
    if os.path.exists(PROVENANCE):
        try:
            old = json.load(open(PROVENANCE))
            if "run" in old:
                prov["run"] = old["run"]
        except (json.JSONDecodeError, OSError):
            pass
    json.dump(prov, open(PROVENANCE, "w"), indent=2)
    return PROVENANCE


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-patients", type=int, default=1)
    ap.add_argument("--out", default=DATA_DIR)
    ap.add_argument("--list-retest", action="store_true",
                    help="report test-retest candidates (>=2 TRACE exams) and exit")
    args = ap.parse_args()

    by_patient = list_trace_exams()
    if args.list_retest:
        rt = {p: ex for p, ex in by_patient.items() if len(ex) >= 2}
        print(f"patients with >=2 TRACE exams (test-retest candidates): {len(rt)}")
        for p, ex in list(rt.items())[:10]:
            print(f"  {p}: {len(ex)} exams  dates={[e['date'] for e in ex]}")
        return 0

    patients = sorted(by_patient)[:args.n_patients]
    records = []
    for p in patients:
        for exam in by_patient[p][:2]:  # up to test+retest per patient
            print(f"[fetch] {p} study ...{exam['study'][-8:]} "
                  f"(mask={'yes' if exam.get('mask_uid') else 'no'})")
            records.append(fetch_exam(exam, args.out))
    prov = write_provenance(records)

    print("\n" + "=" * 76)
    print(f"FETCH SUMMARY -- {COLLECTION} (license {LICENSE}, DOI {DOI})")
    print("=" * 76)
    for r in records:
        print(f"  {r['patient']} | shape {tuple(r['shape'])} | "
              f"b={r['bvals']} | tumor voxels={r.get('n_tumor_voxels','-')}")
    n_rt = sum(1 for ex in by_patient.values() if len(ex) >= 2)
    print(f"  test-retest available in collection: YES ({n_rt} patients with >=2 exams)")
    print(f"  pixel data -> {os.path.relpath(args.out, _ROOT)}/ (git-ignored)")
    print(f"  provenance (committed) -> {os.path.relpath(prov, _ROOT)}")
    print("=" * 76)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
