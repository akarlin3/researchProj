"""Ferry dataset loader — grounds the substrate on a PUBLIC RT dataset, by script.

Dataset (CP0 selection, clinical-leaning venue):
    TCIA **Pancreatic-CT-CBCT-SEG** (Version 2, 2022-08-23)
    DOI 10.7937/TCIA.ESHQ-4D90 — license **CC BY 4.0** (gate-free, no registration).
    See ``LICENSE_DATASET.md`` for the verbatim license + citation.

What is loaded (clean IP — guardrail 2):
  * **NO data blobs are committed.** This module *downloads by script* from the TCIA NBIA
    REST API into a git-ignored cache (``_cache/``), and persists only the tiny derived
    ``G x G`` label + dose grids (also git-ignored). Re-running regenerates them.
  * For one patient it pulls only the **RTSTRUCT** (clinician contours) and **RTDOSE**
    (delivered 3-D dose grid) series — a few tens of MB — never the CT/CBCT image volumes.

Mapping to the twin's anatomy labels:
  * **TUMOR**  = the treatment target VOI (``ROI``)
  * **OAR**    = small bowel + stomach/duodenum (the abdominal organs-at-risk, ``*_planCT``)
  * **NORMAL** = remaining tissue inside the crop

A representative axial slice (max target area) is rasterised onto the RTDOSE in-plane grid,
cropped to the target's bounding box (+margin), and block-downsampled to ``G x G``
(majority-vote labels, mean dose).
"""
from __future__ import annotations

import io
import os
import zipfile

import numpy as np

from ..config import NORMAL, TUMOR, OAR

# ---- dataset identity / provenance (no data, just the public pointers) -----------
COLLECTION = "Pancreatic-CT-CBCT-SEG"
DOI = "10.7937/TCIA.ESHQ-4D90"
LICENSE = "CC BY 4.0"
NBIA_API = "https://services.cancerimagingarchive.net/nbia-api/services/v1"
DEFAULT_PATIENT = "Pancreas-CT-CB_001"

# ROI-name → twin label (planning-CT structure set; case-insensitive substring match)
TARGET_NAMES = ("roi",)                                   # the treatment target VOI
OAR_NAMES = ("bowel_sm", "stomach_duo")                   # abdominal OARs near the target

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(HERE, "_cache")


# --------------------------------------------------------------------- networking --
def _http_get(url: str, timeout: float = 180.0) -> bytes:
    """GET bytes with stdlib only (no extra deps)."""
    import urllib.request
    with urllib.request.urlopen(url, timeout=timeout) as r:    # noqa: S310 (trusted host)
        return r.read()


def _series_for(modality: str, patient: str) -> list:
    import json
    url = f"{NBIA_API}/getSeries?Collection={COLLECTION}&PatientID={patient}&Modality={modality}"
    return json.loads(_http_get(url, timeout=60).decode())


def _download_series_dicoms(series_uid: str):
    """Return the parsed pydicom datasets in a series (downloads a zip into memory)."""
    import pydicom
    raw = _http_get(f"{NBIA_API}/getImage?SeriesInstanceUID={series_uid}")
    out = []
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        for nm in zf.namelist():
            if nm.lower().endswith(".dcm"):
                out.append(pydicom.dcmread(io.BytesIO(zf.read(nm))))
    return out


# ----------------------------------------------------------------- rasterisation --
def _planct_struct(structs):
    """Pick the RTSTRUCT that carries the treatment target (the planning-CT set)."""
    for s in structs:
        names = [r.ROIName.lower() for r in s.StructureSetROISequence]
        if any(any(t in n for t in TARGET_NAMES) for n in names):
            return s
    raise ValueError("no RTSTRUCT in this patient carries the target ROI")


def _roi_contours(struct):
    num2name = {r.ROINumber: r.ROIName for r in struct.StructureSetROISequence}
    out = {}
    for rc in struct.ROIContourSequence:
        nm = num2name.get(rc.ReferencedROINumber, str(rc.ReferencedROINumber))
        out[nm] = [np.asarray(c.ContourData, float).reshape(-1, 3)
                   for c in getattr(rc, "ContourSequence", [])]
    return out


def _match(names_lc, name):
    return any(t in name.lower() for t in names_lc)


def _rasterise(dose_ds, struct, G: int, margin: int = 12):
    """Rasterise target + OAR contours onto the dose grid, crop to target, downsample to GxG.

    Returns ``(labels[G,G], dose_gy[G,G], meta)``.
    """
    from matplotlib.path import Path

    ipp = np.asarray(dose_ds.ImagePositionPatient, float)
    dr, dc = map(float, dose_ds.PixelSpacing)            # row(y), col(x) spacing (mm)
    gfov = np.asarray(dose_ds.GridFrameOffsetVector, float)
    vol = dose_ds.pixel_array.astype(float) * float(dose_ds.DoseGridScaling)  # (nz,ny,nx) Gy
    nz, ny, nx = vol.shape
    zc = ipp[2] + gfov
    x0, y0 = ipp[0], ipp[1]

    cols, rows = np.meshgrid(np.arange(nx), np.arange(ny))
    grid_pts = np.column_stack([cols.ravel(), rows.ravel()])

    def fill(contours, z, tol=1.0):
        mask = np.zeros((ny, nx), bool)
        for poly in contours:
            if abs(poly[0, 2] - z) > tol:
                continue
            cc = (poly[:, 0] - x0) / dc
            rr = (poly[:, 1] - y0) / dr
            inside = Path(np.column_stack([cc, rr])).contains_points(grid_pts).reshape(ny, nx)
            mask ^= inside                                # even-odd fill (handles holes)
        return mask

    contours = _roi_contours(struct)
    target = [c for nm, cs in contours.items() if _match(TARGET_NAMES, nm) for c in cs]
    oar = {nm: cs for nm, cs in contours.items() if _match(OAR_NAMES, nm)}

    # choose the axial slice with the largest target area
    best_z, best_area = zc[0], -1
    for z in zc:
        a = int(fill(target, z).sum())
        if a > best_area:
            best_area, best_z = a, z
    fi = int(np.argmin(np.abs(zc - best_z)))

    tum = fill(target, best_z)
    oar_mask = np.zeros((ny, nx), bool)
    for cs in oar.values():
        oar_mask |= fill(cs, best_z)
    oar_mask &= ~tum
    dose2d = vol[fi]

    ys, xs = np.where(tum)
    if ys.size == 0:
        raise ValueError("target empty on selected slice")
    r0, r1 = max(ys.min() - margin, 0), min(ys.max() + margin + 1, ny)
    c0, c1 = max(xs.min() - margin, 0), min(xs.max() + margin + 1, nx)
    lbl = np.where(tum, TUMOR, np.where(oar_mask, OAR, NORMAL))[r0:r1, c0:c1]
    dcrop = dose2d[r0:r1, c0:c1]

    L = _block_label(lbl, G)
    Dg = _block_mean(dcrop, G)
    meta = dict(slice_frame=fi, slice_z_mm=round(float(best_z), 2),
                crop=[int(r0), int(r1), int(c0), int(c1)],
                target_px_fullres=int(best_area),
                dose_gy_range=[round(float(dcrop.min()), 2), round(float(dcrop.max()), 2)])
    return L, Dg, meta


def _block_label(lbl, G):
    H, W = lbl.shape
    ys = np.linspace(0, H, G + 1).astype(int)
    xs = np.linspace(0, W, G + 1).astype(int)
    out = np.zeros((G, G), int)
    for i in range(G):
        for j in range(G):
            blk = lbl[ys[i]:ys[i + 1], xs[j]:xs[j + 1]]
            if blk.size == 0:
                continue
            vals, cnts = np.unique(blk, return_counts=True)
            out[i, j] = int(vals[np.argmax(cnts)])
    return out


def _block_mean(a, G):
    H, W = a.shape
    ys = np.linspace(0, H, G + 1).astype(int)
    xs = np.linspace(0, W, G + 1).astype(int)
    out = np.zeros((G, G))
    for i in range(G):
        for j in range(G):
            blk = a[ys[i]:ys[i + 1], xs[j]:xs[j + 1]]
            out[i, j] = float(blk.mean()) if blk.size else 0.0
    return out


# ------------------------------------------------------------------- public API ----
def load_substrate(G: int = 32, patient: str = DEFAULT_PATIENT,
                   cache_dir: str = CACHE_DIR, allow_download: bool = True):
    """Load (or build, cached) a :class:`FerrySubstrate` at grid resolution ``G``.

    Tries the git-ignored cache first; otherwise downloads the patient's RTSTRUCT+RTDOSE
    from the TCIA NBIA API, rasterises, and caches the derived ``G x G`` grids. Raises
    ``FerryDataUnavailable`` if there is no cache and download is disabled/unavailable.
    """
    from .substrate import FerrySubstrate

    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, f"substrate_{patient}_G{G}.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return FerrySubstrate(G=int(z["G"]), labels=z["labels"], dose_gy=z["dose_gy"],
                              provenance=dict(z["provenance"].item()))

    if not allow_download:
        raise FerryDataUnavailable(f"no cache at {cache} and download disabled")
    try:
        struct_series = _series_for("RTSTRUCT", patient)
        dose_series = _series_for("RTDOSE", patient)
        if not struct_series or not dose_series:
            raise FerryDataUnavailable(f"{patient}: missing RTSTRUCT or RTDOSE series")
        structs = []
        for s in struct_series:
            structs += _download_series_dicoms(s["SeriesInstanceUID"])
        dose_dcms = _download_series_dicoms(dose_series[0]["SeriesInstanceUID"])
        struct = _planct_struct(structs)
        dose_ds = dose_dcms[0]
    except FerryDataUnavailable:
        raise
    except Exception as e:                                # network / parse failure
        raise FerryDataUnavailable(f"could not fetch {patient} from NBIA: {e!r}") from e

    labels, dose_gy, meta = _rasterise(dose_ds, struct, G)
    provenance = dict(
        dataset=COLLECTION, doi=DOI, license=LICENSE, host="TCIA / NBIA",
        patient=patient, G=G,
        struct_series=struct.SeriesInstanceUID, dose_series=dose_ds.SeriesInstanceUID,
        frame_of_reference=str(getattr(dose_ds, "FrameOfReferenceUID", "")),
        target_roi=[r.ROIName for r in struct.StructureSetROISequence
                    if _match(TARGET_NAMES, r.ROIName)],
        oar_rois=[r.ROIName for r in struct.StructureSetROISequence
                  if _match(OAR_NAMES, r.ROIName)],
        **meta,
    )
    np.savez(cache, G=G, labels=labels, dose_gy=dose_gy,
             provenance=np.array(provenance, dtype=object))
    return FerrySubstrate(G=G, labels=labels, dose_gy=dose_gy, provenance=provenance)


class FerryDataUnavailable(RuntimeError):
    """Raised when the grounded substrate cannot be loaded (no cache and no network)."""
