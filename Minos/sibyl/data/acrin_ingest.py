"""
Real-data ingestion seam: ACRIN-6698 DWI -> UnitTable.

This module is the boundary between the (large, separately-downloaded) ACRIN-6698
imaging and the Tier-2 analysis pipeline. The whole analysis (Arms 1/2, trust,
controls) consumes only a ``UnitTable``; this module is the only place that touches
imaging files, so it can be swapped/validated independently.

What is verified about the public TCIA collection (see VERIFICATION.md):
  * Raw multi-b DWI source series ARE public (not ADC-only). Prefer the derived
    TRACE / MergedMSMB series as the clean per-b signal source.
  * Test/retest "coffee-break" repeat DWI IS public (71 analyzable pairs;
    separate download option in the full collection).
  * b-values are 0,100,600,800. Most sites: one 4-b series (single b0). One site:
    three 2-b series (0/100, 0/600, 0/800) -> the trial builds a MergedMSMB series
    = all non-zero b images + a SINGLE AVERAGED b0. So the canonical 4-b vector is
    [mean(b0s), b100, b600, b800]; this averaging is implemented in reconcile_b0().
  * QA analyzability flags live in the "Full Collection Ancillary Patient
    Information" XLSX (DWI QC ratings), NOT in DICOM. Column header strings must be
    read from the file -- they are parameterized here via ``col_map``.
  * Whole-tumor ROIs ship as DICOM SEG / MR mask objects.

Two entry points:
  * build_units_from_nifti(cases, ...): fully runnable now (NIfTI 4D + bval + mask).
    Use after converting the DICOM to NIfTI, or for testing with synthetic volumes.
  * build_units_from_acrin_dicom(...): walks an NBIA download, matches ACRIN series
    by SeriesDescription, applies the SEG/MR mask. Requires the actual download.

Dependencies: SimpleITK (present) for image I/O; pandas for the QA table.
"""

from pathlib import Path

import numpy as np

from sibyl.data.units import UnitTable
from sibyl.forward_model.ivim import ACRIN_B_SCHEME

# Canonical b-value order for the assembled signal vector.
CANONICAL_BVALS = ACRIN_B_SCHEME.numpy().astype(float)  # [0, 100, 600, 800]

# SeriesDescription substrings used by the ACRIN-6698 derived objects (TCIA data
# descriptor). Matched case-insensitively; first hit wins. TRACE is the cleanest
# per-b source; MergedMSMB is the reconciled 4-b set for the 2-b acquisition site.
ACRIN_DWI_SERIES_PATTERNS = ("dwi trace", "mergedmsmb", "ax dwi")

# Default column names for the ancillary QA spreadsheet. These are PLACEHOLDERS:
# the exact header strings must be confirmed against the downloaded XLSX (the one
# fact that could not be verified from public web text). Override via ``col_map``.
DEFAULT_COL_MAP = {
    "patient_id": "Patient ID",
    "qa_rating": "DWI_QC_rating",   # e.g. poor/moderate/good
    "analyzable": "analyzable",      # e.g. yes/no or pass/fail
}


# ---------------------------------------------------------------------------
# Pure numeric helpers (library-free; unit-tested in tests/)
# ---------------------------------------------------------------------------
def reconcile_b0(b0_values: np.ndarray, axis: int = 0) -> np.ndarray:
    """
    Reconcile multiple b=0 acquisitions into a single averaged b0, matching the
    trial's MergedMSMB construction. Works on scalars (ROI means) or volumes.

    Parameters
    ----------
    b0_values : np.ndarray
        Stack of b0 measurements along ``axis`` (>=1 entries).

    Returns
    -------
    np.ndarray
        The mean b0 (one fewer dimension than the input along ``axis``).
    """
    b0_values = np.asarray(b0_values, dtype=np.float64)
    return b0_values.mean(axis=axis)


def roi_mean_per_b(vol_per_b: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """
    Whole-tumor ROI mean signal at each b-value.

    Parameters
    ----------
    vol_per_b : np.ndarray
        Signal volume, shape [B, ...spatial] (b-value first).
    mask : np.ndarray
        Boolean/0-1 ROI mask over the spatial dims (...spatial).

    Returns
    -------
    np.ndarray
        ROI-mean signal per b-value, shape [B].
    """
    vol_per_b = np.asarray(vol_per_b, dtype=np.float64)
    m = np.asarray(mask).astype(bool)
    assert vol_per_b.shape[1:] == m.shape, f"mask {m.shape} vs volume {vol_per_b.shape[1:]}"
    if m.sum() == 0:
        return np.full(vol_per_b.shape[0], np.nan)
    flat = vol_per_b.reshape(vol_per_b.shape[0], -1)[:, m.reshape(-1)]
    return flat.mean(axis=1)


def assemble_4b_vector(per_b_signal: dict) -> np.ndarray:
    """
    Order a {bvalue: signal} mapping into the canonical [b0, b100, b600, b800]
    vector. Multiple b0 entries (a list/array under key 0) are averaged first.

    Parameters
    ----------
    per_b_signal : dict
        Keys are b-values (int/float). The b=0 entry may be a scalar or an
        array/list of repeated b0 measurements (reconciled by averaging).

    Returns
    -------
    np.ndarray
        Length-4 signal vector in canonical b-order.
    """
    out = []
    for b in CANONICAL_BVALS:
        key = _match_bkey(per_b_signal, b)
        val = per_b_signal[key]
        if b == 0 and np.ndim(val) >= 1:
            val = reconcile_b0(np.asarray(val))
        out.append(float(val))
    return np.array(out, dtype=np.float64)


def _match_bkey(d: dict, b: float, tol: float = 1.0):
    for k in d:
        if abs(float(k) - float(b)) <= tol:
            return k
    raise KeyError(f"b-value {b} not found among {list(d)}")


# ---------------------------------------------------------------------------
# Image I/O (SimpleITK) -- runnable NIfTI path
# ---------------------------------------------------------------------------
def _require_sitk():
    try:
        import SimpleITK as sitk  # noqa
        return sitk
    except Exception as e:  # pragma: no cover
        raise ImportError(
            "SimpleITK is required for ACRIN image ingestion. Install with "
            "`pip install SimpleITK`."
        ) from e


def _read_bvals(bval_path) -> np.ndarray:
    return np.loadtxt(bval_path).ravel().astype(float)


def load_dwi_nifti(dwi_path, bval_path):
    """
    Load a 4D DWI NIfTI + FSL-style .bval file as a {bvalue: [volumes]} mapping,
    grouping repeated acquisitions of the same b-value (so multiple b0s are kept
    for reconcile_b0). Returns (per_b_volumes, sorted_unique_bvals).
    """
    sitk = _require_sitk()
    img = sitk.GetArrayFromImage(sitk.ReadImage(str(dwi_path)))  # [N_vol, Z, Y, X] or [N_vol, Y, X]
    bvals = _read_bvals(bval_path)
    assert img.shape[0] == len(bvals), f"{img.shape[0]} volumes vs {len(bvals)} bvals"
    per_b = {}
    for b in np.unique(np.round(bvals)):
        idx = np.where(np.round(bvals) == b)[0]
        per_b[float(b)] = img[idx]  # [n_rep, ...spatial]
    return per_b, np.unique(np.round(bvals)).astype(float)


def load_mask_nifti(mask_path) -> np.ndarray:
    sitk = _require_sitk()
    return sitk.GetArrayFromImage(sitk.ReadImage(str(mask_path))).astype(bool)


def roi_4b_from_nifti(dwi_path, bval_path, mask_path) -> np.ndarray:
    """
    Full per-scan extraction: NIfTI 4D DWI + bvals + ROI mask -> canonical 4-b
    ROI-mean signal vector [b0, b100, b600, b800], with b0 reconciliation.
    """
    per_b_vols, _ = load_dwi_nifti(dwi_path, bval_path)
    mask = load_mask_nifti(mask_path)
    per_b_signal = {}
    for b, vols in per_b_vols.items():
        # ROI mean for each repeat, then keep all (b0 reconciled later by assemble).
        roi_means = np.array([roi_mean_per_b(v[None], mask)[0] for v in vols])
        per_b_signal[b] = roi_means if (b == 0 and len(roi_means) > 1) else roi_means.mean()
    return assemble_4b_vector(per_b_signal)


def build_units_from_nifti(cases, out_npz=None) -> UnitTable:
    """
    Build a UnitTable from pre-converted NIfTI test/retest scans.

    Parameters
    ----------
    cases : list of dict
        Each: {
          'unit_id': str,
          'dwi_test', 'bval_test', 'mask_test': paths,
          'dwi_retest', 'bval_retest', 'mask_retest': paths,
          'qa_pass': bool,
        }
    out_npz : path or None
        If given, save the UnitTable there.
    """
    ids, st, sr, qa = [], [], [], []
    for c in cases:
        st.append(roi_4b_from_nifti(c["dwi_test"], c["bval_test"], c["mask_test"]))
        sr.append(roi_4b_from_nifti(c["dwi_retest"], c["bval_retest"], c["mask_retest"]))
        ids.append(c["unit_id"])
        qa.append(bool(c.get("qa_pass", True)))
    table = UnitTable(
        unit_id=np.array(ids),
        sig_test=np.array(st),
        sig_retest=np.array(sr),
        qa_pass=np.array(qa, dtype=bool),
        is_synth_id=np.zeros(len(ids), dtype=bool),  # all real, in-vivo
        theta=None,
        snr=None,
    )
    if out_npz is not None:
        table.save(out_npz)
    return table


# ---------------------------------------------------------------------------
# DICOM path (requires the NBIA download) -- documented, keyed to ACRIN series
# ---------------------------------------------------------------------------
def build_units_from_acrin_dicom(*args, **kwargs):
    """
    Walk an NBIA ACRIN-6698 download, match DWI TRACE/MergedMSMB series per subject
    (test & retest), apply the whole-tumor SEG/MR mask, compute the 4-b ROI-mean
    vectors, join QA flags from the ancillary spreadsheet, and write a UnitTable.

    This requires the actual (tens-of-GB) download, which is not present in this
    environment, so the heavy walk is intentionally not executed here. The pure
    numeric core it relies on (reconcile_b0, roi_mean_per_b, assemble_4b_vector,
    roi_4b_from_nifti) is implemented and tested above; wiring the DICOM series
    discovery on top is the remaining step once the data is downloaded.

    Recommended implementation once data is available:
      1. Use SimpleITK ImageSeriesReader to group DICOM files by SeriesInstanceUID.
      2. Select DWI series whose SeriesDescription (0008,103e) matches
         ACRIN_DWI_SERIES_PATTERNS; read b-values from the diffusion tags
         (vendor-specific: 0018,9087 / 0019,100c / 0043,1039).
      3. Identify the test vs retest study for each subject (same-session pair).
      4. Load the whole-tumor SEG (DICOM SEG) or MR mask object; resample onto the
         DWI grid; pass to roi_mean_per_b. Reconcile multiple b0s with reconcile_b0.
      5. Join QA flags from the ancillary XLSX/CSV via DEFAULT_COL_MAP (override
         with the real header strings) and set qa_pass.
    """
    raise NotImplementedError(
        "build_units_from_acrin_dicom requires the NBIA ACRIN-6698 download. "
        "Download the test/retest subset via NBIA Data Retriever, then either "
        "convert to NIfTI and use build_units_from_nifti(), or implement the "
        "DICOM series discovery sketched in this function's docstring on top of "
        "the tested pure-numeric core."
    )


def load_qa_flags(qa_table_path, col_map=None):
    """
    Read QA analyzability flags from the ancillary spreadsheet (CSV or XLSX) into
    a {patient_id: qa_pass(bool)} dict. Column names are parameterized via col_map
    (defaults in DEFAULT_COL_MAP) because the exact ACRIN header strings must be
    confirmed against the downloaded file.
    """
    import pandas as pd
    cm = {**DEFAULT_COL_MAP, **(col_map or {})}
    path = str(qa_table_path)
    df = pd.read_excel(path) if path.lower().endswith((".xlsx", ".xls")) else pd.read_csv(path)

    pid = cm["patient_id"]
    if cm["analyzable"] in df.columns:
        flag = df[cm["analyzable"]].astype(str).str.lower().isin({"yes", "pass", "true", "1", "analyzable"})
    elif cm["qa_rating"] in df.columns:
        flag = ~df[cm["qa_rating"]].astype(str).str.lower().isin({"poor", "low", "fail", "non-analyzable"})
    else:
        raise KeyError(
            f"Neither '{cm['analyzable']}' nor '{cm['qa_rating']}' in columns {list(df.columns)}. "
            "Pass the real header strings via col_map."
        )
    return dict(zip(df[pid].astype(str), flag.astype(bool)))
