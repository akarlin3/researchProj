"""OSIPI open abdomen scan: download-on-demand + provenance.

Target T1 (the 54.7% NLLS D* railing rate) is a REAL-data number. The data is the
in-vivo abdomen acquisition in the OSIPI TF2.4 open Zenodo archive **14605039**
(CC-BY-4.0). Mirroring Lattice's ``osipi.py`` / Lethe's ``fetch_invivo.py``:

* fetched **on demand** into a **gitignored** ``download/`` dir, never redistributed;
* a provenance manifest (record id, file SHA-256) is written;
* importing this module touches no network.

Acquisition (from the archive): 144x144x21 volume, 104 measurements over 12 unique
b-values {0,10,20,30,40,50,75,100,150,250,400,600} s/mm^2 (Philips 3T). The archive
ships ``Data/mask_abdomen_homogeneous.nii.gz`` -- a homogeneous-tissue ROI of **1932**
voxels; Fashion's "1618 high-SNR ROI voxels" is that ROI under an (unstated) SNR cut.
Gnomon documents its own selection: the homogeneous ROI, optionally restricted to a
b0-SNR threshold (b0-SNR = mean / std over the 15 b=0 repeats). The load-bearing
output is the D* boundary-railing **rate**, reported for the full ROI and SNR cuts.
"""
from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ZENODO_RECORD = "14605039"
ZIP_NAME = "OSIPI_TF24_data_phantoms.zip"
ZIP_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/{ZIP_NAME}?download=1"
ZIP_SHA256 = "2a53054d6e6e76335c9fcdae245e6003460db014e88147b6030d9de4dd650b3e"
LICENSE = "CC-BY-4.0 (OSIPI TF2.4 open data)"

_ABDOMEN = "Data/abdomen.nii.gz"
_BVAL = "Data/abdomen.bval"
_MASK = "Data/mask_abdomen_homogeneous.nii.gz"


def default_download_dir():
    return Path(__file__).resolve().parent.parent / "download"


def fetch(dest=None, url=ZIP_URL):
    """Download the OSIPI zip on demand (cached). Returns the local path."""
    dest = Path(dest) if dest else default_download_dir()
    dest.mkdir(parents=True, exist_ok=True)
    zip_path = dest / ZIP_NAME
    if not zip_path.exists():
        import urllib.request
        urllib.request.urlretrieve(url, zip_path)
    return zip_path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance_record(zip_path, verify=True):
    digest = sha256(zip_path)
    return {
        "source": "OSIPI TF2.4 IVIM code collection",
        "zenodo_record": ZENODO_RECORD,
        "zenodo_url": f"https://zenodo.org/records/{ZENODO_RECORD}",
        "license": LICENSE,
        "zip": str(zip_path),
        "sha256": digest,
        "sha256_ok": (digest == ZIP_SHA256) if verify else None,
        "note": "External CC-BY-4.0 data; fetched on demand, not redistributed in-tree.",
    }


def _read_nii_from_zip(zf, name, tmp_dir):
    # nibabel reads .nii.gz transparently from a path; extract to a temp file.
    import nibabel as nib
    out = Path(tmp_dir) / Path(name).name
    out.write_bytes(zf.read(name))
    return np.asarray(nib.load(str(out)).dataobj).astype(float)


@dataclass
class AbdomenROI:
    signals: np.ndarray   # (n, nb) b0-normalized
    bvalues: np.ndarray   # (nb,)
    snr: np.ndarray       # (n,) b0-SNR
    n_roi_total: int      # voxels in the homogeneous ROI before SNR cut
    snr_threshold: float | None


def load_abdomen_roi(zip_path, snr_threshold=None):
    """Load the homogeneous-ROI abdomen voxels, b0-normalized, with per-voxel b0-SNR.

    ``snr_threshold`` (b0-SNR) optionally restricts to high-SNR voxels; ``None`` keeps
    the full homogeneous ROI. Signals are divided by each voxel's mean-b0.
    """
    zip_path = Path(zip_path)
    with zipfile.ZipFile(zip_path) as zf, \
            _temp_dir(zip_path.parent / "_osipi_tmp") as tmp:
        bval = np.array([float(x) for x in zf.read(_BVAL).split()])
        arr = _read_nii_from_zip(zf, _ABDOMEN, tmp)        # (X,Y,Z,nb)
        mask = _read_nii_from_zip(zf, _MASK, tmp) > 0       # (X,Y,Z)
    vox = arr[mask]                                          # (n_roi, nb)
    b0 = bval == 0
    s0 = vox[:, b0].mean(axis=1)
    snr = s0 / (vox[:, b0].std(axis=1) + 1e-9)
    sig = vox / np.clip(s0[:, None], 1e-6, None)
    n_total = sig.shape[0]
    if snr_threshold is not None:
        keep = snr > float(snr_threshold)
        sig, snr = sig[keep], snr[keep]
    return AbdomenROI(signals=sig, bvalues=bval, snr=snr,
                      n_roi_total=n_total, snr_threshold=snr_threshold)


class _temp_dir:
    """Tiny context manager for a scratch dir (no external deps)."""

    def __init__(self, path):
        self.path = Path(path)

    def __enter__(self):
        self.path.mkdir(parents=True, exist_ok=True)
        return self.path

    def __exit__(self, *exc):
        import shutil
        shutil.rmtree(self.path, ignore_errors=True)
        return False
