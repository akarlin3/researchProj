#!/usr/bin/env python
"""Download-on-demand fetch for Echo's real test--retest signals (ACRIN-6698).

POSTURE (identical to Gauge -- this script REUSES ``Gauge/scripts/fetch_invivo.py``, the
sanctioned data-handling template, read-only via importlib):
  * Dataset: TCIA collection **ACRIN-6698 / I-SPY2 Breast DWI**, same-day test--retest arm
    (TrT0/TrT1). License **CC-BY-4.0**, DOI 10.7937/tcia.kk02-6d95.
  * **No pixel data is committed or even stored.** Echo reduces each repeat to a whole-tumor
    **ROI-mean signal** (4 numbers per repeat) on the fly. Only the ROI-mean JSONs (under the
    git-ignored ``data/``) and the committed provenance manifest ``results/invivo_provenance.json``
    persist. Echo redistributes nothing.

Echo's region-level (whole-tumor ROI mean) posture matches Gauge's: the two repeats are
unregistered, so a per-voxel statement is not licensed -- the honest unit is the tumor.

Modes:
  python scripts/fetch_invivo.py --check                 # report data availability (no DL)
  python scripts/fetch_invivo.py --list                  # enumerate test-retest cells (API only)
  python scripts/fetch_invivo.py --fetch [--max-patients N]   # download-on-demand ROI means

Requires: pydicom (DICOM parse), network egress to the NBIA REST API. Gauge must be a
monorepo sibling (read-only) -- it supplies the proven DICOM/NBIA assembly.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from echo_repeat import provenance, _paths  # noqa: E402

DATA_DIR = ROOT / "data" / "invivo_retest"
MANIFEST = ROOT / "results" / "invivo_provenance.json"


def _load_gauge_fetch():
    """Read-only importlib load of Gauge's fetch script as a module (the data template)."""
    gpath = _paths.GAUGE / "scripts" / "fetch_invivo.py"
    if not gpath.exists():
        raise FileNotFoundError(f"Gauge fetch template not found at {gpath}")
    spec = importlib.util.spec_from_file_location("gauge_fetch", gpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _roi_mean_repeat(gf, cell, label):
    """Download one repeat (TrT0/TrT1), reduce to a whole-tumor S0-normalised ROI-mean signal."""
    vol, zlocs, bvals = gf.assemble_trace(gf._get_series_zip(cell["trace"][label]))
    mask = gf.assemble_mask(gf._get_series_zip(cell["mask"][label]), zlocs)
    nvox = int(mask.sum())
    if nvox == 0:
        raise ValueError(f"empty tumor mask for {cell['patient']} {label}")
    roi = vol[mask]                              # (n_vox, 4)
    sig = np.nanmean(roi, axis=0)               # (4,)
    s0 = float(sig[int(np.argmin(bvals))])
    if not np.isfinite(s0) or s0 <= 0:
        raise ValueError(f"bad S0 for {cell['patient']} {label}")
    return (sig / s0), bvals, nvox


def fetch(max_patients=None) -> int:
    gf = _load_gauge_fetch()
    cells = gf.list_test_retest_cells()
    by_patient = {}
    for c in cells:
        by_patient.setdefault(c["patient"], c)
    patients = sorted(by_patient)
    if max_patients:
        patients = patients[:max_patients]
    print(f"[echo-fetch] same-day TrT0/TrT1 cells: {len(cells)} across {len(by_patient)} "
          f"patients; fetching ROI means for {len(patients)}")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    exams, skips = [], []
    for i, pid in enumerate(patients, 1):
        cell = by_patient[pid]
        try:
            sig_a, bvals, nvox_a = _roi_mean_repeat(gf, cell, "TrT0")
            sig_b, _, nvox_b = _roi_mean_repeat(gf, cell, "TrT1")
        except Exception as e:                  # skip+log; never fabricate
            skips.append({"patient": pid, "reason": str(e)[:140]})
            print(f"  [{i}/{len(patients)}] {pid}: SKIP ({str(e)[:80]})")
            continue
        rec = {"tumor": pid, "study": cell["study"], "date": cell["date"],
               "timepoint": cell["timepoint"], "manufacturer": cell["manufacturer"],
               "bvals": bvals.tolist(),
               "signal_a": sig_a.tolist(), "signal_b": sig_b.tolist()}
        (DATA_DIR / f"{pid}.json").write_text(json.dumps(rec))
        exams.append({"tumor": pid, "timepoint": cell["timepoint"],
                      "manufacturer": cell["manufacturer"],
                      "n_tumor_voxels_a": nvox_a, "n_tumor_voxels_b": nvox_b,
                      "sha_a": provenance.sha16(sig_a.tobytes()),
                      "sha_b": provenance.sha16(sig_b.tobytes())})
        print(f"  [{i}/{len(patients)}] {pid}: OK (ROI voxels {nvox_a}/{nvox_b})")
    provenance.write_manifest(
        MANIFEST, provenance.ACRIN6698, exams,
        extra={"reduction": "whole-tumor ROI mean per b per repeat (S0-normalised)",
               "derived_via": "Gauge/scripts/fetch_invivo.py (read-only importlib reuse)",
               "n_pairs": len(exams), "n_skipped": len(skips), "skips": skips})
    print(f"[echo-fetch] wrote {len(exams)} ROI-mean pairs to {DATA_DIR} "
          f"(skipped {len(skips)}); pixel data NOT stored")
    print(f"[echo-fetch] provenance manifest {MANIFEST}")
    return 0 if exams else 2


def list_cells() -> int:
    gf = _load_gauge_fetch()
    cells = gf.list_test_retest_cells()
    by_patient = {}
    for c in cells:
        by_patient.setdefault(c["patient"], c)
    print(json.dumps({"n_cells": len(cells), "n_patients": len(by_patient),
                      "patients": sorted(by_patient)[:5] + ["..."]}, indent=2))
    return 0


def check() -> int:
    have_echo = DATA_DIR.exists() and any(DATA_DIR.glob("*.json"))
    n = len(list(DATA_DIR.glob("*.json"))) if have_echo else 0
    print(json.dumps({
        "echo_roi_signals_present": bool(have_echo),
        "echo_n_pairs": n,
        "echo_data_dir": str(DATA_DIR),
        "manifest": str(MANIFEST) if MANIFEST.exists() else None,
        "gauge_template_present": (_paths.GAUGE / "scripts" / "fetch_invivo.py").exists(),
        "dataset": provenance.ACRIN6698["name"],
        "license": provenance.ACRIN6698["license"],
        "doi": provenance.ACRIN6698["doi"],
    }, indent=2))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="download-on-demand ROI means")
    ap.add_argument("--list", action="store_true", help="enumerate test-retest cells (API only)")
    ap.add_argument("--check", action="store_true", help="report data availability (no DL)")
    ap.add_argument("--max-patients", type=int, default=None)
    args = ap.parse_args()
    if args.check:
        return check()
    if args.list:
        return list_cells()
    if args.fetch:
        return fetch(max_patients=args.max_patients)
    return check()


if __name__ == "__main__":
    raise SystemExit(main())
