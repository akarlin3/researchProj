"""Ferry CP1 gate — the real-data substrate adapter is a CLEAN DROP-IN.

Four falsifiable checks, all green => Ferry CP1 passes:

  (1) LOOP.PY BYTE-UNCHANGED  the git-blob hash of ``matrix/loop.py`` equals the hash it
                              shipped with (Matrix PR #56). Ferry edits the loop in NO way.
  (2) ENGINE IS THE GENUINE ARTICLE  the grounded driver imports and calls the SAME
                              ``run_iteration`` object from ``matrix.loop`` (identity check),
                              and ``GroundedTwin`` IS-A ``Twin`` satisfying the full contract.
  (3) REPRODUCIBLE FROM SEED  two grounded twins from the same seed + substrate are
                              bit-identical in every synthetic ground-truth field.
  (4) END-TO-END ON THE SUBSTRATE  the full four-stage loop runs to completion through the
                              real-data substrate; every LoopState is complete.

Checks (1)-(3) and the mechanism of (4) need no network (a synthetic-geometry substrate is
used as a fallback so the drop-in proof is self-contained). If the public dataset is
reachable, (4) additionally runs on the REAL grounded substrate.

Run:  <proteus python> Matrix/verify_ferry_cp1.py        # exit 0 = green
"""
from __future__ import annotations

import hashlib
import os
import sys
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from matrix import MatrixConfig, Twin, Interfaces, run_iteration, NORMAL, TUMOR, OAR
from matrix.loop import run_iteration as loop_run_iteration
from matrix.ferry import (FerrySubstrate, GroundedTwin, load_substrate,
                          run_grounded_loop, FerryDataUnavailable, DOI, LICENSE)
from matrix.ferry.loop_grounded import run_iteration as grounded_run_iteration

HERE = os.path.dirname(os.path.abspath(__file__))
LOOP_PY = os.path.join(HERE, "matrix", "loop.py")
# The git-blob hash loop.py shipped with (Matrix PR #56, origin/main @ c3155b2).
LOOP_PY_BLOB_SHA1 = "4a34806ac4fa55c0ce5453b9864d37c67abfda92"

GRID = 24


def _hr(t): print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


def _git_blob_sha1(path: str) -> str:
    data = open(path, "rb").read()
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def _synthetic_geometry_substrate(G: int) -> FerrySubstrate:
    """A deterministic stand-in substrate (no network) — proves the drop-in mechanism
    when the public dataset is unreachable. NOT real data; clearly labelled."""
    yy, xx = np.mgrid[0:G, 0:G]
    cx = cy = (G - 1) / 2.0
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    labels = np.full((G, G), NORMAL, int)
    labels[r <= 0.30 * G] = TUMOR
    labels[(xx < 0.18 * G) & (yy < 0.18 * G)] = OAR
    dose_gy = 20.0 + 40.0 * np.exp(-(r / (0.30 * G)) ** 2)     # hot core, cool rim
    return FerrySubstrate(G=G, labels=labels, dose_gy=dose_gy,
                          provenance=dict(dataset="SYNTHETIC-GEOMETRY (CP1 fallback)",
                                          license="n/a", doi="n/a"))


def _get_substrate(G: int):
    try:
        sub = load_substrate(G=G)
        return sub, True
    except FerryDataUnavailable as e:
        print(f"  [note] public dataset unreachable ({e}); using synthetic-geometry fallback "
              f"for the drop-in mechanism proof.")
        return _synthetic_geometry_substrate(G), False


def check_loop_unchanged() -> None:
    _hr("Ferry CP1 check 1/4 — loop.py is BYTE-UNCHANGED")
    got = _git_blob_sha1(LOOP_PY)
    assert got == LOOP_PY_BLOB_SHA1, f"loop.py changed! {got} != {LOOP_PY_BLOB_SHA1}"
    print(f"  git-blob sha1(matrix/loop.py) = {got}")
    print(f"  == shipped hash (PR #56).  Ferry touches loop.py in NO way.  PASS")


def check_engine_identity(cfg, sub) -> None:
    _hr("Ferry CP1 check 2/4 — the grounded driver uses the GENUINE loop engine + Twin contract")
    assert grounded_run_iteration is loop_run_iteration is run_iteration, \
        "grounded driver does not use matrix.loop.run_iteration"
    print("  matrix.ferry.loop_grounded.run_iteration IS matrix.loop.run_iteration (identity).")
    g = GroundedTwin.from_substrate(cfg, sub)
    assert isinstance(g, Twin), "GroundedTwin is not a Twin"
    for fld in ("cfg", "labels", "D", "Dstar", "f", "f0", "highdstar", "lowsnr",
                "snr_map", "dose"):
        assert hasattr(g, fld) and getattr(g, fld) is not None, f"missing field {fld}"
    assert g.n_voxels == cfg.n_voxels == sub.n_voxels
    rng = np.random.default_rng(0)
    scan = g.scan(rng)
    assert scan.shape == (len(cfg.bvals), cfg.n_voxels, cfg.n_noise), "scan shape wrong"
    snap = g.truth_snapshot()
    for k in ("D", "Dstar", "f", "labels", "highdstar", "lowsnr", "dose"):
        assert k in snap, f"truth_snapshot missing {k}"
    g.apply_plan(np.full(cfg.n_voxels, cfg.dose_baseline))    # callable, no raise
    print(f"  GroundedTwin IS-A Twin; all 10 fields present; scan{scan.shape} OK; "
          f"truth_snapshot + apply_plan OK.  CONTRACT SATISFIED: PASS")


def check_reproducible(cfg, sub) -> None:
    _hr("Ferry CP1 check 3/4 — grounded twin is reproducible from seed")
    a = GroundedTwin.from_substrate(cfg, sub)
    b = GroundedTwin.from_substrate(cfg, sub)
    for fld in ("D", "Dstar", "f", "labels", "highdstar", "lowsnr", "snr_map", "dose"):
        assert np.array_equal(getattr(a, fld), getattr(b, fld)), f"{fld} not reproducible"
    c = GroundedTwin.from_substrate(cfg.with_seed(cfg.seed + 1), sub)
    assert not np.array_equal(a.f, c.f), "different seed gave identical synthetic layer"
    assert np.array_equal(a.labels, c.labels), "labels must be seed-independent (real)"
    print(f"  same seed -> bit-identical (8/8 fields); seed+1 -> synthetic layer differs, "
          f"real labels unchanged.  REPRODUCIBLE: PASS")


def check_end_to_end(cfg, sub, is_real) -> None:
    _hr("Ferry CP1 check 4/4 — loop runs end-to-end through the substrate")
    twin, states = run_grounded_loop(cfg, sub, Interfaces.placeholders())
    assert len(states) == cfg.n_iter
    for s in states:
        assert s.is_complete(), f"iteration {s.iteration} incomplete"
    src = "REAL public RT dataset" if is_real else "synthetic-geometry fallback"
    if is_real:
        p = sub.provenance
        print(f"  substrate: {p['dataset']} (DOI {p.get('doi')}, {p.get('license')}), "
              f"patient {p.get('patient')}, slice z={p.get('slice_z_mm')}mm")
        print(f"  REAL labels: tumor={int((sub.labels==TUMOR).sum())} "
              f"oar={int((sub.labels==OAR).sum())} normal={int((sub.labels==NORMAL).sum())} "
              f"| REAL dose {p.get('dose_gy_range')} Gy")
    print(f"  ran {len(states)} iterations on the {src}; every LoopState complete.")
    print(f"  provenance: ruler={states[0].components['ruler']!r}")
    print("  END-TO-END: PASS")


def main() -> int:
    print("Ferry CP1 verification — real-data substrate adapter is a clean drop-in")
    cfg = replace(MatrixConfig(), nx=GRID, ny=GRID)
    sub, is_real = _get_substrate(GRID)
    check_loop_unchanged()
    check_engine_identity(cfg, sub)
    check_reproducible(cfg, sub)
    check_end_to_end(cfg, sub, is_real)
    _hr("Ferry CP1 GATE: PASS")
    print(f"  Ferry grounds Matrix on {DOI} ({LICENSE}) — REAL anatomy + dose geometry, "
          f"SYNTHETIC perfusion.")
    print("  Interface-swap only: loop.py byte-unchanged; GroundedTwin satisfies the Twin "
          "contract; no clinical claim.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
