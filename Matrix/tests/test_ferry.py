"""Ferry — the real-data substrate adapter: drop-in, contract, reproducibility, grounding.

Mechanism tests run OFFLINE on a deterministic synthetic-geometry substrate (no network,
no data blobs). One network-gated test exercises the REAL public RT dataset and is skipped
if TCIA/NBIA is unreachable.
"""
import hashlib
import os
from dataclasses import replace

import numpy as np
import pytest

from matrix import (MatrixConfig, Twin, Interfaces, run_iteration, run_loop,
                    NORMAL, TUMOR, OAR, TREAT)
from matrix.loop import run_iteration as loop_run_iteration
from matrix.ferry import (FerrySubstrate, GroundedTwin, run_grounded_loop,
                          load_substrate, FerryDataUnavailable, rescale_dose_to_band,
                          DOI, LICENSE)
from matrix.ferry.loop_grounded import run_iteration as grounded_run_iteration

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOOP_PY = os.path.join(HERE, "matrix", "loop.py")
LOOP_PY_BLOB_SHA1 = "4a34806ac4fa55c0ce5453b9864d37c67abfda92"
G = 24


def _synthetic_substrate(G=G):
    yy, xx = np.mgrid[0:G, 0:G]
    cx = cy = (G - 1) / 2.0
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    labels = np.full((G, G), NORMAL, int)
    labels[r <= 0.30 * G] = TUMOR
    labels[(xx < 0.18 * G) & (yy < 0.18 * G)] = OAR
    dose_gy = 20.0 + 40.0 * np.exp(-(r / (0.30 * G)) ** 2)
    return FerrySubstrate(G=G, labels=labels, dose_gy=dose_gy,
                          provenance=dict(dataset="SYNTHETIC-GEOMETRY", license="n/a", doi="n/a"))


def _cfg(G=G):
    return replace(MatrixConfig(), nx=G, ny=G)


# ---------------------------------------------------------- drop-in / contract --
def test_loop_py_is_byte_unchanged():
    data = open(LOOP_PY, "rb").read()
    blob_sha1 = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
    assert blob_sha1 == LOOP_PY_BLOB_SHA1, "Ferry must not change loop.py"


def test_grounded_driver_uses_the_genuine_loop_engine():
    # the grounded driver must reuse matrix.loop.run_iteration, not reimplement it.
    assert grounded_run_iteration is loop_run_iteration is run_iteration


def test_grounded_twin_is_a_twin_satisfying_the_contract():
    cfg = _cfg()
    g = GroundedTwin.from_substrate(cfg, _synthetic_substrate())
    assert isinstance(g, Twin)
    for fld in ("cfg", "labels", "D", "Dstar", "f", "f0", "highdstar", "lowsnr",
                "snr_map", "dose"):
        assert getattr(g, fld) is not None
    assert g.n_voxels == cfg.n_voxels
    scan = g.scan(np.random.default_rng(0))
    assert scan.shape == (len(cfg.bvals), cfg.n_voxels, cfg.n_noise)
    snap = g.truth_snapshot()
    assert {"D", "Dstar", "f", "labels", "highdstar", "lowsnr", "dose"} <= set(snap)


def test_grid_mismatch_raises():
    with pytest.raises(ValueError):
        GroundedTwin.from_substrate(replace(MatrixConfig(), nx=8, ny=8), _synthetic_substrate(G))


# --------------------------------------------------------------- reproducibility --
def test_grounded_twin_bit_reproducible_from_seed():
    cfg = _cfg()
    sub = _synthetic_substrate()
    a = GroundedTwin.from_substrate(cfg, sub)
    b = GroundedTwin.from_substrate(cfg, sub)
    for fld in ("D", "Dstar", "f", "labels", "highdstar", "lowsnr", "snr_map", "dose"):
        assert np.array_equal(getattr(a, fld), getattr(b, fld))


def test_synthetic_layer_varies_by_seed_but_labels_are_real_invariant():
    cfg = _cfg()
    sub = _synthetic_substrate()
    a = GroundedTwin.from_substrate(cfg, sub)
    c = GroundedTwin.from_substrate(cfg.with_seed(cfg.seed + 1), sub)
    assert not np.array_equal(a.f, c.f)            # synthetic perfusion reseeds
    assert np.array_equal(a.labels, c.labels)      # real anatomy is seed-independent


def test_labels_are_real_not_synthetic_layout():
    # the grounded twin's labels come from the substrate, not Twin._layout_labels.
    cfg = _cfg()
    sub = _synthetic_substrate()
    g = GroundedTwin.from_substrate(cfg, sub)
    assert np.array_equal(g.labels, sub.flat_labels())


# ------------------------------------------------------------------- grounding ----
def test_real_dose_grounds_into_band_and_preserves_geometry():
    cfg = _cfg()
    sub = _synthetic_substrate()
    band = rescale_dose_to_band(sub.dose_gy, cfg)
    assert band.min() >= cfg.dose_min - 1e-9 and band.max() <= cfg.dose_max + 1e-9
    # geometry preserved: the map is monotonic non-decreasing in real Gy, so sorting voxels
    # by real dose yields a non-decreasing band (higher real dose -> >= band dose).
    order = np.argsort(sub.dose_gy.reshape(-1), kind="mergesort")
    assert np.all(np.diff(band[order]) >= -1e-9)


def test_ground_dose_flag_controls_initial_dose():
    cfg = _cfg()
    sub = _synthetic_substrate()
    ga = GroundedTwin.from_substrate(cfg, sub, ground_dose=False)
    gf = GroundedTwin.from_substrate(cfg, sub, ground_dose=True)
    assert np.allclose(ga.dose, cfg.dose_baseline)      # anatomy-only: flat baseline
    assert gf.dose.std() > 0                             # full: real, non-uniform geometry


# ------------------------------------------------------------- loop still closes --
def test_grounded_loop_closes_end_to_end_and_reproducibly():
    cfg = _cfg()
    sub = _synthetic_substrate()
    twin, states = run_grounded_loop(cfg, sub, Interfaces.placeholders())
    assert len(states) == cfg.n_iter
    assert all(s.is_complete() for s in states)
    twin_b, states_b = run_grounded_loop(cfg, sub, Interfaces.placeholders())
    assert np.array_equal(states[-1].action, states_b[-1].action)
    assert np.array_equal(twin.f, twin_b.f)


def test_grounded_loop_suppresses_action_on_untrusted():
    from matrix import SPARE
    cfg = _cfg()
    sub = _synthetic_substrate()
    _, states = run_grounded_loop(cfg, sub, Interfaces.placeholders())
    # across the whole run, the trust gate forces ESCALATE on untrusted voxels:
    # zero TREAT/SPARE (i.e. zero "actions") survive on untrustworthy voxels.
    gated_actions = sum(int(np.isin(s.action[~s.trustworthy], (TREAT, SPARE)).sum())
                        for s in states)
    assert gated_actions == 0
    # and there WAS something to suppress (the ungated decision acted on some untrusted voxels)
    ungated_actions = sum(int(np.isin(s.action_ungated[~s.trustworthy], (TREAT, SPARE)).sum())
                          for s in states)
    assert ungated_actions > 0


# ---------------------------------------------------------- REAL data (gated) -----
def test_real_substrate_grounds_and_loop_closes():
    try:
        sub = load_substrate(G=32)
    except FerryDataUnavailable as e:
        pytest.skip(f"public RT dataset unavailable: {e}")
    assert sub.provenance["doi"] == DOI
    assert sub.provenance["license"] == LICENSE
    for lab in (NORMAL, TUMOR, OAR):
        assert np.any(sub.labels == lab)
    assert sub.dose_gy.max() > sub.dose_gy.min()         # real, non-uniform dose
    cfg = _cfg(32)
    twin, states = run_grounded_loop(cfg, sub, Interfaces.placeholders())
    assert all(s.is_complete() for s in states)
    # the loop closes: treatment winds down on real geometry
    assert int((states[-1].action == TREAT).sum()) == 0
