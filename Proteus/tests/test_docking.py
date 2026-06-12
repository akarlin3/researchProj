"""Docking orchestration tests (proteus.docking).

AutoDock Vina is a LOCAL-only dependency and isn't loaded here; the Vina call is
injectable, so these tests drive the host-agnostic orchestration — active-site box
placement (real S4 geometry), affinity ranking, top-N selection, seed/box wiring,
and the screen-hit input filter — with a deterministic FAKE scorer.

box_center_from_model uses only S4 (biotite), so the box tests need the control
PDBs but NOT fpocket/vina. Tests skip if the controls aren't fetched.
"""
from __future__ import annotations

import json
import os
import textwrap

import pytest

from proteus import docking
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")

# poly-ALA toy with no Ser/His/Asp -> S4 finds no triad -> not dockable.
_TOY_NO_TRIAD = textwrap.dedent("""\
    ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
    ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
    ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
    ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
    ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00  0.00           C
    TER
    END
""")


def _require(pdb_id):
    p = os.path.join(STRUCT, f"{pdb_id}.pdb")
    if not os.path.exists(p):
        pytest.skip(f"{pdb_id}.pdb not fetched — run controls/fetch_controls.py")
    return p


def test_box_center_is_catalytic_ser_og():
    cfg = load_config()
    pdb = _require("6EQE")
    center, ser = docking.box_center_from_model(pdb, cfg)
    assert ser == 160, "box should centre on IsPETase catalytic Ser160"
    assert center is not None and len(center) == 3
    assert all(abs(float(x)) < 1e4 for x in center)  # finite coordinate
    # it really is the Ser160 OG coordinate
    og = docking._ser_og_coord(pdb, 160)
    assert list(map(float, center)) == pytest.approx(list(map(float, og)))


def test_no_triad_model_is_not_dockable(tmp_path):
    cfg = load_config()
    pdb = tmp_path / "toy.pdb"
    pdb.write_text(_TOY_NO_TRIAD)
    called = []
    rec = docking.dock_model(str(pdb), cfg, lambda **k: called.append(k) or {"affinity": -9.0},
                             cand_id="toy")
    assert rec["docked"] is False
    assert rec["stage_failed"] == "no_catalytic_site"
    assert called == [], "scorer must not be called when there is no catalytic site"


def test_dock_models_ranks_by_affinity_and_keeps_top_n():
    cfg = load_config()
    cfg["docking"]["top_n"] = 2
    inputs = [{"id": i, "pdb": _require(i), "mean_plddt": None}
              for i in ("6EQE", "4EB0", "1TCA")]
    # fake affinities: 1TCA tightest, then 6EQE, then 4EB0
    aff = {"6EQE": -7.5, "4EB0": -6.0, "1TCA": -9.1}

    def fake(**kw):
        # the receptor path tells us which model is being docked
        rid = os.path.splitext(os.path.basename(kw["receptor"]))[0]
        return {"affinity": aff[rid], "n_poses": 9}

    summary = docking.dock_models(inputs, cfg, fake)
    assert summary["n_docked"] == 3 and summary["n_failed"] == 0
    # ranked ascending (most negative = tightest first)
    assert summary["ranking"] == ["1TCA", "6EQE", "4EB0"]
    assert summary["best_affinity"] == -9.1
    # top_n = 2
    assert summary["kept_ids"] == ["1TCA", "6EQE"]
    ranks = {c["id"]: c.get("rank") for c in summary["candidates"]}
    assert ranks == {"1TCA": 1, "6EQE": 2, "4EB0": 3}


def test_dock_passes_seed_and_box_from_config():
    cfg = load_config()
    pdb = _require("6EQE")
    seen = {}

    def fake(**kw):
        seen.update(kw)
        return {"affinity": -8.0, "n_poses": 1}

    rec = docking.dock_model(pdb, cfg, fake, cand_id="6EQE")
    assert rec["docked"] and rec["affinity"] == -8.0
    assert seen["seed"] == cfg["random_seed"] == 1729      # seeded from config, not hardcoded
    assert list(seen["box_size"]) == [20.0, 20.0, 20.0]
    # box centre handed to the scorer is the catalytic Ser OG
    assert list(map(float, seen["center"])) == pytest.approx(rec["box_center"], abs=1e-3)


def test_inputs_from_candidates_docks_only_hits(tmp_path):
    """Wiring from the S4/S5 screen: by default dock only the PETase-like hits."""
    models = tmp_path / "models"
    models.mkdir()
    for cid in ("hitA", "hitB", "missC"):
        (models / f"{cid}.pdb").write_text(_TOY_NO_TRIAD)
    cand = {"candidates": [
        {"id": "hitA", "petase_like_hit": True, "mean_plddt": 90.0},
        {"id": "hitB", "petase_like_hit": True, "mean_plddt": 85.0},
        {"id": "missC", "petase_like_hit": False, "mean_plddt": 88.0},
    ]}
    cj = tmp_path / "s4s5_candidates.json"
    cj.write_text(json.dumps(cand))

    hits = docking._inputs_from_candidates(str(cj), str(models), hits_only=True)
    assert {h["id"] for h in hits} == {"hitA", "hitB"}
    allm = docking._inputs_from_candidates(str(cj), str(models), hits_only=False)
    assert {h["id"] for h in allm} == {"hitA", "hitB", "missC"}


def test_package_is_runnable_as_module():
    """`python -m proteus.docking` must work — the package needs a __main__.py that
    re-exports main(). Guards the documented CLI invocation from regressing."""
    import runpy

    from proteus.docking import __main__ as dunder_main
    assert dunder_main.main is docking.main
    # the module is importable as proteus.docking.__main__ (what `-m` executes)
    spec = runpy.importlib.util.find_spec("proteus.docking.__main__")
    assert spec is not None

