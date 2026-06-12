"""Vast burst launcher tests (proteus.launch).

The launcher only PLANS (and validates) the burst by default — `vastai`/`rsync`/
`ssh` are never invoked unless --execute is passed with real ids/host. So the plan
builder and input validation are pure and fully testable here; no Vast box, no
network.
"""
from __future__ import annotations

import json

import pytest

from proteus import launch
from proteus.utils import load_config


def _manifest(tmp_path, n=3, run_location="vast"):
    p = tmp_path / "s3_job_manifest.json"
    p.write_text(json.dumps({"run_location": run_location, "n_sequences": n,
                             "sequences": [{"id": f"c{i}"} for i in range(n)]}))
    return str(p)


def _shortlist(tmp_path, n=3):
    p = tmp_path / "s2_shortlist.fasta"
    p.write_text("".join(f">c{i}\nMKLV\n" for i in range(n)))
    return str(p)


def test_validate_accepts_coherent_inputs(tmp_path):
    v = launch.validate_inputs(_manifest(tmp_path, 3), _shortlist(tmp_path, 3))
    assert v["ok"] is True and v["errors"] == []
    assert v["manifest"]["n_sequences"] == 3


def test_validate_flags_count_mismatch(tmp_path):
    v = launch.validate_inputs(_manifest(tmp_path, 5), _shortlist(tmp_path, 3))
    assert v["ok"] is False
    assert any("shortlist has 3" in e and "manifest says 5" in e for e in v["errors"])


def test_validate_refuses_non_vast_manifest(tmp_path):
    v = launch.validate_inputs(_manifest(tmp_path, 3, run_location="local"),
                               _shortlist(tmp_path, 3))
    assert v["ok"] is False
    assert any("run_location" in e for e in v["errors"])


def test_validate_missing_files(tmp_path):
    v = launch.validate_inputs(str(tmp_path / "nope.json"), str(tmp_path / "nope.fasta"))
    assert v["ok"] is False
    assert any("manifest not found" in e for e in v["errors"])


def test_build_plan_orders_the_burst_cycle(tmp_path):
    cfg = load_config()
    man, fa = _manifest(tmp_path), _shortlist(tmp_path)
    plan = launch.build_plan(cfg, man, fa, offer_id="12345", instance_id="67890",
                             ssh="root@1.2.3.4", bid="0.20")
    names = [s["name"] for s in plan]
    assert names == ["search_offers", "create_instance", "upload", "fold",
                     "download", "destroy_instance"]
    by = {s["name"]: s["cmd"] for s in plan}

    # create uses the resolved offer + disk; interruptible (config default) => --bid
    assert "12345" in by["create_instance"] and "--bid" in by["create_instance"]
    assert "0.20" in by["create_instance"]
    # upload ships ONLY the shortlist + manifest, to the mounted input volume
    assert man in by["upload"] and fa in by["upload"]
    assert any(a.endswith(":/data/proteus/in/") for a in by["upload"])
    # fold ssh-runs the on-box runner against the uploaded manifest/fasta
    assert by["fold"][:3] == ["ssh", "root@1.2.3.4", "python3"]
    assert "/opt/proteus/run_fold.py" in by["fold"]
    assert "/data/proteus/in/s3_job_manifest.json" in by["fold"]
    assert "/data/proteus/out/" in by["fold"]
    # download pulls results back to the local folded dir; destroy uses the instance id
    assert by["download"][-1].endswith("structures/folded/")
    assert "67890" in by["destroy_instance"]


def test_build_plan_uses_placeholders_when_ids_absent(tmp_path):
    cfg = load_config()
    plan = launch.build_plan(cfg, _manifest(tmp_path), _shortlist(tmp_path))
    by = {s["name"]: s["cmd"] for s in plan}
    assert "<OFFER_ID>" in by["create_instance"]
    assert "<INSTANCE_ID>" in by["destroy_instance"]
    # registry image unset in default config -> a clearly-flagged placeholder
    assert any("REGISTRY_IMAGE" in tok for tok in by["create_instance"])


def test_main_dry_run_does_not_execute(tmp_path, capsys, monkeypatch):
    # if anything tried to actually run, this would raise
    monkeypatch.setattr(launch.subprocess, "run",
                        lambda *a, **k: pytest.fail("dry-run must not execute commands"))
    rc = launch.main(["--manifest", _manifest(tmp_path), "--shortlist", _shortlist(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "DRY-RUN" in out and "vastai search offers" in out


def test_main_returns_error_on_bad_inputs(tmp_path):
    rc = launch.main(["--manifest", _manifest(tmp_path, 5),
                      "--shortlist", _shortlist(tmp_path, 3)])
    assert rc == 2
