"""S3 GCE-runner orchestration tests (gce/run_fold.py).

ESMFold itself runs only on the GCE CUDA box, so it is NEVER loaded here. The
runner is designed with the model backend injectable: these tests drive the
host-agnostic orchestration — manifest integrity, per-sequence checkpoint/resume,
mean-pLDDT keep/drop gating, long-sequence chunk flagging, and the run summary —
with a deterministic FAKE folder. This is the same split the pipeline uses to
keep its GPU step GCE-only (cf. the S3 dry-run manifest test in test_smoke.py).

We import gce/run_fold.py by path (it lives outside src/ so the fold image can
COPY a single self-contained file).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_FOLD = os.path.join(REPO, "gce", "run_fold.py")


def _load_runner():
    spec = importlib.util.spec_from_file_location("proteus_run_fold", RUN_FOLD)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rf = _load_runner()


def _sha(seq: str) -> str:
    return hashlib.sha256(seq.upper().encode()).hexdigest()


# A 2-residue CA-only PDB whose B-factor column carries a chosen pLDDT, so the
# fake folder produces output the real mean_plddt_from_pdb() can also parse.
def _toy_pdb(plddt: float) -> str:
    b = f"{plddt:6.2f}"
    return (
        f"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00{b}           C\n"
        f"ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00{b}           C\n"
        "TER\nEND\n"
    )


def _manifest(entries):
    """entries: list of (id, seq). Build a minimal S3 manifest like s3_fold emits."""
    return {
        "schema": "proteus.s3_fold.job_manifest/v1",
        "run_location": "gce",
        "random_seed": 1729,
        "fold_params": {"model": "esmfold_v1", "device": "cuda",
                        "plddt_min": 70.0, "chunk_size": 400, "max_recycles": 3},
        "sequences": [{"id": i, "length": len(s), "sha256": _sha(s)} for i, s in entries],
    }


def _write_fasta(path, entries):
    with open(path, "w") as fh:
        for i, s in entries:
            fh.write(f">{i}\n{s}\n")


def test_mean_plddt_from_pdb_reads_bfactor():
    assert rf.mean_plddt_from_pdb(_toy_pdb(87.5)) == pytest.approx(87.5)
    assert rf.mean_plddt_from_pdb("HEADER only, no atoms") == 0.0


def test_normalize_plddt_scale_lifts_0_1_to_0_100():
    # transformers writes pLDDT 0–1; normalize_plddt_scale must lift it to 0–100 so the
    # keep/drop gate at plddt_min=70 works and the B-factors are the standard convention.
    pdb01 = _toy_pdb(0.9)
    assert rf.mean_plddt_from_pdb(pdb01) == pytest.approx(0.9)  # raw, would be "dropped"
    fixed = rf.normalize_plddt_scale(pdb01)
    assert rf.mean_plddt_from_pdb(fixed) == pytest.approx(90.0)
    # column widths preserved (still parseable, B-factor in cols 61–66)
    for line in fixed.splitlines():
        if line.startswith("ATOM"):
            assert line[12:16].strip() == "CA"
            assert float(line[60:66]) == pytest.approx(90.0)


def test_normalize_plddt_scale_is_idempotent_on_0_100():
    # already-0–100 (native ESMFold, or a second pass) must be left untouched.
    pdb100 = _toy_pdb(87.5)
    assert rf.normalize_plddt_scale(pdb100) == pdb100
    assert rf.normalize_plddt_scale(rf.normalize_plddt_scale(_toy_pdb(0.9))) \
        == rf.normalize_plddt_scale(_toy_pdb(0.9))


def test_fold_batch_gates_on_mean_plddt(tmp_path):
    entries = [("hi", "MKLV"), ("lo", "GGGG")]
    man = _manifest(entries)
    fasta = tmp_path / "shortlist.fasta"
    _write_fasta(fasta, entries)
    out = tmp_path / "out"

    # high-pLDDT for "hi" (kept), low for "lo" (dropped at plddt_min=70)
    plddt = {"MKLV": 91.0, "GGGG": 42.0}

    def fake_folder(seq, chunk_size, max_recycles):
        return _toy_pdb(plddt[seq]), rf.mean_plddt_from_pdb(_toy_pdb(plddt[seq]))

    summary = rf.fold_batch(man, str(fasta), str(out), folder=fake_folder, device="cpu")

    assert summary["n_folded"] == 2
    assert summary["n_kept"] == 1 and summary["kept_ids"] == ["hi"]
    assert summary["n_dropped"] == 1 and summary["dropped_ids"] == ["lo"]
    # both PDBs + per-seq result jsons + the rolling summary exist
    assert (out / "hi.pdb").exists() and (out / "lo.pdb").exists()
    assert (out / "hi.json").exists() and (out / "s3_results.json").exists()
    rec = json.loads((out / "hi.json").read_text())
    assert rec["kept"] is True and rec["mean_plddt"] == pytest.approx(91.0)
    assert rec["resumed"] is False


def test_fold_batch_resumes_and_skips_completed(tmp_path):
    entries = [("a", "MKLV"), ("b", "WYFW")]
    man = _manifest(entries)
    fasta = tmp_path / "s.fasta"
    _write_fasta(fasta, entries)
    out = tmp_path / "out"

    calls = []

    def folder(seq, chunk_size, max_recycles):
        calls.append(seq)
        return _toy_pdb(88.0), 88.0

    rf.fold_batch(man, str(fasta), str(out), folder=folder, device="cpu")
    assert sorted(calls) == ["MKLV", "WYFW"]  # both folded first pass

    # Second pass: a folder that BLOWS UP if invoked — proves completed ids are
    # skipped (interruptible-safe resume), not refolded.
    def explode(seq, chunk_size, max_recycles):
        raise AssertionError(f"should not refold completed id for seq {seq}")

    summary = rf.fold_batch(man, str(fasta), str(out), folder=explode, device="cpu")
    assert summary["n_folded"] == 2 and summary["n_kept"] == 2
    assert all(r["resumed"] for r in summary["results"]), "resumed ids must be flagged"


def test_fold_batch_flags_sha256_mismatch(tmp_path):
    entries = [("x", "MKLV")]
    man = _manifest(entries)
    # FASTA carries a DIFFERENT sequence than the manifest recorded -> integrity fail
    fasta = tmp_path / "tampered.fasta"
    _write_fasta(fasta, [("x", "MKLVAAAA")])
    out = tmp_path / "out"

    summary = rf.fold_batch(man, str(fasta), str(out), folder=lambda *a: (_toy_pdb(99), 99),
                            device="cpu")
    assert summary["n_folded"] == 0
    assert summary["n_errors"] == 1 and "sha256" in summary["errors"][0]
    assert not (out / "x.pdb").exists()


def test_fold_batch_flags_long_sequence_as_chunked(tmp_path):
    long_seq = "A" * 450  # > chunk_size (400) -> chunked flag set
    entries = [("long", long_seq)]
    man = _manifest(entries)
    fasta = tmp_path / "l.fasta"
    _write_fasta(fasta, entries)
    out = tmp_path / "out"

    seen = {}

    def folder(seq, chunk_size, max_recycles):
        seen["chunk_size"] = chunk_size  # orchestration passes config chunk_size through
        return _toy_pdb(80.0), 80.0

    summary = rf.fold_batch(man, str(fasta), str(out), folder=folder, device="cpu")
    assert seen["chunk_size"] == 400
    assert summary["results"][0]["chunked"] is True
    assert summary["n_kept"] == 1
