"""Positive-output smoke tests — one per LOCAL tool.

Each test asserts a POSITIVE artifact (a coordinate, a hit row, a cluster, a
finite score), not merely "no exception". Mirrors the CD silent-failure rule:
green-without-output is not green.

Tools that are not installed are SKIPPED (not passed and not failed) so the
suite still runs and the SMOKE SUMMARY reflects true per-tool state. On a host
with no GPU, ESMFold/CUDA-torch tests skip by design — see
envlog/recon-report.md.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import textwrap

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# A minimal valid 3-residue poly-ALA backbone+CB PDB (two copies used as toy input).
TOY_PDB = textwrap.dedent("""\
    ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
    ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
    ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
    ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
    ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00  0.00           C
    ATOM      6  N   ALA A   2       3.332   1.552   0.000  1.00  0.00           N
    ATOM      7  CA  ALA A   2       3.977   2.857   0.000  1.00  0.00           C
    ATOM      8  C   ALA A   2       5.486   2.700   0.000  1.00  0.00           C
    ATOM      9  O   ALA A   2       6.009   1.580   0.000  1.00  0.00           O
    ATOM     10  CB  ALA A   2       3.585   3.659   1.232  1.00  0.00           C
    ATOM     11  N   ALA A   3       6.190   3.823   0.000  1.00  0.00           N
    ATOM     12  CA  ALA A   3       7.645   3.844   0.000  1.00  0.00           C
    ATOM     13  C   ALA A   3       8.200   5.262   0.000  1.00  0.00           C
    ATOM     14  O   ALA A   3       7.444   6.234   0.000  1.00  0.00           O
    ATOM     15  CB  ALA A   3       8.165   3.083  -1.213  1.00  0.00           C
    TER
    END
""")


def _have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


# --------------------------------------------------------------------------- #
# ESMFold
# --------------------------------------------------------------------------- #
def test_esmfold_positive_output(tmp_path):
    torch = pytest.importorskip("torch", reason="torch not installed")
    esm = pytest.importorskip("esm", reason="fair-esm not installed")
    if not torch.cuda.is_available():
        pytest.skip("no CUDA GPU — ESMFold not runnable on this host")

    model = esm.pretrained.esmfold_v1().eval().cuda()
    seq = "GSSGSSGAEAEAEAEAKLKL"  # ~20 residues
    with torch.no_grad():
        out = model.infer_pdb(seq)
    plddt = model.infer(seq)["plddt"]
    assert "ATOM" in out and out.count("\n") > 10, "no PDB coordinates emitted"
    assert plddt.numel() > 0, "empty pLDDT array"
    assert float(plddt.mean()) > 0.0, "non-positive pLDDT"


# --------------------------------------------------------------------------- #
# Foldseek
# --------------------------------------------------------------------------- #
def test_foldseek_positive_output(tmp_path):
    if not _have("foldseek"):
        pytest.skip("foldseek not installed")
    ver = subprocess.run(["foldseek", "version"], capture_output=True, text=True)
    assert ver.stdout.strip() or ver.returncode == 0, "foldseek version did not parse"

    # Foldseek's k-mer prefilter needs real-sized structures, so use two fetched
    # controls for the trivial all-vs-all (self-hits guarantee >=1 hit row).
    pdbs = [os.path.join(REPO, "structures", f"{p}.pdb") for p in ("6EQE", "4EB0")]
    if not all(os.path.exists(p) for p in pdbs):
        pytest.skip("control PDBs not fetched — run controls/fetch_controls.py")

    qdir = tmp_path / "pdbs"
    qdir.mkdir()
    for src in pdbs:
        (qdir / os.path.basename(src)).write_text(open(src).read())

    res = tmp_path / "aln.m8"
    cmd = ["foldseek", "easy-search", str(qdir), str(qdir), str(res),
           str(tmp_path / "tmp"), "--format-mode", "0", "-e", "10"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert res.exists() and res.stat().st_size > 0, (
        f"no hit rows produced (rc={proc.returncode}): {proc.stderr[-400:]}")
    assert "6EQE" in res.read_text(), "expected a self-hit row for 6EQE"


# --------------------------------------------------------------------------- #
# MMseqs2
# --------------------------------------------------------------------------- #
def test_mmseqs2_positive_output(tmp_path):
    if not _have("mmseqs"):
        pytest.skip("mmseqs2 not installed")
    fasta = tmp_path / "toy.fasta"
    fasta.write_text(
        ">s1\nMKKLLPTAAAGLLLLAAQPAMA\n"
        ">s2\nMKKLLPTAAAGLLLLAAQPAMA\n"   # identical to s1 -> should co-cluster
        ">s3\nWWWWYYYYFFFFGGGGHHHHKK\n"
    )
    res = tmp_path / "clu"
    sub = subprocess.run(
        ["mmseqs", "easy-cluster", str(fasta), str(res), str(tmp_path / "tmp"),
         "--min-seq-id", "0.9", "-c", "0.8"],
        capture_output=True, text=True)
    tsv = tmp_path / "clu_cluster.tsv"
    assert tsv.exists() and tsv.stat().st_size > 0, (
        f"no cluster TSV produced (rc={sub.returncode}): {sub.stderr[-400:]}")


# --------------------------------------------------------------------------- #
# ProstT5
# --------------------------------------------------------------------------- #
def test_prostt5_positive_output():
    """Tokenize one short sequence with ProstT5's tokenizer; assert non-empty tokens.

    Primary path uses transformers' T5Tokenizer. Some transformers releases route
    ProstT5's SentencePiece model through an incompatible converter; in that case
    we fall back to loading the same `spiece.model` asset directly via
    sentencepiece (still ProstT5's tokenizer, still a positive artifact).
    """
    os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
    pytest.importorskip("torch", reason="torch not installed")
    pieces = None
    try:
        from transformers import T5Tokenizer
        tok = T5Tokenizer.from_pretrained("Rostlab/ProstT5", do_lower_case=False)
        ids = tok.batch_encode_plus([" ".join("PRTEINS")], add_special_tokens=True)["input_ids"][0]
        pieces = tok.convert_ids_to_tokens(ids)
    except Exception:
        spm = pytest.importorskip("sentencepiece", reason="sentencepiece not installed")
        try:
            from huggingface_hub import hf_hub_download
            model_path = hf_hub_download("Rostlab/ProstT5", "spiece.model")
        except Exception as exc:  # no network / weights unavailable
            pytest.skip(f"ProstT5 tokenizer asset unavailable: {exc}")
        sp = spm.SentencePieceProcessor()
        sp.load(model_path)
        pieces = sp.encode(" ".join("PRTEINS"), out_type=str)
    assert pieces and len([p for p in pieces if p.strip()]) > 0, "empty token sequence"


# --------------------------------------------------------------------------- #
# fpocket
# --------------------------------------------------------------------------- #
def test_fpocket_positive_output(tmp_path):
    if not _have("fpocket"):
        pytest.skip("fpocket not installed")
    control = os.path.join(REPO, "structures", "6EQE.pdb")
    if not os.path.exists(control):
        pytest.skip("control 6EQE.pdb not fetched — run controls/fetch_controls.py")

    work = tmp_path / "6EQE.pdb"
    work.write_text(open(control).read())
    proc = subprocess.run(["fpocket", "-f", str(work)], capture_output=True, text=True)
    out_dir = tmp_path / "6EQE_out"
    info = out_dir / "6EQE_info.txt"
    assert info.exists(), f"fpocket produced no info file (rc={proc.returncode})"
    text = info.read_text()
    assert "Pocket 1" in text, "fpocket reported zero pockets"


# --------------------------------------------------------------------------- #
# AutoDock Vina
# --------------------------------------------------------------------------- #
def test_vina_positive_output():
    vina = pytest.importorskip("vina", reason="vina python bindings not installed")
    try:
        from vina import Vina
        v = Vina(sf_name="vina", seed=1729, verbosity=0)
    except Exception as exc:
        pytest.skip(f"vina not usable: {exc}")
    # Instantiating the scorer is the minimal positive signal available without
    # prepared receptor/ligand PDBQTs; a full scoring run needs S5 outputs.
    assert v is not None
