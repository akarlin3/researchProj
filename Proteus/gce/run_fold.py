#!/usr/bin/env python3
"""Proteus S3 — GCE BURST ESMFold runner (runs ON the GCE box, NOT the Mac).

The local pipeline (S0–S2) narrows the corpus and emits, via
`python -m proteus.s3_fold --dry-run`, a job manifest + the S2 shortlist FASTA.
This script is the *contract consumer* on the GCE side (see gce/sync.md): it
loads ESMFold (esmfold_v1) on the configured device (`--device cpu` here, since
this project has no GPU quota; `cuda` when a GPU is available), folds each
shortlisted sequence, writes a PDB + per-residue pLDDT, and keeps only models
whose mean pLDDT >= the manifest's plddt_min. It is intentionally a SINGLE
self-contained file so the fold image only needs to COPY this one runner (plus
ESMFold) — it does NOT import the `proteus` package.

Interruptible-safe: each finished model is written atomically to the output dir
(staged to the GCS bucket), and a completed id is SKIPPED on restart — so a
preempt costs only the single in-flight sequence, never the batch.

The ESMFold model is loaded LAZILY (only when a real fold is requested), so the
orchestration — manifest integrity, resume/checkpoint, pLDDT gating, the run
summary — is exercised by tests on any host with an injected fake folder, exactly
the way the rest of the pipeline keeps its fold step GCE-only.

Usage (inside the GCE instance):
    python3 run_fold.py \
        --manifest /data/proteus/in/s3_job_manifest.json \
        --fasta    /data/proteus/in/s2_shortlist.fasta \
        --out      /data/proteus/out/
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from datetime import datetime, timezone

_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


# --------------------------------------------------------------------------- #
# IO helpers (dependency-free)
# --------------------------------------------------------------------------- #
def parse_fasta(path: str):
    """Yield (id, sequence) pairs. Minimal FASTA reader."""
    rid, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(seq)
                rid = line[1:].split()[0] if len(line) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
    if rid is not None:
        yield rid, "".join(seq)


def _sha256(seq: str) -> str:
    return hashlib.sha256(seq.upper().encode()).hexdigest()


def _atomic_write(path: str, text: str) -> None:
    """Write to a temp file in the same dir then rename — so a reclaim mid-write
    never leaves a half-written file that a restart would mistake for complete."""
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def mean_plddt_from_pdb(pdb_str: str) -> float:
    """Mean pLDDT over CA atoms, read from the PDB B-factor column (ESMFold writes
    pLDDT, 0–100, there). Returns 0.0 if no CA atoms are present."""
    vals = []
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                continue
    return float(sum(vals) / len(vals)) if vals else 0.0


def _max_ca_bfactor(pdb_str: str) -> float:
    """Largest CA B-factor in the PDB (0.0 if none) — used to detect the pLDDT scale."""
    mx = 0.0
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            try:
                mx = max(mx, float(line[60:66]))
            except ValueError:
                continue
    return mx


def _rescale_bfactors(pdb_str: str, factor: float) -> str:
    """Multiply every ATOM/HETATM B-factor (cols 61–66) by `factor`, preserving the
    fixed-width PDB columns. transformers' EsmForProteinFolding.output_to_pdb writes
    pLDDT on a 0–1 scale, but the rest of the toolchain (and the PDB convention ESMFold
    itself uses) expects 0–100, so we lift it back. Only the B-factor field is touched."""
    out = []
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                b = float(line[60:66])
            except ValueError:
                out.append(line)
                continue
            out.append(f"{line[:60]}{b * factor:6.2f}{line[66:]}")
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if pdb_str.endswith("\n") else "")


def normalize_plddt_scale(pdb_str: str) -> str:
    """Ensure the PDB's pLDDT B-factors are on the 0–100 scale. transformers writes
    them 0–1; if the max CA B-factor looks like a probability (<= 1.5) we scale ×100.
    Already-0–100 PDBs (e.g. native ESMFold) are left untouched, so this is idempotent."""
    if 0.0 < _max_ca_bfactor(pdb_str) <= 1.5:
        return _rescale_bfactors(pdb_str, 100.0)
    return pdb_str


# --------------------------------------------------------------------------- #
# ESMFold backend (lazy — imported only when a real fold runs)
# --------------------------------------------------------------------------- #
def load_esmfold(device: str = "cpu", seed: int | None = None):
    """Load ESMFold (facebook/esmfold_v1) via HuggingFace transformers onto `device`
    ('cpu' or 'cuda'). transformers' EsmForProteinFolding is used (not fair-esm) so
    the image builds CPU-only — it needs no OpenFold CUDA kernels. Imports are lazy
    (GCE-only). On CPU it runs fp32: slow + RAM-heavy but fine for the small narrowed
    shortlist; never fold on the Mac's MPS. Returns (model, tokenizer)."""
    import torch  # noqa: PLC0415
    from transformers import AutoTokenizer, EsmForProteinFolding  # noqa: PLC0415
    if seed is not None:
        torch.manual_seed(seed)
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(
            "run_fold.py was asked for --device cuda but this box has no GPU. Use "
            "--device cpu (the default GCE fold here), or launch a GPU VM.")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1",
                                                 low_cpu_mem_usage=True)
    model = model.eval().to(device)
    return model, tokenizer


def make_esmfold_folder(bundle, device: str = "cpu"):
    """Wrap a loaded (model, tokenizer) as a folder callable: (seq, chunk_size,
    max_recycles) -> (pdb_str, mean_plddt). For long sequences we reduce the trunk
    chunk size to fit memory (set_chunk_size trades speed for memory — it does NOT
    split the chain, which would corrupt the fold). transformers writes per-residue
    pLDDT into the PDB B-factor column on a 0–1 scale, so we normalize it back to the
    0–100 convention (normalize_plddt_scale) before mean_plddt_from_pdb reads it."""
    import torch  # noqa: PLC0415
    model, tokenizer = bundle

    def _fold(seq: str, chunk_size: int | None, max_recycles: int | None):
        if chunk_size and len(seq) > int(chunk_size):
            model.trunk.set_chunk_size(int(chunk_size))
        else:
            model.trunk.set_chunk_size(None)
        inputs = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        kw = {} if max_recycles is None else {"num_recycles": int(max_recycles)}
        with torch.no_grad():
            outputs = model(**inputs, **kw)
        pdb_str = normalize_plddt_scale(model.output_to_pdb(outputs)[0])
        return pdb_str, mean_plddt_from_pdb(pdb_str)

    return _fold


# --------------------------------------------------------------------------- #
# Orchestration (host-agnostic; tested with an injected folder)
# --------------------------------------------------------------------------- #
def _result_path(out_dir: str, rid: str) -> str:
    return os.path.join(out_dir, f"{rid}.json")


def _is_complete(out_dir: str, rid: str) -> bool:
    """A sequence is done iff its result JSON AND its PDB both exist."""
    rp = _result_path(out_dir, rid)
    if not os.path.exists(rp):
        return False
    try:
        with open(rp) as fh:
            rec = json.load(fh)
    except (json.JSONDecodeError, OSError):
        return False
    pdb = rec.get("pdb_path")
    return bool(pdb) and os.path.exists(os.path.join(out_dir, os.path.basename(pdb)))


def fold_batch(manifest: dict, fasta_path: str, out_dir: str, folder=None,
               device: str = "cuda", seed: int | None = None) -> dict:
    """Fold every shortlisted sequence with checkpoint/resume + pLDDT gating.

    `folder` is a callable (seq, chunk_size, max_recycles) -> (pdb_str, mean_plddt).
    If None, ESMFold is loaded lazily and wrapped (the real GCE path). Tests pass
    a deterministic fake folder so this logic runs without a GPU.

    Returns a run summary dict (also written to <out_dir>/s3_results.json).
    """
    os.makedirs(out_dir, exist_ok=True)
    fp = manifest.get("fold_params", {})
    plddt_min = float(fp.get("plddt_min", 0.0))
    chunk_size = fp.get("chunk_size")
    max_recycles = fp.get("max_recycles")
    if seed is None:
        seed = manifest.get("random_seed")

    # manifest ids + their recorded sha256 (integrity check against the FASTA)
    man_sha = {e["id"]: e.get("sha256") for e in manifest.get("sequences", [])}
    seqs = dict(parse_fasta(fasta_path))

    built_folder = folder
    results, errors = [], []
    for rid in [e["id"] for e in manifest.get("sequences", [])]:
        if rid not in seqs:
            errors.append(f"{rid}: in manifest but absent from FASTA")
            continue
        seq = seqs[rid].upper()
        if man_sha.get(rid) and _sha256(seq) != man_sha[rid]:
            errors.append(f"{rid}: FASTA sequence sha256 != manifest (input mismatch)")
            continue
        bad = sorted(set(seq) - _AA)
        if bad:
            errors.append(f"{rid}: non-amino-acid chars {bad}")
            continue

        if _is_complete(out_dir, rid):
            with open(_result_path(out_dir, rid)) as fh:
                rec = json.load(fh)
            rec["resumed"] = True
            results.append(rec)
            print(f"[S3][skip] {rid}: already folded (mean_plddt="
                  f"{rec.get('mean_plddt')}, kept={rec.get('kept')})")
            continue

        # Build the real ESMFold folder on first actual need (lazy, GPU-only).
        if built_folder is None:
            print(f"[S3] loading ESMFold on {device} (seed={seed}) …")
            model = load_esmfold(device=device, seed=seed)
            built_folder = make_esmfold_folder(model, device=device)

        chunked = bool(chunk_size and len(seq) > int(chunk_size))
        try:
            pdb_str, mean_plddt = built_folder(seq, chunk_size, max_recycles)
        except Exception as exc:  # noqa: BLE001 — a single fold failure must not kill the batch
            errors.append(f"{rid}: fold failed: {exc}")
            print(f"[S3][FAIL] {rid}: {exc}", file=sys.stderr)
            continue

        kept = mean_plddt >= plddt_min
        pdb_name = f"{rid}.pdb"
        _atomic_write(os.path.join(out_dir, pdb_name), pdb_str)
        rec = {
            "id": rid, "length": len(seq),
            "mean_plddt": round(float(mean_plddt), 3),
            "plddt_min": plddt_min, "kept": bool(kept),
            "chunked": chunked, "num_recycles": max_recycles,
            "pdb_path": pdb_name, "resumed": False,
        }
        _atomic_write(_result_path(out_dir, rid), json.dumps(rec, indent=2) + "\n")
        results.append(rec)
        print(f"[S3][fold] {rid}: mean_plddt={mean_plddt:.1f} "
              f"({'KEEP' if kept else 'drop <'+str(plddt_min)}){' [chunked]' if chunked else ''}")

    kept = [r for r in results if r["kept"]]
    summary = {
        "schema": "proteus.s3_fold.run_results/v1",
        "generated": datetime.now(timezone.utc).isoformat(),
        "stage": "S3_esmfold_batch",
        "device": device,
        "random_seed": seed,
        "plddt_min": plddt_min,
        "n_manifest": len(manifest.get("sequences", [])),
        "n_folded": len(results),
        "n_kept": len(kept),
        "n_dropped": len(results) - len(kept),
        "n_errors": len(errors),
        "kept_ids": [r["id"] for r in kept],
        "dropped_ids": [r["id"] for r in results if not r["kept"]],
        "errors": errors,
        "results": results,
    }
    _atomic_write(os.path.join(out_dir, "s3_results.json"),
                  json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="S3 job manifest (from s3_fold --dry-run)")
    ap.add_argument("--fasta", required=True, help="S2 shortlist FASTA to fold")
    ap.add_argument("--out", required=True, help="output dir (staged to the GCS bucket)")
    ap.add_argument("--device", default=None, help="cpu (default here) | cuda (GPU VM)")
    args = ap.parse_args(argv)

    for p in (args.manifest, args.fasta):
        if not os.path.exists(p):
            print(f"input not found: {p}", file=sys.stderr)
            return 2
    with open(args.manifest) as fh:
        manifest = json.load(fh)

    device = args.device or manifest.get("fold_params", {}).get("device", "cpu")
    try:
        summary = fold_batch(manifest, args.fasta, args.out, device=device)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    print(f"[S3] folded {summary['n_folded']}/{summary['n_manifest']}; "
          f"kept {summary['n_kept']} (mean pLDDT >= {summary['plddt_min']}), "
          f"dropped {summary['n_dropped']}, errors {summary['n_errors']}.")
    print(f"[S3] results -> {os.path.join(args.out, 's3_results.json')}")
    print("[S3] next: rsync the kept PDBs DOWN to structures/folded/ and resume S4/S5 "
          "locally (gce/sync.md).")
    # Non-zero only if nothing got folded at all (e.g. total input mismatch).
    return 0 if summary["n_folded"] > 0 or summary["n_manifest"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
