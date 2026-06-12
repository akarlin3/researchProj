"""Local narrowing pipeline — corpus -> S0 -> S1 -> S2 -> S3 job manifest.

One command for the whole LOCAL (M4) narrowing front: assemble the raw corpus,
dereplicate (S0), tokenise to 3Di (S1), fold-class triage (S2), and emit the S3
fold job manifest for the shortlist. Folding itself (S3) runs on Vast — this
driver stops at the manifest the burst box consumes (see vast/sync.md).

It just chains the existing per-stage functions and reports the NARROWING FUNNEL
(how many sequences survive each gate). Every stage reads thresholds + the seed
from config; nothing here re-implements stage logic.

S2 needs pre-built reference DBs (s2_foldclass_triage.references[].db); pass
`references=` to override (e.g. freshly built ones). The funnel is the deliverable
shipped up to Vast: shortlist FASTA + s3_job_manifest.json.

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.pipeline --out data/interim
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from proteus.corpus import assemble_corpus
from proteus.s0_dereplicate import dereplicate
from proteus.s1_tokenize import tokenize
from proteus.s2_foldclass_triage import (
    _references_from_config,
    triage,
    write_shortlist,
    write_triage_tsv,
)
from proteus.s3_fold import build_manifest, validate_records
from proteus.s3_fold import parse_fasta as _s3_parse
from proteus.utils import DEFAULT_CONFIG, REPO, load_config


def run_local(cfg: dict, references=None, out_dir: str | None = None,
              corpus_glob: str | None = None, prostt5_model: str | None = None) -> dict:
    """Run corpus -> S0 -> S1 -> S2 -> S3 manifest locally. Returns the funnel summary.

    `references`: S2 reference list (see s2_foldclass_triage.triage). Defaults to the
    config's; pass explicit dicts to use freshly built reference DBs.
    """
    out = out_dir or os.path.join(REPO, "data", "interim")
    os.makedirs(out, exist_ok=True)

    # --- corpus front door ------------------------------------------------- #
    corpus_fasta = os.path.join(out, "corpus.fasta")
    cz = assemble_corpus(cfg, corpus_fasta, glob_override=corpus_glob)

    # --- S0 dereplicate ---------------------------------------------------- #
    reps = os.path.join(out, "s0_representatives.fasta")
    clusters = os.path.join(out, "s0_clusters.tsv")
    s0 = dereplicate(corpus_fasta, reps, clusters, cfg)

    # --- S1 tokenize ------------------------------------------------------- #
    s1 = tokenize(reps, os.path.join(out, "s1_3di"), cfg, prostt5_model=prostt5_model)

    # --- S2 fold-class triage --------------------------------------------- #
    s2cfg = cfg["s2_foldclass_triage"]
    refs = references if references is not None else _references_from_config(s2cfg)
    if not refs or all(not r.get("db") for r in refs):
        raise RuntimeError(
            "S2 has no reference DB configured (s2_foldclass_triage.references[].db is "
            "blank) — point it at pre-built curated/broad Foldseek DBs or pass references=")
    shortlist = os.path.join(out, "s2_shortlist.fasta")
    s2 = triage(s1["querydb"], reps, refs, s2cfg)
    write_shortlist(reps, s2["shortlisted"], shortlist)
    write_triage_tsv(s2, os.path.join(out, "s2_triage.tsv"))

    # --- S3 job manifest (Vast hand-off) ---------------------------------- #
    corpus_cfg = cfg.get("corpus", {})
    valid, _errors = validate_records(
        _s3_parse(shortlist), int(corpus_cfg.get("min_length", 0)),
        int(corpus_cfg.get("max_length", 10**9)))
    manifest = build_manifest(shortlist, cfg, valid)
    manifest_path = os.path.join(out, "s3_job_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    funnel = {
        "corpus_read": cz["n_read"],
        "corpus_kept": cz["n_written"],
        "s0_representatives": s0["n_representatives"],
        "s2_shortlisted": s2["n_shortlisted"],
        "s3_manifest_sequences": manifest["n_sequences"],
        "artifacts": {
            "corpus_fasta": corpus_fasta,
            "representatives_fasta": reps,
            "clusters_tsv": clusters,
            "query_db": s1["querydb"],
            "shortlist_fasta": shortlist,
            "s3_manifest": manifest_path,
        },
        "corpus": cz, "s0": {k: s0[k] for k in ("n_input", "n_representatives", "n_clusters")},
        "s1_backend": s1["backend"],
        "s2": {k: s2[k] for k in ("n_input", "n_shortlisted", "dropped")},
    }
    print("[pipeline] NARROWING FUNNEL: "
          f"corpus {funnel['corpus_kept']} (of {funnel['corpus_read']} read) "
          f"-> S0 reps {funnel['s0_representatives']} "
          f"-> S2 shortlist {funnel['s2_shortlisted']} "
          f"-> S3 manifest {funnel['s3_manifest_sequences']} seq")
    print(f"[pipeline] ship to Vast: {os.path.relpath(shortlist, os.getcwd())} + "
          f"{os.path.relpath(manifest_path, os.getcwd())} (see vast/sync.md)")
    return funnel


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(REPO, "data", "interim"),
                    help="output dir for the per-stage artifacts")
    ap.add_argument("--glob", default=None, help="override corpus.fasta_glob")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    try:
        run_local(cfg, out_dir=args.out, corpus_glob=args.glob)
    except (FileNotFoundError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
