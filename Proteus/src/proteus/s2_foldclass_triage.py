"""S2 — Fold-CLASS triage with Foldseek (unseeded).

Match candidates against the alpha/beta-hydrolase fold CLASS, not against
specific PETase templates. The point is to retain anything with the right
*architecture* regardless of sequence-level homology to known enzymes — keeping
the divergent dark tail in play. (The non-PETase serine hydrolases that share the
fold are SUPPOSED to survive S2; they are separated later, at S4/S5.)

BACKEND: Foldseek structural search of the S1 query DB (ProstT5-predicted 3Di)
against pre-built reference 3Di DBs. We search BOTH references and UNION the
survivors (decision recorded in envlog/env-failures.md):

  * CURATED reference (kind: curated) — high-precision alpha/beta-hydrolase
    anchors (e.g. the curated ESTHER/representative set).
  * BROAD reference (kind: broad) — a broad AF-DB/PDB DB that returns hits to
    everything, followed by a fold-class POST-FILTER: only hits whose target is
    on the reference's `fold_class_allowlist` count. This keeps the divergent
    tail (recall) while discarding off-fold matches.

A query representative is SHORTLISTED iff, for at least one reference, it has a
Foldseek hit passing both the e-value and bit-score cutoffs (and, for a broad
reference, to an allow-listed target).

Reads target_fold / evalue / min_bits / references from config
(s2_foldclass_triage). Emits the shortlist FASTA (shipped up to Vast for S3) and
a per-query triage TSV.

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.s2_foldclass_triage \
        --query-db data/interim/s1_3di/querydb \
        --reps     data/interim/s0_representatives.fasta \
        --shortlist data/interim/s2_shortlist.fasta
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_CONFIG = os.path.join(REPO, "config", "proteus.yaml")


def _load_config(path: str) -> dict:
    defaults = {
        "s2_foldclass_triage": {
            "target_fold": "alpha_beta_hydrolase",
            "evalue": 0.01, "min_bits": 50,
            "query_db": os.path.join("data", "interim", "s1_3di", "querydb"),
            "representatives_fasta": os.path.join("data", "interim", "s0_representatives.fasta"),
            "shortlist_fasta": os.path.join("data", "interim", "s2_shortlist.fasta"),
            "triage_tsv": os.path.join("data", "interim", "s2_triage.tsv"),
            "references": [],
        },
    }
    try:
        sys.path.insert(0, os.path.join(REPO, "src"))
        from proteus.utils import load_config  # noqa: PLC0415
        cfg = load_config(path)
        cfg.setdefault("s2_foldclass_triage", {})
        for k, v in defaults["s2_foldclass_triage"].items():
            cfg["s2_foldclass_triage"].setdefault(k, v)
        return cfg
    except Exception:  # noqa: BLE001
        return defaults


def parse_fasta(path: str):
    """Yield (id, sequence) pairs. Minimal dependency-free FASTA reader."""
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


def load_allowlist(path: str) -> set[str]:
    """Read a fold-class allow-list (one target id per line; '#' comments / blanks
    ignored). Empty path -> empty set (caller decides whether that is fatal)."""
    ids: set[str] = set()
    if not path:
        return ids
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                ids.add(line.split()[0])
    return ids


# --------------------------------------------------------------------------- #
# Foldseek helpers
# --------------------------------------------------------------------------- #
def build_reference_db(ref_fasta: str, out_db: str, prostt5_weights: str,
                       foldseek_bin: str = "foldseek") -> str:
    """Build a Foldseek 3Di reference DB from an amino-acid FASTA via ProstT5
    (same CPU path as S1). Production reference DBs are built ONCE; this helper
    prepares them (and the mini-fixture refs in tests). Returns `out_db`."""
    if shutil.which(foldseek_bin) is None:
        raise FileNotFoundError(f"'{foldseek_bin}' not on PATH")
    os.makedirs(os.path.dirname(os.path.abspath(out_db)), exist_ok=True)
    proc = subprocess.run(
        [foldseek_bin, "createdb", ref_fasta, out_db,
         "--prostt5-model", prostt5_weights],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"foldseek createdb (reference) failed "
                           f"(rc={proc.returncode}):\n{proc.stderr[-800:]}")
    return out_db


def _search_reference(query_db: str, ref_db: str, evalue: float,
                      foldseek_bin: str, work: str):
    """Run `foldseek search` + `convertalis`; yield (query, target, evalue, bits)."""
    aln = os.path.join(work, "aln")
    tmp = os.path.join(work, "tmp")
    m8 = os.path.join(work, "aln.m8")
    search = subprocess.run(
        [foldseek_bin, "search", query_db, ref_db, aln, tmp, "-e", str(evalue)],
        capture_output=True, text=True)
    if search.returncode != 0:
        raise RuntimeError(f"foldseek search failed (rc={search.returncode}):\n"
                           f"{search.stderr[-800:]}")
    conv = subprocess.run(
        [foldseek_bin, "convertalis", query_db, ref_db, aln, m8,
         "--format-output", "query,target,evalue,bits"],
        capture_output=True, text=True)
    if conv.returncode != 0:
        raise RuntimeError(f"foldseek convertalis failed (rc={conv.returncode}):\n"
                           f"{conv.stderr[-800:]}")
    if not os.path.exists(m8):
        return
    with open(m8) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            q, t, e, b = line.split("\t")[:4]
            yield q, t, float(e), float(b)


# --------------------------------------------------------------------------- #
# Triage
# --------------------------------------------------------------------------- #
def triage(query_db: str, reps_fasta: str, references: list[dict], cfg_s2: dict,
           foldseek_bin: str = "foldseek", tmp_dir: str | None = None) -> dict:
    """Search the query DB against every reference, threshold + (broad) post-filter
    the hits, and UNION the survivors into a shortlist.

    `references` is a list of dicts: {name, kind: curated|broad, db, allowlist?}.
    Returns a summary with per-query evidence and the shortlisted ids.
    """
    if shutil.which(foldseek_bin) is None:
        raise FileNotFoundError(
            f"'{foldseek_bin}' not on PATH — install Foldseek (see envlog/env-failures.md)")
    if not references:
        raise ValueError("no S2 references configured — set s2_foldclass_triage.references")

    evalue = float(cfg_s2["evalue"])
    min_bits = float(cfg_s2["min_bits"])
    target_fold = cfg_s2.get("target_fold", "alpha_beta_hydrolase")

    rep_ids = [rid for rid, _ in parse_fasta(reps_fasta)]
    # per query: {ref_name: {"best_evalue", "best_bits", "passed": bool}}
    evidence: dict[str, dict] = {rid: {} for rid in rep_ids}

    print(f"[S2] fold-CLASS triage vs '{target_fold}' — UNSEEDED (fold class, not "
          "PETase templates). Non-PETase hydrolases are EXPECTED to pass.")
    print(f"[S2] thresholds: evalue<={evalue} bits>={min_bits}; "
          f"references={[r['name'] for r in references]} (union)")

    cleanup = tmp_dir is None
    base = tmp_dir or tempfile.mkdtemp(prefix="s2_foldseek_")
    try:
        for ref in references:
            kind = ref.get("kind", "curated")
            allow = set()
            if kind == "broad":
                allow = load_allowlist(ref.get("allowlist", "") or "")
                if not allow:
                    print(f"[S2][warn] broad reference '{ref['name']}' has an empty "
                          "fold-class allow-list — its hits cannot be post-filtered "
                          "and are ALL discarded.", file=sys.stderr)
            work = os.path.join(base, ref["name"])
            os.makedirs(work, exist_ok=True)

            # best thresholded (and, for broad, allow-listed) hit per query
            best: dict[str, tuple[float, float]] = {}
            for q, t, e, b in _search_reference(query_db, ref["db"], evalue,
                                                foldseek_bin, work):
                if e > evalue or b < min_bits:
                    continue
                if kind == "broad" and t not in allow:
                    continue  # POST-FILTER: off-fold broad hit -> discard
                cur = best.get(q)
                if cur is None or e < cur[0]:
                    best[q] = (e, b)

            for rid in rep_ids:
                if rid in best:
                    e, b = best[rid]
                    evidence[rid][ref["name"]] = {
                        "best_evalue": e, "best_bits": b, "passed": True}
                else:
                    evidence[rid][ref["name"]] = {
                        "best_evalue": None, "best_bits": None, "passed": False}
            n_pass = sum(1 for rid in rep_ids if evidence[rid][ref["name"]]["passed"])
            print(f"[S2]   {ref['name']} ({kind}): {n_pass}/{len(rep_ids)} pass")
    finally:
        if cleanup:
            shutil.rmtree(base, ignore_errors=True)

    shortlisted = [rid for rid in rep_ids
                   if any(evidence[rid][r["name"]]["passed"] for r in references)]
    dropped = [rid for rid in rep_ids if rid not in shortlisted]
    return {
        "target_fold": target_fold,
        "evalue": evalue, "min_bits": min_bits,
        "references": [r["name"] for r in references],
        "n_input": len(rep_ids),
        "n_shortlisted": len(shortlisted),
        "shortlisted": shortlisted,
        "dropped": dropped,
        "evidence": evidence,
        "rep_ids": rep_ids,
    }


def write_shortlist(reps_fasta: str, shortlisted: list[str], out_fasta: str) -> int:
    """Write the subset of `reps_fasta` whose ids are in `shortlisted`."""
    keep = set(shortlisted)
    os.makedirs(os.path.dirname(os.path.abspath(out_fasta)), exist_ok=True)
    n = 0
    with open(out_fasta, "w") as fh:
        for rid, seq in parse_fasta(reps_fasta):
            if rid in keep:
                fh.write(f">{rid}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
                n += 1
    return n


def write_triage_tsv(summary: dict, out_tsv: str) -> None:
    """Per-query audit: ref-by-ref pass + best e/bits, and the final verdict."""
    refs = summary["references"]
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    header = ["query"]
    for r in refs:
        header += [f"{r}.pass", f"{r}.evalue", f"{r}.bits"]
    header += ["shortlisted"]
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for rid in summary["rep_ids"]:
            row = [rid]
            for r in refs:
                ev = summary["evidence"][rid][r]
                row += [str(ev["passed"]),
                        "" if ev["best_evalue"] is None else f"{ev['best_evalue']:.3e}",
                        "" if ev["best_bits"] is None else f"{ev['best_bits']:.1f}"]
            row.append(str(rid in set(summary["shortlisted"])))
            fh.write("\t".join(row) + "\n")


def _references_from_config(cfg_s2: dict) -> list[dict]:
    """Normalize config references into the triage() shape, resolving paths."""
    out = []
    for ref in cfg_s2.get("references", []) or []:
        entry = {"name": ref.get("name", "ref"), "kind": ref.get("kind", "curated"),
                 "db": ref.get("db", "")}
        pf = ref.get("post_filter") or {}
        entry["allowlist"] = pf.get("fold_class_allowlist", "")
        out.append(entry)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--query-db", default=None, help="S1 Foldseek query DB (overrides config)")
    ap.add_argument("--reps", default=None, help="S0 representative FASTA (overrides config)")
    ap.add_argument("--shortlist", default=None, help="output shortlist FASTA (overrides config)")
    ap.add_argument("--triage-tsv", default=None, help="output triage TSV (overrides config)")
    args = ap.parse_args(argv)

    cfg = _load_config(args.config)
    s2 = cfg["s2_foldclass_triage"]

    def _abs(p):
        return p if os.path.isabs(p) else os.path.join(REPO, p)

    query_db = args.query_db or _abs(s2["query_db"])
    reps = args.reps or _abs(s2["representatives_fasta"])
    shortlist = args.shortlist or _abs(s2["shortlist_fasta"])
    triage_tsv = args.triage_tsv or _abs(s2["triage_tsv"])

    if not os.path.exists(reps):
        print(f"representative FASTA not found: {reps}", file=sys.stderr)
        return 2
    if not (os.path.exists(query_db) or os.path.exists(query_db + ".dbtype")):
        print(f"Foldseek query DB not found: {query_db} (run S1 first)", file=sys.stderr)
        return 2

    references = _references_from_config(s2)
    if not references or all(not r["db"] for r in references):
        print("no S2 reference DB configured (s2_foldclass_triage.references[].db is "
              "blank) — point it at pre-built curated/broad Foldseek DBs.", file=sys.stderr)
        return 2

    try:
        summary = triage(query_db, reps, references, s2)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    n = write_shortlist(reps, summary["shortlisted"], shortlist)
    write_triage_tsv(summary, triage_tsv)
    print(f"[S2] {summary['n_input']} representative(s) -> {summary['n_shortlisted']} "
          f"shortlisted (dropped: {summary['dropped'] or 'none'})")
    print(f"[S2] shortlist FASTA ({n} seq) -> {os.path.relpath(shortlist, os.getcwd())}")
    print(f"[S2] triage audit -> {os.path.relpath(triage_tsv, os.getcwd())}")
    print("[S2] next: ship the shortlist UP to Vast.ai for ESMFold (S3); see vast/sync.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
