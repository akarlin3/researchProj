"""Local resume: screen the structures returned from Vast (S3) through S4 + S5.

After the Vast burst folds the S2 shortlist (S3 -> structures + per-model pLDDT),
the pipeline RESUMES locally (vast/sync.md): each kept model is run through the
S4 catalytic-geometry gate and the S5 cleft filter, then scored RELATIVE to the
positive controls — exactly the calibrated separation, now applied to dark
candidates instead of controls.

A candidate is a PETase-LIKE HIT iff all three gates hold:
  1. S4 geometry  — a Ser-His-Asp triad + oxyanion hole passes (triad_found).
  2. S5 pocket    — a catalytic pocket sits within catalytic_pocket_max_dist of
                    the catalytic Ser OG.
  3. cleft score  — its composite (z-scored against the positive-control anchor,
                    same weights/orientation/peripherality mode as calibration)
                    is >= the calibration operating point (lowest positive control).

Input is an S3 output dir: it reads `s3_results.json` (the runner's kept_ids +
mean pLDDT) when present, else falls back to every `*.pdb` in the dir. The control
anchor + threshold are built from the control structures via the calibration
harness, so screening can never drift from the calibrated operating point.

Local usage, from the repo root (after pulling folded models back):
    PYTHONPATH=src python -m proteus.screen \
        --folded structures/folded --struct-dir structures \
        --out data/processed/s4s5_candidates
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

from proteus.calibrate import analyze_controls, score_analysis
from proteus.s4_geometry import analyze_model
from proteus.s5_cleft_filter import analyze_cleft, composite_from_anchor
from proteus.utils import DEFAULT_CONFIG, REPO, load_config


def build_control_anchor(cfg: dict, struct_dir: str) -> dict:
    """Run the calibration harness on the controls and extract the scoring anchor
    + operating point the screen scores candidates against. Reuses the exact S4/S5
    path so the screen's threshold is the calibrated one, not a fresh guess."""
    analysis = analyze_controls(cfg, struct_dir)
    cal = score_analysis(analysis, cfg)
    op = cal.get("operating_point") or {}
    if "threshold" not in op:
        raise RuntimeError(
            "calibration produced no operating point (need >=1 positive control "
            "with a catalytic pocket) — cannot screen candidates without a threshold")
    return {
        "anchor": cal["anchor"],
        "threshold": op["threshold"],
        "mode": cal["mode"],
        "positive_ids": cal["positive_ids"],
        "separated": cal["verdict"].get("separated"),
        "margin": cal["verdict"].get("margin"),
    }


def _candidate_inputs(folded_dir: str) -> list[dict]:
    """Discover models to screen. Prefer the S3 runner's s3_results.json (so only
    pLDDT-kept models are screened, carrying mean_plddt); else every *.pdb."""
    results = os.path.join(folded_dir, "s3_results.json")
    if os.path.exists(results):
        with open(results) as fh:
            summary = json.load(fh)
        out = []
        for rec in summary.get("results", []):
            if not rec.get("kept", True):
                continue  # dropped at the pLDDT gate in S3 — don't screen it
            pdb = os.path.join(folded_dir, rec.get("pdb_path", f"{rec['id']}.pdb"))
            if os.path.exists(pdb):
                out.append({"id": rec["id"], "pdb": pdb,
                            "mean_plddt": rec.get("mean_plddt")})
        return out
    # Fallback: a plain directory of returned PDBs.
    out = []
    for fn in sorted(os.listdir(folded_dir)):
        if fn.endswith(".pdb"):
            out.append({"id": os.path.splitext(fn)[0],
                        "pdb": os.path.join(folded_dir, fn), "mean_plddt": None})
    return out


def screen_model(pdb_path: str, cfg: dict, anchor: dict, threshold: float,
                 mean_plddt=None, cand_id: str | None = None) -> dict:
    """S4 -> S5 -> anchor-scoring for ONE returned model. Returns the candidate
    record with the gate outcomes and (if scorable) the composite + hit verdict."""
    rec = {"id": cand_id or os.path.basename(pdb_path), "mean_plddt": mean_plddt,
           "triad_found": False, "catalytic_ser": None, "his": None, "acid": None,
           "pocket_ok": False, "composite": None, "above_threshold": None,
           "petase_like_hit": False, "metrics": None}

    s4 = analyze_model(pdb_path, cfg)
    rec["triad_found"] = s4["triad_found"]
    best = s4["best"] if s4["triad_found"] else None
    if best is None:
        rec["stage_failed"] = "S4_geometry"  # no passing catalytic triad
        return rec
    rec["catalytic_ser"] = best["ser"]["res_id"]
    rec["his"] = best["his"]["res_id"]
    rec["acid"] = best["acid"]["res_id"]

    s5 = analyze_cleft(pdb_path, rec["catalytic_ser"], cfg)
    rec["pocket_ok"] = s5["pocket_id"] is not None
    if not rec["pocket_ok"]:
        rec["stage_failed"] = "S5_pocket"  # no catalytic pocket near the Ser OG
        return rec

    scored = composite_from_anchor(s5["metrics"], anchor, cfg)
    rec["metrics"] = s5["metrics"]
    rec["z"] = scored["z"]
    rec["composite"] = scored["composite"]
    rec["above_threshold"] = bool(scored["composite"] >= threshold)
    rec["petase_like_hit"] = rec["above_threshold"]  # all 3 gates passed iff here & above
    return rec


def screen_folded(folded_dir: str, cfg: dict, struct_dir: str) -> dict:
    """Screen every kept S3 model through S4 + S5, scored against the control
    anchor. Returns a ranked summary; raises if no inputs / no anchor."""
    cal = build_control_anchor(cfg, struct_dir)
    inputs = _candidate_inputs(folded_dir)
    if not inputs:
        raise RuntimeError(f"no models to screen in {folded_dir} "
                           "(no s3_results.json and no *.pdb)")

    print(f"[screen] anchor={cal['positive_ids']} mode={cal['mode']} "
          f"threshold={cal['threshold']:.4f} "
          f"(controls separated={cal['separated']}, margin={cal['margin']})")

    candidates = []
    for item in inputs:
        rec = screen_model(item["pdb"], cfg, cal["anchor"], cal["threshold"],
                           mean_plddt=item["mean_plddt"], cand_id=item["id"])
        candidates.append(rec)
        comp = "-" if rec["composite"] is None else f"{rec['composite']:.3f}"
        verdict = ("HIT" if rec["petase_like_hit"] else
                   (f"below-thr({comp})" if rec["pocket_ok"] else
                    rec.get("stage_failed", "rejected")))
        print(f"[screen] {rec['id']}: triad={'Y' if rec['triad_found'] else 'N'} "
              f"pocket={'Y' if rec['pocket_ok'] else 'N'} composite={comp} -> {verdict}")

    # rank scorable candidates by composite (desc); unscorable go last
    scorable = [c for c in candidates if c["composite"] is not None]
    scorable.sort(key=lambda c: c["composite"], reverse=True)
    for rank, c in enumerate(scorable, 1):
        c["rank"] = rank
    hits = [c for c in scorable if c["petase_like_hit"]]
    return {
        "folded_dir": folded_dir,
        "threshold": cal["threshold"],
        "anchor_mode": cal["mode"],
        "positive_ids": cal["positive_ids"],
        "n_screened": len(candidates),
        "n_triad": sum(1 for c in candidates if c["triad_found"]),
        "n_pocket": sum(1 for c in candidates if c["pocket_ok"]),
        "n_hits": len(hits),
        "hit_ids": [c["id"] for c in hits],
        "ranking": [c["id"] for c in scorable],
        "candidates": candidates,
    }


def write_outputs(summary: dict, out_prefix: str) -> tuple[str, str]:
    """Write <prefix>.json + <prefix>.csv. Returns the two paths."""
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    js = out_prefix + ".json"
    with open(js, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    cols = ["rank", "id", "mean_plddt", "triad_found", "catalytic_ser", "his", "acid",
            "pocket_ok", "exposure", "druggability", "depth", "volume", "composite",
            "above_threshold", "petase_like_hit"]
    csv_path = out_prefix + ".csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        order = summary["ranking"] + [c["id"] for c in summary["candidates"]
                                      if c["id"] not in summary["ranking"]]
        by_id = {c["id"]: c for c in summary["candidates"]}
        for cid in order:
            c = by_id[cid]
            m = c.get("metrics") or {}
            w.writerow({
                "rank": c.get("rank", ""), "id": c["id"],
                "mean_plddt": c.get("mean_plddt"),
                "triad_found": c["triad_found"], "catalytic_ser": c["catalytic_ser"],
                "his": c["his"], "acid": c["acid"], "pocket_ok": c["pocket_ok"],
                "exposure": m.get("exposure"), "druggability": m.get("druggability"),
                "depth": m.get("depth"), "volume": m.get("volume"),
                "composite": c["composite"], "above_threshold": c["above_threshold"],
                "petase_like_hit": c["petase_like_hit"],
            })
    return js, csv_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folded", default=os.path.join(REPO, "structures", "folded"),
                    help="dir of S3-returned models (s3_results.json or *.pdb)")
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"),
                    help="dir with the control structures (for the scoring anchor)")
    ap.add_argument("--out", default=os.path.join(REPO, "data", "processed", "s4s5_candidates"),
                    help="output path prefix (.json + .csv)")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if not os.path.isdir(args.folded):
        print(f"folded-model dir not found: {args.folded}", file=sys.stderr)
        return 2
    try:
        summary = screen_folded(args.folded, cfg, args.struct_dir)
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    js, csv_path = write_outputs(summary, args.out)
    print(f"[screen] {summary['n_screened']} model(s): "
          f"{summary['n_triad']} triad, {summary['n_pocket']} pocket, "
          f"{summary['n_hits']} PETase-like hit(s) {summary['hit_ids'] or ''}")
    print(f"[screen] candidates -> {os.path.relpath(csv_path, os.getcwd())} ; "
          f"{os.path.relpath(js, os.getcwd())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
