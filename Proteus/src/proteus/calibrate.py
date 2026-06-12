"""S4 -> S5 calibration & separation report on the control set (Checkpoint 5).

Runs the catalytic-geometry gate (S4) then the cleft filter (S5) across every control
structure, z-scores the cleft metrics against the positive controls, and asks the
make-or-break question: do the known PET hydrolases (IsPETase, LCC-WT) score above
every non-PET serine hydrolase that shares the fold and triad?

Emits:
  - envlog/calibration-report.md   (table + separation verdict + operating point + LOO)
  - data/processed/s5_scores.csv   (ranked composite table)
  - data/processed/s4_triads.json  (per-model triad detail)

N is tiny (2 positive scaffolds, 4 negatives). This establishes FACE VALIDITY and a
PROVISIONAL operating point only — see the honest-stats section of the report.

Usage:
    PYTHONPATH=src python -m proteus.calibrate
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import shutil

from proteus.s4_geometry import analyze_model
from proteus.s5_cleft_filter import (
    PERIPHERALITY_MODES,
    analyze_cleft,
    composite_from_anchor,
    score_controls,
)
from proteus.utils import DEFAULT_CONFIG, REPO, load_config

CONTROLS_CSV = os.path.join(REPO, "controls", "references.csv")
MANIFEST = os.path.join(REPO, "controls", "MANIFEST.json")
TRAP_IDS = {"LCC_ICCG"}  # the S165A inactivated control — expected to have no triad
# Controls that MUST be present to calibrate at all (Checkpoint 0 precondition).
REQUIRED_STRUCTURES = ("6EQE", "4EB0", "6THS", "1TCA")
# Peripherality forms compared in the normalization section (Checkpoint 2), worst->best
# for size-transfer: absolute (size-dependent) -> rg_norm -> percentile (production).
COMPARE_MODES = ("absolute", "rg_norm", "percentile")


# --------------------------------------------------------------------------- #
# Checkpoint 0 — precondition audit
# --------------------------------------------------------------------------- #
def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def preconditions(cfg: dict, struct_dir: str) -> dict:
    """Enumerate every precondition for calibration. Returns {ok: bool, checks: [...]}.

    Calibration cannot run without (a) a working env (torch/biotite/biopython + fpocket
    on PATH), (b) the control structures present with sha256 matching MANIFEST.json, and
    (c) a config with the s4_geometry/s5_cleft_filter blocks."""
    checks = []

    def add(name, ok, detail):
        checks.append({"name": name, "ok": bool(ok), "detail": detail})

    # (a) environment
    for mod in ("torch", "biotite", "Bio"):
        try:
            __import__(mod)
            add(f"import {mod}", True, "ok")
        except Exception as exc:  # noqa: BLE001
            add(f"import {mod}", False, f"{type(exc).__name__}: {exc}")
    fp = shutil.which("fpocket")
    add("fpocket on PATH", fp is not None, fp or "not found — install fpocket")

    # (b) controls present + sha256 match
    manifest_sha = {}
    if os.path.exists(MANIFEST):
        try:
            man = json.load(open(MANIFEST))
            for s in man.get("structures", []):
                if s.get("sha256"):
                    manifest_sha[s["accession"].upper()] = s["sha256"]
            add("MANIFEST.json loads", True, f"{len(manifest_sha)} sha256 records")
        except Exception as exc:  # noqa: BLE001
            add("MANIFEST.json loads", False, str(exc))
    else:
        add("MANIFEST.json present", False, "run controls/fetch_controls.py")

    for pid in REQUIRED_STRUCTURES:
        path = os.path.join(struct_dir, f"{pid}.pdb")
        if not os.path.exists(path):
            add(f"control {pid}.pdb", False, "missing — run controls/fetch_controls.py")
            continue
        want = manifest_sha.get(pid)
        if not want:
            add(f"control {pid}.pdb", False, "no sha256 in MANIFEST — re-run fetch_controls.py")
            continue
        got = _sha256(path)
        add(f"control {pid}.pdb sha256", got == want,
            "match" if got == want else f"MISMATCH got {got[:12]} want {want[:12]}")

    # (c) config blocks
    add("config random_seed", "random_seed" in cfg, str(cfg.get("random_seed")))
    add("config s4_geometry block", "s4_geometry" in cfg,
        "present" if "s4_geometry" in cfg else "MISSING")
    add("config s5_cleft_filter block", "s5_cleft_filter" in cfg,
        "present" if "s5_cleft_filter" in cfg else "MISSING")

    return {"ok": all(c["ok"] for c in checks), "checks": checks}


def print_precondition_report(audit: dict) -> None:
    print("=== PRECONDITION REPORT (Checkpoint 0) ===")
    for c in audit["checks"]:
        print(f"  [{'PASS' if c['ok'] else 'FAIL'}] {c['name']}: {c['detail']}")
    print(f"  -> {'GO' if audit['ok'] else 'STOP — preconditions unmet'}")


def read_structure_controls(path: str = CONTROLS_CSV) -> list:
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            row = {k: (v or "").strip() for k, v in row.items() if k is not None}
            if row.get("type") == "structure":
                rows.append({"id": row["id"], "accession": row["accession"].upper(),
                             "role": row["role"]})
    return rows


def classify(control: dict, positive_ids) -> str:
    if control["id"] in positive_ids:
        return "positive"
    if control["id"] in TRAP_IDS:
        return "trap"
    if control["role"] == "negative":
        return "negative"
    if control["role"] == "recovery":
        return "recovery"  # held-out divergent positive — never enters the anchor
    return "other"


def _cfg_with_mode(cfg: dict, mode: str | None) -> dict:
    """Return cfg unchanged, or a deep copy with s5_cleft_filter.peripherality_mode set.

    Used to run the same calibration under each peripherality form without mutating the
    caller's config. Weights/orientation/positive set are untouched — only the
    normalization of the exposure metric changes."""
    if mode is None:
        return cfg
    if mode not in PERIPHERALITY_MODES:
        raise ValueError(f"unknown peripherality_mode {mode!r}; "
                         f"expected one of {PERIPHERALITY_MODES}")
    cfg = copy.deepcopy(cfg)
    cfg["s5_cleft_filter"]["peripherality_mode"] = mode
    return cfg


def analyze_controls(cfg: dict, struct_dir: str) -> dict:
    """Run S4 (triad) + S5 (cleft) over every control ONCE — fpocket runs a single time
    per structure. The S5 result carries all three peripherality forms, so the same
    analysis can be scored under any `peripherality_mode` without re-running fpocket.
    This keeps the normalization comparison clean: only the exposure metric changes
    between modes, never the (mildly non-deterministic) fpocket pocket metrics."""
    positive_ids = list(cfg["s5_cleft_filter"]["positive_controls"])
    # Recovery controls are HELD OUT of the anchor: they are divergent positives we
    # test the calibration's generalization against, never fit it to. Excluding them
    # here keeps the production anchor/operating point byte-identical. They are scored
    # separately against the finished anchor by recovery_screen().
    controls = [c for c in read_structure_controls() if c["role"] != "recovery"]

    per_control, s4_detail = {}, {}
    for c in controls:
        cls = classify(c, positive_ids)
        pdb = os.path.join(struct_dir, f"{c['accession']}.pdb")
        rec = {"id": c["id"], "accession": c["accession"], "role": c["role"],
               "class": cls, "present": os.path.exists(pdb)}
        if not rec["present"]:
            per_control[c["id"]] = rec
            continue
        s4 = analyze_model(pdb, cfg)
        s4_detail[c["accession"]] = s4
        rec["triad_found"] = s4["triad_found"]
        rec["catalytic_ser"] = s4["best"]["ser"]["res_id"] if s4["best"] else None
        if s4["triad_found"] and rec["catalytic_ser"] is not None:
            s5 = analyze_cleft(pdb, rec["catalytic_ser"], cfg)
            rec["s5"] = s5
            rec["pocket_ok"] = s5["pocket_id"] is not None
        else:
            rec["pocket_ok"] = False
        per_control[c["id"]] = rec

    scored_ids = [cid for cid, r in per_control.items() if r.get("pocket_ok")]
    return {"per_control": per_control, "scored_ids": scored_ids,
            "s4_detail": s4_detail, "positive_ids": positive_ids}


def _exposure_for_mode(s5: dict, mode: str) -> float:
    """The exposure value a given peripherality mode would use, read from the once-
    computed S5 result. Falls back to the stored metric if peripherality is absent
    (e.g. exposure_metric=og_sasa)."""
    periph = s5.get("peripherality")
    if periph is not None and mode in periph:
        return float(periph[mode])
    return float(s5["metrics"]["exposure"])


def score_analysis(analysis: dict, cfg: dict, mode: str | None = None) -> dict:
    """Score the once-analyzed controls under one peripherality `mode` (default: cfg's).

    Only the exposure metric is swapped to the requested mode; every fpocket-derived
    metric is reused verbatim, so cross-mode differences are due solely to normalization.
    Returns a fresh result dict (the input analysis is not mutated)."""
    if mode is None:
        mode = cfg["s5_cleft_filter"].get("peripherality_mode", "percentile")
    if mode not in PERIPHERALITY_MODES:
        raise ValueError(f"unknown peripherality_mode {mode!r}; "
                         f"expected one of {PERIPHERALITY_MODES}")
    positive_ids = list(analysis["positive_ids"])
    scored_ids = list(analysis["scored_ids"])
    per_control = copy.deepcopy(analysis["per_control"])

    # Swap exposure to this mode; build the scoring inputs from the shared metrics.
    metrics_by_id = {}
    for cid in scored_ids:
        s5 = per_control[cid]["s5"]
        s5["metrics"]["exposure"] = _exposure_for_mode(s5, mode)
        s5["peripherality_mode"] = mode
        metrics_by_id[cid] = s5["metrics"]

    scores = score_controls(metrics_by_id, positive_ids, cfg)
    for cid in scored_ids:
        per_control[cid]["composite"] = scores[cid]["composite"]
        per_control[cid]["z"] = scores[cid]["z"]

    ranking = sorted(scored_ids, key=lambda i: per_control[i]["composite"], reverse=True)
    for rank, cid in enumerate(ranking, 1):
        per_control[cid]["rank"] = rank

    pos = [cid for cid in scored_ids if per_control[cid]["class"] == "positive"]
    neg = [cid for cid in scored_ids if per_control[cid]["class"] == "negative"]
    verdict = _separation(per_control, pos, neg)
    op = _operating_point(per_control, pos, neg, verdict)
    loo = _leave_one_out(metrics_by_id, per_control, positive_ids, neg, cfg)

    return {
        "mode": mode,
        "anchor": scores["_anchor"],
        "positive_ids": positive_ids,
        "per_control": per_control,
        "ranking": ranking,
        "verdict": verdict,
        "operating_point": op,
        "loo": loo,
        "s4_detail": analysis["s4_detail"],
        "trap": {cid: r for cid, r in per_control.items() if r.get("class") == "trap"},
    }


def run_calibration(cfg: dict, struct_dir: str, mode: str | None = None) -> dict:
    """Convenience wrapper: analyze the controls once, then score in one mode.

    For the multi-mode report, call analyze_controls() once and score_analysis() per
    mode instead (avoids re-running fpocket and isolates the exposure change)."""
    cfg = _cfg_with_mode(cfg, mode)
    analysis = analyze_controls(cfg, struct_dir)
    return score_analysis(analysis, cfg)


def recovery_screen(cfg: dict, struct_dir: str, cal: dict) -> dict:
    """Score the HELD-OUT divergent-positive recovery controls against the finished
    production anchor (IsPETase/LCC). The anchor is NOT rebuilt — these structures
    test whether the calibrated separation generalizes to sequence-divergent PETases.

    For each recovery structure: run S4 -> S5, score with composite_from_anchor, and
    record whether it clears (a) every negative and (b) the production operating
    point. Then propose a WIDENED operating point that would also keep the recovered
    divergent positives, reporting its precision against the negatives.
    """
    anchor = cal["anchor"]
    mode = cal["mode"]
    prod_thr = cal.get("operating_point", {}).get("threshold")
    per = cal["per_control"]
    neg_comp = {cid: per[cid]["composite"] for cid in per
                if per[cid].get("class") == "negative" and "composite" in per[cid]}
    max_neg = max(neg_comp.values(), default=float("-inf"))
    pos_comp = [per[cid]["composite"] for cid in per
                if per[cid].get("class") == "positive" and "composite" in per[cid]]

    recs = []
    for c in read_structure_controls():
        if c["role"] != "recovery":
            continue
        pdb = os.path.join(struct_dir, f"{c['accession']}.pdb")
        rec = {"id": c["id"], "accession": c["accession"], "present": os.path.exists(pdb),
               "triad_found": None, "catalytic_ser": None, "pocket_ok": None,
               "composite": None, "above_all_negatives": None,
               "above_production_line": None, "status": "missing"}
        if not rec["present"]:
            recs.append(rec)
            continue
        s4 = analyze_model(pdb, cfg)
        rec["triad_found"] = s4["triad_found"]
        rec["catalytic_ser"] = s4["best"]["ser"]["res_id"] if s4["best"] else None
        if not (s4["triad_found"] and rec["catalytic_ser"] is not None):
            rec["status"] = "no_triad"
            recs.append(rec)
            continue
        s5 = analyze_cleft(pdb, rec["catalytic_ser"], cfg)
        rec["pocket_ok"] = s5["pocket_id"] is not None
        if not rec["pocket_ok"]:
            rec["status"] = "no_pocket"
            recs.append(rec)
            continue
        s5["metrics"]["exposure"] = _exposure_for_mode(s5, mode)
        comp = composite_from_anchor(s5["metrics"], anchor, cfg)["composite"]
        rec["composite"] = comp
        rec["above_all_negatives"] = bool(comp > max_neg)
        rec["above_production_line"] = (None if prod_thr is None else bool(comp >= prod_thr))
        rec["status"] = ("RECOVERED" if rec["above_production_line"] else
                         "above_negatives_below_line" if rec["above_all_negatives"] else
                         "MISSED")
        recs.append(rec)

    # Widened operating point: lower the line to also keep every recovered divergent
    # positive that already clears the negatives; report its precision vs the negatives.
    widened = None
    recovered = [r for r in recs if r["composite"] is not None and r["above_all_negatives"]]
    scored = [r for r in recs if r["pocket_ok"]]
    if pos_comp and recovered:
        thr = min(min(pos_comp), min(r["composite"] for r in recovered))
        fp = sorted(cid for cid, v in neg_comp.items() if v >= thr)
        kept_pos = len(pos_comp) + len(recovered)
        widened = {
            "threshold": round(thr, 4),
            "includes_recovery": [r["id"] for r in recovered],
            "false_positives": len(fp),
            "negatives_above_line": fp,
            "precision": round(kept_pos / (kept_pos + len(fp)), 4),
            "divergent_recall": round(len(recovered) / len(scored), 4) if scored else None,
        }
    return {
        "recovery": recs,
        "production_threshold": prod_thr,
        "max_negative": round(max_neg, 4) if neg_comp else None,
        "widened_operating_point": widened,
    }


def _separation(per_control, pos, neg) -> dict:
    if not pos or not neg:
        return {"separated": None, "reason": "need >=1 positive and >=1 negative scored"}
    min_pos = min(per_control[c]["composite"] for c in pos)
    max_neg = max(per_control[c]["composite"] for c in neg)
    lowest_pos = min(pos, key=lambda c: per_control[c]["composite"])
    highest_neg = max(neg, key=lambda c: per_control[c]["composite"])
    return {
        "separated": bool(min_pos > max_neg),
        "min_positive": min_pos, "max_negative": max_neg,
        "margin": round(min_pos - max_neg, 4),
        "lowest_positive": lowest_pos, "highest_negative": highest_neg,
    }


def _operating_point(per_control, pos, neg, verdict) -> dict:
    if not pos:
        return {}
    thr = min(per_control[c]["composite"] for c in pos)  # recall=1.0 on positives
    pos_kept = [c for c in pos if per_control[c]["composite"] >= thr]
    neg_kept = [c for c in neg if per_control[c]["composite"] >= thr]
    tp, fp = len(pos_kept), len(neg_kept)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    return {
        "threshold": thr, "recall_positives": 1.0,
        "precision": round(precision, 4),
        "true_positives": tp, "false_positives": fp,
        "negatives_above_line": sorted(neg_kept),
        "negative_positions": {c: round(per_control[c]["composite"] - thr, 4) for c in neg},
    }


def _leave_one_out(metrics_by_id, per_control, positive_ids, neg, cfg) -> list:
    out = []
    pos_scored = [p for p in positive_ids if p in metrics_by_id]
    for dropped in pos_scored:
        remaining = [p for p in positive_ids if p != dropped]
        re_scores = score_controls(metrics_by_id, remaining, cfg)
        dropped_comp = re_scores[dropped]["composite"]
        max_neg = max((re_scores[c]["composite"] for c in neg), default=float("-inf"))
        out.append({
            "dropped": dropped, "anchored_on": remaining,
            "dropped_composite": dropped_comp,
            "max_negative_composite": round(max_neg, 4),
            "still_above_negatives": bool(dropped_comp > max_neg),
        })
    return out


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_csv(result: dict, path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = ["rank", "id", "accession", "class", "triad_found", "catalytic_ser",
            "peripherality_mode", "exposure", "aromatics", "druggability", "depth",
            "volume", "hydrophobicity", "polarity", "raw_og_sasa", "composite"]
    rows = []
    for cid in result["ranking"]:
        r = result["per_control"][cid]
        m = r["s5"]["metrics"]
        rows.append({
            "rank": r.get("rank"), "id": cid, "accession": r["accession"],
            "class": r["class"], "triad_found": r["triad_found"],
            "catalytic_ser": r["catalytic_ser"],
            "peripherality_mode": r["s5"].get("peripherality_mode"),
            "exposure": m["exposure"],
            "aromatics": int(m["aromatics"]), "druggability": m["druggability"],
            "depth": m["depth"], "volume": m["volume"],
            "hydrophobicity": m["hydrophobicity"], "polarity": m["polarity"],
            "raw_og_sasa": round(r["s5"]["raw_og_sasa"], 3), "composite": r["composite"],
        })
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def _fmt(x, nd=3):
    return f"{x:.{nd}f}" if isinstance(x, (int, float)) else str(x)


def write_report(result: dict, path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    pc = result["per_control"]
    v = result["verdict"]
    op = result["operating_point"]
    L = []
    L.append("# Proteus S4+S5 calibration & separation report")
    L.append("")
    L.append("Generated by `python -m proteus.calibrate`. S4 = catalytic-geometry gate "
             "(Ser-His-Asp triad + oxyanion hole); S5 = cleft filter (fpocket base "
             "metrics + active-site exposure + aromatic subsites), scored relative to "
             "the positive controls.")
    L.append("")
    L.append(f"**Positive anchor set:** {', '.join(result['positive_ids'])}  "
             f"(IsPETase=6EQE, LCC-WT=4EB0)")
    L.append("")
    L.append(f"**Peripherality mode (production):** `{result['mode']}` — a size-invariant "
             "form chosen a priori so the exposure metric transfers across the "
             "size-heterogeneous dark proteome (see the normalization section below). "
             "The `exposure` column and composite below are computed in this mode.")
    L.append("")

    # main table
    L.append("## Per-control results")
    L.append("")
    _expl = {"percentile": "fraction of residues more central than the Ser, 0-1",
             "rg_norm": "OG->centroid / Rg, dimensionless",
             "absolute": "OG->centroid, Angstrom"}.get(result["mode"], result["mode"])
    L.append(f"`exposure` is the catalytic-Ser peripherality in the production "
             f"`{result['mode']}` mode ({_expl}).")
    L.append("")
    L.append("| rank | control | acc | class | triad | cat.Ser | exposure | aromatics "
             "| drugg. | depth | volume | composite |")
    L.append("|---:|---|---|---|:--:|---:|---:|---:|---:|---:|---:|---:|")
    rows = result["ranking"] + [cid for cid in pc
                                if cid not in result["ranking"]]
    for cid in rows:
        r = pc[cid]
        if r.get("pocket_ok"):
            m = r["s5"]["metrics"]
            L.append(f"| {r.get('rank','')} | {cid} | {r['accession']} | {r['class']} "
                     f"| {'Y' if r['triad_found'] else 'N'} | {r['catalytic_ser']} "
                     f"| {_fmt(m['exposure'],2)} | {int(m['aromatics'])} "
                     f"| {_fmt(m['druggability'],3)} | {_fmt(m['depth'],3)} "
                     f"| {_fmt(m['volume'],0)} | {_fmt(r['composite'],3)} |")
        else:
            tf = r.get("triad_found")
            tf = ("Y" if tf else "N") if tf is not None else "-"
            note = "no triad (S165A trap)" if r.get("class") == "trap" else \
                   ("no catalytic pocket" if r.get("present") else "structure absent")
            L.append(f"| - | {cid} | {r['accession']} | {r['class']} | {tf} "
                     f"| {r.get('catalytic_ser','-')} | — | — | — | — | — | {note} |")
    L.append("")

    # raw OG SASA diagnostic + the three peripherality forms side by side
    L.append("### Why exposure is peripherality, not raw OG SASA")
    L.append("")
    L.append("Raw catalytic-Ser OG SASA does **not** separate the classes — it is "
             "dominated by whether the crystal caught the lid open or closed. "
             "Peripherality does, in every form (absolute Angstrom, Rg-normalized, and "
             "residue percentile):")
    L.append("")
    L.append("| control | class | raw OG SASA (A^2) | periph abs (A) | periph rg_norm "
             "| periph percentile |")
    L.append("|---|---|---:|---:|---:|---:|")
    for cid in rows:
        r = pc[cid]
        if r.get("pocket_ok"):
            p = r["s5"].get("peripherality") or {}
            L.append(f"| {cid} | {r['class']} | {_fmt(r['s5']['raw_og_sasa'],2)} "
                     f"| {_fmt(p.get('absolute'),2)} | {_fmt(p.get('rg_norm'),3)} "
                     f"| {_fmt(p.get('percentile'),3)} |")
    L.append("")

    # separation verdict
    L.append("## Separation verdict")
    L.append("")
    if v.get("separated") is None:
        L.append(f"Cannot judge separation: {v.get('reason')}")
    elif v["separated"]:
        L.append(f"**PASS — positives separate from all negatives.** The lowest positive "
                 f"({v['lowest_positive']}, composite {_fmt(v['min_positive'],3)}) scores "
                 f"above the highest negative ({v['highest_negative']}, composite "
                 f"{_fmt(v['max_negative'],3)}).")
        L.append("")
        L.append(f"**Margin = {_fmt(v['margin'],3)}** composite units between the lowest "
                 f"positive and the highest negative.")
    else:
        L.append(f"**FAIL — STOP AND RETHINK.** The lowest positive "
                 f"({v['lowest_positive']}, {_fmt(v['min_positive'],3)}) does NOT score "
                 f"above the highest negative ({v['highest_negative']}, "
                 f"{_fmt(v['max_negative'],3)}). Margin = {_fmt(v['margin'],3)}. The cleft "
                 f"metrics as configured cannot tell PET hydrolases from these decoys. "
                 f"Do not proceed downstream until the filter or the control set is fixed.")
    L.append("")

    # operating point
    if op:
        L.append("## Provisional operating point (recall = 1.0 on positives)")
        L.append("")
        L.append(f"Threshold set at the lowest positive composite = "
                 f"**{_fmt(op['threshold'],4)}** (keeps every known PETase by construction).")
        L.append("")
        L.append(f"- Recall on positives: {op['recall_positives']:.2f} "
                 f"({op['true_positives']}/{op['true_positives']})")
        _kept = op["true_positives"] + op["false_positives"]
        L.append(f"- Precision at this line: **{_fmt(op['precision'],3)}** "
                 f"({op['true_positives']} TP / {_kept} kept)")
        L.append(f"- Negatives above the line (false positives): "
                 f"{op['negatives_above_line'] or 'none'}")
        L.append("")
        L.append("Each negative's composite relative to the line (negative = correctly "
                 "below the threshold):")
        L.append("")
        L.append("| negative | composite - threshold |")
        L.append("|---|---:|")
        for cid, delta in sorted(op["negative_positions"].items(),
                                 key=lambda kv: kv[1], reverse=True):
            L.append(f"| {cid} | {_fmt(delta,3)} |")
        L.append("")

    # leave-one-out
    L.append("## Leave-one-out face validity")
    L.append("")
    L.append("Drop each positive, re-anchor on the remaining positive(s), and check the "
             "dropped positive still scores above every negative:")
    L.append("")
    L.append("| dropped positive | anchored on | dropped composite | max negative | still above? |")
    L.append("|---|---|---:|---:|:--:|")
    for lo in result["loo"]:
        L.append(f"| {lo['dropped']} | {', '.join(lo['anchored_on']) or '(none)'} "
                 f"| {_fmt(lo['dropped_composite'],3)} | {_fmt(lo['max_negative_composite'],3)} "
                 f"| {'YES' if lo['still_above_negatives'] else 'NO'} |")
    L.append("")

    # trap
    L.append("## 6THS S165A trap (expected null)")
    L.append("")
    for cid, r in result["trap"].items():
        tf = r.get("triad_found")
        _verdict = ("EXPECTED NULL OK (catalytic Ser is mutated to Ala; no triad to find)"
                    if tf is False else
                    "UNEXPECTED — a triad was detected in the inactivated mutant")
        L.append(f"- **{cid} ({r['accession']})**: triad_found = {tf} -> {_verdict}")
    L.append("")

    # honest stats
    L.append("## Honest statistics — read this before trusting the threshold")
    L.append("")
    L.append("- **N is tiny:** 2 distinct positive scaffolds (IsPETase, LCC) and 4 "
             "negatives (CalB, AChE, *C. rugosa* lipase, Est2). This is a FACE-VALIDITY "
             "check and a PROVISIONAL operating point — **not** a trustworthy threshold.")
    L.append("- With only 2 positives the per-metric standard deviation is unreliable, so "
             "z-score scales are floored by half the overall control spread "
             "(`_robust_scale`). Treat z-magnitudes as indicative, not calibrated.")
    L.append("- The separation here is carried by **exposure (catalytic-Ser "
             "peripherality)**; druggability and volume separate the big-pocket lipases "
             "but miss the compact esterase (Est2). Raw OG SASA and raw aromatic counts "
             "do **not** separate on this panel (aromatic-gorge esterases inflate the "
             "aromatic count), which is why aromatics is down-weighted.")
    L.append(f"- **Normalization fixes generalization, not N.** Switching exposure to the "
             f"size-invariant `{result['mode']}` form removes the size confound so the "
             "metric can transfer to the size-heterogeneous dark proteome — it does "
             "**not** enlarge the control set. The face-validity caveat (2 positives, 4 "
             "negatives) is unchanged; this is still a provisional operating point, not a "
             "validated threshold. The production mode was chosen a priori on principle "
             "(size-invariance), NOT because it maximised the margin on these 6 controls "
             "(see the normalization section).")
    L.append("- Metric orientations and weights were set from biology a priori, not tuned "
             "to maximise separation on these 6 structures; still, a 6-point fit cannot "
             "be considered validated.")
    L.append("- **Next step:** fold the divergent positives GuaPA and MG8 on the Vast.ai "
             "burst box (S3) to add real, sequence-divergent PETases to the positive set, "
             "expand the decoy panel (more lipase/esterase/cutinase-adjacent folds), then "
             "re-calibrate. Re-check the aromatic-subsite metric on the *folded* models "
             "(rotamer noise) after Chai-1 refinement.")
    L.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


# --------------------------------------------------------------------------- #
# Normalization (size-invariance) comparison — Checkpoints 2 & 3
# --------------------------------------------------------------------------- #
_MODE_DESC = {
    "absolute": "Ser OG -> CA centroid, Angstrom. **Size-DEPENDENT** — kept for "
                "comparison only; cannot transfer across proteins of different size.",
    "rg_norm": "(Ser OG -> CA centroid) / Rg. Dimensionless, size-invariant; the "
               "continuous cross-check.",
    "percentile": "fraction of CA atoms closer to the centroid than the catalytic Ser "
                  "is (0 = core, 1 = surface). Size- AND shape-invariant; the principled "
                  "production choice.",
}


def _mode_summary(result: dict) -> dict:
    """Pull the separation/precision/LOO headline numbers out of a per-mode result."""
    v = result["verdict"]
    op = result.get("operating_point") or {}
    loo = result.get("loo") or []
    return {
        "mode": result["mode"],
        "separated": v.get("separated"),
        "margin": v.get("margin"),
        "precision": op.get("precision"),
        "loo_pass": all(lo["still_above_negatives"] for lo in loo) if loo else None,
    }


def append_normalization_section(results: dict, production_mode: str, path: str) -> None:
    """Append the three-way normalization comparison (Checkpoint 2) + the a-priori mode
    choice and honesty guard (Checkpoint 3) to the report at `path`."""
    L = ["", "## Normalization (size-invariance) comparison", ""]
    L.append("The S5 exposure metric is the catalytic-Ser peripherality. The absolute "
             "Angstrom form separated the controls but is **confounded with protein "
             "size** and cannot transfer to the size-heterogeneous dark proteome. Below "
             "the full S4->S5 calibration is re-run in all three forms; weights, "
             "orientation, and the positive anchor set are identical across them — only "
             "the normalization of exposure changes.")
    L.append("")
    for m in COMPARE_MODES:
        L.append(f"- **{m}**: {_MODE_DESC[m]}")
    L.append("")

    # summary across modes
    L.append("### Summary")
    L.append("")
    L.append("| mode | size-invariant? | positives separate? | margin | "
             "precision @ recall=1.0 | leave-one-out | production |")
    L.append("|---|:--:|:--:|---:|---:|:--:|:--:|")
    for m in COMPARE_MODES:
        s = _mode_summary(results[m])
        inv = "no" if m == "absolute" else "yes"
        sep = "PASS" if s["separated"] else ("FAIL" if s["separated"] is False else "?")
        loo = ("pass" if s["loo_pass"] else "FAIL") if s["loo_pass"] is not None else "-"
        prod = "**<- production**" if m == production_mode else ""
        L.append(f"| `{m}` | {inv} | {sep} | {_fmt(s['margin'],3)} "
                 f"| {_fmt(s['precision'],3)} | {loo} | {prod} |")
    L.append("")

    # per-mode detail tables (same peripherality columns; composite/rank are per mode)
    for m in COMPARE_MODES:
        res = results[m]
        pc = res["per_control"]
        L.append(f"### Mode: `{m}`"
                 + ("  (production)" if m == production_mode else ""))
        L.append("")
        L.append("| control | role | peripherality(abs A) | peripherality(rg_norm) "
                 "| peripherality(percentile) | composite | rank |")
        L.append("|---|---|---:|---:|---:|---:|---:|")
        for cid in res["ranking"]:
            r = pc[cid]
            p = r["s5"].get("peripherality") or {}
            L.append(f"| {cid} | {r['class']} | {_fmt(p.get('absolute'),2)} "
                     f"| {_fmt(p.get('rg_norm'),3)} | {_fmt(p.get('percentile'),3)} "
                     f"| {_fmt(r['composite'],3)} | {r.get('rank','')} |")
        v = res["verdict"]
        op = res.get("operating_point") or {}
        if v.get("separated") is None:
            verdict_line = f"indeterminate ({v.get('reason')})"
        elif v["separated"]:
            verdict_line = (f"**PASS** — lowest positive ({v['lowest_positive']}, "
                            f"{_fmt(v['min_positive'],3)}) > highest negative "
                            f"({v['highest_negative']}, {_fmt(v['max_negative'],3)}); "
                            f"margin = **{_fmt(v['margin'],3)}**, precision@recall=1.0 = "
                            f"{_fmt(op.get('precision'),3)}.")
        else:
            verdict_line = (f"**FAIL** — lowest positive ({v['lowest_positive']}, "
                            f"{_fmt(v['min_positive'],3)}) does NOT exceed highest "
                            f"negative ({v['highest_negative']}, "
                            f"{_fmt(v['max_negative'],3)}); margin = {_fmt(v['margin'],3)}.")
        L.append("")
        L.append(f"Separation: {verdict_line}")
        loo = res.get("loo") or []
        if loo:
            loo_bits = "; ".join(
                f"drop {lo['dropped']} -> {_fmt(lo['dropped_composite'],3)} "
                f"(max neg {_fmt(lo['max_negative_composite'],3)}) "
                f"{'YES' if lo['still_above_negatives'] else 'NO'}"
                for lo in loo)
            L.append("")
            L.append(f"Leave-one-out: {loo_bits}")
        L.append("")

    # a-priori choice + honesty guard (Checkpoint 3)
    inv_modes = [m for m in COMPARE_MODES if m != "absolute"]
    inv_ok = all(results[m]["verdict"].get("separated") for m in inv_modes)
    abs_ok = results["absolute"]["verdict"].get("separated")
    L.append("### Mode choice (a priori) and honesty guard")
    L.append("")
    L.append(f"**Production mode = `{production_mode}`, chosen a priori on principle, "
             "NOT by which best separates the six controls.** The dark proteome spans all "
             "sizes; an absolute Angstrom distance cannot be compared across proteins of "
             "different size, so the production metric MUST be size-invariant. "
             "`percentile` is preferred for shape-awareness (it ranks the Ser against its "
             "own protein's residues); `rg_norm` is the continuous cross-check. Picking "
             "the form with the largest margin would be tuning to N=6.")
    L.append("")
    if inv_ok:
        L.append("**Guard check — separation survives normalization.** Both size-invariant "
                 "forms still rank every positive above every negative, so the original "
                 "separation was carried by burial/peripherality, not by size. Production "
                 "is set to the size-invariant mode (as it would be regardless).")
    elif abs_ok and not inv_ok:
        L.append("**Guard check — STOP AND RETHINK.** The size-invariant forms FAIL to "
                 "separate while `absolute` passes. Per the honesty guard this is the "
                 "finding: the original separation was substantially **size-driven**, not "
                 "burial-driven. Production is set to the size-invariant mode ANYWAY (we "
                 "do not revert to `absolute` to preserve a size-confounded win); S5 needs "
                 "a rethink before production use on the dark proteome.")
    else:
        L.append("**Guard check.** Neither `absolute` nor the size-invariant forms cleanly "
                 "separate on this panel — the discriminator needs rework before "
                 "production use.")
    L.append("")
    with open(path, "a") as fh:
        fh.write("\n".join(L) + "\n")


def write_normalization_csv(results: dict, path: str) -> None:
    """Machine-readable three-way comparison: one row per (mode, scored control)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = ["mode", "control", "role", "periph_abs", "periph_rg_norm",
            "periph_percentile", "exposure_active", "composite", "rank"]
    rows = []
    for m in COMPARE_MODES:
        res = results[m]
        pc = res["per_control"]
        for cid in res["ranking"]:
            r = pc[cid]
            p = r["s5"].get("peripherality") or {}
            rows.append({
                "mode": m, "control": cid, "role": r["class"],
                "periph_abs": p.get("absolute"), "periph_rg_norm": p.get("rg_norm"),
                "periph_percentile": p.get("percentile"),
                "exposure_active": r["s5"]["metrics"]["exposure"],
                "composite": r["composite"], "rank": r.get("rank"),
            })
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def append_recovery_section(recov: dict, path: str) -> None:
    """Append the held-out divergent-positive recovery section to the report."""
    present = [r for r in recov["recovery"] if r["present"]]
    L = ["", "## Divergent-positive recovery (held out of the anchor)", ""]
    L.append("Sequence-divergent PETase **structures** scored against the finished "
             "IsPETase/LCC anchor — they never enter the anchor, so this measures "
             "GENERALIZATION, not fit. (GuaPA and MG8 remain unresolved — no reachable "
             "sequence/structure — so the archaeal PET46/8B4U structure stands in.)")
    L.append("")
    if not present:
        L.append("_No recovery structures present — fetch via controls/fetch_controls.py._")
    else:
        L.append(f"Production operating point = {_fmt(recov['production_threshold'],3)}; "
                 f"highest negative = {_fmt(recov['max_negative'],3)}.")
        L.append("")
        L.append("| id | acc | triad | cat.Ser | composite | > all negatives | "
                 ">= production line | status |")
        L.append("|---|---|:--:|---:|---:|:--:|:--:|---|")
        for r in present:
            comp = "—" if r["composite"] is None else _fmt(r["composite"], 3)
            L.append(f"| {r['id']} | {r['accession']} "
                     f"| {'Y' if r['triad_found'] else 'N'} | {r['catalytic_ser'] or '—'} "
                     f"| {comp} | {'Y' if r['above_all_negatives'] else 'N'} "
                     f"| {'Y' if r['above_production_line'] else 'N'} | {r['status']} |")
        L.append("")
        w = recov["widened_operating_point"]
        if w:
            L.append(f"**Recommended widened operating point = {_fmt(w['threshold'],3)}** — "
                     f"lowers the line to also keep the recovered divergent positive(s) "
                     f"{w['includes_recovery']} while holding precision "
                     f"**{_fmt(w['precision'],3)}** ({w['false_positives']} false positive(s)). "
                     "The production line stays IsPETase/LCC-anchored; this is the "
                     "next-gen threshold once more divergent positives are folded.")
        else:
            L.append("_No widened line proposed (no recovery structure cleared the "
                     "negatives, or no positives scored)._")
    L.append("")
    with open(path, "a") as fh:
        fh.write("\n".join(L) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"))
    ap.add_argument("--report", default=os.path.join(REPO, "envlog", "calibration-report.md"))
    ap.add_argument("--csv", default=os.path.join(REPO, "data", "processed", "s5_scores.csv"))
    _norm_default = os.path.join(REPO, "data", "processed", "s5_normalization_comparison.csv")
    ap.add_argument("--norm-csv", default=_norm_default)
    ap.add_argument("--s4-json", default=os.path.join(REPO, "data", "processed", "s4_triads.json"))
    ap.add_argument("--check-only", action="store_true",
                    help="run only the Checkpoint 0 precondition audit and exit")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    production_mode = cfg["s5_cleft_filter"].get("peripherality_mode", "percentile")

    # Checkpoint 0 — enumerate preconditions; STOP if any fail.
    audit = preconditions(cfg, args.struct_dir)
    print_precondition_report(audit)
    if not audit["ok"]:
        print("Aborting: fix the FAILED preconditions above, then re-run.")
        return 2
    if args.check_only:
        return 0
    print()

    # Checkpoint 2 — analyze the controls ONCE (fpocket per structure runs a single
    # time), then score under every peripherality form. Sharing the analysis isolates
    # the exposure-normalization change from fpocket run-to-run noise. Production output
    # reflects `production_mode`.
    analysis = analyze_controls(cfg, args.struct_dir)
    results = {m: score_analysis(analysis, cfg, mode=m) for m in COMPARE_MODES}
    result = results[production_mode]

    os.makedirs(os.path.dirname(os.path.abspath(args.s4_json)), exist_ok=True)
    with open(args.s4_json, "w") as fh:
        json.dump(result["s4_detail"], fh, indent=2)
        fh.write("\n")
    write_csv(result, args.csv)
    write_report(result, args.report)                       # production body
    append_normalization_section(results, production_mode, args.report)  # Checkpoints 2 & 3
    write_normalization_csv(results, args.norm_csv)

    print(f"=== S4 -> S5 CALIBRATION (production mode = {production_mode}) ===")
    for cid in result["ranking"]:
        r = result["per_control"][cid]
        print(f"  #{r['rank']} {cid:9s} [{r['class']:8s}] composite={r['composite']:+.3f} "
              f"exposure={r['s5']['metrics']['exposure']:.3f}")

    print("--- normalization comparison (size-invariance) ---")
    for m in COMPARE_MODES:
        s = _mode_summary(results[m])
        sep = "PASS" if s["separated"] else ("FAIL" if s["separated"] is False else "?")
        tag = "  <- production" if m == production_mode else ""
        print(f"  {m:11s} separated={sep:4s} margin={_fmt(s['margin'],3):>7s} "
              f"precision={_fmt(s['precision'],3)} "
              f"loo={'pass' if s['loo_pass'] else 'FAIL' if s['loo_pass'] is not None else '-'}"
              f"{tag}")

    # Held-out divergent-positive recovery (e.g. archaeal PET46/8B4U). Scored against
    # the production anchor; never used to fit it. Reports a recommended widened line.
    recov = recovery_screen(cfg, args.struct_dir, result)
    present_recov = [r for r in recov["recovery"] if r["present"]]
    if present_recov:
        print("--- divergent-positive recovery (held out of the anchor) ---")
        for r in present_recov:
            comp = "n/a" if r["composite"] is None else f"{r['composite']:+.3f}"
            print(f"  {r['id']:9s} ({r['accession']}): {r['status']:26s} "
                  f"composite={comp} (production line={_fmt(recov['production_threshold'],3)}, "
                  f"max negative={_fmt(recov['max_negative'],3)})")
        w = recov["widened_operating_point"]
        if w:
            print(f"  WIDENED operating point = {_fmt(w['threshold'],3)} "
                  f"(keeps {w['includes_recovery']}): precision={_fmt(w['precision'],3)}, "
                  f"false_positives={w['false_positives']}")
        append_recovery_section(recov, args.report)
        with open(os.path.join(os.path.dirname(os.path.abspath(args.csv)),
                               "s5_recovery.json"), "w") as fh:
            json.dump(recov, fh, indent=2)
            fh.write("\n")

    inv_ok = all(results[m]["verdict"].get("separated") for m in COMPARE_MODES if m != "absolute")
    abs_ok = results["absolute"]["verdict"].get("separated")
    if inv_ok:
        print("GUARD: separation SURVIVES normalization — size-invariant forms separate. "
              "Burial-driven, not size-driven.")
    elif abs_ok and not inv_ok:
        print("GUARD: STOP AND RETHINK — size-invariant forms FAIL while absolute passes. "
              "Original separation was size-driven. Production stays size-invariant anyway.")
    else:
        print("GUARD: no form cleanly separates — discriminator needs rework.")
    print(f"report   -> {os.path.relpath(args.report, os.getcwd())}")
    print(f"csv      -> {os.path.relpath(args.csv, os.getcwd())}")
    print(f"norm-csv -> {os.path.relpath(args.norm_csv, os.getcwd())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
