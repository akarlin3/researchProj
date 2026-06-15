"""Enlarge the random triad-bearing floor from n=28 toward n >= 200 triad-bearers.

WHY (revision Fix A.1): the central claim — "the PET-branch exposed-cleft rate is
indistinguishable from the random floor (Fisher p=0.30)" — rests on a random
triad-bearing baseline of only n=28 (12/28 = 42.86% above-line). That is a tiny
denominator: its Wilson 95% CI runs 26.5-60.9%, so a *failure to reject* there is
absence of evidence, not equivalence. To turn the null into a powered
equivalence/non-superiority statement we enlarge the baseline.

WHAT: draw ADDITIONAL uniform-random Atlas proteins from the SAME HQ-clust30
universe and pipe them through the *unchanged* screen pipeline
(proteus.floor.{sample_universe, fetch_many, screen_many}) at the SAME pinned
decision line (-1.1587) and the SAME percentile control anchor as the original
floor. The additional draws use a NEW explicit seed (ENLARGE_SEED, distinct from
the original floor seed 1729); block k draws under seed ENLARGE_SEED + k so both
sampling and screening can early-stop and remain bit-for-bit reproducible. New
accessions are de-duplicated against the original 1,500-draw floor and pooled with
it, so the result is a single larger uniform sample. We stop once the pooled
triad-bearer count reaches TARGET_TRIAD_POOLED.

This re-uses the production screen code verbatim (no pipeline change) — only the
sampling driver is new.

Reproduce, from the repo root:
    PYTHONPATH=src python analysis/enlarge_floor.py

Outputs (data/processed/):
    floor_enlarged_new.csv     per-protein rows for the NEW draws (flushed per block)
    floor_enlarged_pooled.csv  original 1,500 + new draws, combined per-protein
    floor_enlarged.json        summary: seeds, per-block funnels, pooled counts,
                               conditional above-line rate + Wilson 95% CI
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import sys

from proteus.atlas_screen import resolve_operating_point
from proteus.floor import (
    fetch_many,
    resolve_lookup,
    sample_universe,
    screen_many,
    wilson_ci,
)
from proteus.utils import REPO, load_config

# --------------------------------------------------------------------------- #
# Pinned constants (printed at run start per the reproducibility guardrail)
# --------------------------------------------------------------------------- #
LINE = -1.1587                 # pinned decision line, identical to the original floor
ORIG_SEED = 1729               # the original floor draw (do NOT reuse for new draws)
ENLARGE_SEED = 20260614        # NEW explicit seed for the additional draws (block k -> +k)
TARGET_TRIAD_POOLED = 200      # stop once pooled triad-bearers reach this
BLOCK_N = 2000                 # new uniform draws per block
MAX_BLOCKS = 10                # hard cap (<= 20,000 new draws) so the run is bounded
FETCH_WORKERS = 16
SCREEN_WORKERS = 8

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)     # worktree root
PROC = os.path.join(WT, "data", "processed")
HITS = os.path.join(WT, "structures", "floor_enlarged_hits")
# Original floor (1,500 draws, seed 1729) lives in the floor-measurement worktree.
ORIG_FLOOR_CSV = os.environ.get(
    "ORIG_FLOOR_CSV",
    os.path.join(os.path.dirname(WT), "floor-measurement",
                 "data", "processed", "floor.csv"),
)

NEW_CSV = os.path.join(PROC, "floor_enlarged_new.csv")
POOLED_CSV = os.path.join(PROC, "floor_enlarged_pooled.csv")
SUMMARY_JSON = os.path.join(PROC, "floor_enlarged.json")

COLS = ["accession", "mean_plddt", "n_res", "triad_found", "catalytic_ser",
        "his", "acid", "pocket_ok", "composite", "above_threshold",
        "petase_like_hit", "source_seed"]


def _tb(v) -> bool:
    return str(v).strip().lower() in ("true", "1")


def load_original(path: str):
    """Return (accessions set, list[row]) for the original 1,500-draw floor."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    accs = {r["accession"] for r in rows}
    return accs, rows


def append_rows(path: str, results: list[dict], seed: int, header: bool):
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        if header:
            w.writeheader()
        for r in results:
            row = {c: r.get(c) for c in COLS}
            row["source_seed"] = seed
            w.writerow(row)


def main() -> int:
    print(f"[enlarge] PINNED line = {LINE} | ORIG_SEED = {ORIG_SEED} | "
          f"ENLARGE_SEED = {ENLARGE_SEED} (block k -> seed {ENLARGE_SEED}+k)")
    print(f"[enlarge] target pooled triad-bearers >= {TARGET_TRIAD_POOLED}")

    cfg = load_config()
    struct_dir = os.path.join(WT, "structures")
    op = resolve_operating_point(cfg, struct_dir)
    anchor = op["anchor"]
    print(f"[enlarge] control anchor={op['positive_ids']} mode={op['mode']} "
          f"separated={op['separated']} margin={op['margin']}; SCREEN PINNED @ {LINE}")

    fc = cfg.get("floor", {})
    lookup_url = fc["lookup_url"]
    timeout = float(cfg.get("atlas", {}).get("request_timeout_s", 60))

    orig_accs, orig_rows = load_original(ORIG_FLOOR_CSV)
    orig_triad = [r for r in orig_rows if _tb(r["triad_found"])]
    orig_above = [r for r in orig_triad if _tb(r["petase_like_hit"])]
    print(f"[enlarge] original floor: {len(orig_rows)} screened, "
          f"{len(orig_triad)} triad+, {len(orig_above)} above-line "
          f"({len(orig_above)}/{len(orig_triad)} = "
          f"{len(orig_above)/len(orig_triad):.4f})")

    final_url, size = resolve_lookup(lookup_url, timeout)

    os.makedirs(PROC, exist_ok=True)
    os.makedirs(HITS, exist_ok=True)

    seen = set(orig_accs)
    new_results: list[dict] = []
    pooled_triad = len(orig_triad)
    block_funnels = []
    wrote_header = False

    for block in range(MAX_BLOCKS):
        if pooled_triad >= TARGET_TRIAD_POOLED:
            break
        seed_k = ENLARGE_SEED + block
        accs, _meta = sample_universe(final_url, size, BLOCK_N, seed_k,
                                      workers=FETCH_WORKERS, timeout=timeout)
        fresh = [a for a in accs if a not in seen]
        seen.update(fresh)
        bdir = os.path.join(HITS, f"block{block}")
        fetched = fetch_many(fresh, bdir, workers=FETCH_WORKERS, timeout=timeout,
                             label=f"blk{block}")
        n_ok = sum(1 for r in fetched.values() if r["pdb"])
        results = screen_many(fetched, cfg, anchor, LINE,
                              workers=SCREEN_WORKERS, label=f"blk{block}")
        b_triad = [r for r in results if r["triad_found"]]
        b_above = [r for r in b_triad if r.get("petase_like_hit")]
        new_results.extend(results)
        pooled_triad += len(b_triad)
        append_rows(NEW_CSV, results, seed_k, header=not wrote_header)
        wrote_header = True

        # keep triad+ PDBs (for active-site-local pLDDT in CP2), drop the rest
        keep = {r["accession"] for r in b_triad}
        for r in fetched.values():
            if r.get("pdb") and r["accession"] not in keep and os.path.exists(r["pdb"]):
                os.remove(r["pdb"])

        funnel = {
            "block": block, "seed": seed_k, "drawn": len(accs),
            "fresh": len(fresh), "fetched_ok": n_ok, "screened": len(results),
            "triad_positive_S4": len(b_triad), "above_line": len(b_above),
        }
        block_funnels.append(funnel)
        print(f"[enlarge] block {block} (seed {seed_k}): fresh={len(fresh)} "
              f"ok={n_ok} triad+={len(b_triad)} above={len(b_above)} | "
              f"POOLED triad+={pooled_triad}")
        # flush summary each block so partial progress is inspectable
        write_summary(orig_rows, orig_triad, orig_above, new_results, block_funnels,
                      final_url, size, done=pooled_triad >= TARGET_TRIAD_POOLED)

    write_pooled(orig_rows, new_results)
    write_summary(orig_rows, orig_triad, orig_above, new_results, block_funnels,
                  final_url, size, done=True)
    print(f"[enlarge] DONE. pooled triad-bearers = {pooled_triad}. "
          f"summary -> {SUMMARY_JSON}")
    return 0


def write_pooled(orig_rows, new_results):
    with open(POOLED_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in orig_rows:
            row = {c: r.get(c) for c in COLS}
            row["source_seed"] = ORIG_SEED
            w.writerow(row)
        for r in new_results:
            row = {c: r.get(c) for c in COLS}
            row["source_seed"] = "new"
            w.writerow(row)


def write_summary(orig_rows, orig_triad, orig_above, new_results, block_funnels,
                  final_url, size, done):
    new_triad = [r for r in new_results if r["triad_found"]]
    new_above = [r for r in new_triad if r.get("petase_like_hit")]
    pooled_triad_n = len(orig_triad) + len(new_triad)
    pooled_above_n = len(orig_above) + len(new_above)
    pooled_screened = len(orig_rows) + len(new_results)
    rate = pooled_above_n / pooled_triad_n if pooled_triad_n else None
    summary = {
        "pinned_line": LINE,
        "orig_seed": ORIG_SEED,
        "enlarge_seed_base": ENLARGE_SEED,
        "seed_scheme": "block k drawn under seed ENLARGE_SEED + k",
        "target_triad_pooled": TARGET_TRIAD_POOLED,
        "complete": done,
        "lookup_final_url": final_url,
        "lookup_bytes": size,
        "blocks": block_funnels,
        "original_floor": {
            "screened": len(orig_rows),
            "triad_positive_S4": len(orig_triad),
            "above_line": len(orig_above),
        },
        "new_draws": {
            "screened": len(new_results),
            "triad_positive_S4": len(new_triad),
            "above_line": len(new_above),
        },
        "pooled_floor": {
            "screened": pooled_screened,
            "triad_positive_S4": pooled_triad_n,
            "above_line": pooled_above_n,
            "conditional_above_given_triad": rate,
            "wilson95": list(wilson_ci(pooled_above_n, pooled_triad_n))
            if pooled_triad_n else None,
        },
    }
    with open(SUMMARY_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    raise SystemExit(main())
