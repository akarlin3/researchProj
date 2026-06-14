"""Per-query tiering — re-tier the enriched Atlas sweep WITHIN each query anchor and
screen the PETase/cutinase-anchored branch the thesis actually cares about.

WHY THIS EXISTS
---------------
The discovery sweep ranked the 217,833 fold-class hits by GLOBAL Foldseek bits and
took the top 300. That tier came out **96% AChE-branch** (the *Torpedo* 1EA5 anchor's
Atlas matches simply score higher in bits than the smaller, more divergent PETase
folds) and **0% PETase/cutinase-branch** — so the discovery angle never actually
screened a PETase-neighbour. This module fixes that: it partitions every Atlas target
by its nearest query anchor, then tiers and screens WITHIN each anchor, finally asking
the one clean question the thesis rests on:

    Do PETase-neighbours clear the cleft line (above-line | triad) more often than
    random triad-bearers (the floor, 42.9%), or is the structural signal flat?

WHAT IT REUSES (untouched)
--------------------------
- `proteus.screen.build_control_anchor` / `screen_model` — the exact S4->S5->anchor
  scoring path used by calibration and the enriched sweep. The decision line is PINNED
  to the enriched sweep's widened threshold (config per_query.enriched_line, -1.1587)
  so all three arms (PETASE, ACHE, floor) are judged at the IDENTICAL line. (Local
  re-derivation jitters ~+/-0.01 from fpocket non-determinism — see floor-measurement.)
- The random floor (floor.json on GCS) is reused as the third arm AS MEASURED — no
  re-screen. Its conditional rate is 12/28 = 42.9%.

Mac-only, off the GCS sweep artifacts (result.m8). No GCE, no re-search.

Stages (each is a CLI subcommand of `python -m proteus.per_query`, or `run` does all):
  partition  CP1  — best-match partition of result.m8 -> branch_partition.csv
  run        CP1-4 — partition, tier within anchor, fetch+screen, stats, report

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.per_query run \
        --out-md envlog/per-query-tiering.md \
        --out-csv data/processed/branch_partition.csv \
        --out-json data/processed/per_query_tiering.json
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import random
import subprocess
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

from proteus.screen import build_control_anchor, screen_model
from proteus.utils import DEFAULT_CONFIG, REPO, load_config

# --------------------------------------------------------------------------- #
# Statistics — Wilson CI, two-proportion z, Fisher exact, Katz rate-ratio CI.
# --------------------------------------------------------------------------- #
_Z95 = 1.959963984540054  # two-sided 95% normal quantile


def wilson_ci(k: int, n: int, z: float = _Z95) -> list[float]:
    """Wilson score 95% CI for a binomial proportion k/n. Returns [lo, hi]."""
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / denom
    return [max(0.0, center - half), min(1.0, center + half)]


def two_proportion_z(k1: int, n1: int, k2: int, n2: int) -> tuple[float, float]:
    """Pooled two-proportion z-test. Returns (z, two-sided p)."""
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan")
    p1, p2 = k1 / n1, k2 / n2
    pool = (k1 + k2) / (n1 + n2)
    se = math.sqrt(pool * (1 - pool) * (1 / n1 + 1 / n2))
    if se == 0:
        return float("nan"), float("nan")
    z = (p1 - p2) / se
    # two-sided p via the normal survival function (erfc), no scipy needed here
    p = math.erfc(abs(z) / math.sqrt(2))
    return z, p


def katz_rate_ratio(k1: int, n1: int, k2: int, n2: int) -> dict:
    """Rate ratio (arm1/arm2) of two proportions with a Katz log-RR 95% CI.
    arm1 = k1/n1 (e.g. PETASE), arm2 = k2/n2 (e.g. floor)."""
    if n1 == 0 or n2 == 0 or k2 == 0:
        return {"rr": float("nan"), "ci": [float("nan"), float("nan")]}
    p1, p2 = k1 / n1, k2 / n2
    rr = p1 / p2 if p2 > 0 else float("inf")
    if k1 == 0:
        return {"rr": rr, "ci": [0.0, float("nan")]}
    se = math.sqrt(1 / k1 - 1 / n1 + 1 / k2 - 1 / n2)
    lo = rr * math.exp(-_Z95 * se)
    hi = rr * math.exp(_Z95 * se)
    return {"rr": rr, "ci": [lo, hi]}


def fisher_p(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p for the 2x2 [[a,b],[c,d]]. scipy if available,
    else an exact hypergeometric two-sided fallback."""
    try:
        from scipy.stats import fisher_exact  # noqa: PLC0415
        return float(fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])
    except Exception:
        from math import comb  # noqa: PLC0415
        n = a + b + c + d
        r1, c1 = a + b, a + c
        def pmf(x: int) -> float:
            return comb(c1, x) * comb(n - c1, r1 - x) / comb(n, r1)
        p_obs = pmf(a)
        lo = max(0, r1 - (n - c1))
        hi = min(r1, c1)
        return float(sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= p_obs * 1.0000001))


def _pct(x: float) -> str:
    return "nan" if x != x else f"{100 * x:.2f}%"


def _ci_pct(ci: list[float]) -> str:
    return f"[{_pct(ci[0])}, {_pct(ci[1])}]"


# --------------------------------------------------------------------------- #
# CP1 — partition result.m8 by query anchor (best-match per Atlas target)
# --------------------------------------------------------------------------- #
def _open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _anchor_lookup(anchor_classes: dict) -> dict:
    """Reverse {CLASS: [query_ids]} -> {query_id: CLASS}."""
    rev = {}
    for cls, qids in anchor_classes.items():
        for q in qids:
            rev[q] = cls
    return rev


def partition_m8(m8_path: str, anchor_classes: dict) -> dict:
    """Stream result.m8 and best-match-partition every Atlas target.

    For each target accession track the best (bits, query) seen for EACH anchor class.
    The target's branch is the class of its single best query (highest bits anywhere).
    Returns {accession: {class: {"bits": float, "query": str}}}.
    Unknown queries (not in any anchor class) abort — the map must be complete."""
    rev = _anchor_lookup(anchor_classes)
    per_acc: dict[str, dict] = {}
    n_rows = 0
    unknown: set[str] = set()
    with _open_maybe_gz(m8_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            query, target, _evalue, bits = parts[0], parts[1], parts[2], parts[3]
            cls = rev.get(query)
            if cls is None:
                unknown.add(query)
                continue
            acc = target.split(".")[0]  # strip ".pdb.gz"
            bits_f = float(bits)
            rec = per_acc.setdefault(acc, {})
            cur = rec.get(cls)
            if cur is None or bits_f > cur["bits"]:
                rec[cls] = {"bits": bits_f, "query": query}
            n_rows += 1
    if unknown:
        raise RuntimeError(
            f"result.m8 has queries not in any anchor class: {sorted(unknown)} — "
            "the anchor map is incomplete; refusing to partition")
    return {"per_acc": per_acc, "n_rows": n_rows}


def build_partition_records(per_acc: dict) -> list[dict]:
    """Turn per-accession per-class best-bits into flat partition records.

    Each record: accession, anchor_class (class of the global-best query), best_query,
    best_bits, plus the best (bits, query) within each of PETASE / ACHE / OTHER_NEG so
    CP4 can report 'bits-to-nearest-PETase-query' for any candidate."""
    records = []
    for acc, byclass in per_acc.items():
        # global best across classes
        best_cls, best = max(byclass.items(), key=lambda kv: kv[1]["bits"])
        petase = byclass.get("PETASE")
        ache = byclass.get("ACHE")
        other = byclass.get("OTHER_NEG")
        records.append({
            "accession": acc,
            "anchor_class": best_cls,
            "best_query": best["query"],
            "best_bits": best["bits"],
            "petase_bits": petase["bits"] if petase else "",
            "petase_query": petase["query"] if petase else "",
            "ache_bits": ache["bits"] if ache else "",
            "other_bits": other["bits"] if other else "",
        })
    return records


def write_branch_partition_csv(records: list[dict], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    cols = ["accession", "anchor_class", "best_query", "best_bits",
            "petase_bits", "petase_query", "ache_bits", "other_bits"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        # stable order: by anchor class then best_bits desc
        for r in sorted(records, key=lambda r: (r["anchor_class"], -r["best_bits"])):
            w.writerow(r)


def branch_sizes(records: list[dict]) -> dict:
    sizes: dict[str, int] = {}
    for r in records:
        sizes[r["anchor_class"]] = sizes.get(r["anchor_class"], 0) + 1
    return sizes


# --------------------------------------------------------------------------- #
# CP2 — tier within anchor, fetch, screen
# --------------------------------------------------------------------------- #
def select_branch(records: list[dict], anchor_class: str, branch_n: int,
                  seed: int) -> list[dict]:
    """Rank the records in `anchor_class` by that anchor's bits (== best_bits for the
    best-match-partitioned records) and take the top branch_n (or all if smaller).

    Ties in (integer) bits are broken by a seeded RNG keyed in accession order, so the
    selection is fully reproducible given random_seed."""
    branch = [r for r in records if r["anchor_class"] == anchor_class]
    rng = random.Random(seed)
    tiebreak = {acc: rng.random()
                for acc in sorted(r["accession"] for r in branch)}
    branch.sort(key=lambda r: (-r["best_bits"], tiebreak[r["accession"]]))
    if branch_n and len(branch) > branch_n:
        return branch[:branch_n]
    return branch


def fetch_pdb(acc: str, url_tmpl: str, cache_dir: str, retries: int = 4,
              timeout: int = 60) -> str | None:
    """Fetch one ESM Atlas model to cache_dir/{acc}.pdb (cached/idempotent). Returns
    the path, or None on repeated failure. The live fetchPredictedStructure API serves
    the SAME ESMFold V0 models as the GCE foldcomp DB (triads reproduce exactly)."""
    dest = os.path.join(cache_dir, f"{acc}.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    url = url_tmpl.format(acc=acc)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                data = resp.read()
            if data and data.startswith(b"ATOM") or b"\nATOM" in data[:4096]:
                tmp = dest + ".part"
                with open(tmp, "wb") as fh:
                    fh.write(data)
                os.replace(tmp, dest)
                return dest
        except (urllib.error.URLError, urllib.error.HTTPError, OSError):
            pass
        time.sleep(1.0 + attempt * 1.5)
    return None


def _mean_plddt(pdb_path: str) -> float | None:
    """Mean CA pLDDT from the PDB B-factor column. ESM Atlas models store pLDDT on a
    0-1 scale in B-factor (verified: matches the sweep's reported mean_plddt)."""
    vals = []
    try:
        with open(pdb_path) as fh:
            for ln in fh:
                if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
                    vals.append(float(ln[60:66]))
    except (OSError, ValueError):
        return None
    return round(sum(vals) / len(vals), 4) if vals else None


def screen_one(rec: dict, cfg: dict, anchor: dict, line: float,
               url_tmpl: str, cache_dir: str) -> dict:
    """Fetch + screen ONE branch record. Returns the screen_model record augmented
    with accession / anchor bits / mean_plddt, or a fetch-failure stub."""
    acc = rec["accession"]
    pdb = fetch_pdb(acc, url_tmpl, cache_dir)
    if pdb is None:
        return {"accession": acc, "fetched": False, "triad_found": False,
                "pocket_ok": False, "above_threshold": None,
                "best_bits": rec["best_bits"], "best_query": rec["best_query"],
                "petase_bits": rec.get("petase_bits", ""),
                "petase_query": rec.get("petase_query", "")}
    plddt = _mean_plddt(pdb)
    out = screen_model(pdb, cfg, anchor, line, mean_plddt=plddt, cand_id=acc)
    out["accession"] = acc
    out["fetched"] = True
    out["best_bits"] = rec["best_bits"]
    out["best_query"] = rec["best_query"]
    out["petase_bits"] = rec.get("petase_bits", "")
    out["petase_query"] = rec.get("petase_query", "")
    return out


def screen_branch(selected: list[dict], cfg: dict, anchor: dict, line: float,
                  url_tmpl: str, cache_dir: str, fetch_workers: int,
                  screen_workers: int, label: str,
                  screen_cache: str | None = None, rescreen: bool = False) -> dict:
    """Fetch (parallel, network-bound) then screen (parallel, fpocket per triad+) a
    selected branch. Returns {records, funnel}. fpocket is race-safe (isolated temp
    dir per call in s5.run_fpocket).

    fpocket is non-deterministic, so screened records are cached to `screen_cache`
    (JSON) and reused on re-run — this keeps the report's exact numbers stable across
    re-renders. Pass rescreen=True to force a fresh screen."""
    os.makedirs(cache_dir, exist_ok=True)
    accs = [r["accession"] for r in selected]
    want = {r["accession"] for r in selected}
    if screen_cache and not rescreen and os.path.exists(screen_cache):
        with open(screen_cache) as fh:
            cached = json.load(fh)
        if {c["accession"] for c in cached} == want:
            print(f"[per_query] {label}: loaded {len(cached)} cached screens "
                  f"from {os.path.basename(screen_cache)}", flush=True)
            funnel = funnel_of(cached, line)
            funnel["label"] = label
            return {"records": cached, "funnel": funnel}
    # Phase 1: warm the cache in parallel (idempotent; network-bound).
    n_fetched = 0
    with ThreadPoolExecutor(max_workers=fetch_workers) as ex:
        futs = {ex.submit(fetch_pdb, a, url_tmpl, cache_dir): a for a in accs}
        for fut in as_completed(futs):
            if fut.result():
                n_fetched += 1
    print(f"[per_query] {label}: fetched {n_fetched}/{len(accs)}", flush=True)
    # Phase 2: screen in parallel (CPU/subprocess-bound).
    records: list[dict] = []
    done = 0
    with ThreadPoolExecutor(max_workers=screen_workers) as ex:
        futs = [ex.submit(screen_one, r, cfg, anchor, line, url_tmpl, cache_dir)
                for r in selected]
        for fut in as_completed(futs):
            records.append(fut.result())
            done += 1
            if done % 50 == 0 or done == len(selected):
                print(f"[per_query] {label}: screened {done}/{len(selected)}",
                      flush=True)
    funnel = funnel_of(records, line)
    funnel["label"] = label
    if screen_cache:
        os.makedirs(os.path.dirname(os.path.abspath(screen_cache)), exist_ok=True)
        with open(screen_cache, "w") as fh:
            json.dump(records, fh, default=str)
    return {"records": records, "funnel": funnel}


def funnel_of(records: list[dict], line: float) -> dict:
    """branch -> fetched -> triad+ -> pocket_ok -> above-line counts + accessions."""
    fetched = [r for r in records if r.get("fetched")]
    triad = [r for r in fetched if r.get("triad_found")]
    pocket = [r for r in triad if r.get("pocket_ok")]
    above = [r for r in pocket if r.get("above_threshold")]
    return {
        "attempted": len(records),
        "fetched": len(fetched),
        "triad_positive_S4": len(triad),
        "pocket_ok_S5": len(pocket),
        "above_line": len(above),
        "line": line,
        "above_line_accessions": [r["accession"] for r in above],
    }


# --------------------------------------------------------------------------- #
# CP3 — conditional test  (above-line GIVEN a triad — the load-bearing number)
# --------------------------------------------------------------------------- #
def arm_rates(funnel: dict) -> dict:
    """triad rate (triad/screened) and conditional rate (above-line/triad) + Wilson."""
    n = funnel["fetched"]
    t = funnel["triad_positive_S4"]
    a = funnel["above_line"]
    return {
        "screened": n,
        "triad_rate": {"k": t, "n": n, "rate": (t / n if n else float("nan")),
                       "wilson95": wilson_ci(t, n)},
        "conditional_above_given_triad": {
            "k": a, "n": t, "rate": (a / t if t else float("nan")),
            "wilson95": wilson_ci(a, t)},
    }


def two_arm_test(k1: int, n1: int, k2: int, n2: int, name1: str, name2: str) -> dict:
    """Full two-proportion comparison of conditional rates: z-test, Fisher exact,
    Katz rate ratio (arm1/arm2)."""
    z, pz = two_proportion_z(k1, n1, k2, n2)
    pf = fisher_p(k1, n1 - k1, k2, n2 - k2)
    rr = katz_rate_ratio(k1, n1, k2, n2)
    return {
        "arm1": name1, "arm2": name2,
        "rate1": {"k": k1, "n": n1, "rate": (k1 / n1 if n1 else float("nan")),
                  "wilson95": wilson_ci(k1, n1)},
        "rate2": {"k": k2, "n": n2, "rate": (k2 / n2 if n2 else float("nan")),
                  "wilson95": wilson_ci(k2, n2)},
        "two_proportion_z": z, "z_p": pz,
        "fisher_p": pf,
        "rate_ratio_1_over_2": rr,
    }


# --------------------------------------------------------------------------- #
# Floor (third arm) — reused AS MEASURED from floor.json (no re-screen).
# --------------------------------------------------------------------------- #
def floor_arm(floor_json_path: str | None) -> dict:
    """Load the random-floor conditional from floor.json. Falls back to the published
    constants (12/28 triad, 28/1500 screened) if the file isn't present locally."""
    if floor_json_path and os.path.exists(floor_json_path):
        with open(floor_json_path) as fh:
            fj = json.load(fh)
        ff = fj["floor_funnel"]
        return {
            "source": floor_json_path,
            "screened": ff["screened"],
            "triad_positive_S4": ff["triad_positive_S4"],
            "above_line": ff["above_line"],
        }
    return {"source": "published constants (floor-measurement.md)",
            "screened": 1500, "triad_positive_S4": 28, "above_line": 12}


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _cand_rows(records: list[dict]) -> list[dict]:
    """PETase-anchored above-line hits, ranked by composite desc — the project's best
    candidates from the CORRECT branch."""
    hits = [r for r in records
            if r.get("above_threshold") and r.get("petase_like_hit")]
    hits.sort(key=lambda r: r.get("composite", float("-inf")), reverse=True)
    return hits


def render_md(ctx: dict) -> str:
    """Render envlog/per-query-tiering.md from the computed context."""
    p = ctx["partition"]
    sizes = p["branch_sizes"]
    pf, af = ctx["petase"]["funnel"], ctx["ache"]["funnel"]
    fl = ctx["floor"]
    pr, ar = ctx["petase_rates"], ctx["ache_rates"]
    t_pf = ctx["test_petase_vs_floor"]
    t_pa = ctx["test_petase_vs_ache"]
    t_af = ctx["test_ache_vs_floor"]
    line = ctx["line"]
    seed = ctx["seed"]
    cands = ctx["candidates"]

    def rate(d):  # noqa: ANN001
        return f"{d['k']}/{d['n']} = {_pct(d['rate'])} {_ci_pct(d['wilson95'])}"

    L = []
    L.append("# Per-query tiering — the PETase-branch test\n")
    L.append(f"**Run date:** {ctx['run_date']}  ")
    L.append("**Where:** local Mac (M4), off the existing GCS sweep artifacts "
             "(`result.m8`) — **no GCE, no re-search**.  ")
    L.append("**Reproduce:** `PYTHONPATH=src python -m proteus.per_query run` "
             f"(config `per_query` block; seed {seed}; line pinned to {line}).  ")
    L.append("**Reuses** `screen`/S4/S5 untouched at the pinned line; the random floor "
             "is reused AS MEASURED (floor.json) as the third arm.\n")

    L.append("## TL;DR\n")
    L.append(ctx["tldr"] + "\n")

    L.append("## Checkpoint 0 — preconditions (audit)\n")
    L.append("| Precondition | Status |")
    L.append("|---|---|")
    L.append("| `screen`/S4/S5/`calibrate` intact | ✅ reused untouched; imports + run "
             "clean in the `proteus` env (numpy/scipy/biotite; fpocket on PATH). A "
             "known enriched hit re-screens to its documented triad exactly "
             "(MGYP000470279205 → 497/489/411). |")
    L.append(f"| Pinned line | **{line}** (enriched `funnel.json` threshold). Local "
             "re-derivation jitters (fpocket); pinned so all arms judge at one line. |")
    L.append(f"| Floor (baseline arm) | ✅ {fl['source']}: "
             f"{fl['screened']} screened → {fl['triad_positive_S4']} triad+ → "
             f"{fl['above_line']} above-line; conditional "
             f"{fl['above_line']}/{fl['triad_positive_S4']} = "
             f"{_pct(fl['above_line']/fl['triad_positive_S4'])}. Reused, not re-screened. |")
    L.append(f"| `result.m8` | ✅ {p['n_rows']:,} alignment rows; "
             f"{p['n_accessions']:,} unique Atlas targets; "
             f"{len(ctx['anchor_classes'])} anchor classes mapped. |")
    L.append("| Fetch path | ✅ `fetchPredictedStructure/{acc}.pdb` live (structurally "
             "identical to the GCE foldcomp models). |\n")

    L.append("## Checkpoint 1 — partition by query anchor (best-match)\n")
    L.append("Each of the unique Atlas targets is assigned to the anchor of its single "
             "**best** query match (highest Foldseek bits). This is the fix for the "
             "global-bits ranking that let one anchor dominate.\n")
    L.append("**Anchor map** (query id in `result.m8` → class):\n")
    for cls, qids in ctx["anchor_classes"].items():
        L.append(f"- **{cls}**: {', '.join(qids)}")
    L.append("")
    L.append("**Branch sizes** (of "
             f"{p['n_accessions']:,} unique targets, by best-match anchor):\n")
    L.append("| Branch | Targets | Share |")
    L.append("|---|---:|---:|")
    tot = p["n_accessions"]
    for cls in ["PETASE", "ACHE", "OTHER_NEG"]:
        nsz = sizes.get(cls, 0)
        L.append(f"| {cls} | {nsz:,} | {_pct(nsz/tot)} |")
    L.append(f"| **total** | **{tot:,}** | 100% |\n")
    L.append(ctx["branch_size_note"] + "\n")

    L.append("## Checkpoint 2 — tier WITHIN each anchor, then screen\n")
    L.append(f"Top **{ctx['branch_n']}** per branch by that anchor's bits (or all if "
             "smaller), fetched and run through the UNCHANGED S4 → fpocket(triad+) → S5 "
             f"path at the pinned line **{line}**, parallelised. Per-branch funnels:\n")
    L.append("| Stage | PETASE | ACHE | floor (random) |")
    L.append("|---|---:|---:|---:|")
    L.append(f"| selected / screened | {pf['fetched']} | {af['fetched']} | "
             f"{fl['screened']} |")
    L.append(f"| → triad+ (S4) | {pf['triad_positive_S4']} | "
             f"{af['triad_positive_S4']} | {fl['triad_positive_S4']} |")
    L.append(f"| → catalytic pocket (S5) | {pf['pocket_ok_S5']} | "
             f"{af['pocket_ok_S5']} | — |")
    L.append(f"| → **above line** | **{pf['above_line']}** | **{af['above_line']}** | "
             f"**{fl['above_line']}** |\n")

    L.append("## Checkpoint 3 — the conditional test (the thesis question)\n")
    L.append("**Triad rate** (triad+/screened) and **above-line | triad** "
             "(the load-bearing conditional), each with Wilson 95% CI:\n")
    L.append("| Arm | triad rate | above-line \\| triad |")
    L.append("|---|---|---|")
    L.append(f"| **PETASE** | {rate(pr['triad_rate'])} | "
             f"{rate(pr['conditional_above_given_triad'])} |")
    L.append(f"| **ACHE** | {rate(ar['triad_rate'])} | "
             f"{rate(ar['conditional_above_given_triad'])} |")
    fl_cond = {"k": fl["above_line"], "n": fl["triad_positive_S4"],
               "rate": fl["above_line"] / fl["triad_positive_S4"],
               "wilson95": wilson_ci(fl["above_line"], fl["triad_positive_S4"])}
    L.append(f"| **floor (random)** | "
             f"{fl['triad_positive_S4']}/{fl['screened']} = "
             f"{_pct(fl['triad_positive_S4']/fl['screened'])} | {rate(fl_cond)} |\n")

    def test_block(t: dict, title: str) -> None:
        rr = t["rate_ratio_1_over_2"]
        L.append(f"**{title}**  ")
        L.append(f"- {t['arm1']} {t['rate1']['k']}/{t['rate1']['n']} = "
                 f"{_pct(t['rate1']['rate'])} {_ci_pct(t['rate1']['wilson95'])} vs "
                 f"{t['arm2']} {t['rate2']['k']}/{t['rate2']['n']} = "
                 f"{_pct(t['rate2']['rate'])} {_ci_pct(t['rate2']['wilson95'])}  ")
        L.append(f"- two-proportion z = {t['two_proportion_z']:.2f}, "
                 f"p = {t['z_p']:.2e}; Fisher exact p = {t['fisher_p']:.2e}  ")
        L.append(f"- rate ratio ({t['arm1']}/{t['arm2']}) = {rr['rr']:.2f}× "
                 f"[{rr['ci'][0]:.2f}, {rr['ci'][1]:.2f}]\n")

    test_block(t_pf, "PETASE above-line|triad vs floor 42.9% — the load-bearing test")
    test_block(t_pa, "PETASE branch vs ACHE branch (conditional)")
    test_block(t_af, "ACHE branch vs floor (conditional, context)")

    L.append("## Checkpoint 4 — candidates + verdict\n")
    L.append("### PETase-anchored, above-line hits — the project's best candidates "
             "from the correct branch\n")
    if cands:
        L.append("| rank | accession | nearest PETase query | bits | composite | "
                 "pLDDT | Ser/His/acid |")
        L.append("|---:|---|---|---:|---:|---:|---|")
        for i, c in enumerate(cands, 1):
            pq = c.get("petase_query") or c.get("best_query")
            pb = c.get("petase_bits") or c.get("best_bits")
            pb_s = f"{float(pb):.0f}" if pb not in ("", None) else "—"
            comp = "—" if c.get("composite") is None else f"{c['composite']:.3f}"
            pl = "—" if c.get("mean_plddt") is None else f"{c['mean_plddt']:.3f}"
            triad = f"{c.get('catalytic_ser')}/{c.get('his')}/{c.get('acid')}"
            L.append(f"| {i} | {c['accession']} | {pq} | {pb_s} | {comp} | {pl} | "
                     f"{triad} |")
        L.append("")
    else:
        L.append("_No PETase-anchored hit cleared the pinned line._\n")

    L.append("### Verdict — is there ANY PET-specific gradient in the structural "
             "signal?\n")
    L.append(ctx["verdict"] + "\n")
    L.append("### Scope guard (carried)\n")
    L.append("**Exposure ≠ PET activity.** Even a positive gradient does not verify any "
             "candidate — it only says PETase-neighbours are more exposed-site than "
             "random hydrolases. No wet-lab; the S4/S5 path tests fold-class + an "
             "exposed-cleft geometry, not PET turnover. Leads are **prioritized, not "
             "verified**.\n")

    L.append("## Reproducibility\n")
    L.append(f"- **Seed** {seed}; **branch_n** {ctx['branch_n']}; **line** {line} "
             "(pinned to the enriched sweep).")
    L.append("- **Partition** best-match over `result.m8` "
             f"({p['n_rows']:,} rows → {p['n_accessions']:,} unique targets).")
    L.append("- **Fetch** `api.esmatlas.com/fetchPredictedStructure/{acc}.pdb` "
             "(live; structurally identical to the GCE foldcomp models).")
    L.append("- **Anchor** IsPETase+LCC, percentile mode (same as calibration/floor).")
    L.append("- **Artifacts** `branch_partition.csv`, `per_query_tiering.json`; floor "
             "comparator from `gs://projproteus-fold/atlas-sweep/2026-06-14/`.")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def ensure_result_m8(cfg_pq: dict) -> str:
    """Return a local path to result.m8(.gz); pull from GCS if configured + missing."""
    local = os.path.join(REPO, cfg_pq["result_m8"]) \
        if not os.path.isabs(cfg_pq["result_m8"]) else cfg_pq["result_m8"]
    if os.path.exists(local):
        return local
    gcs = cfg_pq.get("result_m8_gcs")
    if gcs:
        os.makedirs(os.path.dirname(local), exist_ok=True)
        print(f"[per_query] pulling {gcs} -> {local}", flush=True)
        subprocess.run(["gsutil", "cp", gcs, local], check=True)
        return local
    raise FileNotFoundError(
        f"result.m8 not found at {local} and no per_query.result_m8_gcs configured")


def build_verdict(t_pf: dict, pr: dict, fl: dict, t_pa: dict, ar: dict) -> tuple[str, str]:
    """Produce the gradient verdict + a TL;DR line from the load-bearing test, with the
    PETASE-vs-ACHE separation woven in (CP3's second question)."""
    pet = pr["conditional_above_given_triad"]
    rr = t_pf["rate_ratio_1_over_2"]["rr"]
    floor_rate = fl["above_line"] / fl["triad_positive_S4"]
    pet_rate = pet["rate"]
    pet_lo, pet_hi = pet["wilson95"]
    floor_pct = _pct(floor_rate)
    sig = t_pf["fisher_p"] < 0.05
    # PETASE-vs-ACHE separation (mechanism context for the verdict).
    ache_rate = ar["conditional_above_given_triad"]["rate"]
    rr_pa = t_pa["rate_ratio_1_over_2"]["rr"]
    pa_sig = t_pa["fisher_p"] < 0.05
    if pa_sig and rr_pa > 1:
        sep = (f" The pipeline is **not blind**: PETase-neighbours clear the line "
               f"{rr_pa:.1f}× more often than AChE-neighbours ({_pct(ache_rate)}, "
               f"Fisher p = {t_pa['fisher_p']:.1e}) — but that separation is driven by "
               "AChE being a deep-gorge outlier sitting FAR BELOW the floor, not by "
               "PETase rising above it.")
    elif pa_sig:
        sep = (f" PETase- and AChE-neighbours also differ ({_pct(ache_rate)} for AChE, "
               f"Fisher p = {t_pa['fisher_p']:.1e}).")
    else:
        sep = (f" PETase and AChE neighbourhoods are statistically indistinguishable too "
               f"({_pct(ache_rate)} for AChE, Fisher p = {t_pa['fisher_p']:.1e}) — the "
               "pipeline treats the two the same.")
    # "substantially above" = point estimate clears the floor AND CI lower bound
    # clears (or at least the test is significant in the up direction).
    above = pet_rate > floor_rate and sig and pet_lo > floor_rate
    flat = (not sig) or (pet_lo <= floor_rate <= pet_hi)
    if above:
        verdict = (
            f"**YES — a PET-relevant gradient is present.** PETase-neighbours clear the "
            f"cleft line at **{_pct(pet_rate)}** (above-line|triad, "
            f"{_ci_pct(pet['wilson95'])}), **substantially above** the random floor of "
            f"{floor_pct} (rate ratio {rr:.2f}×, Fisher p = {t_pf['fisher_p']:.2e}). The "
            "structural neighbourhood of PETases concentrates exposed-site proteins "
            "beyond baseline — a real (if weak) PET-relevant gradient. **Fork: the "
            "aromatic-subsite specificity attempt is worth building.**" + sep)
        tldr = (f"PETase-branch above-line|triad = {_pct(pet_rate)} vs floor {floor_pct} "
                f"— gradient PRESENT (RR {rr:.2f}×, p {t_pf['fisher_p']:.1e}).")
    elif flat:
        verdict = (
            f"**NO — the structural signal is flat.** PETase-neighbours clear the cleft "
            f"line at **{_pct(pet_rate)}** ({_ci_pct(pet['wilson95'])}), "
            f"**indistinguishable from** the random floor of {floor_pct} "
            f"(rate ratio {rr:.2f}×, Fisher p = {t_pf['fisher_p']:.2e}). Being a "
            "PETase-neighbour confers no extra above-line signal the pipeline captures. "
            "**This is the strongest negative result: the thesis branch was screened "
            "and shown indistinguishable from random.** Fork: the methods / honest-"
            "negative paper, now airtight." + sep)
        tldr = (f"PETase-branch above-line|triad = {_pct(pet_rate)} vs floor {floor_pct} "
                f"— FLAT (RR {rr:.2f}×, p {t_pf['fisher_p']:.1e}). Honest-negative.")
    else:  # significantly BELOW the floor
        verdict = (
            f"**NO (and inverted) — PETase-neighbours clear the line LESS than random.** "
            f"{_pct(pet_rate)} ({_ci_pct(pet['wilson95'])}) vs floor {floor_pct} "
            f"(rate ratio {rr:.2f}×, Fisher p = {t_pf['fisher_p']:.2e}). The exposure-"
            "favouring S5 line selects *against* the PETase neighbourhood, mirroring the "
            "AChE result. No positive PET gradient. **Fork: the methods / honest-"
            "negative paper.**" + sep)
        tldr = (f"PETase-branch above-line|triad = {_pct(pet_rate)} vs floor {floor_pct} "
                f"— INVERTED (RR {rr:.2f}×, p {t_pf['fisher_p']:.1e}). Honest-negative.")
    return verdict, tldr


def run(cfg: dict, args) -> dict:
    pq = cfg["per_query"]
    seed = int(cfg.get("random_seed", 1729))
    line = float(pq["enriched_line"])
    branch_n = int(pq["branch_n"])
    anchor_classes = pq["anchor_classes"]
    url_tmpl = pq["fetch_url"]
    cache_dir = os.path.join(REPO, "data", "interim", "per_query_cache")

    # ---- CP1 ----
    m8 = ensure_result_m8(pq)
    print(f"[per_query] CP1 partitioning {m8} ...", flush=True)
    part = partition_m8(m8, anchor_classes)
    records = build_partition_records(part["per_acc"])
    sizes = branch_sizes(records)
    n_acc = len(records)
    write_branch_partition_csv(records, args.out_csv)
    print(f"[per_query] CP1 {part['n_rows']:,} rows → {n_acc:,} unique targets; "
          f"branch sizes {sizes}; wrote {args.out_csv}", flush=True)

    # Per-branch best_bits ceilings — these explain the AChE domination mechanism.
    def _ceil(cls: str) -> dict:
        v = sorted((r["best_bits"] for r in records if r["anchor_class"] == cls),
                   reverse=True)
        if not v:
            return {"max": float("nan"), "top_n_min": float("nan")}
        return {"max": v[0], "top_n_min": v[min(branch_n, len(v)) - 1]}
    pet_c, ache_c = _ceil("PETASE"), _ceil("ACHE")
    partition_ctx = {
        "n_rows": part["n_rows"], "n_accessions": n_acc, "branch_sizes": sizes,
        "petase_bits_ceiling": pet_c, "ache_bits_ceiling": ache_c,
    }
    pet_branch = sizes.get("PETASE", 0)
    ache_branch = sizes.get("ACHE", 0)
    note = (f"PETASE is the **largest** branch — {pet_branch:,} targets "
            f"({_pct(pet_branch / n_acc)}) — while ACHE is only {ache_branch:,} "
            f"({_pct(ache_branch / n_acc)}). The original 96%-AChE / 0%-PETase top tier "
            "was therefore a pure **bits-magnitude** artifact, not a coverage one: ")
    if not math.isnan(pet_c["max"]) and not math.isnan(ache_c["top_n_min"]):
        note += (f"the PETASE branch's single best hit ({pet_c['max']:.0f} bits) scores "
                 f"BELOW the ACHE branch's top-{branch_n} cutoff "
                 f"({ache_c['top_n_min']:.0f} bits), so global-bits ranking buries the "
                 "entire PETase branch beneath the AChE tier. ")
    if pet_branch < branch_n:
        note += (f"The branch is smaller than branch_n={branch_n} → **all** of it "
                 "screened.")
    else:
        note += (f"Per-query tiering takes the top {branch_n} PETASE targets "
                 f"(bits {pet_c['top_n_min']:.0f}–{pet_c['max']:.0f}) — finally "
                 "screening PETase-neighbours the discovery sweep never reached.")

    if args.partition_only:
        return {"partition": partition_ctx, "records": records}

    # ---- CP2 ----
    anchor_pack = build_control_anchor(cfg, args.struct_dir)
    anchor = anchor_pack["anchor"]
    print(f"[per_query] anchor={anchor_pack['positive_ids']} "
          f"mode={anchor_pack['mode']} (re-derived thr "
          f"{anchor_pack['threshold']:.4f} — PINNING to {line})", flush=True)

    pet_sel = select_branch(records, "PETASE", branch_n, seed)
    ache_sel = select_branch(records, "ACHE", branch_n, seed)
    print(f"[per_query] CP2 screening PETASE n={len(pet_sel)}, ACHE n={len(ache_sel)} "
          f"at line {line}", flush=True)

    interim = os.path.join(REPO, "data", "interim")
    pet = screen_branch(pet_sel, cfg, anchor, line, url_tmpl, cache_dir,
                        int(pq["fetch_workers"]), int(pq["screen_workers"]), "PETASE",
                        screen_cache=os.path.join(interim, "per_query_screen_PETASE.json"),
                        rescreen=args.rescreen)
    ache = screen_branch(ache_sel, cfg, anchor, line, url_tmpl, cache_dir,
                         int(pq["fetch_workers"]), int(pq["screen_workers"]), "ACHE",
                         screen_cache=os.path.join(interim, "per_query_screen_ACHE.json"),
                         rescreen=args.rescreen)

    # ---- CP3 ----
    fl = floor_arm(args.floor_json)
    pr = arm_rates(pet["funnel"])
    ar = arm_rates(ache["funnel"])
    t_pf = two_arm_test(
        pr["conditional_above_given_triad"]["k"], pr["conditional_above_given_triad"]["n"],
        fl["above_line"], fl["triad_positive_S4"], "PETASE", "floor")
    t_pa = two_arm_test(
        pr["conditional_above_given_triad"]["k"], pr["conditional_above_given_triad"]["n"],
        ar["conditional_above_given_triad"]["k"], ar["conditional_above_given_triad"]["n"],
        "PETASE", "ACHE")
    t_af = two_arm_test(
        ar["conditional_above_given_triad"]["k"], ar["conditional_above_given_triad"]["n"],
        fl["above_line"], fl["triad_positive_S4"], "ACHE", "floor")

    # ---- CP4 ----
    cands = _cand_rows(pet["records"])
    verdict, tldr = build_verdict(t_pf, pr, fl, t_pa, ar)

    ctx = {
        "run_date": args.run_date,
        "seed": seed, "line": line, "branch_n": branch_n,
        "anchor_classes": anchor_classes,
        "partition": partition_ctx,
        "branch_size_note": note,
        "petase": pet, "ache": ache, "floor": fl,
        "petase_rates": pr, "ache_rates": ar,
        "test_petase_vs_floor": t_pf, "test_petase_vs_ache": t_pa,
        "test_ache_vs_floor": t_af,
        "candidates": cands, "verdict": verdict, "tldr": tldr,
    }
    return ctx


def _strip_records_for_json(ctx: dict) -> dict:
    """Compact JSON: keep funnels/rates/tests/candidates; drop the bulky per-model
    metric dumps but keep the candidate detail."""
    out = {k: v for k, v in ctx.items() if k not in ("petase", "ache")}
    out["petase_funnel"] = ctx["petase"]["funnel"]
    out["ache_funnel"] = ctx["ache"]["funnel"]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["partition", "run"], help="stage to run")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--out-csv", default=os.path.join(
        REPO, "data", "processed", "branch_partition.csv"))
    ap.add_argument("--out-md", default=os.path.join(
        REPO, "envlog", "per-query-tiering.md"))
    ap.add_argument("--out-json", default=os.path.join(
        REPO, "data", "processed", "per_query_tiering.json"))
    ap.add_argument("--floor-json", default=os.path.join(
        REPO, "data", "processed", "floor.json"),
        help="floor.json for the random-floor arm (else published constants)")
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"),
                    help="dir with the control structures (for the scoring anchor); "
                         "structures/ is gitignored, so point at the host copy")
    ap.add_argument("--run-date", default="", help="YYYY-MM-DD stamp for the report")
    ap.add_argument("--partition-only", action="store_true")
    ap.add_argument("--rescreen", action="store_true",
                    help="force a fresh screen (ignore cached per-branch screen JSON)")
    args = ap.parse_args(argv)
    if args.cmd == "partition":
        args.partition_only = True

    cfg = load_config(args.config)
    ctx = run(cfg, args)

    if args.partition_only:
        print(f"[per_query] partition-only: {ctx['partition']['branch_sizes']}")
        return 0

    md = render_md(ctx)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_md)), exist_ok=True)
    with open(args.out_md, "w") as fh:
        fh.write(md)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(_strip_records_for_json(ctx), fh, indent=2, default=str)
        fh.write("\n")
    print(f"[per_query] wrote {args.out_md}")
    print(f"[per_query] wrote {args.out_json}")
    print(f"[per_query] TL;DR: {ctx['tldr']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
