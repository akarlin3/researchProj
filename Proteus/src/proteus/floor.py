"""Random-floor measurement + lift decomposition (CP1-CP3 of the floor run).

WHY: the Atlas pilot's floor was 1 event in 60 (1.7%) — a point estimate whose own
95% CI runs ~0-9%, useless as a null. The enriched fold-class sweep
(atlas-sweep/2026-06-14) cleared the widened line at 23/300 = 7.7%. To know whether
that lift is real — and whether it comes from S4 (the fold-class/triad search) or
S5 (the cleft line) — we need a PROPERLY-SAMPLED floor drawn from the SAME
population the enriched set came from, screened through the SAME pipeline at the
SAME line.

WHAT THIS DOES (Mac-only, off the existing GCS artifacts; no GCE, no re-search):

  CP1  Uniformly sample `floor.n` accessions from the HQ clust30 universe — the
       Foldseek DB's own representative list (`highquality_clust30.lookup`,
       36.99 M entries) — by drawing seeded random byte offsets across the whole
       1.03 GB file via HTTP Range. SAME population as the enriched set, NO
       fold-class pre-filter: these are random Atlas proteins, the null. The draw
       is reproducible (seed = config random_seed) and recorded to
       data/processed/floor_sample.txt.

  CP2  Fetch each structure (ESMFold Atlas V0 model) and screen it through the
       UNCHANGED screen pipeline: S4 geometry -> (fpocket only on triad+) -> S5
       cleft -> control-anchored composite, at the enriched run's widened line.
       The enriched 300 are optionally re-screened in the SAME pass (same anchor,
       same line, same fpocket session) for a jitter-immune matched comparison
       (fpocket is non-deterministic run-to-run; S4 geometry is deterministic).

  CP3  Tightened floor (Wilson CI), corrected significance vs the enriched lift
       (Fisher + two-proportion z, rate ratio + CI), and the lift DECOMPOSITION:
         - S4 / triad rate:      triad+ fraction, random vs enriched (~99.3%)
         - S5 / conditional rate: above-line GIVEN triad+, random vs enriched
                                   (23/298 = 7.7%) -- the load-bearing number.
       Writes envlog/floor-measurement.md + data/processed/floor_*.{csv,json}.

Fetch path: the Atlas hosted Foldseek search is down (503) and the prebuilt-DB
fetch on GCE used the local foldcomp DB; on the Mac we fetch the same V0 models by
accession from `fetchPredictedStructure` (verified live and structurally identical
to the GCE foldcomp models — triad residues reproduce exactly). The foldcomp worker
(Range-served DB files) is the documented fallback. See envlog/atlas-endpoints.md.

Usage, from the repo root:
    PYTHONPATH=src python -m proteus.floor \
        --enriched-csv data/interim/atlas_sweep/atlas_candidates.csv \
        --enriched-funnel data/interim/atlas_sweep/funnel.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

from proteus.atlas_screen import resolve_operating_point
from proteus.screen import screen_model
from proteus.utils import DEFAULT_CONFIG, REPO, get_seed, load_config

MGYP_RE = re.compile(rb"MGYP\d{6,}")
_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
       "AppleWebKit/537.36 (KHTML, like Gecko) Proteus/floor")
_FETCH_BASE = "https://api.esmatlas.com/fetchPredictedStructure"


# --------------------------------------------------------------------------- #
# HTTP (stdlib; returns status so callers branch on 403/transient)
# --------------------------------------------------------------------------- #
def _http(url: str, *, rng: tuple[int, int] | None = None, timeout: float = 60,
          max_bytes: int | None = None) -> tuple[int, bytes, str]:
    headers = {"User-Agent": _UA}
    if rng is not None:
        headers["Range"] = f"bytes={rng[0]}-{rng[1]}"
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read(max_bytes) if max_bytes else resp.read()
            return resp.status, body, resp.geturl()
    except urllib.error.HTTPError as exc:
        try:
            body = exc.read()
        except Exception:  # noqa: BLE001
            body = b""
        return exc.code, body, url


# --------------------------------------------------------------------------- #
# CP1 — uniform sample from the HQ clust30 universe (deterministic + parallel)
# --------------------------------------------------------------------------- #
def resolve_lookup(lookup_url: str, timeout: float, log=print) -> tuple[str, int]:
    """Follow the worker redirect once to the Range-served origin and read its size.
    Returns (final_url, content_length)."""
    status, _body, final = _http(lookup_url, rng=(0, 0), timeout=timeout)
    if status not in (200, 206):
        raise RuntimeError(f"lookup HEAD/range -> HTTP {status} for {lookup_url}")
    # content-length of the *range* isn't the file size; re-HEAD the resolved URL.
    req = urllib.request.Request(final, headers={"User-Agent": _UA}, method="HEAD")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        size = int(resp.headers.get("Content-Length", "0"))
    if size <= 0:
        raise RuntimeError(f"lookup size unresolved (Content-Length={size})")
    log(f"[floor:cp1] universe lookup -> {final} ({size:,} bytes)")
    return final, size


def _accession_at(final_url: str, offset: int, win: int, timeout: float,
                  retries: int = 6) -> str | None:
    """Range-read a window at `offset`, skip the (possibly partial) first line, and
    return the first complete MGYP accession found. Retries transient transport
    failures (S3 resets connections under load) with backoff so each offset
    DETERMINISTICALLY resolves to its accession — keeping the seeded sample
    reproducible regardless of network weather."""
    for attempt in range(retries):
        try:
            status, body, _ = _http(final_url, rng=(offset, offset + win - 1),
                                    timeout=timeout)
        except Exception:  # noqa: BLE001 — transport reset/timeout: retry
            time.sleep(0.3 * (attempt + 1))
            continue
        if status not in (200, 206) or not body:
            time.sleep(0.3 * (attempt + 1))
            continue
        nl = body.find(b"\n")
        if nl == -1:
            return None
        m = MGYP_RE.search(body[nl + 1:])
        return m.group(0).decode() if m else None
    return None


def sample_universe(final_url: str, size: int, n: int, seed: int, *,
                    workers: int = 16, win: int = 256, timeout: float = 60,
                    log=print) -> tuple[list[str], dict]:
    """Draw `n` unique accessions uniformly over the lookup by seeded random byte
    offsets. Deterministic regardless of fetch parallelism: a seeded RNG produces
    the offset stream, windows are fetched in parallel, then walked in offset-stream
    order so the resulting set depends only on (seed, file content)."""
    rng = random.Random(seed)
    picked: list[str] = []
    seen: set[str] = set()
    total_offsets = 0
    batch = max(n + n // 4, 64)
    s3_workers = min(workers, 8)  # S3 resets connections under heavier fan-out
    while len(picked) < n:
        offsets = [rng.randrange(0, size - win) for _ in range(batch)]
        total_offsets += len(offsets)
        results: dict[int, str | None] = {}
        with ThreadPoolExecutor(max_workers=s3_workers) as ex:
            futs = {ex.submit(_accession_at, final_url, off, win, timeout): i
                    for i, off in enumerate(offsets)}
            for fut in as_completed(futs):
                results[futs[fut]] = fut.result()
        for i in range(len(offsets)):  # deterministic offset-stream order
            acc = results.get(i)
            if acc and acc not in seen:
                seen.add(acc)
                picked.append(acc)
                if len(picked) == n:
                    break
        log(f"[floor:cp1] sampled {len(picked)}/{n} unique "
            f"(from {total_offsets} offsets)")
        batch = max((n - len(picked)) * 2, 32)
    meta = {"seed": seed, "n": n, "lookup_url": final_url, "lookup_bytes": size,
            "offsets_drawn": total_offsets, "window_bytes": win}
    return picked, meta


def write_sample(accs: list[str], meta: dict, plddt: dict, path: str) -> None:
    """Record the sampled accession list + provenance to floor_sample.txt."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("# Proteus random-floor sample — HQ clust30 uniform draw\n")
        fh.write(f"# seed={meta['seed']}  n={meta['n']}  "
                 f"offsets_drawn={meta['offsets_drawn']}\n")
        fh.write(f"# universe={meta['lookup_url']}\n")
        fh.write(f"# universe_bytes={meta['lookup_bytes']}  "
                 f"window_bytes={meta['window_bytes']}\n")
        fh.write("# Uniform over byte offsets of the Foldseek DB representative "
                 "list (~36.99M entries); reproducible from the seed.\n")
        fh.write("# columns: accession\tmean_plddt(0-1)\tfetch_status\n")
        for acc in accs:
            rec = plddt.get(acc, {})
            mp = rec.get("mean_plddt")
            mp_s = "" if mp is None else f"{mp:.4f}"
            fh.write(f"{acc}\t{mp_s}\t{rec.get('fetch_status', 'unfetched')}\n")


# --------------------------------------------------------------------------- #
# CP2 — fetch (ESMFold V0 model) + screen through the UNCHANGED pipeline
# --------------------------------------------------------------------------- #
def _mean_plddt(pdb_bytes: bytes) -> tuple[float | None, int]:
    vals = []
    for line in pdb_bytes.splitlines():
        if line.startswith(b"ATOM") and line[12:16].strip() == b"CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    if not vals:
        return None, 0
    return round(sum(vals) / len(vals), 4), len(vals)


def fetch_one(acc: str, dest_dir: str, *, timeout: float = 60,
              retries: int = 3) -> dict:
    """Fetch one Atlas V0 model by accession into dest_dir. Returns a record with
    the path, mean pLDDT (0-1), residue count, and fetch_status (200 | <code> |
    error:*). Retries transient failures with backoff."""
    rec = {"accession": acc, "pdb": None, "mean_plddt": None, "n_res": 0,
           "fetch_status": None}
    url = f"{_FETCH_BASE}/{acc}.pdb"
    for attempt in range(retries):
        try:
            status, body, _ = _http(url, timeout=timeout)
        except Exception as exc:  # noqa: BLE001 — transport; retry then give up
            rec["fetch_status"] = f"error:{type(exc).__name__}"
            time.sleep(0.5 * (attempt + 1))
            continue
        rec["fetch_status"] = status
        if status == 200 and (body.startswith(b"HEADER") or b"ATOM" in body[:400]):
            path = os.path.join(dest_dir, f"{acc}.pdb")
            with open(path, "wb") as fh:
                fh.write(body)
            rec["pdb"] = path
            rec["mean_plddt"], rec["n_res"] = _mean_plddt(body)
            return rec
        if status == 404:
            return rec  # not in the served subset — deterministic miss, no retry
        # 403 can be a true "not in served HQ subset" OR transient throttling under
        # load — retry with backoff; only the final attempt's status is recorded.
        time.sleep(0.5 * (attempt + 1))
    return rec


def fetch_many(accs: list[str], dest_dir: str, *, workers: int = 16,
               timeout: float = 60, label: str = "", log=print) -> dict[str, dict]:
    os.makedirs(dest_dir, exist_ok=True)
    out: dict[str, dict] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_one, a, dest_dir, timeout=timeout): a for a in accs}
        for fut in as_completed(futs):
            rec = fut.result()
            out[rec["accession"]] = rec
            done += 1
            if done % 100 == 0 or done == len(accs):
                ok = sum(1 for r in out.values() if r["pdb"])
                log(f"[floor:cp2] fetch{(' ' + label) if label else ''} "
                    f"{done}/{len(accs)} ({ok} ok)")
    return out


def screen_many(fetched: dict[str, dict], cfg: dict, anchor: dict, line: float, *,
                workers: int = 8, label: str = "", log=print) -> list[dict]:
    """Screen every fetched structure through screen_model at `line`, in parallel.
    fpocket runs in an isolated temp dir per call (race-free)."""
    items = [(a, r) for a, r in fetched.items() if r.get("pdb")]

    def _run(acc, rec):
        out = screen_model(rec["pdb"], cfg, anchor, line,
                           mean_plddt=rec.get("mean_plddt"), cand_id=acc)
        out["accession"] = acc
        out["mean_plddt"] = rec.get("mean_plddt")
        out["n_res"] = rec.get("n_res")
        return out

    results: list[dict] = []
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_run, a, r) for a, r in items]
        for fut in as_completed(futs):
            results.append(fut.result())
            done += 1
            if done % 100 == 0 or done == len(items):
                tri = sum(1 for r in results if r["triad_found"])
                hit = sum(1 for r in results if r.get("petase_like_hit"))
                log(f"[floor:cp2] screen{(' ' + label) if label else ''} "
                    f"{done}/{len(items)} (triad+={tri}, above-line={hit})")
    return results


def funnel(results: list[dict], n_attempted: int, n_fetched: int) -> dict:
    triad = [r for r in results if r["triad_found"]]
    pocket = [r for r in triad if r["pocket_ok"]]
    above = [r for r in results if r.get("petase_like_hit")]
    return {
        "attempted": n_attempted,
        "fetched": n_fetched,
        "screened": len(results),
        "triad_positive_S4": len(triad),
        "pocket_ok_S5": len(pocket),
        "above_line": len(above),
        "above_line_accessions": sorted(r["accession"] for r in above),
    }


# --------------------------------------------------------------------------- #
# CP3 — statistics
# --------------------------------------------------------------------------- #
def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score 95% CI for a binomial proportion (robust at small k/extremes)."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    z2 = z * z
    denom = 1 + z2 / n
    centre = (p + z2 / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def two_proportion_z(k1: int, n1: int, k2: int, n2: int) -> dict:
    """Two-proportion z-test (pooled). Returns z and two-sided p (normal approx)."""
    if n1 == 0 or n2 == 0:
        return {"z": None, "p_value": None}
    p1, p2 = k1 / n1, k2 / n2
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    if se == 0:
        return {"z": None, "p_value": None}
    z = (p1 - p2) / se
    p_value = math.erfc(abs(z) / math.sqrt(2))  # two-sided
    return {"z": z, "p_value": p_value}


def rate_ratio_ci(k1: int, n1: int, k2: int, n2: int,
                  z: float = 1.959963984540054) -> dict:
    """Rate ratio (group1/group2) of two proportions with a log (Katz) 95% CI.
    A 0.5 continuity correction is applied if any cell is zero."""
    if n1 == 0 or n2 == 0 or k2 == 0 and k1 == 0:
        return {"ratio": None, "ci": [None, None]}
    a, b = k1, k2
    if a == 0 or b == 0:
        a, b = k1 + 0.5, k2 + 0.5
    p1, p2 = a / n1, b / n2
    rr = p1 / p2
    se_ln = math.sqrt((1 - p1) / a + (1 - p2) / b)
    lo = rr * math.exp(-z * se_ln)
    hi = rr * math.exp(z * se_ln)
    return {"ratio": rr, "ci": [lo, hi]}


def fisher_exact(k1: int, n1: int, k2: int, n2: int) -> float | None:
    """Two-sided Fisher's exact p on the 2x2 [[k1, n1-k1], [k2, n2-k2]]."""
    try:
        from scipy.stats import fisher_exact as _fe  # noqa: PLC0415
    except Exception:  # noqa: BLE001
        return None
    _odds, p = _fe([[k1, n1 - k1], [k2, n2 - k2]])
    return float(p)


def compare(label: str, k1: int, n1: int, k2: int, n2: int) -> dict:
    """A full enriched-vs-floor comparison block for one rate."""
    return {
        "label": label,
        "group1": {"k": k1, "n": n1, "rate": (k1 / n1) if n1 else None,
                   "wilson95": list(wilson_ci(k1, n1))},
        "group2": {"k": k2, "n": n2, "rate": (k2 / n2) if n2 else None,
                   "wilson95": list(wilson_ci(k2, n2))},
        "fisher_p": fisher_exact(k1, n1, k2, n2),
        "two_proportion_z": two_proportion_z(k1, n1, k2, n2),
        "rate_ratio_g1_over_g2": rate_ratio_ci(k1, n1, k2, n2),
    }


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def load_enriched(csv_path: str) -> list[str]:
    with open(csv_path, newline="") as fh:
        return [row["accession"] for row in csv.DictReader(fh) if row.get("accession")]


def run(cfg: dict, struct_dir: str, enriched_csv: str, enriched_funnel: str,
        out_prefix: str, sample_txt: str, hits_dir: str, log=print) -> dict:
    fc = cfg.get("floor", {})
    n = int(fc.get("n", 1500))
    line = float(fc.get("enriched_line", -1.1587))
    fetch_workers = int(fc.get("fetch_workers", 16))
    screen_workers = int(fc.get("screen_workers", 8))
    rescreen = bool(fc.get("rescreen_enriched", True))
    lookup_url = fc.get("lookup_url")
    timeout = float(cfg.get("atlas", {}).get("request_timeout_s", 60))
    seed = get_seed(cfg)

    with open(enriched_funnel) as fh:
        ef = json.load(fh)
    enr_screened = int(ef["screened"])
    enr_triad = int(ef["triad_positive_S4"])
    enr_above = int(ef["above_widened_line"])
    enr_line = ef.get("threshold")
    log(f"[floor] enriched (GCE canonical): screened={enr_screened} "
        f"triad+={enr_triad} above-line={enr_above} @ line={enr_line}")
    if abs(line - float(enr_line)) > 1e-6:
        log(f"[floor][WARN] config line {line} != enriched funnel line {enr_line}; "
            "using config line for the floor (set floor.enriched_line to match).")

    # Anchor: build ONCE from the controls (reuses calibrate untouched). Pin the
    # decision line to the enriched run's value for an apples-to-apples comparison.
    op = resolve_operating_point(cfg, struct_dir)
    anchor = op["anchor"]
    log(f"[floor] control anchor={op['positive_ids']} mode={op['mode']} "
        f"separated={op['separated']} margin={op['margin']}; "
        f"local widened re-derivation={op['threshold']:.4f} (fpocket jitter), "
        f"SCREENING PINNED at {line:.4f}")

    # ---- CP1: sample the universe -----------------------------------------
    final_url, size = resolve_lookup(lookup_url, timeout, log=log)
    accs, meta = sample_universe(final_url, size, n, seed,
                                 workers=fetch_workers, timeout=timeout, log=log)

    # ---- CP2: fetch + screen the floor ------------------------------------
    floor_dir = os.path.join(hits_dir, "floor")
    fetched = fetch_many(accs, floor_dir, workers=fetch_workers, timeout=timeout,
                         label="floor", log=log)
    write_sample(accs, meta, fetched, sample_txt)
    n_fetch_ok = sum(1 for r in fetched.values() if r["pdb"])
    log(f"[floor:cp2] floor fetched {n_fetch_ok}/{len(accs)} "
        f"({len(accs) - n_fetch_ok} miss/err)")
    floor_results = screen_many(fetched, cfg, anchor, line,
                                workers=screen_workers, label="floor", log=log)
    floor_funnel = funnel(floor_results, len(accs), n_fetch_ok)

    # ---- CP2b: re-screen the enriched set locally (jitter-immune match) ----
    enriched_local_funnel = None
    enriched_local_results = []
    if rescreen:
        enr_accs = load_enriched(enriched_csv)
        log(f"[floor:cp2] re-screening enriched set locally: {len(enr_accs)} accs")
        enr_dir = os.path.join(hits_dir, "enriched")
        enr_fetched = fetch_many(enr_accs, enr_dir, workers=fetch_workers,
                                 timeout=timeout, label="enriched", log=log)
        enr_ok = sum(1 for r in enr_fetched.values() if r["pdb"])
        enriched_local_results = screen_many(enr_fetched, cfg, anchor, line,
                                             workers=screen_workers,
                                             label="enriched", log=log)
        enriched_local_funnel = funnel(enriched_local_results, len(enr_accs), enr_ok)

    # ---- CP3: statistics ---------------------------------------------------
    fN = floor_funnel["screened"]
    fTriad = floor_funnel["triad_positive_S4"]
    fAbove = floor_funnel["above_line"]
    fCondAbove = sum(1 for r in floor_results
                     if r["triad_found"] and r.get("petase_like_hit"))

    stats = {
        # 1. tightened floor (overall above-line rate)
        "tightened_floor": {
            "k": fAbove, "n": fN, "rate": (fAbove / fN) if fN else None,
            "wilson95": list(wilson_ci(fAbove, fN)),
            "old_thin_floor": ef.get("random_floor"),
        },
        # 2. corrected significance: enriched overall vs floor overall
        "overall_above_line": compare("enriched vs floor: above-line / screened",
                                      enr_above, enr_screened, fAbove, fN),
        # 3a. S4 / triad rate (fold-class contribution)
        "S4_triad_rate": compare("enriched vs floor: triad+ / screened",
                                 enr_triad, enr_screened, fTriad, fN),
        # 3b. S5 / conditional rate (cleft-line contribution) -- load-bearing
        "S5_conditional_rate": compare(
            "enriched vs floor: above-line / triad+ (GIVEN a triad)",
            enr_above, enr_triad, fCondAbove, fTriad),
    }
    if enriched_local_funnel is not None:
        eN = enriched_local_funnel["screened"]
        eT = enriched_local_funnel["triad_positive_S4"]
        eA = enriched_local_funnel["above_line"]
        eCond = sum(1 for r in enriched_local_results
                    if r["triad_found"] and r.get("petase_like_hit"))
        stats["matched_local"] = {
            "note": "enriched re-screened in the SAME local pass (same anchor, "
                    "line, fpocket session) -> jitter-immune vs the floor.",
            "overall_above_line": compare("enriched(local) vs floor: above/screened",
                                          eA, eN, fAbove, fN),
            "S4_triad_rate": compare("enriched(local) vs floor: triad+/screened",
                                     eT, eN, fTriad, fN),
            "S5_conditional_rate": compare(
                "enriched(local) vs floor: above-line / triad+",
                eCond, eT, fCondAbove, fTriad),
            "enriched_local_funnel": enriched_local_funnel,
            "reproduces_canonical": {"canonical_above": enr_above,
                                     "local_above": eA,
                                     "canonical_triad": enr_triad,
                                     "local_triad": eT},
        }

    summary = {
        "seed": seed, "line": line, "enriched_line": enr_line,
        "local_widened_rederivation": op["threshold"],
        "anchor_mode": op["mode"], "anchor_ids": op["positive_ids"],
        "sample_meta": meta,
        "floor_funnel": floor_funnel,
        "enriched_funnel_canonical": {
            "screened": enr_screened, "triad_positive_S4": enr_triad,
            "above_line": enr_above, "line": enr_line,
            "conditional_above_given_triad": enr_above / enr_triad,
        },
        "floor_conditional_above_given_triad": {
            "k": fCondAbove, "n": fTriad,
            "rate": (fCondAbove / fTriad) if fTriad else None,
            "wilson95": list(wilson_ci(fCondAbove, fTriad)),
        },
        "stats": stats,
        "floor_candidates": sorted(
            [r for r in floor_results if r.get("petase_like_hit")],
            key=lambda r: (r["composite"] is not None, r["composite"] or 0),
            reverse=True),
    }

    # outputs
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    with open(out_prefix + ".json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
        fh.write("\n")
    _write_floor_csv(floor_results, out_prefix + ".csv")
    log(f"[floor] floor funnel: {floor_funnel['screened']} screened -> "
        f"{floor_funnel['triad_positive_S4']} triad+ -> "
        f"{floor_funnel['above_line']} above the {line} line")
    return summary


def _write_floor_csv(results: list[dict], path: str) -> None:
    cols = ["accession", "mean_plddt", "n_res", "triad_found", "catalytic_ser",
            "his", "acid", "pocket_ok", "composite", "above_threshold",
            "petase_like_hit"]
    rows = sorted(results, key=lambda r: (r["composite"] is not None,
                                          r["composite"] or 0), reverse=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--struct-dir", default=os.path.join(REPO, "structures"))
    ap.add_argument("--enriched-csv",
                    default=os.path.join(REPO, "data", "interim", "atlas_sweep",
                                         "atlas_candidates.csv"))
    ap.add_argument("--enriched-funnel",
                    default=os.path.join(REPO, "data", "interim", "atlas_sweep",
                                         "funnel.json"))
    ap.add_argument("--out", default=os.path.join(REPO, "data", "processed", "floor"))
    ap.add_argument("--sample-out",
                    default=os.path.join(REPO, "data", "processed", "floor_sample.txt"))
    ap.add_argument("--hits-dir", default=os.path.join(REPO, "structures", "floor_hits"))
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if "floor" not in cfg:
        print("config has no `floor` block — add it (see config/proteus.yaml).",
              file=sys.stderr)
        return 2
    for p in (args.enriched_csv, args.enriched_funnel):
        if not os.path.exists(p):
            print(f"missing enriched artifact: {p} (pull it from "
                  "gs://projproteus-fold/atlas-sweep/2026-06-14/)", file=sys.stderr)
            return 2
    run(cfg, args.struct_dir, args.enriched_csv, args.enriched_funnel,
        args.out, args.sample_out, args.hits_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
