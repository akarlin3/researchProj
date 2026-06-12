"""Drive the Vast.ai S3 burst from the Mac: up -> fold -> down, in one command.

This is the LOCAL orchestrator that automates the manual steps in vast/sync.md. It
validates the S3 job manifest + shortlist (the contract from `proteus.s3_fold`),
then plans the burst cycle from `compute.vast_burst` in config:

  1. search interruptible CUDA offers (vastai)
  2. create the instance from the pushed fold image
  3. rsync the shortlist + manifest UP to the mounted input volume
  4. ssh-run the on-box runner (vast/run_fold.py) — ESMFold, pLDDT gate, resumable
  5. rsync the folded PDBs + pLDDT DOWN to structures/folded/
  6. destroy the instance (stop billing)

DRY-RUN BY DEFAULT: it prints the exact command plan and never touches Vast unless
`--execute` is passed (and even then the create/offer steps need an --offer and the
fold/sync steps need --ssh). The plan builder + input validation are pure and
unit-tested; the actual `vastai`/`rsync`/`ssh` calls are LOCAL-on-the-Mac.

ESMFold itself never runs here — this just ships the narrowed shortlist up and
pulls the models back so S4/S5 (proteus.screen) can resume locally.

Usage (from the repo root):
    # plan only (safe default):
    PYTHONPATH=src python -m proteus.launch \
        --manifest data/interim/s3_job_manifest.json \
        --shortlist data/interim/s2_shortlist.fasta
    # after picking an offer + the instance is up, run the data round-trip:
    PYTHONPATH=src python -m proteus.launch ... --ssh root@<host> --execute
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys

from proteus.utils import DEFAULT_CONFIG, REPO, load_config


def _count_fasta(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def validate_inputs(manifest_path: str, shortlist_path: str) -> dict:
    """Check the S3 hand-off is coherent before spending any GPU time. Returns
    {ok, errors, manifest} — errors is a list of human-readable problems."""
    errors = []
    if not os.path.exists(shortlist_path):
        errors.append(f"shortlist FASTA not found: {shortlist_path}")
    if not os.path.exists(manifest_path):
        errors.append(f"job manifest not found: {manifest_path}")
        return {"ok": False, "errors": errors, "manifest": None}

    try:
        manifest = json.load(open(manifest_path))
    except (json.JSONDecodeError, OSError) as exc:
        return {"ok": False, "errors": errors + [f"manifest unreadable: {exc}"],
                "manifest": None}

    if manifest.get("run_location") != "vast":
        errors.append("manifest run_location != 'vast' — refusing to launch a burst "
                      "for a non-Vast job")
    n_man = manifest.get("n_sequences")
    if not n_man:
        errors.append("manifest has 0 sequences — nothing to fold")
    if os.path.exists(shortlist_path):
        n_fa = _count_fasta(shortlist_path)
        if n_man is not None and n_fa != n_man:
            errors.append(f"shortlist has {n_fa} sequences but manifest says {n_man} "
                          "— regenerate the manifest from this shortlist")
    return {"ok": not errors, "errors": errors, "manifest": manifest}


def _burst_cfg(cfg: dict) -> dict:
    vb = dict(cfg.get("compute", {}).get("vast_burst", {}))
    vb.setdefault("gpu", "RTX_4090")
    vb.setdefault("num_gpus", 1)
    vb.setdefault("cuda_min", "12.4")
    vb.setdefault("disk_gb", 60)
    vb.setdefault("remote_in", "/data/proteus/in")
    vb.setdefault("remote_out", "/data/proteus/out")
    vb.setdefault("local_out", "structures/folded")
    vb.setdefault("registry_image", "")
    vb.setdefault("instance_pref", "interruptible")
    return vb


def build_plan(cfg: dict, manifest_path: str, shortlist_path: str,
               offer_id: str | None = None, instance_id: str | None = None,
               ssh: str | None = None, bid: str | None = None) -> list[dict]:
    """Build the ordered burst command plan. Unknown runtime ids appear as
    <OFFER_ID>/<INSTANCE_ID> placeholders. Each step: {name, cmd (argv list),
    note}. Pure — no side effects; this is what --dry-run prints and tests assert."""
    vb = _burst_cfg(cfg)
    image = vb["registry_image"] or "<REGISTRY_IMAGE — set compute.vast_burst.registry_image>"
    offer = offer_id or "<OFFER_ID>"
    inst = instance_id or "<INSTANCE_ID>"
    host = ssh or "<root@HOST>"
    remote_in, remote_out, local_out = vb["remote_in"], vb["remote_out"], vb["local_out"]
    interruptible = str(vb["instance_pref"]).lower() == "interruptible"

    man_name = os.path.basename(manifest_path)
    fa_name = os.path.basename(shortlist_path)
    runner = "/opt/proteus/run_fold.py"

    search = ["vastai", "search", "offers",
              f"gpu_name={vb['gpu']} num_gpus={vb['num_gpus']} "
              f"cuda_vers>={vb['cuda_min']} inet_down>200",
              "-o", "dph"]
    if interruptible:
        search.append("--interruptible")

    create = ["vastai", "create", "instance", offer, "--image", image,
              "--disk", str(vb["disk_gb"])]
    if interruptible:
        create += ["--bid", bid or "<BID_PRICE>"]   # interruptible bids; on-demand omits

    up = ["rsync", "-avP", shortlist_path, manifest_path, f"{host}:{remote_in}/"]
    fold = ["ssh", host, "python3", runner,
            "--manifest", f"{remote_in}/{man_name}",
            "--fasta", f"{remote_in}/{fa_name}",
            "--out", f"{remote_out}/"]
    down = ["rsync", "-avP", f"{host}:{remote_out}/", f"{local_out}/"]
    destroy = ["vastai", "destroy", "instance", inst]

    return [
        {"name": "search_offers", "cmd": search,
         "note": "pick a cheap offer id for --offer"},
        {"name": "create_instance", "cmd": create,
         "note": "interruptible: set --bid; on-demand: drop --bid. Note the instance id."},
        {"name": "upload", "cmd": up,
         "note": "ship ONLY the shortlist + manifest up (mounted volume)"},
        {"name": "fold", "cmd": fold,
         "note": "ESMFold on CUDA; resumable + pLDDT-gated (vast/run_fold.py)"},
        {"name": "download", "cmd": down,
         "note": "pull folded PDBs + pLDDT back; then run proteus.screen"},
        {"name": "destroy_instance", "cmd": destroy,
         "note": "stop billing once results are down"},
    ]


def render_plan(plan: list[dict]) -> str:
    lines = ["[launch] Vast burst plan (DRY-RUN — nothing executed):"]
    for i, step in enumerate(plan, 1):
        lines.append(f"  {i}. {step['name']}: {step['note']}")
        lines.append(f"     $ {' '.join(shlex.quote(c) for c in step['cmd'])}")
    return "\n".join(lines)


# Steps that need a real ssh host before they can run.
_NEEDS_SSH = {"upload", "fold", "download"}


def execute_plan(plan: list[dict], ssh: str | None, only: set[str] | None = None) -> int:
    """Execute selected steps (LOCAL, on the Mac). Skips steps whose runtime ids are
    still placeholders, or ssh steps without --ssh. Returns a process exit code."""
    for step in plan:
        if only and step["name"] not in only:
            continue
        if step["name"] in _NEEDS_SSH and not ssh:
            print(f"[launch] skip {step['name']}: needs --ssh root@<host>", file=sys.stderr)
            continue
        if any(tok.startswith("<") and tok.endswith(">") for tok in step["cmd"]):
            print(f"[launch] skip {step['name']}: unresolved placeholder "
                  f"({' '.join(step['cmd'])})", file=sys.stderr)
            continue
        print(f"[launch] run {step['name']}: {' '.join(shlex.quote(c) for c in step['cmd'])}")
        proc = subprocess.run(step["cmd"])
        if proc.returncode != 0:
            print(f"[launch] {step['name']} failed (rc={proc.returncode}) — stopping",
                  file=sys.stderr)
            return proc.returncode
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    _interim = os.path.join(REPO, "data", "interim")
    ap.add_argument("--manifest", default=os.path.join(_interim, "s3_job_manifest.json"))
    ap.add_argument("--shortlist", default=os.path.join(_interim, "s2_shortlist.fasta"))
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--offer", default=None, help="vastai offer id (from search_offers)")
    ap.add_argument("--instance", default=None, help="vastai instance id (for destroy)")
    ap.add_argument("--ssh", default=None, help="ssh target, e.g. root@1.2.3.4")
    ap.add_argument("--bid", default=None, help="interruptible bid price ($/hr)")
    ap.add_argument("--execute", action="store_true",
                    help="actually run the plan (default: dry-run print only)")
    ap.add_argument("--only", default=None,
                    help="comma-separated step names to execute (e.g. upload,fold,download)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    v = validate_inputs(args.manifest, args.shortlist)
    for e in v["errors"]:
        print(f"[launch][invalid] {e}", file=sys.stderr)
    if not v["ok"]:
        print("[launch] fix the manifest/shortlist above before bursting.", file=sys.stderr)
        return 2
    print(f"[launch] inputs OK: {v['manifest']['n_sequences']} sequence(s) to fold "
          f"({os.path.relpath(args.shortlist, os.getcwd())}).")

    plan = build_plan(cfg, args.manifest, args.shortlist, offer_id=args.offer,
                      instance_id=args.instance, ssh=args.ssh, bid=args.bid)
    if not args.execute:
        print(render_plan(plan))
        print("[launch] dry-run only — re-run with --execute "
              "(and --offer/--instance/--ssh) to burst.")
        return 0

    only = set(args.only.split(",")) if args.only else None
    return execute_plan(plan, args.ssh, only=only)


if __name__ == "__main__":
    raise SystemExit(main())
