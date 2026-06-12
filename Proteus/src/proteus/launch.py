"""Drive the GCE S3 burst from the Mac: up -> fold -> down, in one command.

This is the LOCAL orchestrator that automates the manual steps in gce/sync.md. It
validates the S3 job manifest + shortlist (the contract from `proteus.s3_fold`),
then plans the burst cycle from `compute.gce_burst` in config:

  1. gcloud storage cp the shortlist + manifest UP to the staging bucket
  2. gcloud create the fold VM (SPOT or on-demand; GPU Deep-Learning image, or CPU COS)
  3. gcloud ssh -> on the VM: pull inputs from the bucket, `docker run` the fold image
     (gce/run_fold.py — ESMFold, pLDDT gate, resumable), push outputs to the bucket
  4. gcloud storage cp the folded PDBs + pLDDT DOWN to structures/folded/
  5. gcloud delete the VM (stop billing)

Outputs are staged to a GCS bucket, so a preempted SPOT VM doesn't lose finished
models. DRY-RUN BY DEFAULT: it prints the exact command plan and never touches GCP
unless `--execute` is passed. The plan builder + input validation are pure and
unit-tested; the actual `gcloud` calls are LOCAL-on-the-Mac.

ESMFold itself never runs here — this just ships the narrowed shortlist up and
pulls the models back so S4/S5 (proteus.screen) can resume locally.

Usage (from the repo root):
    # plan only (safe default):
    PYTHONPATH=src python -m proteus.launch \
        --manifest data/interim/s3_job_manifest.json \
        --shortlist data/interim/s2_shortlist.fasta
    # once project / bucket / image are set in config, run the burst:
    PYTHONPATH=src python -m proteus.launch ... --execute
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

    if manifest.get("run_location") != "gce":
        errors.append("manifest run_location != 'gce' — refusing to launch a burst "
                      "for a non-GCE job")
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
    vb = dict(cfg.get("compute", {}).get("gce_burst", {}))
    vb.setdefault("project", "")
    vb.setdefault("zone", "us-central1-a")
    vb.setdefault("instance_name", "proteus-fold")
    vb.setdefault("machine_type", "g2-standard-8")
    vb.setdefault("accelerator", "nvidia-l4")
    vb.setdefault("accelerator_count", 1)
    vb.setdefault("image", "")
    vb.setdefault("boot_disk_gb", 100)
    vb.setdefault("spot", True)
    vb.setdefault("bucket", "")
    vb.setdefault("vm_workdir", "/tmp/proteus")  # writable host dir (COS root fs is read-only)
    vb.setdefault("remote_in", "/data/proteus/in")
    vb.setdefault("remote_out", "/data/proteus/out")
    vb.setdefault("local_out", "structures/folded")
    # Host VM the fold container runs on — a Deep Learning image with CUDA + docker.
    vb.setdefault("vm_image_family", "common-cu123")
    vb.setdefault("vm_image_project", "deeplearning-platform-release")
    return vb


def build_plan(cfg: dict, manifest_path: str, shortlist_path: str) -> list[dict]:
    """Build the ordered GCE burst command plan. Unset config (project/bucket/image)
    appears as <PROJECT>/<gs://BUCKET>/<IMAGE> placeholders. Each step: {name, cmd
    (argv list), note}. Pure — no side effects; this is what --dry-run prints and
    tests assert."""
    vb = _burst_cfg(cfg)
    project = vb["project"] or "<PROJECT — set compute.gce_burst.project>"
    bucket = (vb["bucket"] or "<gs://BUCKET — set compute.gce_burst.bucket>").rstrip("/")
    image = vb["image"] or "<IMAGE — set compute.gce_burst.image>"
    name = vb["instance_name"]
    zone = vb["zone"]
    remote_in, remote_out, local_out = vb["remote_in"], vb["remote_out"], vb["local_out"]
    vm_workdir = vb.get("vm_workdir") or "/tmp/proteus"  # host dir mounted to /data/proteus
    spot = bool(vb["spot"])
    prov = "SPOT" if spot else "on-demand"
    accel = (str(vb.get("accelerator") or "")).strip()
    n_accel = int(vb.get("accelerator_count", 0) or 0)
    gpu = bool(accel) and accel.lower() != "none" and n_accel > 0

    man_name = os.path.basename(manifest_path)
    fa_name = os.path.basename(shortlist_path)
    sdk = "google/cloud-sdk:slim"  # COS ships neither gcloud nor gsutil — run them in a container

    # stage_up runs LOCALLY (the Mac has gcloud) -> the bucket's in/ prefix.
    stage_up = ["gcloud", "storage", "cp", shortlist_path, manifest_path, f"{bucket}/in/"]

    create = [
        "gcloud", "compute", "instances", "create", name,
        "--project", project, "--zone", zone,
        "--machine-type", vb["machine_type"],
        "--boot-disk-size", f"{vb['boot_disk_gb']}GB",
        "--scopes", "cloud-platform",  # GCS (gcloud storage) + docker pull (Artifact Registry)
    ]
    if gpu:
        # GPU fold: attach the accelerator on a Deep Learning (CUDA + docker) image.
        create += ["--accelerator", f"type={accel},count={vb['accelerator_count']}",
                   "--maintenance-policy", "TERMINATE",
                   "--image-family", vb.get("vm_image_family", "common-cu123"),
                   "--image-project", vb.get("vm_image_project", "deeplearning-platform-release"),
                   "--metadata", "install-nvidia-driver=True"]
        create_note = f"{prov} GPU VM (Deep Learning image: CUDA + docker)"
        gpu_flag, device = "--gpus all ", "cuda"
    else:
        # CPU fold: Container-Optimized OS ships docker; no GPU, no driver. ESMFold on
        # CPU is slow + RAM-heavy but fine for the small narrowed shortlist.
        create += ["--image-family", "cos-stable", "--image-project", "cos-cloud"]
        create_note = f"{prov} CPU VM (Container-Optimized OS: docker); ESMFold on CPU"
        gpu_flag, device = "", "cpu"
    if spot:
        create += ["--provisioning-model", "SPOT",
                   "--instance-termination-action", "DELETE"]

    # On the VM: pull inputs from the bucket (via the cloud-sdk container, since COS has
    # no gcloud), fold in the container, push outputs back (same container). The VM's
    # service account (--scopes cloud-platform) authenticates these automatically.
    # The on-VM working dir lives on a WRITABLE path (COS mounts / read-only); it is
    # bind-mounted into every container at /data/proteus, so the in-container paths
    # (remote_in/remote_out) stay /data/proteus/{in,out}.
    mount = f"-v {vm_workdir}:/data/proteus"
    # COS docker is not authenticated to Artifact Registry, so pulling the fold image
    # fails with "Unauthenticated request". Log docker in first using an access token
    # from the VM's service account (printed by the cloud-sdk container). Only for
    # Google registries (AR/GCR); a public image needs no login.
    registry = image.split("/")[0] if "/" in image else ""
    login = ""
    if any(host in registry for host in ("pkg.dev", "gcr.io")):
        login = (f'docker login -u oauth2accesstoken '
                 f'-p "$(docker run --rm {sdk} gcloud auth print-access-token)" '
                 f"https://{registry} && ")
    remote = (
        f"mkdir -p {vm_workdir}/in {vm_workdir}/out && "
        f"docker run --rm {mount} {sdk} "
        f"gcloud storage cp {bucket}/in/* {remote_in}/ && "
        f"{login}"
        f"docker run {gpu_flag}{mount} {image} "
        f"--manifest {remote_in}/{man_name} --fasta {remote_in}/{fa_name} "
        f"--out {remote_out}/ --device {device} && "
        f"docker run --rm {mount} {sdk} "
        f"gcloud storage cp --recursive {remote_out}/* {bucket}/out/"
    )
    fold = ["gcloud", "compute", "ssh", name, "--project", project, "--zone", zone,
            "--command", remote]

    # stage_down runs LOCALLY -> pull the folded models from the bucket.
    stage_down = ["gcloud", "storage", "cp", "--recursive", f"{bucket}/out/*", f"{local_out}/"]
    delete = ["gcloud", "compute", "instances", "delete", name,
              "--project", project, "--zone", zone, "--quiet"]

    return [
        {"name": "stage_up", "cmd": stage_up,
         "note": "ship ONLY the shortlist + manifest up to the GCS staging bucket"},
        {"name": "create_instance", "cmd": create, "note": create_note},
        {"name": "fold", "cmd": fold,
         "note": "on the VM: pull inputs, docker-run ESMFold (resumable + pLDDT-gated), push out"},
        {"name": "stage_down", "cmd": stage_down,
         "note": "pull folded PDBs + pLDDT down; then run proteus.screen"},
        {"name": "delete_instance", "cmd": delete,
         "note": "stop billing once results are down"},
    ]


def render_plan(plan: list[dict]) -> str:
    lines = ["[launch] GCE burst plan (DRY-RUN — nothing executed):"]
    for i, step in enumerate(plan, 1):
        lines.append(f"  {i}. {step['name']}: {step['note']}")
        lines.append(f"     $ {' '.join(shlex.quote(c) for c in step['cmd'])}")
    return "\n".join(lines)


def execute_plan(plan: list[dict], only: set[str] | None = None) -> int:
    """Execute selected steps (LOCAL, on the Mac). Skips steps that still carry an
    unresolved <PLACEHOLDER> (project/bucket/image unset). Returns a process exit code."""
    for step in plan:
        if only and step["name"] not in only:
            continue
        if any("<" in tok and ">" in tok for tok in step["cmd"]):
            print(f"[launch] skip {step['name']}: unresolved placeholder — set "
                  "compute.gce_burst.{project,bucket,image}", file=sys.stderr)
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
    ap.add_argument("--execute", action="store_true",
                    help="actually run the plan (default: dry-run print only)")
    ap.add_argument("--only", default=None,
                    help="comma-separated step names to execute (e.g. stage_up,create_instance)")
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

    plan = build_plan(cfg, args.manifest, args.shortlist)
    if not args.execute:
        print(render_plan(plan))
        print("[launch] dry-run only — set compute.gce_burst.{project,bucket,image} and "
              "re-run with --execute to burst.")
        return 0

    # Fail fast: refuse to --execute with unresolved config, so we never create a VM
    # only to skip the fold (and leave a half-run). build_plan leaves <…> placeholders
    # when project/bucket/image are unset.
    vb = _burst_cfg(cfg)
    missing = [k for k in ("project", "bucket", "image") if not str(vb.get(k) or "").strip()]
    if missing:
        print(f"[launch] refusing to --execute: set compute.gce_burst.{{{', '.join(missing)}}} "
              "in config first (they show as <…> placeholders in the dry-run above).",
              file=sys.stderr)
        return 2

    # stage_down does `gcloud storage cp … <local_out>/`, which errors ("Destination URL
    # must name an existing directory") if local_out doesn't exist yet. Create it up front.
    os.makedirs(os.path.join(os.getcwd(), vb["local_out"]), exist_ok=True)

    only = set(args.only.split(",")) if args.only else None
    return execute_plan(plan, only=only)


if __name__ == "__main__":
    raise SystemExit(main())
