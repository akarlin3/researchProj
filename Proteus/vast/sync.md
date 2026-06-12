# Vast.ai burst — ship S2 shortlist up, fold (S3), pull structures back

ESMFold (S3) and any Chai-1 refinement / GPU docking are GPU-heavy and run as a
**burst** on a Vast.ai Linux+CUDA box, **not** on the M4 MacBook Air. The local
pipeline (S0–S2) narrows the corpus first, so only the **S2 shortlist** is shipped
up. Nothing here installs into the local conda env.

## What moves, in which direction

```
   LOCAL (M4 Air)                         VAST.AI burst (Linux + CUDA)
   ──────────────                         ────────────────────────────
   S0 dereplicate ─┐
   S1 ProstT5 3Di  ├─ narrow locally
   S2 foldclass   ─┘
        │  (1) UP: shortlist FASTA + S3 job manifest
        ▼ ───────────────────────────────►  S3 ESMFold batch fold
                                             (+ optional Chai-1 refine)
        ◄───────────────────────────────  (2) DOWN: folded PDBs + pLDDT
   S4 geometry  ◄─┐
   S5 cleft     ◄─┘  resume locally on returned structures
```

## 0. Build & push the image (once)

```bash
# On a Linux/CUDA box or CI — NOT on the Mac:
docker build -t <registry>/proteus-fold:cu124 -f vast/Dockerfile.fold .
docker push <registry>/proteus-fold:cu124        # Docker Hub / GHCR Vast can pull
```

## 1. Emit the job manifest locally (dry-run, on the Mac)

```bash
PYTHONPATH=src python -m proteus.s3_fold --dry-run \
  --fasta data/interim/s2_shortlist.fasta \
  --out   data/interim/s3_job_manifest.json
```

The manifest lists each sequence (id, length, sha256) + the fold params resolved
from `config/proteus.yaml` (`plddt_min`, `chunk_size`, `max_recycles`,
`random_seed`). It is the contract the burst runner consumes.

## 2. Launch a Vast instance and ship the inputs UP

```bash
# Pick a cheap interruptible CUDA offer (see "interruptible vs on-demand" below)
vastai search offers 'gpu_name=RTX_4090 num_gpus=1 cuda_vers>=12.4 inet_down>200' \
  --interruptible -o dph

vastai create instance <OFFER_ID> \
  --image <registry>/proteus-fold:cu124 \
  --disk 60 --bid <PRICE>                 # interruptible: bid; on-demand: drop --bid

# Push only the shortlist + manifest (NOT the whole corpus)
rsync -avP data/interim/s2_shortlist.fasta data/interim/s3_job_manifest.json \
  vast:/data/proteus/in/
```

## 3. Fold on the box

```bash
# Inside the instance: fold to a MOUNTED/persistent volume (see checkpointing)
python3 /opt/proteus/run_fold.py \
  --manifest /data/proteus/in/s3_job_manifest.json \
  --fasta    /data/proteus/in/s2_shortlist.fasta \
  --out      /data/proteus/out/            # PDBs + per-model pLDDT here
```

## 4. Pull structures back DOWN to the Mac

```bash
rsync -avP vast:/data/proteus/out/ structures/folded/
# resume locally — screen the returned models through S4 (triad geometry) + S5
# (cleft), scored against the positive-control anchor at the calibrated operating
# point. Reads the runner's s3_results.json (only pLDDT-kept models) automatically:
PYTHONPATH=src python -m proteus.screen \
  --folded structures/folded --struct-dir structures \
  --out data/processed/s4s5_candidates        # ranked PETase-like hits (.csv + .json)
vastai destroy instance <INSTANCE_ID>      # stop billing once results are down
```

`proteus.screen` chains the two gates and the control-anchored cleft score: a
returned model is a **PETase-like hit** only if it (1) passes the S4 catalytic
triad + oxyanion-hole gate, (2) has a catalytic pocket (S5), and (3) scores at or
above the calibration operating point (lowest positive control). Non-PETase serine
hydrolases clear (1)+(2) but fall below (3) — exactly as the controls calibrate.

## Interruptible vs on-demand (and why checkpointing matters)

- **Interruptible** (default; `compute.vast_burst.instance_pref: interruptible` in
  `config/proteus.yaml`) — cheapest, but the instance **can be reclaimed mid-job**.
  The burst runner MUST write each folded model to the mounted/persistent volume as
  soon as it finishes and skip already-done ids on restart, so a reclaim costs only
  the single in-flight sequence, not the whole batch.
- **On-demand** — guaranteed, pricier. Use for a final time-boxed pass, or when a
  single long fold must not be interrupted.

> Persistence rule: never write fold outputs only to the container's ephemeral
> root. Mount a Vast volume (or sync to a bucket) at `/data/proteus/out` so a
> reclaim/restart resumes instead of refolding.
