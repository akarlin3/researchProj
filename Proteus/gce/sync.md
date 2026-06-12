# GCE burst — ship S2 shortlist up, fold (S3), pull structures back

ESMFold (S3) runs as a **burst** on a Google Compute Engine VM, **not** on the M4
MacBook Air. This project has **no GPU quota**, so we fold on **CPU** (a C4 VM) —
ESMFold on CPU is slow + RAM-heavy but the narrowed shortlist is small. The local
pipeline (S0–S2) narrows the corpus first, so only the **S2 shortlist** is shipped
up. (To move to GPU once quota is granted: set `compute.gce_burst.accelerator` and
rebuild the image on a CUDA base — `proteus.launch` then builds the GPU plan.)

> **One-command option:** `proteus.launch` automates the steps below — it validates
> the manifest/shortlist, then plans/runs **stage_up → create → fold → stage_down →
> delete** from `compute.gce_burst`. It is **dry-run by default** (just prints the
> plan) and only touches GCP with `--execute`:
> ```bash
> PYTHONPATH=src python -m proteus.launch \
>   --manifest data/interim/s3_job_manifest.json \
>   --shortlist data/interim/s2_shortlist.fasta        # prints the burst plan
> ```
> Set `compute.gce_burst.{project,bucket,image}` first. The manual steps below are
> the same cycle, spelled out.

## What moves, in which direction

```
   LOCAL (M4 Air)                         GCE burst (C4 CPU VM, SPOT)
   ──────────────                         ──────────────────────────────
   S0 dereplicate ─┐
   S1 ProstT5 3Di  ├─ narrow locally
   S2 foldclass   ─┘
        │  (1) UP: shortlist FASTA + S3 job manifest  ──► gs://BUCKET/in/
        ▼ ───────────────────────────────►  S3 ESMFold batch fold on the VM
                                             (writes gs://BUCKET/out/, survives preempt)
        ◄───────────────────────────────  (2) DOWN: folded PDBs + pLDDT
   S4 geometry  ◄─┐
   S5 cleft     ◄─┘  resume locally on returned structures
```

A GCS staging bucket decouples the data from the VM, so a preempted **SPOT** VM
never loses finished models.

## 0. Build & push the image to Artifact Registry (once)

```bash
gcloud artifacts repositories create proteus --repository-format=docker --location=us-central1
# Build on Cloud Build (x86, so the image runs on the x86 GCE VM even though you submit
# from an arm64 Mac — no local Docker/buildx). The image tag is set in cloudbuild.yaml:
gcloud builds submit --config cloudbuild.yaml .
# -> us-central1-docker.pkg.dev/projproteus/proteus/proteus-fold:cpu (put this in
#    compute.gce_burst.image). The CPU torch + transformers build is pinned in
#    gce/Dockerfile.fold; the ~2.5 GB esmfold_v1 weights pull from the Hub on first fold.
```

> Do NOT `docker build` on the arm64 Mac (the VM is x86). Cloud Build above avoids that;
> or cross-build: `docker buildx build --platform linux/amd64 -t $IMAGE -f gce/Dockerfile.fold --push .`

## 1. Emit the job manifest locally (dry-run, on the Mac)

```bash
PYTHONPATH=src python -m proteus.s3_fold --dry-run \
  --fasta data/interim/s2_shortlist.fasta \
  --out   data/interim/s3_job_manifest.json
```

The manifest lists each sequence (id, length, sha256) + the fold params resolved
from `config/proteus.yaml` (`plddt_min`, `chunk_size`, `max_recycles`,
`random_seed`). It is the contract the burst runner consumes (`run_location: gce`).

## 2. Stage the inputs UP to the bucket

```bash
BUCKET=gs://<BUCKET>
gcloud storage cp data/interim/s2_shortlist.fasta data/interim/s3_job_manifest.json $BUCKET/in/
```

## 3. Create the SPOT CPU VM and fold on it

```bash
gcloud compute instances create proteus-fold \
  --project projproteus --zone us-central1-a \
  --machine-type c4-highmem-8 \
  --image-family cos-stable --image-project cos-cloud \
  --boot-disk-size 100GB --scopes cloud-platform \
  --provisioning-model SPOT --instance-termination-action DELETE

# On the VM (Container-Optimized OS has docker but NO gcloud/gsutil), so the GCS copies
# run in the google/cloud-sdk container; the fold runs in our image. The VM service
# account (--scopes cloud-platform) authenticates both automatically.
gcloud compute ssh proteus-fold --project projproteus --zone us-central1-a --command "
  mkdir -p /tmp/proteus/in /tmp/proteus/out &&
  docker run --rm -v /tmp/proteus:/data/proteus google/cloud-sdk:slim \
    gcloud storage cp $BUCKET/in/* /data/proteus/in/ &&
  docker login -u oauth2accesstoken -p \"\$(docker run --rm google/cloud-sdk:slim gcloud auth print-access-token)\" https://us-central1-docker.pkg.dev &&
  docker run -v /tmp/proteus:/data/proteus $IMAGE \
    --manifest /data/proteus/in/s3_job_manifest.json \
    --fasta    /data/proteus/in/s2_shortlist.fasta \
    --out      /data/proteus/out/ --device cpu &&
  docker run --rm -v /tmp/proteus:/data/proteus google/cloud-sdk:slim \
    gcloud storage cp --recursive /data/proteus/out/* $BUCKET/out/"
```

The runner checkpoints per sequence to `/data/proteus/out` and pushes to the bucket,
so a preempt costs only the single in-flight sequence (skipped-already-done on restart).

## 4. Pull structures back DOWN and resume locally

```bash
gcloud storage cp --recursive $BUCKET/out/* structures/folded/
gcloud compute instances delete proteus-fold --project projproteus --zone us-central1-a --quiet

# resume locally — screen the returned models through S4 (triad geometry) + S5
# (cleft), scored against the positive-control anchor at the calibrated operating
# point. Reads the runner's s3_results.json (only pLDDT-kept models) automatically:
PYTHONPATH=src python -m proteus.screen \
  --folded structures/folded --struct-dir structures \
  --out data/processed/s4s5_candidates        # ranked PETase-like hits (.csv + .json)
```

`proteus.screen` chains the two gates and the control-anchored cleft score: a
returned model is a **PETase-like hit** only if it (1) passes the S4 catalytic
triad + oxyanion-hole gate, (2) has a catalytic pocket (S5), and (3) scores at or
above the calibration operating point (lowest positive control). Non-PETase serine
hydrolases clear (1)+(2) but fall below (3) — exactly as the controls calibrate.

## SPOT vs on-demand (and why bucket staging matters)

- **SPOT** (default; `compute.gce_burst.spot: true`) — cheapest, but the VM **can be
  preempted mid-job**. Because the runner writes each model to `/data/proteus/out`
  and pushes to the bucket as soon as it finishes (and skips already-done ids on
  restart), a preempt costs only the single in-flight sequence, not the whole batch.
- **On-demand** (`spot: false`) — guaranteed, pricier. Use for a final time-boxed
  pass, or when a single long fold must not be interrupted.

> Persistence rule: never keep fold outputs only on the VM's boot disk. Stage them
> to the GCS bucket (`gs://BUCKET/out`) so a preempt/restart resumes instead of
> refolding. The launcher's `stage_down` then pulls from the bucket, not the VM.
