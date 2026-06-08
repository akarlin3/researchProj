# Vast.ai burst — Chai-1 + GNINA

GPU-heavy stages (Chai-1 cofolding, GNINA / DiffDock docking) run as **burst**
jobs on Vast.ai rather than on the local RTX 4080. This dir holds the burst
image and launch notes. Nothing here installs into the local conda env.

## Image

`Dockerfile.chai1` — pinned CUDA base (`nvidia/cuda:12.4.1-...`). Chai-1 and
GNINA versions are left as `TODO(P#)` pins: set them on the first real build and
record the resolved versions next to the local lockfiles in `envlog/`.

```bash
docker build -t proteus-burst:chai1 -f vast/Dockerfile.chai1 .
# push to a registry Vast can pull from (Docker Hub / GHCR)
docker push <registry>/proteus-burst:chai1
```

## Interruptible vs on-demand

- **Interruptible** (default; `compute.vast_burst.instance_pref: interruptible`
  in `config/proteus.yaml`) — cheapest, can be reclaimed mid-job. Use for large
  embarrassingly-parallel batches where partial results + resume are fine.
  Checkpoint per-candidate outputs so a reclaim costs only the in-flight item.
- **On-demand** — guaranteed, pricier. Use for a final, time-boxed pass or when
  a single long job must not be interrupted.

```bash
# pick the cheapest interruptible P4-class offer
vastai search offers 'gpu_name=RTX_4090 num_gpus=1 inet_down>200' --interruptible -o dph

# launch (bid for interruptible)
vastai create instance <OFFER_ID> \
  --image <registry>/proteus-burst:chai1 \
  --disk 60 --bid <PRICE>

# on-demand instead: drop --bid and add --price for an on-demand offer
```

## Data flow

1. Local pipeline (S0–S5) produces ranked candidate PDBs + ligands.
2. Stage that batch to the burst instance (rsync/scp or a bucket).
3. Run cofold + dock on GPU; write scores back.
4. Pull results into `data/processed/` locally for final ranking.
