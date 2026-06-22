# Matrix — the synthetic-twin closed loop (Keystone's no-scanner mode)

Matrix is a **working closed-loop harness on a synthetic digital twin**: it closes the loop

```
scan → posterior → trust gate → action gate → dose replan → re-scan
```

on a simulated abdominal target with known ground-truth IVIM `(D, D*, f)` and a dose/response
model — with **no scanner** and **no real patient data** anywhere in the tree or history.

Matrix is **Keystone's no-scanner run-mode**, built standalone and early to de-risk lab access.
Keystone (the capstone) runs the same loop **real-time** (Forge + scanner), **offline**, or as
**Matrix** (synthetic twin). Building Matrix now means the capstone exists even if real MR-Linac
access, the Forge dose engine, or IRB approval never lands.

> **Scope ceiling (read this).** Matrix is a *working closed-loop harness on synthetic data*, **not
> a validated clinical loop**. Every result means **"the loop closes and behaves sensibly on a
> synthetic twin"** — never a clinical claim. See [`ASSUMPTIONS.md`](ASSUMPTIONS.md).

## Status

| checkpoint | gate | status |
|---|---|---|
| **CP1** — synthetic twin + harness | twin bit-reproducible from seed; harness runs end-to-end (stubs) | **PASS** (`verify_cp1.py`) |
| **CP2** — Fashion ruler + Minos trust/action gates | gates fire correctly; provisional flags in place | **PASS** (`verify_cp2.py`) |
| **CP3** — Forge-shaped dose-replan stage | stage consumes the decision, returns a replan; Forge-drop-in-ready | **PASS** (`verify_cp3.py`) |
| **CP4** — closed-loop run + honest scope | loop closes reproducibly; scope honest; no clinical claim | **PASS** (`verify_cp4.py`) |

Headline (synthetic twin, seed `20260622`, 12×12 voxels, 6 iterations; 95% bootstrap CIs):

- **Ruler** calibrates the perfusion-fraction error bar: 95% empirical coverage **0.94**, ECE_f **0.007**.
- **Trust gate** flags the untrustworthy (low-SNR) zone: AUROC **1.00**; fires **1.00** inside, **0.02** outside.
- **Action gate + trust gate** suppress action on untrustworthy voxels: action rate **0.00 [0.00, 0.00]**
  gated vs **0.42 [0.36, 0.49]** ungated, while trustworthy voxels stay free to act (**0.57 [0.54, 0.61]**).
- **Closed loop converges:** trusted tumour perfusion drops **0.176 [0.167, 0.185]** under treatment
  (CI excludes 0) while untrusted tumour is **held** (0.000); treatment winds down (TREAT actions → 0).

## The three consumed components (each stubbed behind a clean interface)

None is final, so each sits behind a documented interface with a clearly-labelled placeholder; the
real component drops in **without touching `loop.py`** (see [`PROMOTION.md`](PROMOTION.md)).

| role | component | interface | placeholder | real component status |
|---|---|---|---|---|
| calibrated error bars | **Fashion** | `interfaces/ruler.py :: Ruler` | `PlaceholderRuler` (NOT-Fashion) | in review @ *NMR in Biomedicine* |
| trust + action gates | **Minos** | `interfaces/gates.py :: TrustGate`, `ActionGate` | `Placeholder{Trust,Action}Gate` (NOT-Minos) | applied half provisional (PR #49) |
| dose engine | **Forge** | `interfaces/dose.py :: DoseEngine` | `PlaceholderDoseEngine` (NOT-Forge) | **deferred to 2027 — not built** |

## Layout

```
Matrix/
  README.md  ASSUMPTIONS.md  PROMOTION.md  reproduce.sh
  verify_cp1.py … verify_cp4.py     per-checkpoint falsifiable gate scripts
  matrix/
    config.py        MatrixConfig — every knob, one seeded object
    twin.py          synthetic digital twin: ground-truth (D,D*,f) + dose/response
    forward.py       the twin's synthetic IVIM scanner (biexponential + Rician noise)
    fit.py           posterior stage — segmented IVIM fit → raw (mu, sigma)
    state.py         LoopState — the shared object threaded through the four stages
    loop.py          the four-stage harness + Interfaces bundle
    evaluate.py      bootstrap CIs + AUROC for the load-bearing numbers
    interfaces/
      ruler.py       Fashion-shaped   (calibrated error bars)
      gates.py       Minos-shaped     (trust gate + action gate)
      dose.py        Forge-shaped     (dose engine; placeholder only)
  tests/             pytest suite
  results/           seeded run printouts (RESULTS_CP4.md)
```

## Run it

Uses the `proteus` conda env (numpy / scipy):

```bash
PROT=/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python
bash Matrix/reproduce.sh        # ONE command: pytest → CP1 → CP2 → CP3 → CP4

# or individually:
PYTHONPATH=Matrix $PROT -m pytest Matrix/tests -q
$PROT Matrix/verify_cp4.py      # the closed-loop run + honest scope
```

## What "the loop closes" means (the falsifiable check)

On the seeded synthetic twin, `scan → posterior → trust → action → replan → re-scan` must satisfy:
1. **Reproducible** — same seed → bit-identical twin and identical final state.
2. **End-to-end** — every stage writes its output every iteration; nothing raises.
3. **Sensible** (each with a bootstrap CI): the trust gate suppresses action on untrustworthy
   voxels; the action gate changes dose only where warranted; the loop converges (trusted tumour
   perfusion falls under treatment; untrusted tumour is held).

All three are properties of the **synthetic** twin + **placeholder** components — they validate the
*harness*, not any clinical effect.
