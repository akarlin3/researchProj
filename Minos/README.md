# projSybil — Sibyl

Out-of-distribution detection and trustworthiness flagging for quantitative
diffusion-MRI microstructure estimation (IVIM via µGUIDE).

The core idea: a per-input OOD flag should track *when the estimator's output stops
being trustworthy*. Two detector families are studied — **Family 1** (summary-space:
Mahalanobis/MMD on the µGUIDE embedding) and **Family 2** (signal-space: prediction
residual / density-ratio / conformal).

## Tiers

- **Tier 1** (`sibyl/experiments/tier1.py`) — fully synthetic. Inject controlled
  shifts (noise model, SNR, corruption) on a dense 10-b scheme and show the OOD score
  couples to calibration failure of the µGUIDE posterior.
- **Tier 2** (`sibyl/experiments/tier2.py`) — in-vivo breast DWI at the ACRIN-6698
  4-b scheme. Show the detector flags the real synthetic→in-vivo shift and that the
  flag tracks loss of **ADC test-retest** trustworthiness. See
  [`docs/tier2.md`](docs/tier2.md) and the pre-coding [`VERIFICATION.md`](VERIFICATION.md).
  *(D\*/f are not validated in vivo here — that waits for the liver data.)*
- **Tier 3** (`sibyl/experiments/tier3.py`) — glioma / liver cohort (stub).

## Layout

```
sibyl/
  forward_model/ivim.py     IVIM biexponential, ACRIN scheme, mono-exp ADC fit
  data/
    synthetic.py            breast IVIM priors + dense ID dataset
    shift.py                Rician/Gaussian noise, corruption, Tier-1 shift axes
    acrin_reference.py      synthetic 4-b ID reference + b0 normalization (Arm 1)
    imputation.py           segmented IVIM fit → dense-grid imputation (Arm 2)
    units.py                UnitTable intermediate + synthetic test/retest generators
    acrin_ingest.py         real ACRIN DWI → UnitTable (NIfTI runnable; DICOM documented)
  detectors/
    family1.py              Mahalanobis on embedding
    family2.py              residual conformal (dense)
    signal_space.py         estimator-free kNN density-ratio / conformal (Arm 1)
  estimator/wrapper.py      µGUIDE train/inference + torch-2.x compat shims
  metrics/eval.py           detection AUROC, ADC repeatability, coupling, controls
  experiments/              tier1 / tier2 / tier3
tests/                      pytest suite (pure-numeric, estimator-free, µGUIDE integ.)
```

## Setup

```bash
# µGUIDE is vendored under uGUIDE/ and imported via a user-site .pth (or `pip install ./uGUIDE`).
pip install pyro-ppl SimpleITK            # runtime deps not in the base env
export PYTHONPATH=$(pwd)
python3 -m pytest tests/ -q
python3 -m sibyl.experiments.tier2        # synthetic validation harness
```
