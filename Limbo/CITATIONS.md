# CITATIONS.md — Limbo's verified-claim ledger

> **No phantom citations, no drift (the Ouroboros / Augur lesson).** Every entry in
> `limbo.bib` appears below with (a) a resolvable identifier — DOI, arXiv eprint, or stable
> proceedings/JMLR URL for venues that mint no DOI — and (b) a one-line *verified claim* that the
> survey may cite, traceable to the source. Each identifier was **resolved against its primary
> source during CP1 (checked 2026-06-22)**: title + authors + year matched the resolved record.
> `verify_citations.py` enforces the bib ↔ ledger correspondence and fails the build on any orphan
> on either side. The script re-runs (with `--online` DOI/arXiv HEAD checks) at CP3.
>
> **Quote status.** Most rows carry a verbatim abstract sentence (held in the per-bucket subagent
> verification logs). Four rows are marked **[thesis-level]**: the claim is the paper's central,
> title/abstract-evident thesis, but a load-bearing *number* from these is not quoted here and must
> be re-pulled from the PDF before any verbatim use — `koo2016`, `blandaltman1986` (publisher pages
> blocked direct abstract fetch; DOI resolves and bibliographic metadata confirmed via Crossref),
> and `ling2000`, `vanhoudt2021qib` (abstract sentence surfaced only as a fragment).

## A — Parameter estimation for quantitative/diffusion body MRI and its uncertainty

| key | identifier | verified claim |
|---|---|---|
| `zhang2013` | doi:10.1109/EMBC.2013.6609549 | Cramér–Rao lower bounds quantify the SNR-dependent estimation uncertainty of the IVIM diffusion coefficient and perfusion fraction (e.g. ~3.9% on D and ~11.7% on f at SNR 40). |
| `jalnefjord2019` | doi:10.1002/mrm.27826 | Analytical CRLBs for D and f from segmented IVIM fitting can be minimised to derive SNR-efficient b-value sampling schemes. |
| `barbieri2016` | doi:10.1002/mrm.25765 | IVIM parameter estimates differ significantly across fitting algorithms in upper-abdominal organs; a Bayesian fit gives the lowest inter-reader and inter-subject variability. |
| `while2017` | doi:10.1002/mrm.26598 | In simulation, Bayesian IVIM fitting consistently outperforms least-squares with lower relative error and deviation and cleaner parameter maps. |
| `gurneychampion2018` | doi:10.1371/journal.pone.0194590 | Among six IVIM fit algorithms applied to pancreatic-cancer DWI, a Bayesian fit gives the best trade-off between tumour contrast and precision for D and f. |
| `bertleff2017` | doi:10.1002/nbm.3833 | A neural network estimating a combined IVIM-kurtosis model yields parameter estimates in better agreement with the literature than least-squares, which biases f, D* and K. |
| `barbieri2020` | doi:10.1002/mrm.27910 | A deep network (including an unsupervised physics-informed IVIM-NET) fits the IVIM model with lower error than least-squares and Bayesian approaches. |
| `kaandorp2021` | doi:10.1002/mrm.28852 | An improved unsupervised physics-informed network gives substantially better test–retest repeatability for D and f than state-of-the-art IVIM fitting methods. |
| `vasylechko2021` | doi:10.1002/mrm.28989 | A self-supervised, physics-based forward-model network estimates IVIM parameters with substantially improved robustness to noise relative to conventional fitting. |
| `casali2025` | arXiv:2508.04588 | A probabilistic deep-learning framework for IVIM MRI estimates total predictive uncertainty and decomposes it into aleatoric and epistemic components. |
| `schmid2006` | doi:10.1109/TMI.2006.884210 | A Bayesian spatial (GMRF) approach to DCE-MRI tracer-kinetic estimation reduces parameter variability in tumour regions while preserving tissue boundaries. |

## B — Trust machinery: calibration, conformal prediction, reliability standards

| key | identifier | verified claim |
|---|---|---|
| `angelopoulos2021` | arXiv:2107.07511 | Conformal prediction is a model-agnostic paradigm that produces statistically rigorous, distribution-free uncertainty sets/intervals for any black-box predictor. |
| `shafer2008` | url:jmlr.org/papers/v9/shafer08a | Conformal prediction outputs, at a user-chosen error rate ε, a label set that contains the true label with probability 1−ε. |
| `kuleshov2018` | url:proceedings.mlr.press/v80/kuleshov18a | A simple post-hoc recalibration procedure can make any regression algorithm produce calibrated uncertainty estimates given enough data. |
| `guo2017` | url:proceedings.mlr.press/v70/guo17a | Modern deep neural networks are poorly calibrated, and temperature scaling is a simple, effective post-hoc remedy. |
| `lakshminarayanan2017` | arXiv:1612.01474 | Deep ensembles are a simple, scalable, non-Bayesian method that yields high-quality predictive uncertainty estimates. |
| `gal2016` | arXiv:1506.02142 | Dropout training can be cast as approximate Bayesian inference, enabling model-uncertainty estimates from standard dropout networks. |
| `amini2020` | arXiv:1910.02600 | Deep evidential regression learns aleatoric and epistemic uncertainty in a single forward pass without sampling or out-of-distribution training data. |
| `ovadia2019` | arXiv:1906.02530 | Predictive uncertainty quality degrades under dataset shift, so trustworthy uncertainty must be evaluated against distribution shift, not just in-distribution. |
| `raunig2015` | doi:10.1177/0962280214537344 | QIBA technical-performance assessment of imaging biomarkers rests on three metrology areas: linearity/bias, repeatability, and reproducibility. |
| `obuchowski2015` | doi:10.1177/0962280214537392 | Bias, precision, and agreement are the appropriate statistical targets for assessing and comparing quantitative imaging biomarker algorithms. |
| `sullivan2015` | doi:10.1148/radiol.2015142202 | Consistent metrology terminology and standards are needed to define and assess the technical performance of quantitative imaging biomarkers. |
| `koo2016` | doi:10.1016/j.jcm.2016.02.012 | **[thesis-level]** Provides the standard guideline for selecting, reporting, and interpreting intraclass correlation coefficients (ICCs) in reliability research. |
| `blandaltman1986` | doi:10.1016/S0140-6736(86)90837-8 | **[thesis-level]** The limits-of-agreement method, not correlation, is the appropriate way to assess agreement between two clinical measurement methods. |

## C — Diffusion/perfusion body-MRI biomarkers and their measured reliability

| key | identifier | verified claim |
|---|---|---|
| `lebihan1986` | doi:10.1148/radiology.161.2.3763909 | The original IVIM paper introduced imaging of intravoxel incoherent motion to separate diffusion and perfusion effects. |
| `lebihan1988` | doi:10.1148/radiology.168.2.3393671 | IVIM was developed to separately measure molecular diffusion and capillary-network microcirculation within a voxel, rather than only the combined ADC. |
| `lebihan2019` | doi:10.1016/j.neuroimage.2017.12.062 | The IVIM concept estimates tissue perfusion by treating randomly oriented capillary blood flow as a pseudo-diffusion process. |
| `tofts1991` | doi:10.1002/mrm.1910170208 | Established the foundational compartmental model for measuring the DCE transfer constant (permeability–surface-area product per unit tissue volume) from dynamic MR imaging. |
| `tofts1999` | doi:10.1002/(SICI)1522-2586(199909)10:3<223::AID-JMRI2>3.0.CO;2-S | The consensus paper standardising DCE-MRI quantities, defining Ktrans, ve and kep with kep = Ktrans/ve. |
| `koh2007` | doi:10.2214/AJR.06.1403 | Body DWI provides qualitative and quantitative (ADC) information reflecting cellular-level change for tumour evaluation, but clinical implementation faces technical obstacles. |
| `padhani2009` | doi:10.1593/neo.81328 | The NCI consensus on DW-MRI as a cancer biomarker recommends that patient-level reproducibility studies be built into study designs, acknowledging measurement variability. |
| `miquel2016` | doi:10.4329/wjr.v8.i1.21 | Across the body-DWI literature, in-vivo ADC reproducibility is substantially worse than phantom reproducibility, with wide cross-study variation. |
| `sun2017` | doi:10.1097/MD.0000000000006866 | In rectal cancer at 3T, IVIM perfusion parameters f and D* have very poor short-term test–retest repeatability (repeatability coefficients ~126% and ~197%) versus far better ADC and D. |
| `sun2019` | doi:10.1016/j.acra.2018.08.012 | In rectal cancer, IVIM pseudo-diffusion D* correlates only weakly with DCE Ktrans (r=0.389), while the f·D* product correlates moderately (r=0.533). |
| `yang2019` | doi:10.1177/0284185118791201 | In rectal cancer, IVIM D* shows only moderate reproducibility (ICC=0.55, CV≈20%) and no significant correlation with DCE Ktrans. |
| `simchick2024` | doi:10.1002/mrm.30237 | Multi-scanner liver IVIM reproducibility varies markedly by parameter, with within-subject CVs from ~5% (D) to ~14% (perfusion-fraction components). |

## D — Value of information: decision-theoretic use of trusted uncertainty

| key | identifier | verified claim |
|---|---|---|
| `vickers2006` | doi:10.1177/0272989X06295361 | Decision curve analysis evaluates prediction models by plotting clinical net benefit against the decision threshold probability. |
| `vickers2016` | doi:10.1136/bmj.i6 | Net benefit places harms and benefits on a common scale via an exchange rate to judge whether a model, marker, or test does more good than harm. |
| `vickers2019` | doi:10.1186/s41512-019-0064-7 | In decision curve analysis the chosen threshold probability encodes the clinician's relative weighting of a missed disease against an unnecessary intervention. |
| `steyerberg2010` | doi:10.1097/EDE.0b013e3181c30fb2 | Decision-analytic measures such as net benefit should be reported when a prediction model is intended to support clinical decisions, beyond discrimination and calibration. |
| `baker2009` | doi:10.1111/j.1467-985X.2009.00592.x | Relative utility curves quantify how much of the achievable clinical utility a risk-prediction model captures relative to perfect prediction. |
| `zhao2021` | arXiv:2107.05719 | Decision calibration requires a model's predicted distribution to be indistinguishable from the true distribution to a set of downstream decision-makers. |
| `lacostejulien2011` | url:proceedings.mlr.press/v15/lacoste_julien11a | Standard approximate inference ignores the decision task and its losses, motivating loss-calibrated approximations tailored to the decision at hand. |
| `vadera2021` | arXiv:2106.06997 | Post-hoc loss-calibration learns posterior approximations that favour high-utility decisions under the relevant loss. |
| `felli1998` | doi:10.1177/0272989X9801800117 | Expected value of perfect information is a superior sensitivity measure because it combines the probability a decision changes with the marginal benefit of that change. |
| `claxton1999` | doi:10.1016/S0167-6296(98)00039-3 | Clinical/economic decisions should be made on expected (mean) net benefit; conventional rules of statistical inference are irrelevant to the choice. |

## E — MR-guided adaptive radiotherapy and biomarker-/uncertainty-driven dose

| key | identifier | verified claim |
|---|---|---|
| `lagendijk2008` | doi:10.1016/j.radonc.2007.10.034 | Integrating MRI with a radiotherapy accelerator makes MRI's soft-tissue and functional imaging available for high-precision, real-time image-guided radiotherapy. |
| `raaymakers2009` | doi:10.1088/0031-9155/54/12/N01 | A prototype integrating a 1.5 T MRI scanner with a 6 MV accelerator demonstrated simultaneous, unhampered MR imaging during irradiation. |
| `raaymakers2017` | doi:10.1088/1361-6560/aa9517 | First-patient treatments on a 1.5 T MRI-Linac established that high-precision, high-field MRI-guided radiotherapy is clinically feasible. |
| `lagendijk2014` | doi:10.1016/j.semradonc.2014.02.009 | The MRI-linac integrates a fully functional 1.5 T diagnostic-quality MRI with a 6 MV linear accelerator in one system. |
| `ling2000` | doi:10.1016/S0360-3016(00)00467-3 | **[thesis-level]** Multidimensional radiotherapy proposes deriving a biological target volume from functional/biological images to enable dose painting and biological conformality. |
| `bentzen2005` | doi:10.1016/S1470-2045(05)01737-7 | Theragnostic imaging uses molecular/functional imaging to prescribe radiation dose in four dimensions — the basis of dose-painting by numbers. |
| `bentzen2011` | doi:10.1016/j.semradonc.2010.10.001 | Dose painting prescribes a non-uniform dose to the target based on functional or molecular images indicating local risk of relapse. |
| `unkelbach2018` | doi:10.1088/1361-6560/aae659 | Robust and probabilistic optimisation methods directly incorporate motion and uncertainty into RT plan optimisation, rather than relying on PTV margins. |
| `otazo2021` | doi:10.1148/radiol.2020202747 | Several quantitative MRI-derived biomarkers reflect radiotherapy response earlier than anatomic imaging, supporting MRI-guided adaptive radiation oncology. |
| `vanhoudt2021` | doi:10.3389/fonc.2020.615643 | Frequent quantitative MRI monitoring of volume, shape, and biological characteristics lets the treatment plan be updated to accommodate observed tumour response. |
| `vanhoudt2021qib` | doi:10.1016/j.ejca.2021.04.041 | **[thesis-level]** Quantitative imaging biomarkers on MR-linacs require technical validation (repeatability/reproducibility) before they can drive clinical decisions, because the integrated scanner differs from conventional diagnostic MRI. |
| `leibfarth2018` | doi:10.1016/j.ctro.2018.09.002 | Diffusion-weighted MRI shows promise as a radiotherapy biomarker for outcome prediction, response assessment, and tumour delineation, but faces geometric-accuracy and quantification challenges. |
| `hall2022` | doi:10.3322/caac.21707 | MR guidance presents a strong opportunity to bring adaptive radiation therapy into routine radiation-oncology clinical practice. |

## In-portfolio cross-references (NOT external citations; not part of the survey body)

These are the author's own program, cited — if at all — only as a *minority* of survey entries on
equal footing with external work, never as the organising centre (the distinctness-from-Augur rule;
see `ASSUMPTIONS.md`). They are tracked here for provenance, not asserted as published literature.

| anchor | status | location |
|---|---|---|
| Fashion | PROVISIONAL (in review, NMR in Biomedicine) | `Fashion/paper_retool/` |
| Minos / Plumbline | PROVISIONAL | `Minos/theory/plumbline.md` |
| Gauge | PROVISIONAL | `Gauge/` |
| Lethe | PROVISIONAL | `Lethe/` |
