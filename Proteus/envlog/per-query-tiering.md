# Per-query tiering — the PETase-branch test

**Run date:** 2026-06-14  
**Where:** local Mac (M4), off the existing GCS sweep artifacts (`result.m8`) — **no GCE, no re-search**.  
**Reproduce:** `PYTHONPATH=src python -m proteus.per_query run` (config `per_query` block; seed 1729; line pinned to -1.1587).  
**Reuses** `screen`/S4/S5 untouched at the pinned line; the random floor is reused AS MEASURED (floor.json) as the third arm.

## TL;DR

PETase top-300 above-line|triad = 32.65% vs floor 42.86% — FLAT overall, but a real tail gradient (top 25 = 80.00%, p 0.01) confined to near-homolog bits. Honest-negative for DISCOVERY.

## Checkpoint 0 — preconditions (audit)

| Precondition | Status |
|---|---|
| `screen`/S4/S5/`calibrate` intact | ✅ reused untouched; imports + run clean in the `proteus` env (numpy/scipy/biotite; fpocket on PATH). A known enriched hit re-screens to its documented triad exactly (MGYP000470279205 → 497/489/411). |
| Pinned line | **-1.1587** (enriched `funnel.json` threshold). Local re-derivation jitters (fpocket); pinned so all arms judge at one line. |
| Floor (baseline arm) | ✅ /Users/averykarlin/projProteus/.claude/worktrees/per-query-tiering/data/processed/floor.json: 1500 screened → 28 triad+ → 12 above-line; conditional 12/28 = 42.86%. Reused, not re-screened. |
| `result.m8` | ✅ 1,081,416 alignment rows; 217,833 unique Atlas targets; 3 anchor classes mapped. |
| Fetch path | ✅ `fetchPredictedStructure/{acc}.pdb` live (structurally identical to the GCE foldcomp models). |

## Checkpoint 1 — partition by query anchor (best-match)

Each of the unique Atlas targets is assigned to the anchor of its single **best** query match (highest Foldseek bits). This is the fix for the global-bits ranking that let one anchor dominate.

**Anchor map** (query id in `result.m8` → class):

- **PETASE**: 6EQE, 4EB0, 8B4U, 4WFI, 4CG1
- **ACHE**: 1EA5
- **OTHER_NEG**: 1TCA_A, 1CRL_A, 1EVQ

**Branch sizes** (of 217,833 unique targets, by best-match anchor):

| Branch | Targets | Share |
|---|---:|---:|
| PETASE | 169,487 | 77.81% |
| ACHE | 6,649 | 3.05% |
| OTHER_NEG | 41,697 | 19.14% |
| **total** | **217,833** | 100% |

PETASE is the **largest** branch — 169,487 targets (77.81%) — while ACHE is only 6,649 (3.05%). The original 96%-AChE / 0%-PETase top tier was therefore a pure **bits-magnitude** artifact, not a coverage one: the PETASE branch's single best hit (1717 bits) scores BELOW the ACHE branch's top-300 cutoff (1757 bits), so global-bits ranking buries the entire PETase branch beneath the AChE tier. Per-query tiering takes the top 300 PETASE targets (bits 962–1717) — finally screening PETase-neighbours the discovery sweep never reached.

## Checkpoint 2 — tier WITHIN each anchor, then screen

Top **300** per branch by that anchor's bits (or all if smaller), fetched and run through the UNCHANGED S4 → fpocket(triad+) → S5 path at the pinned line **-1.1587**, parallelised. Per-branch funnels:

| Stage | PETASE | ACHE | floor (random) |
|---|---:|---:|---:|
| selected / screened | 300 | 300 | 1500 |
| → triad+ (S4) | 294 | 297 | 28 |
| → catalytic pocket (S5) | 288 | 294 | — |
| → **above line** | **96** | **26** | **12** |

## Checkpoint 3 — the conditional test (the thesis question)

**Triad rate** (triad+/screened) and **above-line | triad** (the load-bearing conditional), each with Wilson 95% CI:

| Arm | triad rate | above-line \| triad |
|---|---|---|
| **PETASE** | 294/300 = 98.00% [95.71%, 99.08%] | 96/294 = 32.65% [27.55%, 38.21%] |
| **ACHE** | 297/300 = 99.00% [97.10%, 99.66%] | 26/297 = 8.75% [6.04%, 12.52%] |
| **floor (random)** | 28/1500 = 1.87% | 12/28 = 42.86% [26.51%, 60.93%] |

**PETASE above-line|triad vs floor 42.9% — the load-bearing test**  
- PETASE 96/294 = 32.65% [27.55%, 38.21%] vs floor 12/28 = 42.86% [26.51%, 60.93%]  
- two-proportion z = -1.09, p = 2.74e-01; Fisher exact p = 2.98e-01  
- rate ratio (PETASE/floor) = 0.76× [0.48, 1.20]

**PETASE branch vs ACHE branch (conditional)**  
- PETASE 96/294 = 32.65% [27.55%, 38.21%] vs ACHE 26/297 = 8.75% [6.04%, 12.52%]  
- two-proportion z = 7.18, p = 7.11e-13; Fisher exact p = 3.27e-13  
- rate ratio (PETASE/ACHE) = 3.73× [2.49, 5.58]

**ACHE branch vs floor (conditional, context)**  
- ACHE 26/297 = 8.75% [6.04%, 12.52%] vs floor 12/28 = 42.86% [26.51%, 60.93%]  
- two-proportion z = -5.37, p = 7.94e-08; Fisher exact p = 9.50e-06  
- rate ratio (ACHE/floor) = 0.20× [0.12, 0.36]

### Gradient WITHIN the PETASE branch — above-line | triad by closeness to a PETase query

The load-bearing test averages over the whole top-300 tier. Stratifying the PETASE triad-bearers by Foldseek bits (closeness to the nearest PETase query) shows where any signal actually lives:

| tier (by bits) | min bits | above-line \| triad | vs floor 42.9% (Fisher p, RR) |
|---|---:|---|---|
| top 25 | 1251 | 20/25 = 80.00% [60.87%, 91.14%] | p = 0.011, RR 1.87× |
| top 50 | 1180 | 32/50 = 64.00% [50.14%, 75.86%] | p = 0.096, RR 1.49× |
| top 100 | 1090 | 44/100 = 44.00% [34.67%, 53.77%] | p = 1.000, RR 1.03× |
| rank 101–294 | 962 | 52/194 = 26.80% [21.07%, 33.44%] | p = 0.116, RR 0.63× |

The gradient decays monotonically with bits and crosses the floor (~43%) in the mid-tier: the closest PETase-neighbours clear well above random, the bulk at or below it. Crucially, the above-floor tail sits at **near-homolog bits** — where a plain sequence search already reaches — so it does **not** rescue the divergent-dark-tail discovery the project targets.

## Checkpoint 4 — candidates + verdict

### PETase-anchored, above-line hits — the project's best candidates from the correct branch

| rank | accession | nearest PETase query | bits | composite | pLDDT | Ser/His/acid |
|---:|---|---|---:|---:|---:|---|
| 1 | MGYP000644355661 | 6EQE | 1026 | 2.758 | 0.947 | 154/236/204 |
| 2 | MGYP003569112580 | 4EB0 | 1032 | 2.529 | 0.918 | 153/230/198 |
| 3 | MGYP001248801471 | 4EB0 | 1035 | 2.127 | 0.901 | 149/229/197 |
| 4 | MGYP001305161758 | 6EQE | 1048 | 1.988 | 0.897 | 157/238/206 |
| 5 | MGYP003526118276 | 4CG1 | 1508 | 1.940 | 0.949 | 164/244/212 |
| 6 | MGYP001285521820 | 6EQE | 1021 | 1.921 | 0.903 | 137/217/185 |
| 7 | MGYP001262149556 | 4CG1 | 1642 | 1.568 | 0.961 | 155/233/201 |
| 8 | MGYP000617819104 | 4EB0 | 1046 | 1.492 | 0.921 | 150/232/200 |
| 9 | MGYP001400559074 | 4EB0 | 1057 | 1.481 | 0.907 | 171/248/216 |
| 10 | MGYP001320156145 | 4CG1 | 1567 | 1.396 | 0.921 | 161/239/207 |
| 11 | MGYP001162562203 | 6EQE | 1064 | 1.337 | 0.933 | 150/230/198 |
| 12 | MGYP001295327930 | 6EQE | 989 | 1.329 | 0.939 | 137/216/184 |
| 13 | MGYP003109166174 | 6EQE | 991 | 1.299 | 0.890 | 150/232/200 |
| 14 | MGYP000300759853 | 4CG1 | 978 | 1.233 | 0.869 | 152/228/196 |
| 15 | MGYP001201181085 | 4CG1 | 1655 | 1.199 | 0.922 | 165/243/211 |
| 16 | MGYP000368171655 | 4CG1 | 1693 | 1.181 | 0.963 | 169/247/215 |
| 17 | MGYP000601513684 | 4EB0 | 1056 | 1.165 | 0.957 | 36/113/81 |
| 18 | MGYP001235833815 | 6EQE | 978 | 1.145 | 0.932 | 146/225/193 |
| 19 | MGYP000544200538 | 4CG1 | 1316 | 1.144 | 0.962 | 140/219/187 |
| 20 | MGYP000931547452 | 6EQE | 1347 | 1.071 | 0.942 | 153/231/199 |
| 21 | MGYP003575547131 | 4CG1 | 1437 | 1.052 | 0.970 | 155/234/202 |
| 22 | MGYP001410896016 | 4CG1 | 1703 | 0.978 | 0.950 | 169/247/215 |
| 23 | MGYP002634375525 | 6EQE | 1063 | 0.965 | 0.928 | 158/240/208 |
| 24 | MGYP001489421514 | 4CG1 | 1717 | 0.801 | 0.951 | 174/252/220 |
| 25 | MGYP000414108240 | 4CG1 | 978 | 0.775 | 0.898 | 151/231/199 |
| 26 | MGYP001166430017 | 6EQE | 1221 | 0.764 | 0.926 | 148/227/195 |
| 27 | MGYP001286498480 | 4CG1 | 1428 | 0.704 | 0.955 | 150/228/196 |
| 28 | MGYP000060581275 | 4CG1 | 1326 | 0.698 | 0.947 | 152/231/199 |
| 29 | MGYP001821672909 | 4CG1 | 1087 | 0.660 | 0.886 | 135/211/181 |
| 30 | MGYP003974971481 | 4EB0 | 1014 | 0.617 | 0.915 | 140/222/190 |
| 31 | MGYP003471375709 | 6EQE | 1059 | 0.466 | 0.916 | 138/213/183 |
| 32 | MGYP003790300779 | 4CG1 | 990 | 0.452 | 0.928 | 100/177/146 |
| 33 | MGYP001149455765 | 4CG1 | 1515 | 0.388 | 0.928 | 148/227/195 |
| 34 | MGYP000564582922 | 4EB0 | 1069 | 0.261 | 0.912 | 172/249/217 |
| 35 | MGYP000740062793 | 6EQE | 1511 | 0.050 | 0.896 | 173/252/220 |
| 36 | MGYP000726735156 | 8B4U | 1152 | 0.020 | 0.979 | 113/225/194 |
| 37 | MGYP003700266689 | 6EQE | 1157 | -0.002 | 0.949 | 147/226/194 |
| 38 | MGYP003450003198 | 6EQE | 1140 | -0.011 | 0.880 | 185/260/230 |
| 39 | MGYP000629608983 | 6EQE | 1046 | -0.012 | 0.892 | 177/259/227 |
| 40 | MGYP003407001727 | 4EB0 | 1077 | -0.024 | 0.897 | 171/254/221 |
| 41 | MGYP001247230757 | 4CG1 | 1374 | -0.031 | 0.925 | 183/260/229 |
| 42 | MGYP001368910707 | 4WFI | 966 | -0.144 | 0.924 | 59/137/105 |
| 43 | MGYP001368387099 | 6EQE | 1272 | -0.174 | 0.951 | 109/186/154 |
| 44 | MGYP001477053272 | 4EB0 | 1078 | -0.194 | 0.898 | 136/218/186 |
| 45 | MGYP003476882473 | 6EQE | 1244 | -0.196 | 0.942 | 93/170/139 |
| 46 | MGYP002714884928 | 6EQE | 1032 | -0.203 | 0.964 | 91/169/137 |
| 47 | MGYP000210912620 | 6EQE | 1209 | -0.214 | 0.937 | 115/193/161 |
| 48 | MGYP002719446038 | 4CG1 | 1044 | -0.218 | 0.919 | 170/246/214 |
| 49 | MGYP000173184570 | 4CG1 | 1306 | -0.242 | 0.930 | 174/252/220 |
| 50 | MGYP000389168146 | 6EQE | 1185 | -0.368 | 0.948 | 104/182/150 |
| 51 | MGYP001151488865 | 6EQE | 1036 | -0.369 | 0.932 | 117/202/170 |
| 52 | MGYP002817680355 | 6EQE | 1095 | -0.378 | 0.915 | 194/271/239 |
| 53 | MGYP001245784375 | 8B4U | 1180 | -0.389 | 0.971 | 112/231/198 |
| 54 | MGYP001318806059 | 4CG1 | 1066 | -0.393 | 0.960 | 78/156/124 |
| 55 | MGYP001176807024 | 4CG1 | 1140 | -0.427 | 0.952 | 81/159/127 |
| 56 | MGYP001323382937 | 6EQE | 1080 | -0.432 | 0.889 | 160/237/205 |
| 57 | MGYP001466411950 | 8B4U | 1184 | -0.452 | 0.980 | 112/234/206 |
| 58 | MGYP003575394538 | 4CG1 | 1052 | -0.471 | 0.958 | 91/170/138 |
| 59 | MGYP000388388267 | 4CG1 | 1101 | -0.475 | 0.940 | 99/177/145 |
| 60 | MGYP000153285637 | 8B4U | 1044 | -0.482 | 0.906 | 128/254/220 |
| 61 | MGYP000975344424 | 8B4U | 967 | -0.484 | 0.963 | 137/259/231 |
| 62 | MGYP000647201728 | 8B4U | 983 | -0.498 | 0.972 | 116/229/200 |
| 63 | MGYP002716753343 | 4CG1 | 1100 | -0.534 | 0.962 | 165/240/209 |
| 64 | MGYP002784815584 | 6EQE | 1258 | -0.574 | 0.951 | 110/189/157 |
| 65 | MGYP000456884793 | 8B4U | 1205 | -0.699 | 0.967 | 114/236/204 |
| 66 | MGYP002395463326 | 8B4U | 1006 | -0.736 | 0.902 | 119/229/200 |
| 67 | MGYP000855197346 | 8B4U | 1206 | -0.752 | 0.968 | 113/232/203 |
| 68 | MGYP001815130755 | 8B4U | 963 | -0.757 | 0.979 | 107/228/199 |
| 69 | MGYP001583139274 | 8B4U | 970 | -0.763 | 0.931 | 114/228/200 |
| 70 | MGYP001251303412 | 8B4U | 1277 | -0.770 | 0.958 | 114/234/202 |
| 71 | MGYP000686585230 | 4CG1 | 1135 | -0.779 | 0.933 | 96/174/142 |
| 72 | MGYP000403492051 | 4CG1 | 1557 | -0.828 | 0.909 | 122/200/168 |
| 73 | MGYP002625649965 | 8B4U | 987 | -0.830 | 0.931 | 132/253/224 |
| 74 | MGYP002514292131 | 8B4U | 963 | -0.845 | 0.855 | 173/297/267 |
| 75 | MGYP001413951859 | 8B4U | 1015 | -0.892 | 0.960 | 119/241/208 |
| 76 | MGYP000280442270 | 8B4U | 1236 | -0.893 | 0.954 | 114/233/202 |
| 77 | MGYP001376702680 | 4EB0 | 1032 | -0.894 | 0.926 | 177/254/222 |
| 78 | MGYP001490773307 | 8B4U | 1053 | -0.895 | 0.967 | 116/232/203 |
| 79 | MGYP000415728319 | 8B4U | 1237 | -0.938 | 0.967 | 120/240/208 |
| 80 | MGYP001401963747 | 8B4U | 973 | -0.944 | 0.965 | 105/228/199 |
| 81 | MGYP001351152434 | 8B4U | 1035 | -0.956 | 0.965 | 113/229/201 |
| 82 | MGYP000858230106 | 8B4U | 1093 | -0.976 | 0.977 | 117/239/206 |
| 83 | MGYP001186384904 | 8B4U | 1035 | -0.984 | 0.964 | 112/235/203 |
| 84 | MGYP001229627836 | 8B4U | 1156 | -0.984 | 0.967 | 116/238/206 |
| 85 | MGYP000870037404 | 8B4U | 1066 | -1.007 | 0.981 | 115/237/204 |
| 86 | MGYP002521931278 | 8B4U | 1004 | -1.016 | 0.950 | 119/244/214 |
| 87 | MGYP003608607993 | 8B4U | 1203 | -1.016 | 0.979 | 114/231/203 |
| 88 | MGYP001205384655 | 8B4U | 1090 | -1.053 | 0.967 | 110/231/202 |
| 89 | MGYP001809843667 | 8B4U | 1017 | -1.057 | 0.944 | 125/245/213 |
| 90 | MGYP000471455435 | 8B4U | 976 | -1.060 | 0.974 | 111/226/198 |
| 91 | MGYP000989336048 | 8B4U | 1048 | -1.062 | 0.957 | 107/225/196 |
| 92 | MGYP003609754630 | 8B4U | 1234 | -1.067 | 0.977 | 113/234/202 |
| 93 | MGYP000162028012 | 8B4U | 1011 | -1.092 | 0.950 | 118/227/200 |
| 94 | MGYP000087980178 | 4CG1 | 1122 | -1.125 | 0.860 | 160/238/206 |
| 95 | MGYP000938081187 | 8B4U | 1156 | -1.139 | 0.974 | 113/232/203 |
| 96 | MGYP001817575622 | 8B4U | 1085 | -1.155 | 0.952 | 117/238/206 |

### Verdict — is there ANY PET-specific gradient in the structural signal?

**NO — the structural signal is flat.** PETase-neighbours clear the cleft line at **32.65%** ([27.55%, 38.21%]), **indistinguishable from** the random floor of 42.86% (rate ratio 0.76×, Fisher p = 2.98e-01). Being a PETase-neighbour confers no extra above-line signal the pipeline captures. **At the discovery tier (top-300) the thesis branch was screened and shown indistinguishable from random.** The pipeline is **not blind**: PETase-neighbours clear the line 3.7× more often than AChE-neighbours (8.75%, Fisher p = 3.3e-13) — but that separation is driven by AChE being a deep-gorge outlier sitting FAR BELOW the floor, not by PETase rising above it. **But a bits-stratified analysis finds a real gradient confined to the closest PETase-neighbours:** the top 25 by bits (≥1251) clear at **80.00%** — significantly above the floor (Fisher p = 0.011, 1.87×) — decaying monotonically through the floor in the bulk. The gradient lives only at **near-homolog distances** (high bits), exactly where a sequence search already reaches; in the divergent dark tail the project targets, the structural signal is flat. So the gradient is real but does NOT generalise to the discovery tier. **Fork:** the honest-negative / methods paper stands for DISCOVERY (the divergent dark tail is flat); any tail gradient is too thin and too near-homolog to justify the aromatic-subsite specificity build on its own.

### Scope guard (carried)

**Exposure ≠ PET activity.** Even a positive gradient does not verify any candidate — it only says PETase-neighbours are more exposed-site than random hydrolases. No wet-lab; the S4/S5 path tests fold-class + an exposed-cleft geometry, not PET turnover. Leads are **prioritized, not verified**.

## Caveats

- **Floor triad+ n = 28 (< 50)** — the floor conditional CI is wide [26.5%, 60.9%]. The FLAT discovery-tier verdict is robust regardless: PETASE top-300 sits AT/BELOW the floor point estimate, so a tighter floor cannot manufacture a positive gradient (PETASE would have to rise ABOVE it). The closest-neighbour tail result (top-25 vs floor) already clears significance (p = 0.011) despite the small floor.
- **fpocket is non-deterministic** run-to-run (composites swing ~±1.0; the re-derived line wanders ~±0.01). Handled by PINNING the decision line to -1.1587 and caching the screened records: the headline conditional (96/294) reproduced EXACTLY across two independent screen passes.
- **The above-floor signal is near-homolog.** The tail that beats the floor sits at high Foldseek bits — sequence search already reaches there. The divergent dark tail (lower bits) that motivates a *structure-first* discovery is exactly where the signal is flat.

## Reproducibility

- **Seed** 1729; **branch_n** 300; **line** -1.1587 (pinned to the enriched sweep).
- **Partition** best-match over `result.m8` (1,081,416 rows → 217,833 unique targets).
- **Fetch** `api.esmatlas.com/fetchPredictedStructure/{acc}.pdb` (live; structurally identical to the GCE foldcomp models).
- **Anchor** IsPETase+LCC, percentile mode (same as calibration/floor).
- **Artifacts** `branch_partition.csv`, `per_query_tiering.json`; floor comparator from `gs://projproteus-fold/atlas-sweep/2026-06-14/`.
