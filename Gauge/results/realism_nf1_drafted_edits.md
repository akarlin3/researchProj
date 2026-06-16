# Drafted manuscript edits — realism cohort friction pass (NF1–NF4)

**STATUS: DRAFTED FOR SIGN-OFF — NOT committed to `gauge_v3_revised.tex`.**
All numbers below trace to `results/realism_report.txt` / `results/realism_multiseed.json`
(16-seed bands) and `results/realism_provenance.json`. The NF1 §4.12 reframe is
thesis-adjacent and is held for your explicit decision.

Bib keys used (all already present in `refs.bib`, all abdominal IVIM patient cohorts):
`gurneychampion2018` (pancreatic), `kaandorp2021` (pancreatic), `barbieri2020` (abdominal DL-IVIM).

---

## NF1 — §4.12 reframe (THESIS-ADJACENT — STOP for decision)

**Verdict from the count-matched control: PARTIAL — both rare-bin calibration sparsity
AND a genuine residual wall contribute.**

> NUMBERS BELOW ARE LOCKED to the gated 16-seed bands (`results/realism_multiseed.json`,
> n_seeds=16; GATE 3 PASS). seed-0 byte-identity preserved (native 0.4514...).

Native fixed-edge hi-D\* recal coverage **0.533**\,[0.480, 0.594] → count-matched
**0.661**\,[0.580, 0.721] (uniform-prior anchor **0.807**\,[0.774, 0.831]); **47%**\,[26%, 64%]
of the native→uniform gap closed by count-matching alone. hi-D\* calibration count
**73 → 504**. Fidelity gate PASSES on all 16 seeds (up-sampled voxels KS-closer to native
realistic than to uniform; differ from uniform). Decisively, count-matched coverage stays
**≤0.72 in every seed** — never reaching the 0.807 anchor or 0.90 nominal. The count-matched
easy bins over-cover (lo 0.997, mid 0.922), the expected signature of a single global
conformal offset widened to serve the rare tail.

**Location:** `gauge_v3_revised.tex:919–922` (the rarity-attribution sentence).

**Current:**
> The severest under-coverage concentrates in a clinically \emph{rare}
> high-$\Dstar$ tail (fixed high-$\Dstar$ prevalence $0.049$\band{0.044}{0.056}):
> under a liver-like prior the ill-posed regime is sparsely populated -- a scoping
> caveat we report as such, not as a weakened wall.

**Proposed replacement (PARTIAL variant; numbers LOCKED to the gated 16-seed bands):**
> The severest under-coverage concentrates in a clinically \emph{rare} high-$\Dstar$
> tail (fixed high-$\Dstar$ prevalence $0.049$\band{0.044}{0.056}), which entangles
> two effects: the rarity of the bin leaves the recalibration split few high-$\Dstar$
> examples for the single global conformal offset, and the identifiability wall itself.
> A calibration-size control separates them. Up-sampling the fixed high-$\Dstar$
> \emph{calibration} count to the uniform-prior level ($\sim\!73\to504$ voxels), drawn
> from the realistic generator at fixed within-bin distribution -- a fidelity gate
> confirms the up-sampled voxels match the realistic, not the uniform, high-$\Dstar$
> law -- with train, test and prior held fixed, lifts fixed-edge high-$\Dstar$ coverage
> from $0.533$ only to $0.661$\band{0.580}{0.721}, closing $\approx47\%$\band{26\%}{64\%}
> of the way to the uniform-prior anchor $0.807$\band{0.774}{0.831} and never reaching it
> (count-matched coverage stays $\le0.72$ in every seed; Fig.~\ref{fig:realism}D). Rare-bin
> calibration sparsity thus accounts for roughly half of the fixed-edge deficit; the
> remainder is a genuine identifiability wall that survives count-matching. The rarity is
> a scoping caveat on \emph{prevalence} and a partial contributor to the fixed-edge
> number, not an explanation of the wall.

**Figure caption addendum** (`gauge_v3_revised.tex:949`, end of caption, before the
seed line): add panel D —
> (D)~NF1 calibration-size control: fixed-edge high-$\Dstar$ recalibrated coverage
> under native vs.\ count-matched calibration (hi-$\Dstar$ calibration count
> $73\to504$), against the uniform-prior anchor; only $\approx47\%$ of the
> native$\to$uniform gap is closed by count-matching, isolating a residual wall.

---

## NF1 — cross-reference (sub-unity CRLB ratio + coverage gap → §4.8 dissociation)

**Location:** `gauge_v3_revised.tex:917–919` (immediately after the CRLB-ratio clause).

**Current:**
> ... with the fixed-edge CRLB-to-bin-width ratio reaching $0.895$\band{0.849}{0.961}
> at high $\Dstar$ (Fig.~\ref{fig:realism}A,B).

**Proposed — append one sentence:**
> ... with the fixed-edge CRLB-to-bin-width ratio reaching $0.895$\band{0.849}{0.961}
> at high $\Dstar$ (Fig.~\ref{fig:realism}A,B). As in the shrinkage dissociation
> (Sec.~\ref{sec:dissociation-shrinkage}), a sub-unity CRLB-to-width ratio coexisting
> with a large conditional-coverage gap is expected, not contradictory: the ratio
> bounds \emph{point} precision, whereas the gap is a \emph{conditional-coverage}
> failure on a latent axis -- the same point-precision-vs-coverage dissociation,
> here with the global conformal offset set by a calibration split that is
> overwhelmingly low/mid-$\Dstar$.

---

## NF2 — rarity does not lower clinical stakes (one line)

**Location:** `gauge_v3_revised.tex:922`, appended to the rarity sentence (or as the
next sentence in the §4.12 paragraph).

**Proposed — add:**
> That the regime is rare does not make it low-stakes: a high-$\Dstar$ voxel is often
> the salient finding (brisk perfusion in a lesion), so under-coverage concentrated in
> the rare tail falls exactly where a wide, honest interval matters most.

---

## NF3 — prior provenance sentence + SI block (provenance fix)

### (a) Provenance sentence
**Location:** `gauge_v3_revised.tex:895–898`.

**Current:**
> a clinically realistic \emph{joint} prior -- correlated, right-skewed marginals
> (log-normal $D,\Dstar$; logit-normal $f$) with a within-subject $\Dstar$
> coefficient of variation in the published $50$--$110\%$ band, shaped to
> representative abdominal IVIM distributions -- and a measurement-level nuisance

**Proposed replacement:**
> a clinically realistic \emph{joint} prior -- correlated, right-skewed marginals
> (log-normal $D,\Dstar$; logit-normal $f$) with a within-subject $\Dstar$
> coefficient of variation in the published $50$--$110\%$ band, with hyperparameters
> \emph{shaped to approximate} (not fitted to) representative abdominal IVIM patient
> cohorts~\citep{gurneychampion2018,kaandorp2021,barbieri2020}; the exact joint
> specification is documented design constants (Table~\ref{tab:realism-prior}, Supporting
> Information) and deliberately avoids any circular plug-in source -- and a
> measurement-level nuisance

### (b) SI table (add to the Supporting Information section, `gauge_v3_revised.tex:1127+`)
```latex
\begin{table}[h]
\centering
\caption{Realistic joint-prior specification (Arm~A / NF1 control). Hyperparameters are
documented design constants \emph{shaped to approximate} representative abdominal IVIM
patient cohorts~\citep{gurneychampion2018,kaandorp2021,barbieri2020}, not fitted; the
joint is a Gaussian copula on latent normals with the marginals below. Verbatim in
\texttt{results/realism\_report.txt}.}
\label{tab:realism-prior}
\begin{tabular}{lll}
\hline
Parameter & Marginal & Hyperparameters \\
\hline
$D$        & log-normal   & mean $1.1\times10^{-3}$ mm$^2$/s, CV $0.30$ \\
$\Dstar$   & log-normal   & mean $28.0\times10^{-3}$ mm$^2$/s, CV $0.80$ (target band $0.50$--$1.10$) \\
$f$        & logit-normal & logit-mean $\ln(0.20/0.80)$, logit-sd $0.55$ (median $f\approx0.20$) \\
SNR$(b{=}0)$ & log-normal & mean $35.0$, CV $0.40$, clipped to $[8.0, 120.0]$ \\
\hline
\multicolumn{3}{l}{latent-normal correlation $(D,\Dstar,f)$:
$\begin{psmallmatrix}1.00 & 0.00 & -0.20\\ 0.00 & 1.00 & 0.35\\ -0.20 & 0.35 & 1.00\end{psmallmatrix}$}\\
\hline
\end{tabular}
\end{table}
```
(`\begin{psmallmatrix}` needs `amsmath`+`mathtools`; if `mathtools` is unavailable, use a
plain `\left(\begin{array}{rrr}...\end{array}\right)`.)

---

## NF4 — Arm B D\* stability = interval width, not robust identification (one line)

**Location:** `gauge_v3_revised.tex:930–932` (the "$\Dstar$ interval stays near nominal"
clause).

**Current:**
> while the already-wide $\Dstar$ interval stays near nominal
> ($0.899$\band{0.888}{0.911}), and recalibrating within the nuisance-laden
> distribution restores it ($\Dstar$ $0.899$\band{0.886}{0.914} at $\eta=1$).

**Proposed — append:**
> ... at $\eta=1$). This $\Dstar$ stability reflects interval \emph{width}, not robust
> identification: the high-$\Dstar$ band is already wide enough to absorb the added
> nuisance, so flat $\Dstar$ coverage is the characterized-but-unresolved wall showing
> through, not evidence that $\Dstar$ is well-determined under nuisance.

---

## §7 reproducibility clause (optional, no new entry point)

**Location:** `gauge_v3_revised.tex:1044–1051` — the control runs inside the existing
`\texttt{gauge.realism}` / `\texttt{gauge.realism sweep 16}` commands; no new entry point.
Optionally add, after the realism command mention:
> (the realism run also emits the NF1 calibration-size control and the verbatim
> prior-provenance table).
