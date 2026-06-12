# POSITIONING — what Minos is *not*

Stub. No results claims here — those live in `RESULTS.md`. This file fixes the
must-distinguish neighbours so the novelty (pricing the **error bar** via VoC and VoTG) is
not silently reframed as one of them.

Minos prices **the per-voxel error bar itself** — its *width* (Value of Calibration) and its
*trustworthiness under shift* (Value of the Trust-Gate) — on a treat / spare / escalate
decision. The object being priced is the uncertainty, not the point estimate, and not a
population parameter.

| Neighbour | What it prices | How Minos differs |
|---|---|---|
| **Vickers decision-curve analysis / net benefit** | The decision value of a **marker's point prediction** across threshold probabilities; "is acting on this classifier better than treat-all / treat-none?" | Minos holds the point estimate fixed and prices the **calibration and trustworthiness of the error bar around it** (VoC, VoTG). Net benefit has no notion of the uncertainty being mis-scaled or untrustworthy; VoC is identically invisible to it. |
| **ISPOR VoI canon (EVPI / EVPPI / EVSI)** | The population value of **resolving or reducing parameter uncertainty** by collecting more data — what a perfect/partial/sample study is worth. | Minos prices a **per-voxel, already-reported** error bar at decision time: VoC is the loss from the reported bar being mis-*scaled*; VoTG is the gain from *detecting* that the bar is untrustworthy. No data-collection or parameter-learning step is being valued. (Minos does carry an EVPI-*analog* as a degenerate-limit sanity anchor, explicitly labelled as such — it is not the headline.) |
| **ARCliDS adaptive-RT CDSS** | An end-to-end **clinical decision-support / RL policy** that recommends adaptive radiotherapy actions from patient state. | Minos is not a policy-learning system or a CDSS; it is a **formal accounting of the decision-value of an uncertainty estimate**. It asks "what is this error bar worth, and when should it be distrusted?", independent of how any policy or estimator was trained. |

**One-line guard.** If an output of Minos is ever described as "net benefit," "decision-curve
analysis," or "EVPI/EVPPI of a parameter," it has been mis-framed: the priced object is the
error bar (VoC / VoTG), not the point estimate and not a population parameter.
