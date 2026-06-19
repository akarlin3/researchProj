# The ruler as a standard — and Datum's role in the Casali differentiation

This note states, for Fashion's revision, what Datum *is* with respect to the
calibration **ruler-as-standard**, and how that differentiates Fashion from a
method-paper like Casali (2026). It is a positioning document, **not** a
novel-result claim.

## Two different kinds of object

- A **method** answers "here is a new way to produce IVIM uncertainty" (e.g. Casali
  2026: a deep-ensemble mixture-density network that yields aleatoric/epistemic
  uncertainty; it documents D\* overconfidence as a property of that method).
- A **standard** answers "here is the fixed yardstick by which *any* method's
  uncertainty is judged honest." Fashion's contribution is of this second kind: a
  calibration ruler — coverage / ECE / sharpness at frozen nominal levels — that is
  method-agnostic.

A method can only report its own calibration. A standard lets you put *every*
method on the same axis. That is the difference Fashion draws against Casali:
Fashion is not proposing yet another UQ method; it is fixing the ruler.

## Datum operationalizes the standard

Datum is the artifact that makes "Fashion's ruler is a standard" concrete and
checkable:

1. It **freezes the ruler into a benchmark** — a fixed task, a curated baseline
   panel, reference numbers, and a submission interface (`datum.submit`).
2. It scores **diverse methods on one axis**: parametric Gaussian error bars, a
   segmented reference, conformal corrections (split / CQR / Mondrian), a learned
   normalizing-flow posterior, and *any* submitted method — all judged by the same
   ruler. A Casali-style MDN method would be scored the same way, by submitting its
   quantiles; no special-casing.
3. It reproduces, as benchmark numbers, the substantive findings the standard was
   built to expose: Gaussian/segmented/learned error bars **under-cover D\***;
   conformal restores marginal coverage; and the **high-D\* identifiability wall**
   persists for marginal conformal and is bought back by CQR/Mondrian only by
   inflating width. (See [`../results/REFERENCE.md`](../results/REFERENCE.md).)

So the differentiation is demonstrated, not just asserted: "ruler-as-standard"
means *there exists a benchmark on which any method, including a Casali-style one,
is scored against Fashion's yardstick* — and Datum is that benchmark.

## Honest scope

- Datum makes **no claim** that any particular method (Fashion's, Gauge's, or a
  submission's) is best; it reports what each scores.
- Every number is **PROVISIONAL** because the ruler is in review. If Fashion's
  ruler changes in revision, the standard changes and Datum re-validates
  (`python -m datum.run`). The differentiation argument (standard vs method) is
  structural and survives that; the specific numbers do not until the ruler locks.
- The external **OSIPI DRO** validation shows the calibration story is not an
  artifact of the authors' own forward model — it holds on an independent synthetic
  phantom.
