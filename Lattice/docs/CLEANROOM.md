# Clean-room / IP statement

Lattice packages synthetic IVIM machinery of the kind also present in `Gauge/`.
It was built **clean-room** to keep a strict one-way dependency and avoid
entangling any non-synthetic data.

## What "clean-room" means here

The five forward models are **re-implementations of the published IVIM signal
equations** (physics, not borrowed code): the bi-exponential model, the Laplace
transform of a Gamma pseudo-diffusivity, a Gauss-Hermite quadrature over a
log-normal pseudo-diffusivity, a stretched-exponential perfusion envelope, and a
three-compartment sum. These equations are standard literature; the
implementation in `lattice/generators.py` is original.

## Synthetic-only, in tree

- **No clinical / patient / scanner data** is present in this repository — no
  DICOM, no NIfTI, no `pancData3`, no MSK data, no ACRIN/in-vivo arrays.
- All cohorts are generated from seeds at runtime; nothing real is committed.
- The **optional** OSIPI integration (`lattice/osipi.py`,
  `scripts/fetch_osipi.py`) fetches external CC-BY-4.0 reference data
  *on demand* and writes a provenance manifest; it is **never redistributed**
  in-tree, mirroring Gauge's posture. Importing the module touches no network.

Notably, the only non-synthetic asset anywhere in Gauge — optional ACRIN-6698
in-vivo breast DWI via download-on-demand — was **deliberately not extracted**
into Lattice.

## One-way dependency

```
Caliper, papers (Fashion/Gauge/Minos), Vernier, Echo  ──consume──►  Lattice
Lattice  ──imports──►  (nothing back)
```

The `lattice` core package imports only numpy (scipy/requests are optional
extras for self-check / external data). It imports **nothing** from Caliper or
any paper project. The single place Caliper appears is
`examples/evaluate_with_caliper.py`, an optional example that demonstrates the
consumption direction; it degrades gracefully if Caliper is not importable.
