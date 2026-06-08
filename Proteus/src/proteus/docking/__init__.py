"""Docking — AutoDock Vina wrappers (CPU, P4).

Local CPU docking of PET-mimic ligands into ranked candidate clefts from S5.
GPU docking (GNINA / DiffDock) and Chai-1 cofolding are Vast.ai burst targets,
not part of this local module — see vast/.
"""

# TODO(P4): receptor prep (Meeko/ADFR), define box from S5 cleft centroid, run Vina,
#           parse affinities; seed any stochastic search from config.random_seed.
