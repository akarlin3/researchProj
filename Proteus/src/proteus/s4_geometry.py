"""S4 — Catalytic geometry gate: Ser-His-Asp triad + oxyanion hole.

Given folded models from S3, test for a serine-hydrolase catalytic machine:
the Ser-His-Asp(/Glu) triad in hydrogen-bonding geometry plus a backbone
oxyanion hole. Geometry distance windows come from config (s4_geometry).

This is a positive structural assertion — a candidate passes only if the
triad and oxyanion hole are present within tolerance.
"""

# TODO(P3): parse PDBs with biotite/biopython, enumerate SER/HIS/(ASP|GLU) triples,
#           score Ser-OG..His-NE2, His-ND1..acid, and oxyanion-hole donor distances;
#           emit passing residue triads with measured geometry.
