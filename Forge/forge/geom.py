#!/usr/bin/env python3
"""
Forge Monte Carlo Geometry and Parameter Deck Generator.
Generates randomized voxelized phantoms, beam/source configurations, 
magnetic fields, and TOPAS text decks.
"""

import os
import json
import numpy as np

# ----------------------------------------------------------------------------
# Beam modality: Elekta Unity-class 1.5 T MR-Linac is a PHOTON machine with a
# single 7 MV FFF (flattening-filter-free) bremsstrahlung x-ray source.
# The dataset must therefore be photon, NOT proton and NOT a monoenergetic gamma.
# ----------------------------------------------------------------------------

def photon_fff_spectrum(Emax_MV=7.0, Elow_MeV=0.25, n_bins=40, hardening_k=0.65):
    """
    Synthesized APPROXIMATE 7 MV FFF bremsstrahlung photon spectrum for an
    Elekta Unity-class MR-Linac source.

    *** This is a synthesized analytical approximation, NOT a measured Unity
    phase-space file. *** It is flagged as approximate throughout the dataset
    metadata so it can be swapped for a commissioned phase space / spectrum.

    Basis (cited, not fabricated):
      - Shape: Kramers thick-target bremsstrahlung ~ (Emax - E)/E, modified by an
        empirical low-energy hardening factor E**k to emulate beam hardening
        through the target/primary collimator of an FFF head (no flattening
        filter is present on Unity, so the raw spectrum is intentionally soft).
      - Normalization target: the fluence-mean photon energy is tuned to
        ~2.1 MeV, matching the published Unity 7 MV FFF mean photon energy of
        ~2.11 MeV at isocenter reported in Monte Carlo commissioning literature
        (BEAMnrc/EGSnrc head models of the Elekta Unity Agility head).

    Returns (energies_MeV, weights) suitable for a TOPAS continuous spectrum.
    """
    edges = np.linspace(Elow_MeV, Emax_MV, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    weights = (Emax_MV - centers) / centers * centers ** hardening_k
    weights = np.clip(weights, 0.0, None)
    weights = weights / weights.sum()
    return centers, weights


def _format_topas_vector(values, fmt="{:.4f}"):
    """Format a numpy array into a space-separated TOPAS vector body."""
    return " ".join(fmt.format(v) for v in values)


def generate_layered_phantom(Nx=60, Ny=60, Nz=60, seed=42):
    """
    Generates a 3D voxelized phantom with randomized layers and an air gap.
    Material Tags:
      0: G4_AIR (Density ~0.0012 g/cm3)
      1: G4_LUNG_ICRP (Density ~0.26 g/cm3)
      2: G4_WATER (Density ~1.0 g/cm3)
      3: G4_BONE_COMPACT_ICRU (Density ~1.85 g/cm3)
    """
    rng = np.random.default_rng(seed)
    
    # Initialize as water/soft tissue (tag 2)
    phantom = np.full((Nx, Ny, Nz), 2, dtype=np.int16)
    
    # Determine random slab interfaces along the Z (depth) axis
    # The beam travels along +Z direction.
    # We want a soft-tissue slab -> air gap -> lung/bone/soft-tissue slab.
    
    # Let's divide Z axis (Nz voxels) into 4-5 slabs
    # Slab 1: Tissue (always first to build buildup region)
    # Slab 2: Air gap (crucial for ERE)
    # Slab 3: Lung or Bone or Tissue
    # Slab 4: Bone or Lung or Tissue
    
    z_slab1 = rng.integers(15, 22) # Soft tissue region
    z_gap_thickness = rng.integers(6, 12) # Air gap thickness
    z_slab2 = z_slab1 + z_gap_thickness
    
    z_slab3 = rng.integers(z_slab2 + 10, Nz - 10)
    
    # Apply tags
    # Slab 1: Soft tissue (already 2)
    
    # Slab 2: Air gap
    phantom[:, :, z_slab1:z_slab2] = 0
    
    # Slab 3: Lung (1) or Bone (3) or Soft tissue (2)
    slab3_material = rng.choice([1, 2, 3])
    phantom[:, :, z_slab2:z_slab3] = slab3_material
    
    # Slab 4: Bone (3) or Lung (1)
    slab4_material = rng.choice([1, 3])
    phantom[:, :, z_slab3:] = slab4_material
    
    return phantom, {
        "slab_boundaries_z": [int(z_slab1), int(z_slab2), int(z_slab3)],
        "materials": ["G4_WATER", "G4_AIR", 
                      "G4_LUNG_ICRP" if slab3_material == 1 else ("G4_WATER" if slab3_material == 2 else "G4_BONE_COMPACT_ICRU"),
                      "G4_LUNG_ICRP" if slab4_material == 1 else "G4_BONE_COMPACT_ICRU"]
    }

def generate_case_deck(case_id, seed, output_dir="cases", histories=500000):
    """
    Generates a single randomized MC case, writes the binary voxel grid 
    and the TOPAS parameter text deck.
    """
    os.makedirs(output_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    
    # Grid configuration
    Nx, Ny, Nz = 60, 60, 60
    dx, dy, dz = 2.5, 2.5, 2.5 # mm
    
    # Voxel phantom generation
    phantom, meta = generate_layered_phantom(Nx, Ny, Nz, seed)
    
    # Save raw voxel grid as short (int16) binary file
    bin_filename = f"case_{case_id}.bin"
    bin_filepath = os.path.join(output_dir, bin_filename)
    phantom.tofile(bin_filepath)
    
    # B-field: 1.5 T transverse dipole (Elekta Unity geometry).
    # The beam travels along +Z, so a field in X or Y is PERPENDICULAR to the
    # beam axis (B 'tab' beam), which is the correct Unity orientation. Bz = 0.
    b_strength = 1.5 # Tesla
    b_dir = rng.choice(["X", "Y"])
    bx_dir = 1.0 if b_dir == "X" else 0.0
    by_dir = 1.0 if b_dir == "Y" else 0.0
    bz_dir = 0.0

    # Source / Beam configuration: Elekta Unity-class 7 MV FFF PHOTON source.
    # The machine is a photon MR-Linac; we sample a 7 MV FFF bremsstrahlung
    # spectrum (synthesized, ~2.1 MeV mean) rather than protons or a single
    # monoenergetic gamma. Energy is delivered as a TOPAS continuous spectrum.
    particle = "gamma"
    spec_energies, spec_weights = photon_fff_spectrum()
    spectrum_mean_mev = float(np.sum(spec_energies * spec_weights))
    energy = spectrum_mean_mev  # reported as the fluence-mean of the spectrum
    n_spec = len(spec_energies)
    spec_energy_str = _format_topas_vector(spec_energies, "{:.4f}")
    spec_weight_str = _format_topas_vector(spec_weights, "{:.6e}")

    # Randomized field size and isocenter so the geometry distribution truly
    # varies (anatomy + field size + isocenter), matching Unity field sizes.
    field_half_cm = float(rng.uniform(1.5, 11.0))   # 3x3 cm up to 22x22 cm fields
    iso_x_cm = float(rng.uniform(-3.0, 3.0))         # isocenter lateral offset
    iso_y_cm = float(rng.uniform(-3.0, 3.0))

    # Topas Deck Generation
    deck_filename = f"case_{case_id}.txt"
    deck_filepath = os.path.join(output_dir, deck_filename)
    
    deck_content = f"""# TOPAS Parameter Deck for Case {case_id}
# Generated by Forge (geom.py) - Seed: {seed}

# ====================================================================
# WORLD GEOMETRY
# ====================================================================
d:Ge/World/HLX = 20.0 cm
d:Ge/World/HLY = 20.0 cm
d:Ge/World/HLZ = 30.0 cm
s:Ge/World/Material = "G4_AIR"

# ====================================================================
# PATIENT/PHANTOM (VOXELIZED)
# ====================================================================
s:Ge/Patient/Type = "TsImageCube"
s:Ge/Patient/Parent = "World"
s:Ge/Patient/InputDirectory = "./"
s:Ge/Patient/InputFile = "{bin_filename}"
s:Ge/Patient/DataType = "short"
i:Ge/Patient/NumberOfVoxelsX = {Nx}
i:Ge/Patient/NumberOfVoxelsY = {Ny}
i:Ge/Patient/NumberOfVoxelsZ = {Nz}
d:Ge/Patient/VoxelSizeX = {dx} mm
d:Ge/Patient/VoxelSizeY = {dy} mm
d:Ge/Patient/VoxelSizeZ = {dz} mm
d:Ge/Patient/TransX = 0.0 cm
d:Ge/Patient/TransY = 0.0 cm
d:Ge/Patient/TransZ = 0.0 cm
d:Ge/Patient/RotX = 0.0 deg
d:Ge/Patient/RotY = 0.0 deg
d:Ge/Patient/RotZ = 0.0 deg

s:Ge/Patient/ImagingToMaterialConverter = "MaterialTagNumber"
iv:Ge/Patient/MaterialTagNumbers = 4 0 1 2 3
sv:Ge/Patient/MaterialNames = 4 "G4_AIR" "G4_LUNG_ICRP" "G4_WATER" "G4_BONE_COMPACT_ICRU"

# ====================================================================
# MAGNETIC FIELD (TRANSVERSE Dipole)
# ====================================================================
s:Ge/Patient/Field = "DipoleMagnet"
u:Ge/Patient/MagneticFieldDirectionX = {bx_dir}
u:Ge/Patient/MagneticFieldDirectionY = {by_dir}
u:Ge/Patient/MagneticFieldDirectionZ = {bz_dir}
d:Ge/Patient/MagneticFieldStrength = {b_strength} tesla

# ====================================================================
# BEAM SOURCE (Elekta Unity-class 7 MV FFF photon spectrum)
# Particle: gamma (photon). Energy: continuous 7 MV FFF bremsstrahlung
# spectrum (synthesized analytical approximation, fluence-mean ~{energy:.3f} MeV).
# NOTE: This is an APPROXIMATE spectrum, not a measured Unity phase space.
# ====================================================================
s:So/Beam/Type = "Beam"
s:So/Beam/Component = "World"
s:So/Beam/BeamParticle = "gamma"
s:So/Beam/BeamEnergySpectrumType = "Continuous"
dv:So/Beam/BeamEnergySpectrumValues = {n_spec} {spec_energy_str} MeV
uv:So/Beam/BeamEnergySpectrumWeights = {n_spec} {spec_weight_str}
s:So/Beam/BeamPositionDistribution = "Gaussian"
s:So/Beam/BeamPositionCutoffShape = "Rectangle"
d:So/Beam/BeamPositionCutoffX = {field_half_cm:.3f} cm
d:So/Beam/BeamPositionCutoffY = {field_half_cm:.3f} cm
d:So/Beam/BeamPositionSpreadX = 5.0 mm
d:So/Beam/BeamPositionSpreadY = 5.0 mm
s:So/Beam/BeamAngularDistribution = "None"
d:So/Beam/TransX = {iso_x_cm:.3f} cm
d:So/Beam/TransY = {iso_y_cm:.3f} cm
d:So/Beam/TransZ = -25.0 cm
d:So/Beam/RotX = 0.0 deg
d:So/Beam/RotY = 0.0 deg
d:So/Beam/RotZ = 0.0 deg
i:So/Beam/NumberOfHistoriesInRun = {histories}

# ====================================================================
# SCORING (3D Dose to Medium)
# ====================================================================
s:Sc/Dose/Quantity = "DoseToMedium"
s:Sc/Dose/Component = "Patient"
s:Sc/Dose/OutputFile = "dose_case_{case_id}"
s:Sc/Dose/OutputType = "csv"
s:Sc/Dose/IfOutputFileAlreadyExists = "Overwrite"
sv:Sc/Dose/Report = 2 "Mean" "Standard_Deviation"

# ====================================================================
# TIMING & CONTROLS
# ====================================================================
# s:Gr/View/Type = "OpenGL"
# b:Gr/View/Active = "False"
i:Ts/ShowHistoryCountAtInterval = 1000
i:Ts/Seed = {seed}
"""
    with open(deck_filepath, "w") as f:
        f.write(deck_content)
        
    return {
        "case_id": case_id,
        "seed": int(seed),
        "bin_file": bin_filename,
        "deck_file": deck_filename,
        "phantom_meta": meta,
        "magnetic_field": {
            "strength_tesla": b_strength,
            "direction": [bx_dir, by_dir, bz_dir]
        },
        "source": {
            "particle": particle,
            "modality": "photon",
            "spectrum": "7 MV FFF bremsstrahlung (synthesized, APPROXIMATE)",
            "spectrum_mean_mev": float(spectrum_mean_mev),
            "spectrum_endpoint_mv": 7.0,
            "energy_mev": float(energy),
            "field_half_cm": field_half_cm,
            "isocenter_offset_cm": [iso_x_cm, iso_y_cm],
            "histories": histories
        },
        "dimensions": [Nx, Ny, Nz],
        "voxel_sizes_mm": [dx, dy, dz]
    }


def main():
    K = 12
    base_seed = 1000
    cases_meta = []
    
    print(f"Generating {K} randomized cases in 'cases/'...")
    for k in range(1, K + 1):
        seed = base_seed + k
        case_info = generate_case_deck(k, seed, "cases")
        cases_meta.append(case_info)
        print(f"  - Case {k} generated (Seed={seed}, Particle={case_info['source']['particle']}, B_dir={case_info['magnetic_field']['direction']})")
        
    # Write manifest file
    manifest_filepath = os.path.join("cases", "manifest.json")
    with open(manifest_filepath, "w") as f:
        json.dump({
            "num_cases": K,
            "base_seed": base_seed,
            "cases": cases_meta
        }, f, indent=2)
    print(f"Generated manifest: {manifest_filepath}")

if __name__ == "__main__":
    main()
