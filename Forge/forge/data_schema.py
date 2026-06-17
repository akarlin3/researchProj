#!/usr/bin/env python3
"""
Forge FNO HDF5 Data Schema Contract.
Defines the on-disk data format shared by the Monte Carlo run and FNO trainer.
"""

import os
import h5py
import numpy as np

def write_case_to_h5(
    file_path,
    case_id,
    density_grid,
    magnetic_field_grid,
    dose_grid,
    metadata
):
    """
    Writes a single case's inputs, target, and metadata to the HDF5 file.
    
    Args:
        file_path (str): Path to the HDF5 file.
        case_id (int): Unique identifier for the case.
        density_grid (np.ndarray): Shape (Nx, Ny, Nz), float32. Voxel densities.
        magnetic_field_grid (np.ndarray): Shape (Nx, Ny, Nz, 3), float32. B-field vector at each voxel.
        dose_grid (np.ndarray): Shape (Nx, Ny, Nz), float32. Target 3D dose.
        metadata (dict): Key-value metadata for the case.
    """
    # Open file in append mode (creates if not exists)
    with h5py.File(file_path, "a") as f:
        group_name = f"case_{case_id:05d}"
        
        # If group already exists, delete it to overwrite cleanly
        if group_name in f:
            del f[group_name]
            
        g = f.create_group(group_name)
        
        # Create datasets
        g.create_dataset("density", data=density_grid.astype(np.float32), compression="gzip", compression_opts=4)
        g.create_dataset("magnetic_field", data=magnetic_field_grid.astype(np.float32), compression="gzip", compression_opts=4)
        g.create_dataset("dose", data=dose_grid.astype(np.float32), compression="gzip", compression_opts=4)
        
        # Write metadata as attributes
        for key, val in metadata.items():
            if isinstance(val, str):
                g.attrs[key] = bytes(val, 'utf-8') # HDF5 handles ASCII strings best as bytes
            else:
                g.attrs[key] = val

def read_case_from_h5(file_path, case_id):
    """
    Reads a single case from the HDF5 file and returns datasets and metadata.
    
    Args:
        file_path (str): Path to the HDF5 file.
        case_id (int): Case number to read.
        
    Returns:
        tuple: (density_grid, magnetic_field_grid, dose_grid, metadata)
    """
    with h5py.File(file_path, "r") as f:
        group_name = f"case_{case_id:05d}"
        if group_name not in f:
            raise KeyError(f"Case {case_id} not found in {file_path}")
            
        g = f[group_name]
        
        # Read datasets into memory
        density = g["density"][:]
        magnetic_field = g["magnetic_field"][:]
        dose = g["dose"][:]
        
        # Read attributes
        metadata = {}
        for key, val in g.attrs.items():
            if isinstance(val, bytes):
                metadata[key] = val.decode("utf-8")
            else:
                metadata[key] = val
                
        return density, magnetic_field, dose, metadata

def list_cases_in_h5(file_path):
    """Lists all cases present in the HDF5 file."""
    if not os.path.exists(file_path):
        return []
    with h5py.File(file_path, "r") as f:
        return sorted(list(f.keys()))

def validate_schema():
    """Runs a round-trip validation by writing and reading 2 dummy cases."""
    test_h5 = "test_fno_contract.h5"
    if os.path.exists(test_h5):
        os.remove(test_h5)
        
    print(f"Validating FNO HDF5 Schema Contract on '{test_h5}'...")
    
    Nx, Ny, Nz = 60, 60, 60
    
    # Create dummy data for Case 1
    rng1 = np.random.default_rng(42)
    density_1 = rng1.uniform(0.0012, 1.85, (Nx, Ny, Nz))
    magnetic_1 = rng1.uniform(-1.5, 1.5, (Nx, Ny, Nz, 3))
    dose_1 = rng1.exponential(1.0, (Nx, Ny, Nz))
    meta_1 = {
        "case_id": 1,
        "particle": "gamma",
        "beam_energy_mev": 2.09,  # fluence-mean of 7 MV FFF spectrum
        "seed": 1001,
        "histories": 5000,
        "core_hours": 0.45
    }
    
    # Create dummy data for Case 2
    rng2 = np.random.default_rng(43)
    density_2 = rng2.uniform(0.0012, 1.85, (Nx, Ny, Nz))
    magnetic_2 = rng2.uniform(-1.5, 1.5, (Nx, Ny, Nz, 3))
    dose_2 = rng2.exponential(1.0, (Nx, Ny, Nz))
    meta_2 = {
        "case_id": 2,
        "particle": "gamma",
        "beam_energy_mev": 6.0,
        "seed": 1002,
        "histories": 10000,
        "core_hours": 0.90
    }
    
    # Write to HDF5
    write_case_to_h5(test_h5, 1, density_1, magnetic_1, dose_1, meta_1)
    write_case_to_h5(test_h5, 2, density_2, magnetic_2, dose_2, meta_2)
    
    print("  - Successfully wrote 2 dummy cases to HDF5.")
    
    # Read back and validate
    cases = list_cases_in_h5(test_h5)
    print(f"  - Detected cases in file: {cases}")
    assert len(cases) == 2, f"Expected 2 cases, got {len(cases)}"
    
    # Validate Case 1
    d1, m1, ds1, mt1 = read_case_from_h5(test_h5, 1)
    np.testing.assert_array_almost_equal(d1, density_1)
    np.testing.assert_array_almost_equal(m1, magnetic_1)
    np.testing.assert_array_almost_equal(ds1, dose_1)
    assert mt1["particle"] == "gamma"
    assert mt1["beam_energy_mev"] == 2.09
    assert mt1["core_hours"] == 0.45
    print("  - Case 1 readback validation: PASS")
    
    # Validate Case 2
    d2, m2, ds2, mt2 = read_case_from_h5(test_h5, 2)
    np.testing.assert_array_almost_equal(d2, density_2)
    np.testing.assert_array_almost_equal(m2, magnetic_2)
    np.testing.assert_array_almost_equal(ds2, dose_2)
    assert mt2["particle"] == "gamma"
    assert mt2["beam_energy_mev"] == 6.0
    assert mt2["core_hours"] == 0.90
    print("  - Case 2 readback validation: PASS")
    
    # Clean up test file
    if os.path.exists(test_h5):
        os.remove(test_h5)
    print("HDF5 Schema Contract validation complete: SUCCESS!")

if __name__ == "__main__":
    validate_schema()
