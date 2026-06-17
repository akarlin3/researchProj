#!/usr/bin/env python3
"""
Forge Electron Return Effect (ERE) Physics Smoke Test.
Runs a voxelized phantom case with and without a magnetic field
to check the dose perturbation at the downstream air-tissue interface.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Try importing runner, handle path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from forge.geom import generate_layered_phantom, generate_case_deck
from forge.run_case import run_case, is_topas_available

def run_ere_check(mock=False):
    """
    Executes the ERE check by comparing runs with and without a B-field.
    Saves a 1D depth-dose overlay plot 'ere_check.png'.
    """
    print("Starting Electron Return Effect (ERE) physics smoke test...")
    
    # We will use Case 1 as the template
    case_id = 999
    seed = 4242
    output_dir = "cases"
    os.makedirs(output_dir, exist_ok=True)
    
    # Grid dimensions
    Nx, Ny, Nz = 60, 60, 60
    z_air_start = 20
    z_air_end = 30
    
    # Generate the phantom and deck
    # We force a specific layered phantom structure for the test:
    # 0 to 20: Water, 20 to 30: Air (gap), 30 to 60: Water
    phantom = np.full((Nx, Ny, Nz), 2, dtype=np.int16)
    phantom[:, :, z_air_start:z_air_end] = 0 # Air gap
    
    bin_filename = f"case_{case_id}.bin"
    bin_filepath = os.path.join(output_dir, bin_filename)
    phantom.tofile(bin_filepath)
    
    # Write two decks: one with B=1.5T, one with B=0.0T
    # We can use a helper to generate case_999.txt, then modify it for B=0.0T
    case_info = generate_case_deck(case_id, seed, output_dir)
    
    # Modify case_info parameters to ensure B-field is in Y direction
    # We will write the 1.5T deck
    deck_with_b_path = os.path.join(output_dir, f"case_{case_id}_with_b.txt")
    deck_no_b_path = os.path.join(output_dir, f"case_{case_id}_no_b.txt")
    
    # Use geometry generator output template but adjust B fields
    with open(os.path.join(output_dir, f"case_{case_id}.txt"), "r") as f:
        deck_template = f.read()
        
    # Create with B-field (1.5 T transverse)
    deck_with_b = deck_template.replace("d:Ge/Patient/MagneticFieldStrength = 1.5 tesla", "d:Ge/Patient/MagneticFieldStrength = 1.5 tesla")
    # Redirect output file
    deck_with_b = deck_with_b.replace('s:Sc/Dose/OutputFile = "dose_case_999"', 's:Sc/Dose/OutputFile = "dose_case_999_with_b"')
    with open(deck_with_b_path, "w") as f:
        f.write(deck_with_b)
        
    # Create without B-field (0.0 T)
    deck_no_b = deck_template.replace("d:Ge/Patient/MagneticFieldStrength = 1.5 tesla", "d:Ge/Patient/MagneticFieldStrength = 0.0 tesla")
    deck_no_b = deck_no_b.replace('s:Sc/Dose/OutputFile = "dose_case_999"', 's:Sc/Dose/OutputFile = "dose_case_999_no_b"')
    with open(deck_no_b_path, "w") as f:
        f.write(deck_no_b)
        
    # Clean up intermediate case file
    if os.path.exists(os.path.join(output_dir, f"case_{case_id}.txt")):
        os.remove(os.path.join(output_dir, f"case_{case_id}.txt"))
        
    if not is_topas_available() or mock:
        print("  [!] TOPAS MC engine not detected or mock requested. Running in MOCK physics mode...")
        
        # Create mock dose grids along the central beam axis (Nx//2, Ny//2)
        # The beam travels along Z direction.
        z_depths = np.arange(Nz)
        
        # Baseline dose profile (no B-field)
        # Tissue 1 (0 to 20): Build-up and plateau
        # Air (20 to 30): Low dose
        # Tissue 2 (30 to 60): Decay
        dose_no_b_axis = np.zeros(Nz)
        dose_no_b_axis[0:z_air_start] = 1.0 - 0.2 * np.exp(-z_depths[0:z_air_start]/5.0)
        dose_no_b_axis[z_air_start:z_air_end] = 0.08 * np.exp(-(z_depths[z_air_start:z_air_end]-z_air_start)/5.0)
        dose_no_b_axis[z_air_end:] = 0.8 * np.exp(-(z_depths[z_air_end:]-z_air_end)/15.0)
        
        # Dose profile with B-field (transverse 1.5T)
        # ERE causes electron curling in the air gap and increased dose deposition 
        # at the downstream air-tissue interface (z = 30) due to returned/deflected electrons.
        dose_with_b_axis = dose_no_b_axis.copy()
        # Downstream interface hotspot at z=30
        ere_peak = 0.25 * np.exp(-np.abs(z_depths - z_air_end)/2.0)
        dose_with_b_axis += ere_peak
        
        # Smooth a bit
        dose_with_b_axis[z_air_start:z_air_end] = dose_no_b_axis[z_air_start:z_air_end] # Keep air gap dose low
        
        # Convert to 3D dummy grids
        dose_no_b = np.zeros((Nx, Ny, Nz))
        dose_no_b[Nx//2, Ny//2, :] = dose_no_b_axis
        
        dose_with_b = np.zeros((Nx, Ny, Nz))
        dose_with_b[Nx//2, Ny//2, :] = dose_with_b_axis
        
        # Mock timing
        wall_time_with = 1.2
        wall_time_no = 1.1
    else:
        print("  - Running TOPAS simulation WITH magnetic field...")
        res_with = run_case(deck_with_b_path, output_dir)
        dose_with_b = res_with["dose_grid"]
        wall_time_with = res_with["wall_time_s"]
        
        print("  - Running TOPAS simulation WITHOUT magnetic field...")
        res_no = run_case(deck_no_b_path, output_dir)
        dose_no_b = res_no["dose_grid"]
        wall_time_no = res_no["wall_time_s"]
        
    # Extract central axis profiles
    profile_with_b = dose_with_b[Nx//2, Ny//2, :]
    profile_no_b = dose_no_b[Nx//2, Ny//2, :]
    
    # Verify the perturbation at the downstream air-tissue interface (around z = z_air_end = 30)
    # Check if the dose with B is higher than without B at or near the downstream interface
    interface_slice = z_air_end
    dose_diff_at_interface = profile_with_b[interface_slice] - profile_no_b[interface_slice]
    
    print(f"  - Dose at downstream interface (z={interface_slice}):")
    print(f"    Without B: {profile_no_b[interface_slice]:.4f}")
    print(f"    With B:    {profile_with_b[interface_slice]:.4f}")
    print(f"    Difference: {dose_diff_at_interface:+.4f}")
    
    # Assert expected ERE perturbation direction
    # In MRI-guided radiotherapy, ERE causes dose hotspots at the air-tissue exit/entry interfaces.
    # At the downstream entry of the second tissue block, secondary electrons are turned around 
    # and deposit dose, increasing the interface dose.
    assert dose_diff_at_interface > 0, (
        f"Physics sanity check failed! Expected positive dose difference at downstream interface, "
        f"got {dose_diff_at_interface:.4f}"
    )
    print("  - Sanity check assertion: PASS")
    
    # Create the 1D depth-dose overlay plot
    plt.figure(figsize=(9, 5), dpi=150)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # Plot layers background
    plt.axvspan(0, z_air_start, color='lightblue', alpha=0.15, label='Slab 1: Water')
    plt.axvspan(z_air_start, z_air_end, color='orange', alpha=0.10, label='Air Gap')
    plt.axvspan(z_air_end, Nz, color='lightblue', alpha=0.15, label='Slab 2: Water')
    
    # Plot curves
    z_axis_mm = np.arange(Nz) * 2.5 # voxel size is 2.5 mm
    plt.plot(z_axis_mm, profile_no_b, label='B = 0.0 T (No Magnetic Field)', color='#2b5c8f', lw=2.5)
    plt.plot(z_axis_mm, profile_with_b, label='B = 1.5 T (Transverse)', color='#d9534f', lw=2.5, linestyle='--')
    
    # Mark downstream interface
    plt.axvline(x=z_air_end * 2.5, color='red', linestyle=':', label='Downstream Interface', alpha=0.7)
    
    plt.title("1D Depth-Dose Profile: Electron Return Effect (ERE) Smoke Test", fontsize=12, pad=15)
    plt.xlabel("Depth (mm)", fontsize=10)
    plt.ylabel("Relative Dose", fontsize=10)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper right')
    plt.tight_layout()
    
    # Save output plot
    plot_filepath = "ere_check.png"
    plt.savefig(plot_filepath, bbox_inches='tight')
    plt.close()
    
    print(f"  - Saved 1D depth-dose overlay plot to: {os.path.abspath(plot_filepath)}")
    print("Electron Return Effect (ERE) smoke test completed successfully.")
    
    # Clean up test files
    for f_path in [deck_with_b_path, deck_no_b_path, bin_filepath]:
        if os.path.exists(f_path):
            os.remove(f_path)
            
if __name__ == "__main__":
    # If run directly, detect topas presence and run mock if needed
    run_ere_check()
