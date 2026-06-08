#!/usr/bin/env python3
"""
Forge Single-Case Monte Carlo Runner.
Wrapper that executes a single TOPAS deck, captures execution timings,
and parses the resulting 3D binned dose grid.
"""

import os
import time
import subprocess
import shutil
import numpy as np
import resource

def is_topas_available():
    """Checks if topas executable is available on the system path."""
    return shutil.which("topas") is not None

def parse_topas_csv(csv_path):
    """
    Parses a binned CSV output from TOPAS.
    The file format consists of comment header lines starting with '#' 
    followed by the voxel index and binned value columns:
    BinX, BinY, BinZ, Value, ...
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"TOPAS output file not found: {csv_path}")
        
    dims = [60, 60, 60] # Default fallback
    
    # Try to extract actual dimensions from the comment header
    with open(csv_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                line_lower = line.lower()
                # TOPAS header example: "# X in 60 bins of 2.5 mm"
                if "bins" in line_lower:
                    parts = line_lower.replace(":", " ").replace(",", " ").split()
                    for idx, char in enumerate(["x", "y", "z"]):
                        if char in parts:
                            try:
                                # Find the first integer after 'x'/'y'/'z' or before 'bins'
                                c_idx = parts.index(char)
                                for p in parts[c_idx:]:
                                    if p.isdigit():
                                        dims[idx] = int(p)
                                        break
                            except Exception:
                                pass
            else:
                break
                
    try:
        # Load CSV using numpy, ignoring header lines starting with '#'
        data = np.loadtxt(csv_path, delimiter=',', comments='#')
        
        # If TOPAS output has indices and values (X, Y, Z, Value, ...)
        if len(data.shape) == 2 and data.shape[1] >= 4:
            ix = data[:, 0].astype(int)
            iy = data[:, 1].astype(int)
            iz = data[:, 2].astype(int)
            val = data[:, 3]
            
            Nx = max(ix.max() + 1, dims[0])
            Ny = max(iy.max() + 1, dims[1])
            Nz = max(iz.max() + 1, dims[2])
            
            dose_grid = np.zeros((Nx, Ny, Nz), dtype=np.float32)
            dose_grid[ix, iy, iz] = val
            return dose_grid
            
        elif len(data.shape) == 1:
            # Flat list of values - reshape using the detected/fallback dimensions
            return data.reshape(dims).astype(np.float32)
        else:
            raise ValueError(f"Unexpected data shape: {data.shape}")
            
    except Exception as e:
        raise ValueError(f"Failed parsing TOPAS CSV output at '{csv_path}': {e}")

def run_case(deck_path, output_dir="cases"):
    """
    Executes the TOPAS Monte Carlo engine on the specified deck file, 
    measuring performance timings and reading the 3D dose grid.
    """
    if not is_topas_available():
        raise FileNotFoundError(
            "TOPAS Monte Carlo engine was not found on your system PATH.\n"
            "Please follow the installation instructions in RESULTS_mc_floor.md to set it up."
        )
        
    deck_name = os.path.basename(deck_path)
    case_name = os.path.splitext(deck_name)[0]
    
    # Command to run TOPAS
    cmd = ["topas", deck_name]
    
    t0_wall = time.perf_counter()
    usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0_cpu_children = usage_start.ru_utime + usage_start.ru_stime
    
    # Run the simulation
    # TOPAS output files are written to the current working directory (output_dir)
    result = subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True)
    
    t1_wall = time.perf_counter()
    usage_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    t1_cpu_children = usage_end.ru_utime + usage_end.ru_stime
    
    if result.returncode != 0:
        raise RuntimeError(
            f"TOPAS simulation failed with return code {result.returncode}.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        
    wall_time = t1_wall - t0_wall
    cpu_time = t1_cpu_children - t0_cpu_children
    
    # Locate dose output CSV
    csv_filename = f"dose_{case_name}.csv"
    csv_filepath = os.path.join(output_dir, csv_filename)
    
    dose_grid = parse_topas_csv(csv_filepath)
    
    # Extract number of threads from the deck if possible, default to 1
    threads = 1
    try:
        with open(deck_path, "r") as f:
            for line in f:
                if "NumberOfThreads" in line and not line.strip().startswith("#"):
                    parts = line.split("=")
                    threads = int(parts[-1].strip())
    except Exception:
        pass
        
    return {
        "dose_grid": dose_grid,
        "wall_time_s": wall_time,
        "cpu_time_s": cpu_time,
        "threads": threads
    }

if __name__ == "__main__":
    # Self-test checks if TOPAS is present and alerts the user
    if is_topas_available():
        print("TOPAS Monte Carlo engine is available on PATH.")
    else:
        print("TOPAS Monte Carlo engine is NOT available on PATH. (This is expected for this workspace)")
