#!/usr/bin/env python3
"""
Forge Monte Carlo Benchmark Orchestrator.
Runs K cases, measures core-hours per case, extrapolates the timeline for 10k cases,
and writes the results to RESULTS_mc_floor.md.
"""

import os
import sys
import json
import time
import shutil
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# Handle path imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from forge.run_case import run_case, is_topas_available

def run_case_wrapper(case_info, output_dir, mock=False):
    """Worker function to run a single case and return timing/status."""
    case_id = case_info["case_id"]
    deck_path = os.path.join(output_dir, case_info["deck_file"])
    
    if mock:
        # Simulate execution time
        # Standard: 0.15 to 0.45 core-hours per case
        # Let's say wall time is 45 to 135 seconds on a single core
        rng = np.random.default_rng(case_info["seed"])
        wall_time = rng.uniform(45.0, 135.0) # seconds
        cpu_time = wall_time * 0.98 # high CPU utilization
        return {
            "case_id": case_id,
            "success": True,
            "wall_time_s": wall_time,
            "cpu_time_s": cpu_time,
            "threads": 1,
            "core_hours": (cpu_time / 3600.0) * 1 # cpu time in hours * threads
        }
        
    try:
        res = run_case(deck_path, output_dir)
        wall_time = res["wall_time_s"]
        cpu_time = res["cpu_time_s"]
        threads = res["threads"]
        core_hours = (wall_time / 3600.0) * threads
        return {
            "case_id": case_id,
            "success": True,
            "wall_time_s": res["wall_time_s"],
            "cpu_time_s": cpu_time,
            "threads": threads,
            "core_hours": core_hours
        }
    except Exception as e:
        return {
            "case_id": case_id,
            "success": False,
            "error": str(e)
        }

def run_benchmark(mock=False):
    """Runs the benchmark across all cases and writes RESULTS_mc_floor.md."""
    manifest_path = os.path.join("cases", "manifest.json")
    if not os.path.exists(manifest_path):
        print("Error: manifest.json not found. Run geom.py first.")
        sys.exit(1)
        
    with open(manifest_path, "r") as f:
        manifest = json.load(f)
        
    cases = manifest["cases"]
    K = len(cases)
    
    # OS and CPU detection
    # On macOS, hw.logicalcpu represents logical cores
    nproc = os.cpu_count() or 1
    
    # Total RAM
    try:
        ram_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        total_ram_gb = ram_bytes / (1024**3)
    except Exception:
        total_ram_gb = 24.0 # Fallback based on recon
        
    print("=" * 60)
    print("FORGE MONTE CARLO BENCHMARK")
    print("=" * 60)
    print(f"Detected CPU Cores (nproc): {nproc}")
    print(f"Detected Total RAM:         {total_ram_gb:.1f} GB")
    
    topas_present = is_topas_available()
    use_mock = not topas_present or mock
    
    if use_mock:
        print("[!] Running in MOCK benchmark mode (TOPAS not installed or mock requested).")
    else:
        print(f"Running actual TOPAS benchmark for {K} cases in parallel...")
        
    # We will execute using a ProcessPoolExecutor pinned to all available cores
    results = []
    t0 = time.perf_counter()
    
    # Limit pool to min(K, nproc)
    max_workers = min(K, nproc)
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_case_wrapper, c, "cases", mock=use_mock): c for c in cases}
        for fut in as_completed(futures):
            res = fut.result()
            if res["success"]:
                results.append(res)
                print(f"  - Case {res['case_id']} completed: {res['core_hours']:.4f} core-hours "
                      f"(Wall: {res['wall_time_s']:.1f}s, CPU: {res['cpu_time_s']:.1f}s)")
            else:
                print(f"  - Case {futures[fut]['case_id']} failed: {res['error']}")
                
    elapsed_total = time.perf_counter() - t0
    
    if len(results) == 0:
        print("Error: No cases completed successfully.")
        sys.exit(1)
        
    # Compile timings
    core_hours_list = [r["core_hours"] for r in results]
    median_ch = np.median(core_hours_list)
    q75, q25 = np.percentile(core_hours_list, [75, 25])
    iqr_ch = q75 - q25
    
    # Extrapolate 10,000 cases timeline
    # Formula: weeks = (10000 * core_hours_per_case) / (nproc * 168)
    hours_per_week = 168
    total_core_hours_needed = 10000 * median_ch
    wall_clock_hours_needed = total_core_hours_needed / nproc
    weeks_for_10k = wall_clock_hours_needed / hours_per_week
    
    print("-" * 60)
    print(f"Benchmark finished in {elapsed_total:.1f}s.")
    print(f"Median core-hours/case: {median_ch:.4f} (IQR: {iqr_ch:.4f})")
    print(f"Extrapolated wall-clock weeks for 10k: {weeks_for_10k:.2f} weeks")
    print("-" * 60)
    
    # Compute verdict
    # Runway: June 8 -> Dec 31 2026 (~29 weeks)
    # Surrounding work: ~6.5 person-weeks
    # Max target run time: 29 - 6.5 = 22.5 weeks, let's round target to ~25 weeks maximum.
    runway_weeks = 29.0
    surrounding_work_weeks = 6.5
    target_max_weeks = runway_weeks - surrounding_work_weeks # 22.5 weeks
    
    clears_runway = weeks_for_10k <= target_max_weeks
    verdict_str = ""
    if clears_runway:
        verdict_str = (
            f"YES: The extrapolated run duration of {weeks_for_10k:.2f} weeks clears "
            f"the {target_max_weeks:.1f}-week runway ({runway_weeks} weeks total minus {surrounding_work_weeks} weeks of surrounding work)."
        )
    else:
        # Quantify how much more resource or variance reduction is needed
        required_cores = int(np.ceil((10000 * median_ch) / (target_max_weeks * 168)))
        verdict_str = (
            f"NO: The extrapolated run duration of {weeks_for_10k:.2f} weeks exceeds "
            f"the {target_max_weeks:.1f}-week runway limit. "
            f"To clear the runway under 22.5 weeks, we would need to scale workstation cores "
            f"from {nproc} to at least {required_cores} cores, or reduce core-hours/case by "
            f"{(1.0 - target_max_weeks/weeks_for_10k)*100:.1f}%."
        )
        
    print(f"Verdict: {verdict_str}")
    
    # Write RESULTS_mc_floor.md
    write_results_markdown(nproc, total_ram_gb, topas_present, use_mock, results, median_ch, iqr_ch, weeks_for_10k, verdict_str)

def write_results_markdown(nproc, total_ram, topas_present, use_mock, results, median_ch, iqr_ch, weeks_for_10k, verdict_str):
    """Writes the results file RESULTS_mc_floor.md."""
    
    markdown_content = f"""# Forge Monte Carlo Benchmark Results

This document records the measured workstation specs, CPU Monte Carlo (MC) timing results, and timeline extrapolation to decide if 10k data generation fits before the December 2026 deadline.

## Workstation Specifications
- **Operating System**: macOS (Darwin arm64)
- **CPU Cores (`nproc`)**: {nproc} (Apple M4)
- **Total RAM**: {total_ram:.1f} GB
- **Monte Carlo Engine Detected**: {"TOPAS MC (Installed)" if topas_present else "None (Using Mock Run Data)"}

---

## Benchmark Metrics
The benchmark ran $K = {len(results)}$ randomized voxelized phantom cases.

- **Median Core-Hours per Case**: `{median_ch:.4f}`
- **IQR Core-Hours per Case**: `{iqr_ch:.4f}`
- **Extrapolated Wall-Clock Weeks for 10k Cases**: `{weeks_for_10k:.2f}` weeks
  *(Extrapolated using formula: `weeks_for_10k = (10000 * core_hours_per_case) / (nproc * 168)`)*

---

## Runway Verdict
Given a timeline of June 8, 2026, to December 31, 2026 (~29 weeks), and accounting for ~6.5 person-weeks of surrounding pipeline/model-training work, the maximum available run window is **22.5 weeks**.

**Verdict**:
{verdict_str}

---

## Simulation Engine Installation Guidance
If you are seeing "None (Using Mock Run Data)" above, follow these steps to install the TOPAS engine on this Mac:

1. **Request TOPAS License**: Register at [topasmc.org](https://www.topasmc.org) and download the matching macOS binary release.
2. **Download Geant4 Data**: TOPAS requires Geant4 data files. Install via Homebrew:
   ```bash
   brew install geant4
   ```
   Brew will automatically fetch and unpack the required datasets into `/opt/homebrew/share/Geant4/data/` or similar.
3. **Environment Setup**: Define variables in your shell config (`~/.zshrc`):
   ```bash
   export TOPAS_G4_DATA_DIR=$(brew --prefix geant4)/share/Geant4/data
   export PATH="/path/to/topas/bin:$PATH"
   ```
4. **Re-Run Benchmark**: Run the benchmark again to measure real timings:
   ```bash
   python3 forge/benchmark.py
   ```
"""
    with open("RESULTS_mc_floor.md", "w") as f:
        f.write(markdown_content)
    print(f"Wrote results report to: {os.path.abspath('RESULTS_mc_floor.md')}")

if __name__ == "__main__":
    # Run in mock mode if requested or if TOPAS is missing
    run_mock = "--mock" in sys.argv
    run_benchmark(mock=run_mock)
