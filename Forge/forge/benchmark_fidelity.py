#!/usr/bin/env python3
"""
Forge Monte Carlo Training-Fidelity Benchmark Orchestrator.
Runs N=8 cases from the randomized geometry distribution in parallel,
measures timing and achieved statistical uncertainty,
scales to 1% high-dose region uncertainty (clinical/training grade),
runs the ERE physics gate check, and writes RESULTS_mc_floor_fidelity.md.
"""

import os
import sys
import json
import time
import shutil
import platform
import numpy as np
import resource
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from forge.geom import generate_case_deck
from forge.run_case import run_case, is_topas_available

def run_case_wrapper(case_info, output_dir, histories):
    case_id = case_info["case_id"]
    seed = case_info["seed"]
    particle = case_info["source"]["particle"]
    energy = case_info["source"]["energy_mev"]
    
    # Re-generate the parameter deck with custom histories and Mean+StdDev reporting
    case_meta = generate_case_deck(case_id, seed, output_dir, histories=histories)
    deck_path = os.path.join(output_dir, f"case_{case_id}.txt")
    dims = case_meta["dimensions"]
    voxel_mm = case_meta["voxel_sizes_mm"]
    field_half_cm = case_meta["source"].get("field_half_cm")
    
    try:
        t0 = time.perf_counter()
        res = run_case(deck_path, output_dir)
        t1 = time.perf_counter()
        
        wall_time_s = res["wall_time_s"]
        cpu_time_s = res["cpu_time_s"]
        threads = res["threads"]
        
        dose_grid = res["dose_grid"]
        std_grid = res["std_grid"]
        
        if std_grid is None:
            raise ValueError(f"Standard deviation grid is missing for case {case_id}.")
            
        max_dose = np.max(dose_grid)
        if max_dose <= 0:
            raise ValueError(f"Max dose is zero or negative for case {case_id}.")
            
        threshold = 0.5 * max_dose
        high_dose_idx = np.where(dose_grid >= threshold)
        
        dose_hd = dose_grid[high_dose_idx]
        std_hd = std_grid[high_dose_idx]
        
        # Calculate relative uncertainty of the mean (SEM / Mean)
        # SEM = std_hd / sqrt(histories)
        rel_unc = std_hd / (dose_hd * np.sqrt(histories))
        
        median_unc = np.median(rel_unc)
        max_unc = np.max(rel_unc)
        p90_unc = np.percentile(rel_unc, 90)
        
        # Target histories for 1% (0.01) uncertainty using median high-dose uncertainty
        target_histories = int(np.ceil(histories * (median_unc / 0.01)**2))
        # Ensure we don't end up with 0 or negative
        target_histories = max(target_histories, 1000)
        
        # Scale core hours (CPU time) to target histories
        # topas runs single-threaded, so core-hours is CPU time in hours
        measured_core_hours_cpu = (cpu_time_s / 3600.0)
        projected_core_hours_cpu = measured_core_hours_cpu * (target_histories / histories)
        
        # Also compute wall-time based core hours
        measured_core_hours_wall = (wall_time_s / 3600.0) * threads
        projected_core_hours_wall = measured_core_hours_wall * (target_histories / histories)
        
        return {
            "case_id": case_id,
            "success": True,
            "particle": particle,
            "energy_mev": energy,
            "grid_dims": dims,
            "voxel_mm": voxel_mm,
            "field_half_cm": field_half_cm,
            "primaries_run": histories,
            "achieved_median_sigma": median_unc,
            "achieved_p90_sigma": p90_unc,
            "achieved_max_sigma": max_unc,
            "target_primaries": target_histories,
            "measured_wall_s": wall_time_s,
            "measured_cpu_s": cpu_time_s,
            "measured_ch_cpu": measured_core_hours_cpu,
            "projected_ch_cpu": projected_core_hours_cpu,
            "projected_ch_wall": projected_core_hours_wall
        }
        
    except Exception as e:
        return {
            "case_id": case_id,
            "success": False,
            "error": str(e)
        }

def run_ere_fidelity(histories=500000):
    print("Starting ERE physics gate check (500k histories)...")
    case_id = 999
    seed = 4242
    output_dir = "cases"
    
    Nx, Ny, Nz = 60, 60, 60
    z_air_start = 20
    z_air_end = 30
    
    # Generate the phantom
    phantom = np.full((Nx, Ny, Nz), 2, dtype=np.int16)
    phantom[:, :, z_air_start:z_air_end] = 0 # Air gap
    bin_filename = f"case_{case_id}.bin"
    bin_filepath = os.path.join(output_dir, bin_filename)
    phantom.tofile(bin_filepath)
    
    # Write case deck
    generate_case_deck(case_id, seed, output_dir, histories=histories)
    
    deck_with_b_path = os.path.join(output_dir, f"case_{case_id}_with_b.txt")
    deck_no_b_path = os.path.join(output_dir, f"case_{case_id}_no_b.txt")
    
    with open(os.path.join(output_dir, f"case_{case_id}.txt"), "r") as f:
        deck_template = f.read()
        
    # With B-field (1.5T transverse)
    deck_with_b = deck_template.replace('s:Sc/Dose/OutputFile = "dose_case_999"', 's:Sc/Dose/OutputFile = "dose_case_999_with_b"')
    with open(deck_with_b_path, "w") as f:
        f.write(deck_with_b)
        
    # Without B-field (0.0T)
    deck_no_b = deck_template.replace("d:Ge/Patient/MagneticFieldStrength = 1.5 tesla", "d:Ge/Patient/MagneticFieldStrength = 0.0 tesla")
    deck_no_b = deck_no_b.replace('s:Sc/Dose/OutputFile = "dose_case_999"', 's:Sc/Dose/OutputFile = "dose_case_999_no_b"')
    with open(deck_no_b_path, "w") as f:
        f.write(deck_no_b)
        
    if os.path.exists(os.path.join(output_dir, f"case_{case_id}.txt")):
        os.remove(os.path.join(output_dir, f"case_{case_id}.txt"))
        
    print("  - Running simulation WITH magnetic field...")
    res_with = run_case(deck_with_b_path, output_dir)
    dose_with_b = res_with["dose_grid"]
    
    print("  - Running simulation WITHOUT magnetic field...")
    res_no = run_case(deck_no_b_path, output_dir)
    dose_no_b = res_no["dose_grid"]
    
    # Clean up test files
    for f_path in [deck_with_b_path, deck_no_b_path, bin_filepath]:
        if os.path.exists(f_path):
            os.remove(f_path)
            
    profile_with_b = dose_with_b[Nx//2, Ny//2, :]
    profile_no_b = dose_no_b[Nx//2, Ny//2, :]

    # ------------------------------------------------------------------
    # Reference ERE magnitude (CITED, not fabricated):
    # Geant4 Monte Carlo of a 1.5 T transverse field reports an exit-dose
    # enhancement at the air-tissue interface due to the Electron Return
    # Effect of ~15.4% (10x10 cm field) to ~17.9% (22x22 cm field).
    # (Magnetic-field effect on dose in low-density regions, Geant4, 1.5 T.)
    # NOTE: this is a literature magnitude, NOT an in-repo validated reference
    # dose profile. Per the gate spec, PASS requires comparison against an
    # in-repo reference; absent one the gate is reported UNCALIBRATED.
    # ------------------------------------------------------------------
    REFERENCE_ERE_PCT = (15.4, 17.9)   # (10x10, 22x22) Geant4, 1.5 T transverse
    TOLERANCE_PCT = 3.0                 # absolute dose-difference tolerance (pp)
    HAS_INREPO_REFERENCE = False        # no commissioned reference profile in repo

    # Scan a small window around the distal (downstream) interface for the peak
    # enhancement rather than a single voxel (the ERE peak is localized).
    win = slice(max(z_air_end - 1, 0), min(z_air_end + 4, Nz))
    max_dose_no_b = profile_no_b.max()
    local_diff = profile_with_b[win] - profile_no_b[win]
    peak_diff = float(np.max(local_diff)) if local_diff.size else 0.0
    interface_slice = z_air_end
    dose_diff = float(profile_with_b[interface_slice] - profile_no_b[interface_slice])
    rel_diff_pct = (dose_diff / max_dose_no_b) * 100 if max_dose_no_b > 0 else 0.0
    peak_rel_diff_pct = (peak_diff / max_dose_no_b) * 100 if max_dose_no_b > 0 else 0.0

    # Agreement vs the cited reference band (using the peak enhancement).
    ref_lo, ref_hi = REFERENCE_ERE_PCT
    within_ref = (ref_lo - TOLERANCE_PCT) <= peak_rel_diff_pct <= (ref_hi + TOLERANCE_PCT)

    if not HAS_INREPO_REFERENCE:
        # No in-repo commissioned reference -> cannot certify PASS.
        verdict = "UNCALIBRATED"
        is_pass = False
    elif within_ref:
        verdict = "PASS"
        is_pass = True
    else:
        verdict = "FAIL"
        is_pass = False
    
    # Generate overlay plot for ERE
    import matplotlib.pyplot as plt
    plt.figure(figsize=(9, 5), dpi=150)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    plt.axvspan(0, z_air_start * 2.5, color='lightblue', alpha=0.15, label='Slab 1: Water')
    plt.axvspan(z_air_start * 2.5, z_air_end * 2.5, color='orange', alpha=0.10, label='Air Gap')
    plt.axvspan(z_air_end * 2.5, Nz * 2.5, color='lightblue', alpha=0.15, label='Slab 2: Water')
    
    z_axis_mm = np.arange(Nz) * 2.5
    plt.plot(z_axis_mm, profile_no_b, label='B = 0.0 T (No Magnetic Field)', color='#2b5c8f', lw=2.5)
    plt.plot(z_axis_mm, profile_with_b, label='B = 1.5 T (Transverse)', color='#d9534f', lw=2.5, linestyle='--')
    plt.axvline(x=z_air_end * 2.5, color='red', linestyle=':', label='Downstream Interface', alpha=0.7)
    
    plt.title(f"1D Depth-Dose Profile: Electron Return Effect (ERE) Physics Gate ({histories:,} histories)", fontsize=11, pad=15)
    plt.xlabel("Depth (mm)", fontsize=10)
    plt.ylabel("Relative Dose", fontsize=10)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper right')
    plt.tight_layout()
    plt.savefig("ere_check.png", bbox_inches='tight')
    plt.close()
    
    print(f"ERE Check complete: interface diff = {dose_diff:+.4e} ({rel_diff_pct:+.3f}% of max), "
          f"peak enhancement = {peak_rel_diff_pct:+.3f}% of max, "
          f"reference band = {REFERENCE_ERE_PCT} +/- {TOLERANCE_PCT}pp, Verdict = {verdict}")
    return {
        "is_pass": is_pass,
        "dose_no_b": float(profile_no_b[interface_slice]),
        "dose_with_b": float(profile_with_b[interface_slice]),
        "dose_diff": float(dose_diff),
        "rel_diff_pct": float(rel_diff_pct),
        "peak_rel_diff_pct": float(peak_rel_diff_pct),
        "reference_ere_pct": REFERENCE_ERE_PCT,
        "tolerance_pct": TOLERANCE_PCT,
        "has_inrepo_reference": HAS_INREPO_REFERENCE,
        "within_reference": bool(within_ref),
        "verdict": verdict
    }

def main():
    if not is_topas_available():
        print("Error: TOPAS MC is not available. Please make sure it is installed.")
        sys.exit(1)
        
    manifest_path = "cases/manifest.json"
    if not os.path.exists(manifest_path):
        print("Error: manifest.json not found. Run geom.py first.")
        sys.exit(1)
        
    with open(manifest_path, "r") as f:
        manifest = json.load(f)
        
    cases = manifest["cases"]
    # We run N=8 representative cases
    N = 8
    target_cases = cases[:N]
    
    histories_per_case = 500000 # 500k histories
    
    print("=" * 70)
    print(f"FORGE TRAINING-FIDELITY BENCHMARK (N={N} cases, {histories_per_case:,} histories/case)")
    print("=" * 70)
    
    results = []
    nproc = os.cpu_count() or 1
    max_workers = min(N, nproc)
    
    print(f"Running {N} cases in parallel using {max_workers} processes...")
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_case_wrapper, c, "cases", histories_per_case): c for c in target_cases}
        for fut in as_completed(futures):
            res = fut.result()
            if res["success"]:
                results.append(res)
                print(f"  - Case {res['case_id']} ({res['particle']}): Achieved high-dose 1σ uncertainty = {res['achieved_median_sigma']*100:.3f}% "
                      f"-> Target histories for 1%: {res['target_primaries']:,} (Projected core-hours: {res['projected_ch_cpu']:.4f})")
            else:
                print(f"  - Case {futures[fut]['case_id']} failed: {res['error']}")
                
    elapsed = time.perf_counter() - t0
    print(f"Benchmark run finished in {elapsed:.1f}s.")
    
    if len(results) < N:
        print(f"Warning: Only {len(results)} out of {N} cases completed successfully.")
        
    # Run ERE Check
    ere_res = run_ere_fidelity(histories=histories_per_case)
    
    # Calculate statistics on projected core hours
    projected_ch = [r["projected_ch_cpu"] for r in results]
    median_ch = np.median(projected_ch)
    q75, q25 = np.percentile(projected_ch, [75, 25])
    iqr_ch = q75 - q25
    min_ch = np.min(projected_ch)
    max_ch = np.max(projected_ch)
    
    # Extrapolate for 10k cases on 10 cores
    # weeks = (10000 * ch_per_case) / (10 * 168)
    weeks_median = (10000 * median_ch) / (10 * 168)
    weeks_min = (10000 * min_ch) / (10 * 168)
    weeks_max = (10000 * max_ch) / (10 * 168)
    
    # Compare with the two prior figures: 0.0465 (smoke test) and 0.326
    # (proton-contaminated training-fidelity median).
    smoke_ch = 0.0465
    proton_contaminated_ch = 0.326
    mult_vs_smoke = median_ch / smoke_ch
    ratio_vs_proton = median_ch / proton_contaminated_ch

    total_core_hours = 10000 * median_ch
    total_hours_10cores = total_core_hours / 10
    if total_hours_10cores <= 48:
        one_liner_verdict = "About a weekend"
    elif weeks_median <= 3:
        one_liner_verdict = "A few weeks"
    else:
        one_liner_verdict = "Longer than a few weeks (cluster-scale)"

    # Voxel grid (uniform across the randomized distribution)
    gd = results[0]["grid_dims"]
    vm = results[0]["voxel_mm"]
    grid_extent_mm = [gd[i] * vm[i] for i in range(3)]

    # ERE gate text
    ref_lo, ref_hi = ere_res["reference_ere_pct"]
    ere_gate_line = (
        f"{ere_res['verdict']} (peak enhancement {ere_res['peak_rel_diff_pct']:+.3f}% of max "
        f"vs reference band {ref_lo:.1f}-{ref_hi:.1f}% +/- {ere_res['tolerance_pct']:.1f}pp)"
    )

    # Write RESULTS_mc_floor_photon.md
    md_content = f"""# Forge Monte Carlo Training-Fidelity Benchmark Report (PHOTON)

Per-case simulation costs at **training-grade fidelity** (1% 1-sigma uncertainty in the
high-dose region) for the corrected **Elekta Unity-class 7 MV FFF photon** beam at 1.5 T,
plus the **Electron Return Effect (ERE)** physics gate against a cited reference magnitude.

## Workstation Specifications
- **Operating System**: {platform.system()} ({platform.machine()})
- **CPU Cores (`os.cpu_count`)**: {nproc}
- **Monte Carlo Engine**: OpenTOPAS v4.2.3 + Geant4 11.3.2

## Voxel Grid (REPORTED)
- **Dimensions**: {gd[0]} x {gd[1]} x {gd[2]} voxels
- **Resolution**: {vm[0]} x {vm[1]} x {vm[2]} mm
- **Physical extent**: {grid_extent_mm[0]:.0f} x {grid_extent_mm[1]:.0f} x {grid_extent_mm[2]:.0f} mm
- Note: 2.5 mm is coarse for penumbra/ERE structure; 1% sigma is comparatively cheap on
  this grid. The grid was NOT coarsened for speed; it is the documented dataset grid.

---

## 1. Re-Timing Benchmark Results (N = {len(results)} Photon Cases)
Cases run at {histories_per_case:,} histories; achieved 1-sigma uncertainty in the
high-dose region (>=50% isodose) scaled to the 1% target. Core-hours are true CPU time.

| Case | Particle | Field (cm) | Grid | Primaries | Achieved sigma (med) | sigma (max) | Target prim. for 1% | CPU (s) | Target core-h |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
"""
    for r in results:
        field_cm = (r.get("field_half_cm") or 0) * 2.0
        grid_s = f"{r['grid_dims'][0]}x{r['grid_dims'][1]}x{r['grid_dims'][2]}@{r['voxel_mm'][0]}mm"
        md_content += (f"| {r['case_id']} | {r['particle']} | {field_cm:.1f} | {grid_s} | {r['primaries_run']:,} "
                       f"| {r['achieved_median_sigma']*100:.3f}% | {r['achieved_max_sigma']*100:.3f}% | {r['target_primaries']:,} "
                       f"| {r['measured_cpu_s']:.2f} | `{r['projected_ch_cpu']:.4f}` |\n")

    md_content += f"""
### Core-Hours per Case Distribution
- **Median**: `{median_ch:.4f}` core-hours/case
- **IQR**: `{iqr_ch:.4f}` core-hours/case
- **Min / Max**: `{min_ch:.4f}` / `{max_ch:.4f}` core-hours/case

---

## 2. Electron Return Effect (ERE) Gate vs Reference (Photon, 1.5 T)
Water-air-water interface, B-field 1.5 T transverse on vs off, {histories_per_case:,} histories.
Peak distal-interface enhancement compared to a **cited Geant4 reference band**.

- **Reference ERE magnitude**: {ref_lo:.1f}% (10x10 cm) - {ref_hi:.1f}% (22x22 cm) exit-dose
  enhancement at 1.5 T transverse (Geant4 MC, magnetic-field effect in low-density regions).
  This is a *literature magnitude*, not an in-repo validated reference profile.
- **In-repo reference profile available**: {ere_res['has_inrepo_reference']}
- **Tolerance**: +/- {ere_res['tolerance_pct']:.1f} percentage points
- **Measured peak enhancement**: `{ere_res['peak_rel_diff_pct']:+.3f}%` of max dose
- **Single-voxel interface diff**: `{ere_res['rel_diff_pct']:+.3f}%` of max dose
- **ERE Gate Verdict**: **{ere_gate_line}**

Because there is no in-repo commissioned reference profile, the gate is reported
**UNCALIBRATED** rather than PASS even if the magnitude lands in-band.

---

## 3. Extrapolation to 10k (Photon)
`weeks = (10000 * core_h_per_case) / (10 * 168)`

- **Median**: `{weeks_median:.2f}` weeks  |  **Min-Max**: `{weeks_min:.2f}` - `{weeks_max:.2f}` weeks

### Reconciliation with prior figures
- **0.0465** core-h/case (smoke test, 50k histories, no sigma analysis) -> photon is `{mult_vs_smoke:.1f}x` this.
- **0.326** core-h/case (prior "training-fidelity" median) was **proton-contaminated**: 7/8
  timed cases were 120-180 MeV protons (~1100 s CPU each) which do not belong in a photon
  MR-Linac dataset. The corrected photon median is `{median_ch:.4f}` core-h/case
  (`{ratio_vs_proton:.2f}x` the contaminated figure).

### One-Line Verdict
**{one_liner_verdict}** - corrected photon estimate is **{weeks_median:.2f} weeks** (median) for the 10k set.
"""

    with open("RESULTS_mc_floor_photon.md", "w") as f:
        f.write(md_content)
    print("Wrote results to RESULTS_mc_floor_photon.md")

if __name__ == "__main__":
    main()
