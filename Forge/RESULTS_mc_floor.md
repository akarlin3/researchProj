# Forge Monte Carlo Benchmark Results

This document records the measured workstation specs, CPU Monte Carlo (MC) timing results, and timeline extrapolation to decide if 10k data generation fits before the December 2026 deadline.

## Workstation Specifications
- **Operating System**: macOS (Darwin arm64)
- **CPU Cores (`nproc`)**: 10 (Apple M4)
- **Total RAM**: 24.0 GB
- **Monte Carlo Engine Detected**: TOPAS MC (Installed)

---

## Benchmark Metrics
The benchmark ran $K = 12$ randomized voxelized phantom cases.

- **Median Core-Hours per Case**: `0.0465`
- **IQR Core-Hours per Case**: `0.0044`
- **Extrapolated Wall-Clock Weeks for 10k Cases**: `0.28` weeks
  *(Extrapolated using formula: `weeks_for_10k = (10000 * core_hours_per_case) / (nproc * 168)`)*

---

## Runway Verdict
Given a timeline of June 8, 2026, to December 31, 2026 (~29 weeks), and accounting for ~6.5 person-weeks of surrounding pipeline/model-training work, the maximum available run window is **22.5 weeks**.

**Verdict**:
YES: The extrapolated run duration of 0.28 weeks clears the 22.5-week runway (29.0 weeks total minus 6.5 weeks of surrounding work).

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
