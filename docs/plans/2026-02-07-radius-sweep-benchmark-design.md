# Radius Sweep Benchmark Design

**Date:** 2026-02-07
**Status:** Approved

## Problem Statement

The CITEgeist Module 3 deconvolution uses a spatial neighborhood radius parameter that was calibrated for simulated data (spacing ~1.6 units) but is being applied to Xenium pseudo-Visium data (spacing 100µm). The current default `radius=4.0` gives ~18 neighbors on simulated data but 0 neighbors on Xenium data.

This design:
1. Fixes the default radius for Xenium benchmarks
2. Creates a parameter sweep to determine optimal neighborhood size

## Coordinate System Analysis

**Simulated data:** Nearest neighbor distance ~1.6 units
- radius=4 → ~18 neighbors (2 rings)

**Xenium pseudo-Visium data:** Nearest neighbor distance = 100µm
- radius=4 → 0 neighbors (broken)

### Correct Radius-to-Ring Mapping for Xenium

| Rings | Theoretical Neighbors | Radius | Actual Median Neighbors |
|-------|----------------------|--------|------------------------|
| 0 | 0 (center only) | 50 | 0 |
| 1 | 6 | 105 | 6 |
| 2 | 18 | 205 | 18 |
| 3 | 36 | 305 | 36 |

## Scope

### 1. Fix Current Benchmark

Update default radius from 4.0 to 205.0 (2 rings, ~18 neighbors) to match original intent.

**Files to modify:**
- `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`
  - `run_manual_benchmark()` signature: `radius: float = 4.0` → `radius: float = 205.0`
  - `run_hierarchical_benchmark()` signature: `radius: float = 4.0` → `radius: float = 205.0`
  - argparse default: `--radius` default 4.0 → 205.0

### 2. Radius Sweep Experiment

Test 0, 1, 2, 3 ring neighborhoods to find optimal spatial context.

**Parameters:**
- Radii: [50, 105, 205, 305]
- Regions: 0-4 (5 total)
- Mode: Manual (achievable-7 profiles)
- Total runs: 4 radii × 5 regions = 20

**Metrics:**
- JSD (Jensen-Shannon Divergence) - lower is better
- Pearson correlation - higher is better
- RMSE - lower is better

## Implementation

### New Files

#### `run_radius_sweep.py`

```python
RADII = [50, 105, 205, 305]
RING_NAMES = {50: "0_rings", 105: "1_ring", 205: "2_rings", 305: "3_rings"}

# For each region (passed via --region-id):
for radius in RADII:
    output_dir = f"output/radius_sweep/radius_{radius}"
    run_manual_benchmark(
        region_id=region_id,
        radius=radius,
        ...
    )
```

#### `run_radius_sweep.sh`

SLURM array job: `--array=0-4` (one task per region)
Each task runs all 4 radii sequentially.

#### `analyze_radius_sweep.py`

Aggregates results across radii and generates:
1. Per-region CSV with metrics across radii
2. Summary visualization (line plots with error bars)
3. Recommendation output

### Output Structure

```
output/radius_sweep/
├── radius_50/      # 0 rings
│   ├── Xenium_region_0/
│   ├── Xenium_region_1/
│   ├── Xenium_region_2/
│   ├── Xenium_region_3/
│   └── Xenium_region_4/
├── radius_105/     # 1 ring
│   └── ...
├── radius_205/     # 2 rings
│   └── ...
└── radius_305/     # 3 rings
    └── ...
```

## Evaluation

Reuses existing evaluation framework:
- `evaluate_benchmark.py` computes JSD, Pearson, RMSE against ground truth
- New `analyze_radius_sweep.py` aggregates across radii

### Expected Output

**Per-region CSV:**
```
region,radius,rings,jsd_mean,pearson_mean,rmse_mean
0,50,0,0.42,0.65,0.18
0,105,1,0.35,0.72,0.15
...
```

**Summary plots:** Line plots showing each metric vs ring count

## Success Criteria

1. All 20 benchmark runs complete successfully
2. Clear visualization comparing metrics across radii
3. Identification of optimal radius (or documented trade-offs)
