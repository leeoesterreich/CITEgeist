# Cellpose-Compatible Synthetic Images for Simulation Benchmarking

**Date**: 2026-02-12
**Status**: Approved
**Branch**: `dev`

## Overview

Generate synthetic morphology images from scCube simulation data (high_seg and mixed conditions) to enable discrete cell assignment benchmarking using Cellpose segmentation.

## Motivation

The discrete cell assignment pipeline (IQP solver + Cellpose nuclei counts) is currently benchmarked only on Xenium pseudo-Visium data. To validate performance on controlled simulations with known ground truth, we need Cellpose-compatible images generated from scCube single-cell positions.

## Design

### Dual-Mode Image Generation

Two image styles will be generated to match different Cellpose models:

| Mode | Background | Nuclei | Cytoplasm | Cellpose Model |
|------|------------|--------|-----------|----------------|
| `dapi` | Black (0,0,0) | White/gray Gaussian (σ=3.5) | None | `nuclei` |
| `h_and_e` | White (250,250,250) | Purple Gaussian (140,90,130) | Pink rings (240,220,220) | `cyto2` |

**Color rationale** (derived from actual data):
- DAPI mode matches Xenium morphology images (grayscale, R=G=B)
- H&E mode matches Visium tissue images (hematoxylin purple nuclei, eosin pink cytoplasm)

### Image Specifications

- **Size**: 8000×8000 pixels (default)
- **Pixel density**: ~12,000 pixels/cell (matches Xenium hires crops)
- **Nucleus rendering**: Gaussian blob with σ=3.5 pixels (~20-30px effective diameter)
- **Cell positions**: Original scCube coordinates preserved (no jitter)

### Data Sources

**Input** (Brent's scCube simulations):
```
/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/scCube_12k/replicates/
├── high_seg/ST_sim/Wu_ST_{0-4}_index.csv
└── mixed/ST_sim/Wu_ST_{0-4}_index.csv
```

Each `index.csv` contains:
- `Cell`: Cell ID
- `Cell_type`: Cell type label
- `point_x`, `point_y`: Cell coordinates (0-50 range)
- `spot`: Assigned spot ID

**Simulation characteristics**:
- 5,000 cells per replicate
- 951 spots per replicate
- ~5.3 cells per spot (sparser than Xenium's ~17/spot)
- Square coordinate space (50×50 units)

### Output Structure

```
Benchmarking/simulation_benchmarking/CITEgeist/{condition}/
├── images/
│   ├── dapi/
│   │   ├── Wu_rep_0_cellpose.png
│   │   └── ...
│   └── h_and_e/
│       ├── Wu_rep_0_cellpose.png
│       └── ...
└── nuclei_counts/
    ├── Wu_rep_0_nuclei_counts.csv
    └── ...
```

**Nuclei counts derivation**:
```python
df = pd.read_csv(f"Wu_ST_{id}_index.csv")
nuclei_counts = df.groupby('spot').size()
```

## Implementation

### New Files

| File | Purpose |
|------|---------|
| `simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py` | Generate synthetic images + nuclei counts |
| `simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py` | Run discrete benchmark on simulation data |
| `simulation_benchmarking/CITEgeist/slurm/sbatch_generate_images.sh` | SLURM job for image generation |
| `simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_simulation.sh` | SLURM job for benchmark pipeline |

### Script: generate_cellpose_images.py

**Arguments**:
```
--replicate-id INT      Replicate number (0-4)
--condition STR         high_seg | mixed
--mode STR              dapi | h_and_e | both (default: both)
--output-dir PATH       Output directory
--image-size INT        Image size in pixels (default: 8000)
--nucleus-sigma FLOAT   Gaussian sigma for nuclei (default: 3.5)
--input-dir PATH        Path to scCube replicates directory
```

**Algorithm**:
1. Load `ST_sim/Wu_ST_{id}_index.csv`
2. Extract cell coordinates and scale to image pixels
3. For each cell, draw:
   - DAPI mode: White Gaussian blob at cell position
   - H&E mode: Pink circle (cytoplasm) + purple Gaussian (nucleus)
4. Save image to appropriate mode folder
5. Compute and save ground truth nuclei counts per spot

### Script: benchmark_discrete_simulation.py

**Pipeline**:
1. Load synthetic image from `images/{mode}/`
2. Run Cellpose with appropriate model:
   - `dapi` → `model_type="nuclei"`
   - `h_and_e` → `model_type="cyto2"`
3. Assign detected nuclei to spots
4. Compare predicted vs ground truth nuclei counts
5. Load CITE and GEX h5ad files
6. Run discrete cell assignment (IQP)
7. Evaluate proportions against `Wu_ST_{id}_prop.csv`
8. Run GEX deconvolution and evaluate

### Cellpose Model Selection

The segmentation module (`CITEgeist/model/segmentation.py`) currently hardcodes `model_type="nuclei"`. The benchmark script will need to either:
- Pass model type as parameter to segmentation functions
- Use Cellpose directly with the appropriate model

**Recommended**: Add `model_type` parameter to `run_cellpose_nuclei_segmentation()`.

## SLURM Configuration

**Image generation** (CPU, quick):
- Partition: `htc`
- Time: 30 minutes
- Memory: 16GB

**Benchmark pipeline** (GPU for Cellpose):
- Partition: `gpu` or `a100`
- Time: 2 hours per replicate
- Memory: 32GB
- GPUs: 1

## Success Criteria

1. Images generated for all 5 replicates × 2 conditions × 2 modes = 20 images
2. Ground truth nuclei counts match cell counts from index files
3. Cellpose successfully segments both image types
4. Benchmark metrics computed and comparable to Xenium results
5. Results tracked in version control

## Dependencies

- PIL/Pillow (image generation)
- NumPy (Gaussian rendering)
- pandas (data loading)
- cellpose (segmentation)
- Existing CITEgeist discrete assignment pipeline

## References

- Existing scResolve image generation: `simulation_benchmarking/scResolve/src/generate_synthetic_images.py`
- Xenium discrete benchmark: `xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py`
- scCube data location: Brent's folder (see Data Sources above)
