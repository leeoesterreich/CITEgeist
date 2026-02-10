# CARD Benchmarking Design

**Date**: 2026-02-10
**Status**: Ready for implementation
**Priority**: Before submission (critical path)

## Motivation

CITEgeist claims spatial data incorporation as a key advantage. CARD (Ma & Zhou, Nature Biotechnology 2022) also incorporates spatial correlation in its deconvolution model, making it a relevant comparison point. Additionally, CARD claims reference-free capability which we want to assess.

## Scope

### Datasets

| Dataset | Conditions | Ground Truth |
|---------|------------|--------------|
| Simulation | mixed + high_seg | 10 replicates each |
| Xenium | 5 regions | protein_gt only (10 cell types) |

### Modes

- **Reference mode**: Standard CARD with scRNA-seq reference (same reference as RCTD/Cell2Location)
- **Reference-free mode**: CARD without external reference (tests their claim)

### Total Runs

- Xenium: 5 regions × 2 modes = 10 runs
- Simulation: 10 reps × 2 conditions × 2 modes = 40 runs

## Directory Structure

```
Benchmarking/
├── environments/
│   └── card_env.yml
│
├── simulation_benchmarking/
│   └── CARD/
│       ├── mixed/
│       │   ├── CARD_pipeline_mixed.R
│       │   └── CARD_pipeline_mixed.sh
│       └── high_seg/
│           ├── CARD_pipeline_high_seg.R
│           └── CARD_pipeline_high_seg.sh
│
└── xenium_benchmarking/
    └── CARD/
        ├── src/
        │   ├── convert_h5ad.py          # Reuse from RCTD if needed
        │   └── run_benchmark.R          # Unified CARD wrapper
        ├── slurm/
        │   ├── run_all_regions.sh       # Reference mode array job
        │   └── run_all_regions_reffree.sh  # Reference-free array job
        ├── output_protein_gt/           # Reference mode results
        └── output_protein_gt_reffree/   # Reference-free results
```

## Environment Setup

### card_env.yml

```yaml
name: CARD_env
channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base>=4.0
  - r-devtools
  - r-matrix
  - r-ggplot2
  - r-rcolorbrewer
  - r-optparse
  - r-jsonlite
  - bioconductor-singlecellexperiment
  - bioconductor-summarizedexperiment
```

### Installation

```bash
conda env create -f card_env.yml
conda activate CARD_env
Rscript -e 'devtools::install_github("YMa-lab/CARD")'
```

## Core Script Design

### run_benchmark.R

Unified script supporting both modes via `--mode` flag:

```r
# Command line arguments
--mode            # "reference" or "reference_free"
--region-id       # 0-4 for Xenium, replicate ID for simulation
--ref-counts      # scRNA-seq reference counts (ignored if mode=reference_free)
--ref-cell-types  # scRNA-seq cell type labels
--spatial-counts  # Spatial transcriptomics counts
--spatial-coords  # Spatial coordinates
--output-dir      # Output directory

# Mode-specific workflow
if (mode == "reference") {
  CARD_obj <- createCARDObject(
    sc_count = ref_counts,
    sc_meta = ref_meta,
    spatial_count = spatial_counts,
    spatial_location = coords,
    ...
  )
  CARD_obj <- CARD_deconvolution(CARD_obj)

} else {
  # Reference-free mode
  CARD_obj <- createCARDfreeObject(
    markerList = marker_genes,
    spatial_count = spatial_counts,
    spatial_location = coords,
    ...
  )
  CARD_obj <- CARD_refFree(CARD_obj)
}

# Standardized output format
write.csv(proportions, "{sample}_cell_prop_finetuned_results.csv")
```

## SLURM Job Structure

### Xenium Array Jobs

```bash
#!/bin/bash
#SBATCH --job-name=CARD_xenium
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Separate array jobs for reference and reference_free modes
# Can run in parallel since independent
```

### Simulation Array Jobs

```bash
#!/bin/bash
#SBATCH --job-name=CARD_sim
#SBATCH --array=1-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Separate jobs for:
# - high_seg reference
# - high_seg reference_free
# - mixed reference
# - mixed reference_free
```

### Data Reuse

CSV files generated for RCTD can be reused for CARD - same input format. No regeneration needed.

## Evaluation Integration

### Output Format

Matches existing methods:
```
CARD/output_protein_gt/Xenium_region_0/
  └── Xenium_region_0_cell_prop_finetuned_results.csv

CARD/output_protein_gt_reffree/Xenium_region_0/
  └── Xenium_region_0_cell_prop_finetuned_results.csv
```

### Method Labels

- Reference mode: `CARD`
- Reference-free mode: `CARD (ref-free)`

### Metrics

No changes needed - existing metrics (JSD, RMSE, MAE, Pearson) are method-agnostic.

### Evaluation Script Updates

Add `"CARD"` and `"CARD (ref-free)"` to method lists in:
- `xenium_benchmarking/evaluation/src/compare_*.py`
- `simulation_benchmarking/src/benchmarking_spot_deconv.py`

## Implementation Order

1. Create conda environment + verify CARD installs
2. Write unified `run_benchmark.R`
3. Test on single Xenium region (both modes)
4. Scale to all Xenium regions via array jobs
5. Adapt scripts for simulation benchmarks
6. Run full simulation benchmarks
7. Update evaluation scripts, run evaluation
8. Integrate results into manuscript figures

## References

- CARD GitHub: https://github.com/YMa-lab/CARD
- Ma, Y., Zhou, X. (2022). Spatially informed cell-type deconvolution for spatial transcriptomics. Nature Biotechnology, 40, 1349–1359.
