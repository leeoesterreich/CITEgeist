# Xenium Benchmarking Framework

This directory contains benchmarking methods for evaluating deconvolution algorithms on **real spatial multi-omic data** from 10x Xenium. Pseudo-Visium data and ground truth are generated in `../xenium_pseudovisium/`.

## Overview

The framework benchmarks deconvolution methods against ground truth cell proportions derived from Xenium single-cell resolution data. Each method has its own subdirectory with consistent structure.

## Directory Structure

```
xenium_benchmarking/
├── CITEgeist/                        # CITEgeist method (implemented)
│   ├── src/
│   │   ├── run_benchmark.py         # CITEgeist wrapper
│   │   └── run_full_pipeline.py     # Full pipeline (Modules 1-5)
│   ├── slurm/
│   │   ├── xenium_benchmark.sh      # Cell proportion benchmark
│   │   ├── xenium_benchmark_rna_gt.sh  # RNA-based GT benchmark
│   │   ├── run_full_pipeline.sh     # Full pipeline job
│   │   └── run_integration.sh       # Module 5 integration
│   ├── output/                      # Protein-based GT results
│   └── output_rna_gt/               # RNA-based GT results
│
├── Cell2Location/                    # Cell2Location (placeholder)
│   ├── src/
│   ├── slurm/
│   └── output/
│
├── RCTD/                             # RCTD (placeholder)
│   ├── src/
│   ├── slurm/
│   └── output/
│
├── Tangram/                          # Tangram (placeholder)
│   ├── src/
│   ├── slurm/
│   └── output/
│
├── evaluation/                       # Shared evaluation code
│   ├── src/
│   │   ├── evaluate_metrics.py      # JSD, RMSE, MAE, Pearson
│   │   └── evaluate_benchmark.py    # Evaluation utilities
│   └── slurm/
│       ├── evaluate_all.sh          # Evaluate all regions
│       └── evaluate_rna_gt.sh       # Evaluate RNA-based GT
│
├── metrics/                          # Aggregated metrics from all methods
├── figures/                          # Comparison figures
├── results/                          # Summary results
└── README.md
```

## Implemented Methods

### CITEgeist (Complete)

CITEgeist is fully implemented with two benchmarking modes:

1. **Protein-based GT** (`xenium_benchmark.sh`)
   - Uses protein markers for ground truth cell types
   - 10 cell types

2. **RNA-based GT** (`xenium_benchmark_rna_gt.sh`) - **Recommended**
   - Uses RNA clustering for ground truth
   - 6 cell types
   - Avoids circular logic

## Usage

### 1. Run CITEgeist Benchmark

```bash
cd Benchmarking/xenium_benchmarking/CITEgeist/slurm

# Protein-based ground truth
sbatch xenium_benchmark.sh

# RNA-based ground truth (recommended)
sbatch xenium_benchmark_rna_gt.sh
```

### 2. Evaluate Results

After all regions complete:

```bash
cd Benchmarking/xenium_benchmarking/evaluation/slurm
sbatch evaluate_all.sh
```

## Metrics

Benchmark metrics (same as simulation benchmarks):

- **JSD** (Jensen-Shannon Divergence): Distribution similarity [0-1], lower is better
- **RMSE** (Root Mean Square Error): Per cell type and global
- **MAE** (Mean Absolute Error): Per cell type and global
- **Pearson r**: Correlation between predicted and true proportions

## Adding New Methods

To add a new deconvolution method:

1. Create directory: `{Method}/src/`, `{Method}/slurm/`, `{Method}/output/`
2. Implement `run_benchmark.py` that:
   - Loads data from `../../xenium_pseudovisium/data/` or `data_rna_gt/`
   - Runs deconvolution
   - Saves predictions in standard format
3. Update evaluation scripts to include new method

## Data Location

Ground truth data is in `../xenium_pseudovisium/`:
- `data/` - Protein-based ground truth (10 cell types)
- `data_rna_gt/` - RNA-based ground truth (6 cell types, recommended)

## Files Generated

```
CITEgeist/output/
└── Xenium_region_{0-4}/
    ├── cell_proportions.csv          # CITEgeist predictions
    ├── Xenium_region_X_deconv_predictions.csv
    └── run_stats.json

metrics/
├── metrics_region_{0-4}.json         # Per-region metrics
├── summary_table.csv                 # All regions summary
├── aggregated_metrics.json           # Aggregated results
└── final_benchmark_report.json       # Complete report
```

## References

- CITEgeist: Chang et al. 2025, bioRxiv
- Cell2Location: Kleshchevnikov et al. 2022, Nature Biotechnology
- RCTD: Cable et al. 2022, Nature Biotechnology
- Tangram: Biancalani et al. 2021, Nature Methods
