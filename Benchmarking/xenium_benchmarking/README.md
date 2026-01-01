# Xenium Benchmarking Framework for CITEgeist

This directory contains a benchmarking framework for evaluating CITEgeist on **real spatial multi-omic data** from 10x Xenium. Unlike simulation-based benchmarks, this approach uses single-cell resolution Xenium data to create pseudo-Visium spots, providing ground truth cell proportions from actual tissue.

## Overview

The Xenium platform provides single-cell resolution spatial transcriptomics with protein detection. By aggregating single cells into pseudo-Visium spots (55µm diameter), we can:

1. Create pseudo-Visium expression profiles (GEX + protein)
2. Use single-cell protein expression to define ground truth cell types
3. Calculate true cell type proportions per spot
4. Benchmark CITEgeist deconvolution against this ground truth

## Data Source

**10x Xenium Renal Cell Carcinoma Dataset:**
- Location: `/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/`
- **465,534 cells** with single-cell resolution
- **405 genes** (Gene Expression)
- **27 proteins** including immune markers (CD4, CD8A, CD68, CD3E, CD20, etc.)

## Directory Structure

```
xenium_benchmarking/
├── src/                          # Python modules
│   ├── load_xenium.py           # Load Xenium h5/parquet data
│   ├── create_pseudo_spots.py   # Hexagonal grid + cell aggregation
│   ├── define_cell_types.py     # Protein-based cell typing
│   ├── generate_ground_truth.py # Ground truth proportions
│   ├── split_regions.py         # Region splitting
│   ├── run_benchmark.py         # CITEgeist wrapper
│   └── evaluate_metrics.py      # Benchmark metrics
├── data/
│   ├── h5ad_objects/            # Pseudo-Visium h5ad files
│   │   ├── Xenium_region_0_GEX.h5ad
│   │   ├── Xenium_region_0_CITE.h5ad
│   │   └── ...
│   └── ground_truth/            # Ground truth proportions
│       ├── Xenium_region_0_prop.csv
│       └── ...
├── output/                       # CITEgeist outputs
├── metrics/                      # Benchmark results
├── slurm/                        # SLURM job scripts
│   ├── xenium_benchmark.sh      # Main benchmark array job
│   └── evaluate_all.sh          # Evaluation script
└── README.md
```

## Cell Types (10 defined + Unassigned)

Based on available Xenium proteins:

| Cell Type | Major Markers | Minor Markers |
|-----------|---------------|---------------|
| CD4+ T cells | CD3E, CD4 | CD45 |
| CD8+ T cells | CD3E, CD8A | CD45, GranzymeB |
| B cells | CD20 | CD45 |
| Plasma cells | CD138 | - |
| Macrophages | CD68 | CD163, CD45, HLA-DR |
| Dendritic cells | CD11c, HLA-DR | CD45 |
| NK cells | CD16, GranzymeB | CD45 |
| Epithelial | PanCK | E-Cadherin |
| Endothelial | CD31 | - |
| Fibroblasts | alphaSMA | Vimentin |

## Usage

### Step 1: Generate Pseudo-Visium Data (Already Done)

The data has been pre-generated. To regenerate:

```bash
python -c "
import sys
sys.path.insert(0, 'Benchmarking/xenium_benchmarking/src')
from split_regions import create_all_region_datasets

create_all_region_datasets(
    data_dir='/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma',
    output_dir='Benchmarking/xenium_benchmarking/data',
    n_regions=5,
)
"
```

### Step 2: Run CITEgeist Benchmark

Submit the SLURM array job:

```bash
cd Benchmarking/xenium_benchmarking/slurm
mkdir -p slurm_log
sbatch xenium_benchmark.sh
```

This runs CITEgeist on all 5 regions in parallel.

### Step 3: Evaluate Results

After all regions complete:

```bash
sbatch evaluate_all.sh
```

Or run interactively:

```bash
python -c "
import sys
sys.path.insert(0, 'Benchmarking/xenium_benchmarking/src')
from evaluate_metrics import evaluate_all_regions

results = evaluate_all_regions(
    input_dir='Benchmarking/xenium_benchmarking/data',
    output_dir='Benchmarking/xenium_benchmarking/metrics',
    predictions_dir='Benchmarking/xenium_benchmarking/output',
)
print(f'JSD: {results[\"JSD_median_mean\"]:.4f} ± {results[\"JSD_median_std\"]:.4f}')
print(f'RMSE: {results[\"RMSE_global_mean\"]:.4f} ± {results[\"RMSE_global_std\"]:.4f}')
print(f'Pearson: {results[\"Pearson_r_mean\"]:.4f} ± {results[\"Pearson_r_std\"]:.4f}')
"
```

## Metrics

Benchmark metrics (same as simulation benchmarks):

- **JSD** (Jensen-Shannon Divergence): Distribution similarity [0-1], lower is better
- **RMSE** (Root Mean Square Error): Per cell type and global
- **MAE** (Mean Absolute Error): Per cell type and global
- **Pearson r**: Correlation between predicted and true proportions

## Dataset Statistics

| Metric | Value |
|--------|-------|
| Total cells | 465,534 |
| Total pseudo-spots | 7,054 |
| Cells per spot (mean) | 18.1 |
| Cells per spot (median) | 16.0 |
| Regions | 5 |
| Spots per region | ~1,350-1,465 |
| Genes | 405 |
| Proteins | 27 |

## Comparison to Simulation Benchmarks

| Aspect | Simulation | Xenium |
|--------|------------|--------|
| Ground truth source | Synthetic cell mixtures | Real tissue composition |
| Cell types | Simulated profiles | Protein-based classification |
| Spatial patterns | scCube simulation | Actual tissue architecture |
| Sample size | 5 replicates × ~2k spots | 5 regions × ~1.4k spots |
| Genes | Wu et al. 12k | Xenium panel (405) |
| Proteins | Simulated (18) | Measured (27) |

## Key Differences from Simulation

1. **Real tissue heterogeneity**: Xenium captures actual tumor microenvironment
2. **Imperfect cell typing**: Protein-based classification has uncertainty
3. **Limited gene panel**: Xenium uses targeted panel, not whole transcriptome
4. **Different protein markers**: Xenium proteins may differ from CITE-seq panels

## Files Generated

After running the benchmark:

```
data/
├── h5ad_objects/
│   ├── Xenium_region_{0-4}_GEX.h5ad    # Gene expression per region
│   └── Xenium_region_{0-4}_CITE.h5ad   # Protein expression per region
├── ground_truth/
│   └── Xenium_region_{0-4}_prop.csv    # Ground truth proportions
└── dataset_summary.json

output/
└── Xenium_region_{0-4}/
    ├── cell_proportions.csv             # CITEgeist predictions
    └── run_stats.json                   # Runtime statistics

metrics/
├── metrics_region_{0-4}.json            # Per-region metrics
├── summary_table.csv                    # All regions summary
└── aggregated_metrics.json              # Aggregated results
```

## References

- 10x Xenium: https://www.10xgenomics.com/platforms/xenium
- CITEgeist: Chang et al. 2025, bioRxiv
- Pseudo-Visium approach: Weber et al. 2025, Genome Biology
