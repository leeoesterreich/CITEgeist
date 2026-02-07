# Xenium Pseudo-Visium Data Generation

This directory contains scripts for creating pseudo-Visium spots from 10x Xenium single-cell resolution data. The generated datasets are used by benchmarking methods in `../xenium_benchmarking/`.

## Overview

The Xenium platform provides single-cell resolution spatial transcriptomics with protein detection. By aggregating single cells into pseudo-Visium spots (55µm diameter), we can:

1. Create pseudo-Visium expression profiles (GEX + protein)
2. Use single-cell protein/RNA expression to define ground truth cell types
3. Calculate true cell type proportions per spot
4. Provide ground truth for benchmarking deconvolution methods

## Data Source

**10x Xenium Renal Cell Carcinoma Dataset:**
- Location: `/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/`
- **465,534 cells** with single-cell resolution
- **405 genes** (Gene Expression)
- **27 proteins** including immune markers (CD4, CD8A, CD68, CD3E, CD20, etc.)

## Directory Structure

```
xenium_pseudovisium/
├── src/                              # Python modules
│   ├── load_xenium.py               # Load Xenium h5/parquet data
│   ├── create_pseudo_spots.py       # Hexagonal grid + cell aggregation
│   ├── split_regions.py             # Region splitting for cross-validation
│   ├── define_cell_types.py         # Protein-based cell typing (10 types)
│   ├── rna_cell_types.py            # RNA-based cell typing (6 types)
│   ├── generate_ground_truth.py     # Cell proportion ground truth
│   ├── generate_expression_weighted_gt.py  # Expression-weighted proportions
│   └── generate_gex_ground_truth.py # Per-cell-type GEX ground truth
├── slurm/
│   └── generate_gex_gt.sh           # SLURM script for GEX ground truth
├── data/                             # Protein-based ground truth
│   ├── h5ad_objects/                # Pseudo-Visium h5ad files
│   ├── ground_truth/                # Cell proportion CSVs
│   └── ground_truth_gex/            # Per-cell-type GEX
├── data_rna_gt/                      # RNA-based ground truth (recommended)
│   ├── h5ad_objects/
│   ├── ground_truth/
│   └── ground_truth_gex/
└── _archive/                         # Archived scripts
```

## Two Ground Truth Approaches

### 1. Protein-based (`data/`)
- Uses 27 Xenium protein markers for cell classification
- Defines **10 cell types**: CD4+ T, CD8+ T, B cells, Plasma, Macrophages, Dendritic, NK, Epithelial, Endothelial, Fibroblasts
- Cell type markers defined in `define_cell_types.py`

### 2. RNA-based (`data_rna_gt/`) - **Recommended**
- Uses Xenium's built-in RNA clustering for cell classification
- Defines **6 cell types**: T cells, B cells, Macrophages, Fibroblasts, Epithelial, Endothelial
- Avoids circular logic (doesn't use protein markers for both GT and deconvolution)
- Reference: Zhao et al. (2025) BMC Bioinformatics

## Usage

### Generate Pseudo-Visium Data (Already Done)

The data has been pre-generated. To regenerate:

```bash
python -c "
import sys
sys.path.insert(0, 'Benchmarking/xenium_pseudovisium/src')
from split_regions import create_all_region_datasets

create_all_region_datasets(
    data_dir='/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma',
    output_dir='Benchmarking/xenium_pseudovisium/data',
    n_regions=5,
)
"
```

### Generate GEX Ground Truth

Submit SLURM job:

```bash
cd Benchmarking/xenium_pseudovisium/slurm
sbatch generate_gex_gt.sh
```

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

## Files Generated

```
data/
├── h5ad_objects/
│   ├── Xenium_region_{0-4}_GEX.h5ad    # Gene expression per region
│   └── Xenium_region_{0-4}_CITE.h5ad   # Protein expression per region
├── ground_truth/
│   └── Xenium_region_{0-4}_prop.csv    # Ground truth proportions
├── ground_truth_gex/
│   └── Xenium_region_{0-4}/            # Per-cell-type GEX
│       └── {CellType}_GT.csv
├── cell_to_spot_mapping.csv
└── dataset_summary.json
```

## References

- 10x Xenium: https://www.10xgenomics.com/platforms/xenium
- Zhao et al. (2025) Benchmarking cell type annotation methods for 10x Xenium. BMC Bioinformatics 26(1):25
