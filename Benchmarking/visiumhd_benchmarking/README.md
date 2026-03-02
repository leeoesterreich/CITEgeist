# Visium HD H&E Morphology Benchmarking

Benchmarking framework for H&E morphology-based single-cell type assignment using pre-trained Vision Transformer (ViT) features with Multiple Instance Learning (MIL) aggregation.

## Overview

This benchmark evaluates morphology-based cell type prediction on Visium HD pILC data:
- **Input**: H&E whole slide images + Visium HD single-cell annotations
- **Method**: ViT (UNI) feature extraction + Proportion-guided MIL
- **Evaluation**: Single-cell assignment accuracy with ground truth from ScType

## Data Sources

| Sample | Cells | Key Cell Types |
|--------|-------|----------------|
| TP08-2202 | 464,680 | Fibroblast (170K), Plasma (45K), Macrophage (37K), Cancer (32K) |
| TP12-880 | 120,577 | Fibroblast (58K), Cancer (16K), Macrophage (3K) |
| TP15-M509 | 440,535 | Fibroblast (175K), Macrophage (58K), Cancer (45K), T_Cell (38K) |

**Total:** ~1M cells with ScType ground truth (excluding Unknown)

**Location:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project/`

## Directory Structure

```
visiumhd_benchmarking/
├── src/                    # Source code
│   ├── create_pseudo_visium.py    # Generate pseudo-Visium spots
│   ├── run_cellpose_he.py         # Cellpose segmentation
│   ├── extract_patches_he.py      # H&E patch extraction
│   ├── vit_extractor.py           # ViT feature extraction
│   ├── proportion_mil.py          # MIL aggregation module
│   ├── train_mil.py               # Training pipeline
│   ├── evaluate_single_cell.py    # Evaluation pipeline
│   └── run_benchmark.py           # Main benchmark script
├── data/                   # Generated data files
│   └── {sample}/
│       ├── pseudo_visium_spots.csv
│       ├── cell_to_spot_mapping.csv
│       ├── ground_truth_proportions.parquet
│       └── cellpose_masks.npy
├── patches/                # Extracted H&E patches
│   └── {sample}/spot_{id}_patches.npy
├── models/                 # Trained model checkpoints
│   ├── mil_head.pt
│   └── vae_baseline.pt
├── results/                # Evaluation results
│   ├── feature_comparison/
│   └── evaluation/
├── slurm/                  # SLURM submission scripts
│   └── sbatch_benchmark.sh
└── logs/                   # SLURM job logs
```

## Pipeline

1. **Pseudo-Visium Creation**: Generate 55μm spots from Visium HD data
2. **Cellpose Segmentation**: Re-segment H&E images for consistency
3. **Patch Extraction**: Crop 224×224 RGB patches per nucleus
4. **ViT Feature Extraction**: Extract 768-dim embeddings using UNI model
5. **MIL Training**: Train attention-based aggregation with proportion supervision
6. **Single-Cell Assignment**: Hungarian assignment with predicted proportions
7. **Evaluation**: Compare against ground truth ScType annotations

## Usage

```bash
# Run full benchmark pipeline
sbatch slurm/sbatch_benchmark.sh

# Or run individual steps
python src/create_pseudo_visium.py --sample TP08-2202
python src/run_cellpose_he.py --sample TP08-2202
python src/extract_patches_he.py --sample TP08-2202
python src/train_mil.py --train_samples TP08-2202 TP12-880
python src/evaluate_single_cell.py --sample TP15-M509
```

## Dependencies

- `timm` - PyTorch Image Models (for ViT)
- `cellpose` - Nuclear segmentation
- `torch`, `scipy`, `numpy`, `pandas` - Standard scientific stack
- UNI model weights (requires Hugging Face access)

## Related Documents

- Design: `docs/plans/2026-03-02-he-morphology-vit-mil-design.md`
- Implementation Plan: `docs/plans/2026-03-02-he-morphology-vit-mil-implementation.md`
