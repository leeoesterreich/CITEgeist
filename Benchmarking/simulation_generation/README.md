# Spatial Transcriptomics Simulation Framework

This folder contains the simulation framework for generating synthetic spatial transcriptomics datasets to benchmark deconvolution methods using [scCube](https://github.com/ZJUFanLab/scCube).

## Overview

The simulation framework generates synthetic spatial transcriptomics data with known ground truth cell type proportions and gene expression profiles. Two types of simulations are included:

1. **Highly Segmented**: Simulates distinct spatial regions with relatively homogeneous cell type compositions
2. **Mixed**: Simulates more complex tissue architectures with heterogeneous cell type mixing

## Data Sources

The simulations are based on:
- Single-cell RNA sequencing reference data from Wu et al. 2021 breast cancer dataset
- Spatial coordinates from 10x Visium spatial transcriptomics platform, generated using scCube's hex-based Visium profile

## Simulation Parameters

- Number of spots: ~5000 per simulation
- Cell types: 9 major breast cancer cell types
  - Cancer Epithelial
  - Normal Epithelial  
  - T-cells
  - B-cells
  - Myeloid cells
  - PVL
  - Plasmablasts
  - Endothelial cells
  - CAFs

## Notebooks

### `scCube_simulation_12k.ipynb`
Main simulation notebook. Uses scCube to generate ~5,000 pseudo-Visium spots per condition from Wu et al. 2021 breast cancer scRNA-seq reference.

### `build_replicate_objects.ipynb`
Assembles 3 replicates per condition into h5ad objects ready for benchmarking methods.

### `ST_GEX_GT_sim.ipynb`
Generates per-cell-type ground truth gene expression layers for GEX benchmarking.

## Outputs

Generated data is consumed by `../simulation_benchmarking/` methods:
- `{condition}/replicate_{1-3}_GEX.h5ad` — simulated gene expression
- `{condition}/replicate_{1-3}_CITE.h5ad` — simulated protein data
- `{condition}/replicate_{1-3}_prop.csv` — ground truth proportions
- `{condition}/replicate_{1-3}/` — per-cell-type GEX ground truth layers

Conditions: `high_seg` (highly segmented) and `mixed` (heterogeneous mixing).
