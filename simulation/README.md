# Spatial Transcriptomics Simulation Framework

This folder contains the simulation framework for generating synthetic spatial transcriptomics datasets to benchmark deconvolution methods using [scCube](https://github.com/ZJUFanLab/scCube).

## Overview

The simulation framework generates synthetic spatial transcriptomics data with known ground truth cell type proportions and gene expression profiles. Two types of simulations are included:

1. **Highly Segmented**: Simulates distinct spatial regions with relatively homogeneous cell type compositions
2. **Mixed**: Simulates more complex tissue architectures with heterogeneous cell type mixing

## Data Sources

The simulations are based on:
- Single-cell RNA sequencing reference data from Wu et al. 2021 breast cancer dataset
- Spatial coordinates from 10x Visium spatial transcriptomics platform

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
