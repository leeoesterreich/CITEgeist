# Benchmarking Documentation

This document describes the benchmarking framework used to evaluate CITEgeist's performance against other state-of-the-art spatial transcriptomics deconvolution methods.

## Overview

The benchmarking system evaluates two main aspects:
1. Cell Type Deconvolution Performance
2. Gene Expression Prediction Accuracy

## Metrics

### Cell Type Deconvolution Metrics

The following metrics are used to evaluate cell type proportion predictions:

1. **Jensen-Shannon Divergence (JSD)**
   - Measures the similarity between predicted and ground truth cell type distributions
   - Scale: 0 (identical) to 1 (completely different)
   - Calculated per spot and reported as median value

2. **Root Mean Square Error (RMSE)**
   - Calculated both per cell type and globally
   - Measures the square root of the average squared differences between predictions and ground truth
   - Lower values indicate better performance

3. **Mean Absolute Error (MAE)**
   - Calculated both per cell type and globally
   - Measures the average absolute differences between predictions and ground truth
   - Lower values indicate better performance

4. **Pearson Correlation**
   - Measures the linear correlation between predicted and ground truth proportions
   - Scale: -1 to 1 (1 indicates perfect positive correlation)

### Gene Expression Metrics

For gene expression prediction evaluation:

1. **RMSE (Root Mean Square Error)**
   - Calculated after log1p transformation of expression values
   - Computed per cell type and averaged across all cell types
   - Lower values indicate better performance

2. **NRMSE (Normalized RMSE)**
   - RMSE normalized by either:
     - Range of ground truth values ('range' normalization)
     - Mean of ground truth values ('mean' normalization)
   - Allows for comparison across different scales
   - Lower values indicate better performance

3. **MAE (Mean Absolute Error)**
   - Average absolute difference between predicted and ground truth expression
   - Calculated after log1p transformation
   - Less sensitive to outliers than RMSE
   - Lower values indicate better performance

## Benchmarking Process

### Cell Type Deconvolution

1. Input Requirements:
   - Ground truth cell type proportions matrix
   - Predicted cell type proportions matrix
   - Matching spot IDs and cell type names between matrices

2. Process:
   - Calculates JSD for each spot
   - Computes RMSE and MAE for each cell type
   - Calculates global metrics (Sum RMSE, Sum MAE)
   - Determines overall correlation

### Gene Expression Prediction

1. Input Requirements:
   - Directory of ground truth expression files (one per cell type)
   - Directory of predicted expression files
   - Files must share common gene names and spot IDs

2. Process:
   - Log1p transformation of expression values
   - Calculation of metrics per cell type
   - Computation of overall statistics (average/median)
   - Optional normalization for NRMSE calculation

## Visualization and Analysis

The benchmarking results are visualized using R scripts that generate:
- Box plots comparing methods across metrics
- Statistical analysis including ANOVA and Tukey's HSD tests
- Separate analyses for different experimental conditions (e.g., highly segmented vs. mixed patterns)

## Runtime Analysis

Runtime performance is tracked for each method:
- Execution time is measured for each replicate
- Results are aggregated and compared across methods
- Statistical significance of runtime differences is assessed

## File Structure

```
benchmarking/
├── src/
│   ├── benchmarking_gex.py       # Gene expression benchmarking
│   ├── benchmarking_spot_deconv.py # Cell type deconvolution benchmarking
│   └── citegeist_bench_wrapper.py # Wrapper for running benchmarks
├── Figures/
│   ├── CITEgeist_GEX_figures.Rmd # Gene expression visualization
│   └── CITEgeist_prop_figures.Rmd # Proportion prediction visualization
└── Runtime/
    ├── high_seg/                 # Runtime analysis for highly segmented patterns
    └── mixed/                    # Runtime analysis for mixed patterns
```
