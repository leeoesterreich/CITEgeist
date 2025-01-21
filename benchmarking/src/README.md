Benchmarking scripts
# Benchmarking Scripts

This folder contains scripts for benchmarking spatial transcriptomics deconvolution methods.

## Core Benchmarking Scripts

### Main Benchmarking Modules
- `benchmarking_gex.py` - Core module for benchmarking gene expression predictions:
  - Calculates RMSE, NRMSE and MAE metrics between predicted and ground truth gene expression layers
  - Supports different normalization methods ('range' or 'mean')
  - Handles per-cell-type and overall statistics

- `benchmarking_spot_deconv.py` - Core module for benchmarking cell type deconvolution:
  - Evaluates accuracy of predicted cell type proportions against ground truth
  - Calculates correlation metrics, RMSE, MAE, and Jensen-Shannon divergence
  - Provides both spot-level and cell-type-level analysis

### Method-Specific Wrappers
- `cell2location_bench_wrapper.py` - Wrapper for benchmarking Cell2location results
- `tangram_bench_wrapper.py` - Wrapper for benchmarking Tangram results  
- `RCTD_bench_wrapper.py` - Wrapper for benchmarking RCTD results
- `seurat_bench_wrapper.py` - Wrapper for benchmarking Seurat results

## Usage

The core benchmarking modules can be used directly or through the method-specific wrappers.
Please refer to the 'help' message of each script for more information on how to use them.
Each script can be executed from the command line with the appropriate arguments.