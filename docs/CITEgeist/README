# CITEgeist Documentation

## Overview
CITEgeist (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends) is a computational method for deconvolving spatial transcriptomics data using spatially-resolved CITE-seq measurements. The pipeline performs both cell-type proportion estimation and gene expression deconvolution in a two-pass approach, leveraging both protein and RNA measurements from the same spatial locations.

## Requirements
- Python 3.7+
- Required packages:
  - numpy
  - pandas
  - scanpy
  - pyarrow
  - scipy
  - Gurobi (for optimization)

## Input Data
The pipeline requires a single AnnData object containing both gene expression and protein measurements from spatial CITE-seq data:

1. **Spatial CITE-seq Data** (`.h5ad` format):
   - Gene expression counts matrix
   - Antibody capture counts matrix
   - Spatial coordinates for each spot
   - Both measurements must be from the same spatial locations

Note: Unlike other reference-based spatial deconvolution methods, CITEgeist does not require a single-cell RNA sequencing reference dataset. Instead, it leverages the protein measurements from the spatial CITE-seq data directly for deconvolution.

## Pipeline Parameters

### Key Parameters:
- `radius`: Radius for neighbor detection in spatial analysis
- `lambda_reg`: Regularization strength for elastic net (default: 0.001)
- `alpha_elastic`: Elastic net mixing parameter (0=L2, 1=L1)
- `max_y_change`: Maximum allowed change in Y values between iterations (0-1)
- `global_enrichment_weight`: Weight for global enrichment in gene expression optimization (default: 0.5)
- `local_enrichment_weight`: Weight for local enrichment in gene expression optimization (default: 0.5)

### Optional Parameters:
- `profiling_only`: If set, only compute cell-type proportions without gene expression deconvolution
- `skip_pass2`: If set, skip pass 2 and only run pass 1
- `max_workers`: Number of parallel workers for computation
- `checkpoint_interval`: Interval for saving checkpoints during processing

## Pipeline Workflow

1. **Initialization**:
   - Load and validate spatial CITE-seq data
   - Set up logging and output directories
   - Initialize model parameters

2. **Pass 1 - Initial Deconvolution**:
   - Process antibody capture data to identify cell-type profiles
   - Optimize cell-type proportions using protein measurements
   - Perform initial gene expression deconvolution
   - Export results for each cell type

3. **Pass 2 (Optional) - Refinement**:
   - Use Pass 1 results to compute global priors
   - Refine cell-type proportions using both protein and gene expression information
   - Optimize gene expression with prior information
   - Generate final deconvolution results

## Output Files

The pipeline generates several output files:

1. **Cell Type Proportions**:
   - Format: CSV files
   - Contents: Estimated proportions of each cell type per spot
   - Location: `{output_dir}/{sample_name}_cellprop.csv`

2. **Gene Expression Profiles**:
   - Format: Parquet files
   - Contents: Deconvolved gene expression for each cell type
   - Location: `{output_dir}/{sample_name}_gene_expression_pass{1,2}.parquet`

3. **Layer Files**:
   - Format: CSV files
   - Contents: Cell-type specific expression layers
   - Location: `{output_dir}/{sample_name}_pass{1,2}/layers/`

4. **Logs**:
   - Format: Text files
   - Contents: Processing logs and metrics
   - Location: `{output_dir}/{sample_name}.log`

## Usage Example

```bash
python run_citegeist.py \
    --radius 4 \
    --lambda_reg 1.0 \
    --alpha_elastic 0.7 \
    --max_y_change 0.2 \
    --input_folder /path/to/input \
    --output_folder /path/to/output \
    --sample_prefix "sample_name"
```

## Benchmarking

The pipeline includes benchmarking capabilities for:
- Cell type proportion accuracy (RMSE, MAE)
- Gene expression deconvolution accuracy
- Performance comparison with other methods (Cell2Location, RCTD, etc.)

## Troubleshooting

Common issues and solutions:

1. **Memory Issues**:
   - Increase available memory (recommended: 64GB+)
   - Reduce batch size or use checkpointing

2. **Runtime Optimization**:
   - Adjust `max_workers` for parallel processing
   - Use GPU acceleration when available
   - Enable checkpointing for long runs

3. **Data Compatibility**:
   - Ensure input data is properly normalized
   - Verify gene name consistency between reference and spatial data
   - Check for required columns in metadata

## References

For more information about the method and its applications, please refer to:
- CITEgeist GitHub repository
- Associated publications
- Documentation for dependent packages (Scanpy, Gurobi)
