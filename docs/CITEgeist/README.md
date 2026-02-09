# CITEgeist

CITEgeist (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends) is a computational method for deconvolving spatial transcriptomics data using spatially-resolved CITE-seq measurements. The pipeline performs both cell-type proportion estimation and gene expression deconvolution in a two-pass approach, leveraging both protein and RNA measurements from the same spatial locations.

## System Requirements

### Software Dependencies
- Operating System: 
  - Linux (tested on Ubuntu 20.04 LTS, CentOS 7)
- Python 3.10
- Required Python packages:
  - numpy==1.23.0
  - pandas==2.2.2
  - scanpy==1.10.1
  - pyarrow==17.0.0
  - scipy==1.10.1
  - anndata==0.10.7
  - scikit-learn==1.5.0
  - matplotlib==3.9.4
  - seaborn==0.13.2
  - tqdm==4.66.4
  - h5py==3.11.0
  - networkx==3.2.1
  - Gurobi >= 9.5.0 (requires license)
  - squidpy==1.6.1 (for spatial analysis)
  - spatialdata==0.2.1 (for spatial data handling)
  - dask==2024.7.1 (for parallel computing)
  - zarr==2.18.2 (for data storage)
  - python-igraph==0.11.5 (for graph analysis)
  - leidenalg==0.10.2 (for clustering)

### Hardware Requirements
- RAM: Minimum 16GB, Recommended 64GB+
- Storage: 10GB minimum for installation and basic analysis
- CPU: Multi-core processor recommended (8+ cores for optimal performance)
- GPU: Optional, but recommended for large datasets

## Installation Guide

### Step 1: Set up Python Environment
```bash
# Create and activate a new conda environment
conda create -n citegeist python=3.10
conda activate citegeist
```

### Step 2: Install Dependencies
```bash
# Install main dependencies
pip install numpy==1.21.0 pandas==1.3.0 scanpy==1.8.0 pyarrow==6.0.0 scipy==1.7.0

# Install Gurobi
conda install -c gurobi gurobi
```

### Step 3: Install CITEgeist
```bash
git clone https://github.com/yourusername/CITEgeist.git
cd CITEgeist
pip install -e .
```

Typical installation time: 10-15 minutes on a standard desktop computer

## Demo

### Running the Demo
1. Download the demo dataset:
```bash
python scripts/download_demo_data.py
```

2. Run the demo analysis:
```bash
python run_citegeist.py \
    --radius 4 \
    --lambda_reg 1.0 \
    --alpha_elastic 0.7 \
    --max_y_change 0.2 \
    --input_folder demo/data \
    --output_folder demo/results \
    --sample_prefix "demo_sample"
```

### Expected Output
The demo will generate the following files in the `demo/results` directory:
- `demo_sample_cellprop.csv`: Cell type proportions per spot
- `demo_sample_gene_expression_pass1.parquet`: Initial gene expression profiles
- `demo_sample.log`: Processing logs and metrics

Expected runtime for demo: 30-45 minutes on a standard desktop computer (16GB RAM, 4 cores)

## Instructions for Use

### Running CITEgeist on Your Data

1. Prepare your input data:
   - Format your CITE-seq data as an AnnData object (`.h5ad` format)
   - Ensure both gene expression and protein measurements are included
   - Include spatial coordinates for each spot

2. Run the analysis:
```bash
python run_citegeist.py \
    --radius <radius_value> \
    --lambda_reg <regularization_strength> \
    --alpha_elastic <elastic_net_mixing> \
    --max_y_change <max_change_value> \
    --input_folder /path/to/your/data \
    --output_folder /path/to/output \
    --sample_prefix "your_sample_name"
```

Key Parameters:
- `radius`: Radius for neighbor detection (default: 4)
- `lambda_reg`: Regularization strength (default: 0.001)
- `alpha_elastic`: Elastic net mixing parameter (default: 0.7)
- `max_y_change`: Maximum allowed change in Y values (default: 0.2)

Optional Parameters:
- `profiling_only`: Set for cell-type proportions only
- `skip_pass2`: Skip refinement pass
- `max_workers`: Number of parallel workers
- `checkpoint_interval`: Checkpoint saving interval

### Performance Optimization
- For large datasets (>100k spots):
  - Increase available RAM to 64GB+
  - Use `max_workers` parameter for parallel processing
  - Enable checkpointing with `checkpoint_interval`
- For faster processing:
  - Use SSD storage for input/output
  - Enable GPU acceleration if available
  - Adjust batch sizes based on available memory

## Reproduction Instructions

Canonical manuscript reproduction (v5) now lives under `repro/`:

```bash
python -m repro.cli.validate_env --config repro/config/example.paths.yaml
python -m repro.cli.repro --target v5_figures --config repro/config/example.paths.yaml
```

For full instructions:

- `repro/README.md`
- `repro/runbooks/reviewer_quickstart.md`
- `repro/runbooks/full_reproduction.md`
