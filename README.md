# CITEgeist
Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends

## System Requirements

### Software Dependencies
- Operating System: 
  - Linux (tested on Ubuntu 20.04 LTS, CentOS 7)
  - macOS (tested on Monterey 12.0+)
  - Windows 10 with WSL2
- Python 3.7+ (tested on Python 3.7.9, 3.8.10, and 3.9.7)
- Required Python packages:
  - numpy >= 1.21.0
  - pandas >= 1.3.0
  - scanpy >= 1.8.0
  - pyarrow >= 6.0.0
  - scipy >= 1.7.0
  - Gurobi >= 9.5.0 (requires license)

### Hardware Requirements
- RAM: Minimum 16GB, Recommended 64GB+
- Storage: 10GB minimum for installation and basic analysis
- CPU: Multi-core processor recommended (8+ cores for optimal performance)
- GPU: Optional, but recommended for large datasets

## Installation Guide

### Step 1: Set up Python Environment
```bash
# Create and activate a new conda environment
conda create -n citegeist python=3.9
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
- `demo_sample_gene_expression_pass2.parquet`: Refined gene expression profiles
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

To reproduce the results from our manuscript:
1. Download the benchmark datasets:
```bash
python scripts/download_benchmark_data.py
```

2. Run the benchmark analysis:
```bash
python scripts/run_benchmarks.py \
    --dataset_dir benchmark/data \
    --output_dir benchmark/results
```

The benchmark scripts will generate figures and tables matching those in the manuscript. Expected runtime: 4-6 hours on a standard desktop computer.

For detailed methodology and additional analysis scripts, please refer to our [benchmarking documentation](docs/Benchmarking/README.md).

Copyright (c) 2024 Lee Oesterreich Lab
