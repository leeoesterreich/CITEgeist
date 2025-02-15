
# CITEgeist: Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends

CITEgeist is a computational method for deconvolving spatial transcriptomics data using spatially-resolved CITE-seq measurements. The pipeline performs both cell-type proportion estimation and gene expression deconvolution in a two-pass approach, leveraging both protein and RNA measurements from the same spatial locations.

## Table of Contents
1. [System Requirements](#system-requirements)
    - [Software Dependencies](#software-dependencies)
    - [Hardware Requirements](#hardware-requirements)
2. [Getting Started](#getting-started)
    - [1. Download the Code and Data](#1-download-the-code-and-data)
    - [2. Set Up the Environment](#2-set-up-the-environment)
    - [3. Obtain Gurobi License](#3-obtain-gurobi-license)
    - [4. Running CITEgeist](#4-running-citegeist)
3. [Benchmarking and Reproducibility](#benchmarking-and-reproducibility)
4. [Run the Analysis](#run-the-analysis)

## System Requirements

### Software Dependencies

- **Operating System**:  
  - Linux  
  - macOS  
  - Windows 10 with WSL2  

- **Python**: 3.10
- **Gurobi** [version > 3.9](https://www.gurobi.com/downloads/gurobi-software/)
#### Key Python Dependencies
- scanpy==1.10.4
- anndata==0.11.3
- numpy==1.26.4
- pandas==2.2.3
- scipy==1.13.1
- scikit-learn==1.6.1
- gurobipy==11.0.2 (requires license)
- matplotlib==3.10.0
- seaborn==0.13.2
- h5py==3.12.1
- squidpy==1.6.2
- spatialdata==0.2.5.post0

It is recommended to install the dependencies in the `CITEgeist_env.yml` file for running the notebooks.

### Hardware Requirements
- **RAM**: Minimum 16GB, Recommended 64GB+
- **Storage**: 16GB minimum for installation and basic analysis
- **CPU**: Multi-core processor recommended (8+ cores for optimal performance)

---

## Getting Started

### 1. Download the Code and Data (Instructions for Peer Reviewers)

1. Download the code from Figshare: [https://figshare.com/s/34e456fd7786e5211acc](https://figshare.com/s/34e456fd7786e5211acc)
2. Unzip the downloaded file to your preferred location.
3. Download the data from GEO (Reviewers: see Data Availability section in the manuscript for the GEO link and Private Access Token).
4. Run the following code to strip unique identifiers required by GEO:

```bash
# go to the 'data' directory
cd data

# untar the raw files
mkdir -pv ./GEO_data
tar -xvf GEO_data_RAW.tar -C ./GEO_data

# run the py preprocessing script
## round 1) aggregate the files by sample
python3 ./delete_all_but_essential.py --folder GEO_data # select option: 1

## round 2) remove the prefix from necessary files
python3 ./delete_all_but_essential.py --folder GEO_data # select option: 2
```

Note: When prompted, select Option 1 or 2 and type 'Yes' to confirm.

### 2. Set Up the Environment

- Install dependencies using the provided environment file:

```bash
conda env create -f CITEgeist_env.yml
```

- Activate the environment and set up a Jupyter kernel:

```bash
conda activate CITEgeist_env
```

### 3. Obtain Gurobi License

CITEgeist requires a Gurobi license (free for academic use):

1. Sign up for an academic license at: [https://www.gurobi.com/downloads/end-user-license-agreement-academic/](https://www.gurobi.com/downloads/end-user-license-agreement-academic/)
2. Follow the instructions to download and install your license.
3. Update the license file path in the notebooks to match your local license location.

### 4. Running CITEgeist

You can run CITEgeist in two ways:

#### A. Using Jupyter Notebooks
- Update data paths in the top of the notebooks to match your local directory structure.
- Expected runtime on a standard computer (16 threads, 32GB RAM):
  - Vignette 1: ~2 hours
  - Vignette 2: ~2 hours
  - Vignette 3: ~10 hours

#### Key Parameters:
- `radius`: Radius for neighbor detection (default: 4)
- `lambda_reg`: Regularization strength for cell proportion estimation (default: 0.001)
- `alpha_elastic`: Elastic net mixing parameter for cell proportion estimation (default: 0.7)
- `max_y_change`: Maximum allowed change in Y values (default: 0.2)

#### Optional Parameters:
- `profiling_only`: Set for cell-type proportions only.
- `max_workers`: Number of parallel workers.
- `checkpoint_interval`: Checkpoint saving interval.

#### B. Using SLURM Distribution
For large-scale analyses, you can use the provided `CITEgeist/examples/sbatch_sample.sh` script for distributed computing.

---

## Benchmarking and Reproducibility

For specific reproduction of benchmarking tests and detailed methodology, please refer to the 'examples' and 'benchmarking' section in the documentation.

---

## Run the Analysis

You can either:

#### A. Run the notebooks directly:
- Update data paths in the notebooks to match your local directory structure.
- Expected runtime: ~2 hours on a standard computer (16 threads, 32GB RAM).

#### B. Use SLURM distribution:
- Use the provided `examples/sbatch_sample.sh` script for distributed computing.

---

## Additional System Requirements
- **RAM**: 32GB (minimum)
- **CPU**: 16 threads (recommended)
- **Storage**: Sufficient space for the GEO dataset
- **Operating System**: Linux/Unix recommended (Windows users may need additional configuration)
