# CITEgeist: Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends

CITEgeist is a computational method for deconvolving spatial transcriptomics data using spatially-resolved CITE-seq measurements. The pipeline performs both cell-type proportion estimation and gene expression deconvolution in a two-pass approach, leveraging both protein and RNA measurements from the same spatial locations.

## Quick Installation

You can now also install CITEgeist using pip:

```bash
pip install citegeist
```

## Table of Contents
1. [System Requirements](#system-requirements)
    - [Software Dependencies](#software-dependencies)
    - [Hardware Requirements](#hardware-requirements)
2. [Getting Started](#getting-started)
    - [1. Installation](#1-installation)
    - [2. Set Up the Environment](#2-set-up-the-environment)
    - [3. Running CITEgeist](#3-running-citegeist)
3. [Benchmarking and Reproducibility](#benchmarking-and-reproducibility)

## System Requirements

### Software Dependencies

- **Operating System**:
  - Linux
  - macOS
  - Windows 10 with WSL2

- **Python**: 3.10
- **NVIDIA GPU** with 8GB+ VRAM (required for cuOPT QP solver)
#### Key Python Dependencies

See `requirements.txt` for exact version pins. Core dependencies include:
- scanpy, anndata, numpy, pandas, scipy, scikit-learn
- cuopt (GPU-accelerated QP; installed via conda — see `CITEgeist_env.yml`)
- matplotlib, seaborn, squidpy, spatialdata

It is recommended to install the dependencies in the `CITEgeist_env.yml` file for running the notebooks.

### Hardware Requirements
- **RAM**: Minimum 16GB, Recommended 64GB+
- **GPU**: NVIDIA GPU required for Module 3 (cuOPT QP); 8GB+ VRAM, 16GB+ recommended for large datasets
- **Storage**: 16GB minimum for installation; 100GB+ for full multi-patient runs
- **CPU**: Multi-core processor recommended (8+ cores for optimal performance)

---

## Getting Started

### 1. Installation

Install CITEgeist using pip:

```bash
pip install citegeist
```

For development installation:

```bash
git clone https://github.com/leeoesterreich/CITEgeist.git
cd CITEgeist
pip install -e .[dev]
```

### 2. Set Up the Environment

Create and activate the CITEgeist conda environment:

```bash
conda env create -f CITEgeist_env.yml
conda activate CITEgeist_env
```

### 3. Running CITEgeist

#### A. Python API

```python
from CITEgeist.model.citegeist_model import CitegeistModel

model = CitegeistModel(
    sample_name="my_sample",
    adata=adata,                  # AnnData from sq.read.visium()
    output_folder="output/my_sample",
    simulation=False,             # True for simulated data
    resolution="spot",            # "spot" (Visium) or "cell" (VisiumHD/Xenium)
)
model.split_adata()               # Required for real patient data (splits GEX + antibody)
```

Key constructor parameters:
- `sample_name` (str): Sample identifier used in output filenames
- `adata` (AnnData): Spatial data from `sq.read.visium()` — call `split_adata()` before preprocessing
- `output_folder` (str): Directory for results; created automatically
- `simulation` (bool): `True` for simulated data — pass `gene_expression_adata` + `antibody_capture_adata` instead of `adata`
- `resolution` (str): `"spot"` for standard Visium, `"cell"` for single-cell resolution (VisiumHD/Xenium)
- `resolution_overrides` (dict, optional): Override individual resolution parameters (e.g. `{"lambda_sparse": 0.2}`)

**Note:** Module 3 (cuOPT QP) requires a GPU node. Runs on CPU nodes will fail silently.

#### B. Using SLURM (HPC)

For large-scale analyses on HPC clusters, use the provided sbatch scripts in `examples/`:

```bash
sbatch examples/sbatch_patient_phase1.sh      # StarDist nuclei segmentation
sbatch examples/sbatch_patient_phase2.sh      # M3: cuOPT QP proportions (GPU)
sbatch examples/sbatch_patient_phase3.sh      # ViT feature extraction
sbatch examples/sbatch_patient_phase4.sh      # MIL training + cell assignment
sbatch examples/sbatch_patient_phase5_validate.sh  # Validation
sbatch examples/sbatch_bivariate_morans_v3.sh # M4b: bivariate Moran's I
```

Adapt paths and `--gres=gpu:1` directives for your cluster.

---

## Benchmarking and Reproducibility

For specific reproduction of benchmarking tests and detailed methodology, please refer to the 'examples' and 'benchmarking' section in the documentation.

---

