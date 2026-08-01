# CITEgeist

CITEgeist (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends) is a computational method for deconvolving spatial transcriptomics data using spatially-resolved CITE-seq measurements. The pipeline performs both cell-type proportion estimation and gene expression deconvolution in a two-pass approach, leveraging both protein and RNA measurements from the same spatial locations.

## System Requirements

### Software Dependencies
- Operating System:
  - Linux (tested on Ubuntu 20.04 LTS, CentOS 7)
- Python 3.10
- NVIDIA GPU with 8GB+ VRAM (required for cuOPT QP solver, Module 3)
- Required Python packages: see `requirements.txt` for exact versions
- cuOPT: **not on PyPI and not in `CITEgeist_env.yml`** — installed separately from the
  NVIDIA index or a RAPIDS container. There is no CPU fallback for Module 3.

### Hardware Requirements
- RAM: Minimum 16GB, Recommended 64GB+
- Storage: 10GB minimum for installation and basic analysis
- CPU: Multi-core processor recommended (8+ cores for optimal performance)
- GPU: NVIDIA GPU with 8GB+ VRAM required for Module 3 (cuOPT QP)

## Installation Guide

### Quick install (PyPI)

```bash
pip install citegeist
```

The `cuOPT` GPU solver (Module 3) is **not on PyPI and is not included in any environment
file in this repository** — install it separately via NVIDIA RAPIDS/NGC
(`pip install cuopt-cu12 --extra-index-url https://pypi.nvidia.com`, or use a pre-built
RAPIDS container). CPU-only environments **cannot run** the cell-type proportion solver.

### From source

#### Step 1: Set up Python Environment

```bash
# Install the analysis environment (scientific-Python stack; cuOPT NOT included)
conda env create -f CITEgeist_env.yml
conda activate CITEgeist_env

# Or install core dependencies manually (see requirements.txt for versions)
pip install -r requirements.txt

# Then add the GPU solver separately (required for Module 3):
pip install cuopt-cu12 --extra-index-url https://pypi.nvidia.com
```

#### Step 2: Install CITEgeist
```bash
git clone https://github.com/leeoesterreich/CITEgeist.git
cd CITEgeist
pip install -e .
```

Typical installation time: 10-15 minutes on a standard desktop computer

## Demo

See `examples/scripts` for runnable end-to-end pipeline examples demonstrating:
1. Biopsy heterogeneity analysis
2. ESR1 D538G mutant spatial profiling
3. Macrophage treatment response

### Python API

```python
from CITEgeist.model.citegeist_model import CitegeistModel

model = CitegeistModel(
    sample_name="my_sample",
    adata=adata,
    output_folder="output/my_sample",
    simulation=False,
)
model.split_adata()  # Required for real patient data
```

See `CITEgeist/README.md` for full API documentation and SLURM examples.

## Instructions for Use

For detailed usage instructions, API documentation, and SLURM batch processing examples, see `CITEgeist/README.md`.

## Reproduction Instructions

End-to-end pipeline examples are in [`examples/scripts`](../../examples/scripts).
For a full walkthrough on real Visium data, see [`docs/quickstart_real_visium.md`](../../docs/quickstart_real_visium.md).
