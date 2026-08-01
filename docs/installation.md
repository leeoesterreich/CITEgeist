# Installation

## Quick install (PyPI)

```bash
pip install citegeist
```

The `cuOPT` GPU solver (Module 3) is **not on PyPI and is not included in any environment
file in this repository** — install it separately via NVIDIA RAPIDS/NGC
(`pip install cuopt-cu12 --extra-index-url https://pypi.nvidia.com`, or use a pre-built
RAPIDS container). There is **no CPU fallback**: CPU-only environments cannot run the
cell-type proportion solver.

The optional StarDist nucleus segmentation used by single-cell assignment (`stardist`,
`csbdeep`) is likewise not in `CITEgeist_env.yml` and must be installed separately.

## From source (HPC / dev environment)

### Prerequisites

- Conda or Mamba
- NVIDIA GPU (required for Module 3 cuOPT QP solver)

### Environment Setup

Create the analysis environment (note: name is case-sensitive). This environment
provides the scientific-Python stack only; the GPU stack (cuOPT) and the optional
segmentation stack (StarDist) are installed separately, as above.

```bash
conda env create -f CITEgeist_env.yml -n CITEgeist_env
conda activate CITEgeist_env
pip install -e .
```

## HPC / SLURM

CITEgeist is designed for HPC clusters using SLURM. Key requirements:

- **GPU nodes required for Module 3 (QP deconvolution).** The cuOPT solver
  fails silently on CPU nodes. Submit with `--gres=gpu:1`.
- Module 3b (StarDist segmentation) requires `HDF5_USE_FILE_LOCKING=FALSE`
  for array jobs on NFS.

## Verifying Installation

```bash
conda activate CITEgeist_env
pytest tests/ -m "not slow and not requires_data and not requires_cuopt" -q
```

All tests should pass. Tests marked `requires_cuopt` need a GPU node.

## Data Requirements

- Patient data: SpaceRanger output directory (GEX + antibody TIFF)
- No reference scRNA-seq dataset is required — cell-type profiles are defined from the
  antibody panel (see `docs/quickstart_real_visium.md`)
- See `examples/scripts/run_module3_unified.py` for a complete single-sample workflow
