# Installation

## Prerequisites

- Conda or Mamba
- NVIDIA GPU (required for Module 3 cuOPT QP solver)

## Environment Setup

Create the production environment (note: name is case-sensitive):

```bash
conda env create -f requirements.txt -n CITEgeist_env
conda activate CITEgeist_env
pip install -e .
```

## HPC / SLURM

CITEgeist is designed for HPC clusters using SLURM. Key requirements:

- **GPU nodes required for Module 3 (QP deconvolution).** The cuOPT solver
  fails silently on CPU nodes. Submit with `--gres=gpu:1`.
- On PITT CRC, use `gpu-race-submit.sh <script.sbatch>` to race across
  `l40s`, `a100`, and `a100_nvlink` partitions.
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
- Reference scRNA-seq: AnnData `.h5ad` with cell type annotations
- See `examples/scripts/run_module3_unified.py` for a complete single-sample workflow
