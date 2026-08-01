# CITEgeist — Production Example Scripts

These scripts illustrate the full CITEgeist analysis pipeline on Visium + CITE-seq data.
They are illustrative entry points; you must supply your own SpaceRanger output and CITE-seq
antibody count matrix. The QP deconvolution step (Module 3 cuOPT) requires an NVIDIA GPU
with the cuOPT library installed.

| Script | Module / Step | Description |
|--------|--------------|-------------|
| `run_module12_discovery.py` | Modules 1 & 2 | Marker interest detection (M1) and spatial colocalization → cell-type profile discovery (M2a/b/c) |
| `run_module3_unified.py` | Module 3 | Unified deconvolution: cuOPT QP cell-type proportions + SACE GEX allocation |
| `run_cuopt_qp_patient.py` | Module 3 (QP only) | Standalone cuOPT quadratic-programming proportion solver for a single patient sample |
| `run_module4_discovery.py` | Module 4 | Anchored NMF gene-program discovery per cell type |
| `run_module4b_bivariate.py` | Module 4b | Bivariate Moran's I between program pairs (spatial co-variation) |
| `run_module5_integration.py` | Module 5 | Cross-sample program integration and conservation scoring |
| `run_single_cell_assignment.py` | Module 3-post | Per-cell type assignment from spot proportions (Bayesian or Hungarian) |

## Requirements

- SpaceRanger output directory (Visium GEX) and matched CITE-seq antibody counts (`.h5ad`)
- Python environment with CITEgeist installed (`pip install citegeist`, or from source —
  `git clone` the repo, then `pip install -e .`)
- **GPU required** for `run_cuopt_qp_patient.py` and `run_module3_unified.py` (cuOPT QP solver)
- See each script's header for dataset-specific parameters to configure before running
