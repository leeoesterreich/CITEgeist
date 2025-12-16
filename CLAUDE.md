# CLAUDE.md - AI Assistant Guide for CITEgeist

## Project Overview

**CITEgeist** (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends) is a computational framework for analyzing spatial multi-omic data, with a focus on integrating CITE-seq and spatial transcriptomics.

### Purpose
CITEgeist performs deconvolution of spatial transcriptomics data using spatially-resolved CITE-seq measurements, enabling:
- Cell-type proportion estimation
- Gene expression deconvolution
- Spatial pattern analysis in tissue architecture
- Analysis of tumor microenvironments

### Key Context
This is an active research project with an associated paper currently in preparation (bioRxiv preprint available). The codebase is under active development for publication.

## Repository Structure

```
CITEgeist/
├── CITEgeist/                    # Main package directory
│   ├── model/                    # Core model implementation
│   │   ├── citegeist_model.py   # Main CitegeistModel class
│   │   ├── gurobi_impl.py       # Gurobi optimization routines
│   │   ├── utils.py             # Utility functions
│   │   ├── checkpoints.py       # Checkpoint management
│   │   └── analysis_functions.py # Analysis helper functions
│   ├── examples/                 # Example scripts and vignettes
│   │   ├── vignette_1_biopsy_heterogeneity.ipynb
│   │   ├── vignette_2_surgical_d538g.ipynb
│   │   ├── vignette_3_responder_macrophages.ipynb
│   │   └── sbatch_sample.sh     # SLURM batch script template
│   └── tests/                    # Test scripts and validation
│       ├── test_citegeist_simulated.py
│       ├── visualize_prop_gridsearch.py
│       └── visualize_simulated_benchmarks.py
├── Benchmarking/                 # Benchmarking framework
│   ├── benchmarking/             # Benchmark implementations
│   │   ├── src/                  # Benchmark wrapper scripts
│   │   ├── CITEgeist/           # CITEgeist benchmark results
│   │   ├── Cell2Location/       # Cell2Location comparisons
│   │   ├── RCTD/                # RCTD comparisons
│   │   └── Tangram/             # Tangram comparisons
│   └── simulation/               # Simulation framework
├── .github/                      # GitHub configuration
│   └── workflows/
│       └── pylint.yml           # Python linting workflow
├── CITEgeist_env.yml            # Conda environment specification
├── README.md                     # Main project documentation
└── LICENSE                       # BSD 3-Clause License
```

## Tech Stack

### Core Language & Version
- **Python 3.10** (strictly required)

### Key Dependencies

#### Mathematical Optimization
- **gurobipy 11.0.2** - Core optimization engine (requires academic/commercial license)

#### Scientific Computing
- **numpy 1.26.4** - Numerical computations
- **scipy 1.13.1** - Scientific algorithms
- **pandas 2.2.3** - Data manipulation

#### Bioinformatics
- **scanpy 1.10.4** - Single-cell analysis
- **anndata 0.11.3** - Annotated data structures
- **squidpy 1.6.2** - Spatial transcriptomics
- **spatialdata 0.2.5.post0** - Spatial data structures

#### Machine Learning
- **scikit-learn 1.6.1** - ML utilities

#### Visualization
- **matplotlib 3.10.0** - Plotting
- **seaborn 0.13.2** - Statistical visualization

#### Additional Tools
- **h5py 3.12.1** - HDF5 file handling
- **pyarrow** - Parquet file support
- **networkx** - Graph algorithms

### Environment Setup

```bash
# Create conda environment from specification
conda env create -f CITEgeist_env.yml
conda activate CITEgeist_env

# Or install via pip
pip install citegeist

# Development installation
git clone https://github.com/leeoesterreich/CITEgeist.git
cd CITEgeist
pip install -e .[dev]
```

### Gurobi License Requirement
CITEgeist requires a valid Gurobi license (free for academic use):
- Academic license: https://www.gurobi.com/downloads/end-user-license-agreement-academic/
- License file location must be configured in code or environment

## Development Workflow

### Git Branching Strategy

**CRITICAL: Branch Protection Rules**

1. **`main` branch is PROTECTED**
   - Main branch should remain UNTOUCHED until the paper is accepted and published
   - NO direct commits to main
   - NO merges to main until publication approval

2. **`dev` branch is the active development branch**
   - ALL feature branches should branch from `dev`
   - ALL pull requests should target `dev`
   - This is where active development happens

3. **Feature branches**
   - Branch naming: `feature/<descriptive-name>` or `claude/<session-id>`
   - Always branch from `dev`, not `main`
   - Keep branches focused and short-lived

### Typical Workflow

```bash
# Start new feature
git checkout dev
git pull origin dev
git checkout -b feature/my-feature

# Make changes and commit
git add .
git commit -m "descriptive message"

# Push and create PR to dev
git push -u origin feature/my-feature
# Create PR targeting dev branch
```

### Commit Message Conventions
- Use clear, descriptive commit messages
- Prefix with type: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`
- Example: `feat: add spatial neighbor detection with fixed radius`
- Example: `fix: correct cell proportion normalization in Pass 1`

## Code Conventions

### Python Style
- Follow PEP 8 style guidelines
- Pylint is configured via GitHub Actions (runs on push)
- Target Python versions in CI: 3.8, 3.9, 3.10

### Import Organization
```python
# Standard library imports
import os
import logging
from typing import Dict, Any, Optional

# Third-party imports
import numpy as np
import pandas as pd
import scanpy as sc

# Local imports
from .gurobi_impl import optimize_cell_proportions
from .utils import validate_cell_profile_dict
```

### Docstrings
Use Google-style docstrings:
```python
def function_name(arg1, arg2):
    """
    Brief description.

    Args:
        arg1 (type): Description.
        arg2 (type): Description.

    Returns:
        type: Description.
    """
```

### Type Hints
Use type hints where applicable:
```python
def process_data(adata: sc.AnnData, radius: int = 4) -> Dict[str, np.ndarray]:
    ...
```

### Memory Management
- Use `cleanup_memory()` utility for large data processing
- Be mindful of AnnData object memory footprint
- Consider checkpointing for long-running computations

## Key Model Components

### CitegeistModel Class
Primary interface in `CITEgeist/model/citegeist_model.py`

**Initialization Modes:**
1. **Standard mode**: Single AnnData object with both GEX and antibody data
2. **Simulation mode**: Separate gene_expression_adata and antibody_capture_adata

**Key Methods:**
- `preprocess_gene_expression()` - Normalizes and filters gene expression data
- `preprocess_antibody()` - Processes antibody capture data
- `map_antibodies_to_cell_profiles()` - Maps antibodies to cell type profiles
- `optimize_proportions()` - Pass 1: Cell type proportion optimization
- `optimize_gene_expression()` - Pass 2: Gene expression deconvolution
- `export_results()` - Saves results to output folder

### Two-Pass Approach
1. **Pass 1**: Cell-type proportion estimation using antibody data
2. **Pass 2**: Gene expression deconvolution using proportions from Pass 1

### Gurobi Implementation
Located in `gurobi_impl.py`:
- `optimize_cell_proportions()` - Quadratic programming for proportions
- `optimize_gene_expression()` - Gene expression optimization
- `compute_global_prior()` - Computes spatial priors for regularization

## Testing & Validation

### Test Structure
Located in `CITEgeist/tests/`:
- `test_citegeist_simulated.py` - Simulated data validation
- `visualize_prop_gridsearch.py` - Hyperparameter visualization
- `visualize_simulated_benchmarks.py` - Benchmark result visualization

### Running Tests
```bash
cd CITEgeist/tests
python test_citegeist_simulated.py
```

### Benchmarking
Benchmarking framework in `Benchmarking/` directory:
- Compares CITEgeist against Cell2Location, RCTD, Tangram, Seurat
- Metrics: JSD, RMSE, MAE, Pearson correlation
- Separate evaluations for cell proportions and gene expression

## Working with Data

### Input Data Requirements
- **Spatial transcriptomics data**: 10x Visium or similar
- **CITE-seq/antibody capture**: Spatially resolved protein measurements
- **Data format**: AnnData objects (.h5ad files)

### Real Patient Data vs Simulated Data Loading

**CRITICAL DISTINCTION**: Real patient data and simulated data are loaded differently.

#### Real Patient Data (from 10x SpaceRanger)

Uses a **combined AnnData object** loaded from SpaceRanger outputs, then split:

```python
import squidpy as sq
from model import CitegeistModel

# Load from SpaceRanger outs/ folder (contains both GEX and antibody data)
adata = sq.read.visium(
    path_to_visium_folder,  # e.g., "/path/to/sample/outs/"
    counts_file='filtered_feature_bc_matrix.h5',
    load_images=True,
    gex_only=False  # CRITICAL: Must be False to load antibody data too
)

# Initialize model with combined adata
model = CitegeistModel(sample_name="HCC22-088-P1-S2", adata=adata, output_folder='output')

# REQUIRED: Split into separate GEX and antibody AnnData objects
model.split_adata()  # Splits by feature_types column

# Then continue with preprocessing...
model.preprocess_gex()
model.preprocess_antibody()
```

**Key notes for real patient data:**
- Antibody names have "-1" suffix (e.g., "EPCAM-1", "CD68-1", "CD3E-1")
- Data directory: `/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/`
- `min_counts=100` for high-quality biopsy samples, `min_counts=25` for surgical samples

#### Simulated Data (pre-separated files)

Uses **pre-separated** gene expression and antibody AnnData objects:

```python
import scanpy as sc
from model import CitegeistModel

# Load two separate h5ad files
adata_gex = sc.read_h5ad("h5ad_objects/Wu_rep_1_GEX.h5ad")
adata_cite = sc.read_h5ad("h5ad_objects/Wu_rep_1_CITE.h5ad")

# Initialize with simulation=True flag
model = CitegeistModel(
    sample_name="Wu_rep_1",
    output_folder='output',
    simulation=True,  # CRITICAL: Must be True for simulated data
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite
)

# NO split_adata() needed - data is already separate
model.preprocess_gex()
model.preprocess_antibody()
```

**Key notes for simulated data:**
- Location: `replicates/high_seg/h5ad_objects/Wu_rep_X_CITE.h5ad` and `Wu_rep_X_GEX.h5ad`
- Marker names don't have "-1" suffix (e.g., "B-cells_Protein_1", "T-cells_Protein_1")
- Ground truth: 18 cell-type specific proteins vs 82 "Nonspecific_Protein_X" markers

#### Quick Reference Table

| Aspect | Real Patient Data | Simulated Data |
|--------|------------------|----------------|
| **Data Source** | `sq.read.visium()` from SpaceRanger | Separate `_GEX.h5ad` and `_CITE.h5ad` |
| **gex_only** | `False` | N/A |
| **split_adata()** | **Required** | Not needed |
| **simulation flag** | `False` (default) | `True` |
| **Antibody suffix** | "-1" (e.g., "CD68-1") | None (e.g., "B-cells_Protein_1") |

### AnnData Structure Expectations
```python
adata.X                    # Gene expression matrix
adata.obs                  # Observation metadata (spots)
adata.var                  # Variable metadata (genes)
adata.obsm['spatial']      # Spatial coordinates
adata.layers['antibody']   # Antibody capture data (if combined)
```

### Output Files
CITEgeist generates:
- Cell proportion matrices
- Deconvolved gene expression per cell type
- Spatial neighbor graphs
- Checkpoint files for resume capability
- Log files for debugging

## Common Tasks for AI Assistants

### Adding New Features
1. Branch from `dev`, NOT `main`
2. Review existing model structure in `citegeist_model.py`
3. Check if Gurobi optimization changes are needed
4. Update utilities if adding preprocessing steps
5. Add tests for new functionality
6. Update relevant docstrings

### Fixing Bugs
1. Check logs in output folder for error traces
2. Review checkpoint files if optimization failed
3. Verify Gurobi license status if optimization errors occur
4. Check memory usage for large datasets
5. Validate input data structure matches expectations

### Optimizing Performance
1. Profile with memory and time checkpoints
2. Consider parallelization opportunities (`max_workers` parameter)
3. Review Gurobi solver parameters
4. Check checkpoint interval for long runs
5. Consider data subsetting for testing

### Updating Documentation
1. Update relevant README files
2. Keep docstrings synchronized with code
3. Update vignettes if API changes
4. Maintain this CLAUDE.md file for significant structural changes

## Important Files & Their Purposes

### Core Model Files
- **citegeist_model.py**: Main model class, orchestrates workflow
- **gurobi_impl.py**: All Gurobi optimization routines
- **utils.py**: Shared utilities (logging, validation, neighbor detection)
- **checkpoints.py**: Save/load model state for long runs
- **analysis_functions.py**: Post-processing and analysis helpers

### Configuration Files
- **CITEgeist_env.yml**: Complete conda environment specification
- **.gitignore**: Excludes .h5ad, .h5seurat, model files, images
- **LICENSE**: BSD 3-Clause (all derivatives must be open source)

### CI/CD
- **.github/workflows/pylint.yml**: Automated linting on push

## Debugging Tips

### Common Issues

1. **Gurobi License Errors**
   - Verify license file path
   - Check license expiration
   - Ensure license matches machine/user

2. **Memory Issues**
   - Use checkpoints to resume
   - Reduce `max_workers` parameter
   - Process smaller spatial regions
   - Clear intermediate results

3. **Optimization Failures**
   - Check input data for NaN/Inf values
   - Verify cell profile dictionary is valid
   - Review lambda_reg and alpha_elastic parameters
   - Check neighbor detection radius

4. **Import Errors**
   - Ensure correct Python 3.10 environment
   - Verify all dependencies from CITEgeist_env.yml
   - Check PYTHONPATH if using development installation

## Performance Considerations

### Runtime Expectations
- Small dataset: ~2 hours (16 threads, 32GB RAM)
- Medium dataset: ~4 hours
- Large dataset: ~10 hours

### Hardware Recommendations
- **Minimum**: 16GB RAM, 4 cores
- **Recommended**: 64GB+ RAM, 8+ cores
- **Storage**: 16GB minimum for installation + analysis

### Scalability
- Use SLURM distribution for HPC clusters (see `examples/sbatch_sample.sh`)
- Parallelize across samples, not within samples
- Consider checkpoint_interval for very large analyses

## API Stability Notes

This is an active research codebase:
- API may change before paper publication
- Breaking changes possible in dev branch
- Pin versions for reproducible analyses
- Check commit history for recent changes

## License & Citation

### License
BSD 3-Clause License - All modifications and derivatives must be open source

### Citation
If using CITEgeist in research, cite the bioRxiv preprint:
```bibtex
@article{ChangSchlegelCITEgeistCellularIndexing2025,
  title = {{{CITEgeist}}: {{Cellular Indexing}} of {{Transcriptomes}} and {{Epitopes}}
           for {{Guided Exploration}} of {{Intrinsic Spatial Trends}}},
  author = {Chang, Alexander Chih-Chieh and Schlegel, Brent T. and
            Carleton, Neil and McAulife, Priscilla F. and
            Oesterreich, Steffi and Schwartz, Russell and Lee, Adrian V.},
  date = {2025-02-17},
  eprinttype = {bioRxiv},
  doi = {10.1101/2025.02.15.638331}
}
```

## Contact & Support

- **Lab Website**: https://leeoesterreich.org/
- **Issues**: https://github.com/leeoesterreich/CITEgeist/issues
- **Email**: alc376@pitt.edu

## AI Assistant Guidelines

### When Working on This Codebase:

1. **Always verify you're on the correct branch** (`dev`, not `main`)
2. **Understand the two-pass algorithm** before modifying optimization code
3. **Test with simulated data** before running on real datasets
4. **Check Gurobi compatibility** when modifying optimization routines
5. **Maintain backward compatibility** where possible (paper reproducibility)
6. **Document parameter changes** that affect results
7. **Consider memory implications** of changes to data structures
8. **Update checkpointing** if adding new model state
9. **Preserve scientific accuracy** - this is research code with biological interpretation
10. **Respect the publication timeline** - main branch is frozen until publication

### Before Committing:

- [ ] Code follows PEP 8 style
- [ ] Docstrings are complete and accurate
- [ ] No hardcoded paths or credentials
- [ ] Memory cleanup is appropriate
- [ ] Error handling is robust
- [ ] Logging is informative
- [ ] Changes are tested (at least manually)
- [ ] Targeting `dev` branch, not `main`

### Asking for Clarification:

If uncertain about:
- Scientific/biological interpretation → Ask user
- Expected algorithm behavior → Check paper/vignettes
- Parameter sensitivity → Check benchmarking results
- Data structure requirements → Check test_citegeist_simulated.py

---

**Last Updated**: 2025-11-14
**Repository**: https://github.com/leeoesterreich/CITEgeist
**Status**: Active Development (Pre-Publication)
