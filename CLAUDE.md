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

**CRITICAL: Before running any Gurobi-dependent code (tests, optimization, benchmarks), you MUST load the Gurobi module:**

```bash
module load gurobi/12.0.3
```

If the version is unavailable, check available versions with:
```bash
module avail gurobi
```

This applies to:
- Running tests (`pytest`)
- Running CITEgeist optimization (Module 3)
- Any script that imports `gurobipy`
- SLURM batch scripts (add to script header)

For SLURM scripts, add this line after the `#SBATCH` directives:
```bash
module load gurobi/12.0.3
```

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

## CITEgeist Pipeline: Module Definitions (1-5)

CITEgeist is organized into 5 modules that form a complete spatial multi-omics analysis pipeline:

### Module 1: Marker Interest Detection (Auto Protein Marker Prediction)

**File**: `CITEgeist/model/marker_interest.py`

**Purpose**: Identifies spatially-variable protein markers from CITE-seq antibody data worth analyzing further.

**Method**: Uses three statistical tests to score marker "interestingness":
- **Kurtosis gating**: Peaked distributions (kurtosis > 2.0) indicate real signal vs flat noise
- **Gaussian Mixture Model (GMM)**: 2-component GMM separates signal from background, calculates SNR
- **Moran's I**: Spatial autocorrelation on raw data to detect spatially clustered markers

**Output**: `MarkerInterestResult` with ranked markers and decision logic
- `interesting_markers`: Pass EITHER kurtosis OR Moran's I gate (plus GMM filter)
- `boring_markers`: Fail both gates or GMM filter

**Test**: `tests/test_module_one_marker_detection.py`

---

### Module 2: Profile Assembly (Cell Type Marker Profile Discovery)

**File**: `CITEgeist/model/spatial_colocalization.py`

**Purpose**: Automatically discovers cell type protein marker profiles from spatial colocalization patterns.

#### Module 2a: Marker Colocalization Analysis
- Analyzes spatial relationships: same-spot co-occurrence, expression correlation, adjacent-spot enrichment, bivariate Moran's I
- Output: `ColocalizationResult` with scored marker pairs

#### Module 2b: Profile Discovery
- Builds significance-filtered graph from colocalization p-values
- Finds connected components (separate lineages)
- Hierarchical clustering + dynamic tree cutting to extract profiles
- Output: `ProfileDiscoveryResult` with discovered profiles and dendrogram

#### Module 2c: Profile Selection
- Function: `select_profiles_by_reconstruction()`
- Uses reconstruction accuracy to determine optimal number of profiles
- Ensures selected profiles have good coverage and non-redundancy

**Tests**: `tests/test_module2_profile_discovery.py`, `tests/test_module2c_profile_selection.py`

---

### Module 3: Cell Type Proportion & Gene Expression Deconvolution

**Files**: `CITEgeist/model/gurobi_impl.py` (optimization), `CITEgeist/model/citegeist_model.py` (orchestration)

**Purpose**: Two-pass optimization for deconvoluting spatial transcriptomics using CITE-seq antibody references.

#### Pass 1: Cell-type Proportion Estimation (`run_cell_proportion_model()`)
- Maps antibodies to cell type profiles (from Module 2)
- Solves quadratic programming: minimize reconstruction error
- Features: per-marker beta optimization, Laplacian spatial smoothing, validation thresholds
- Finetuning: Optimizes proportions locally using neighbors
- Output: `cell_prop_global_results.csv`, `cell_prop_finetuned_results.csv`

#### Pass 2: Gene Expression Deconvolution (`run_cell_expression_pass1()`)
- Uses proportions from Pass 1 to deconvolve gene expression
- Minimizes reconstruction error with enrichment weights (global + local)
- Produces deconvolved gene expression layers (one per cell type per spot)
- Output: `gene_expression_pass1.parquet` + AnnData layers

**Key Functions**:
- `optimize_cell_proportions()` / `optimize_cell_proportions_per_marker()`: Pass 1
- `optimize_gene_expression()`: Pass 2
- `compute_global_prior()`: Creates expression prior

**Test**: `tests/test_citegeist_simulated.py`

#### Module 3 Alternative: Discrete Cell Assignment (IQP)

**Files**: `CITEgeist/model/gurobi_impl.py`, `CITEgeist/model/citegeist_model.py`, `CITEgeist/model/segmentation.py`

**Purpose**: When nuclei counts are available (from StarDist segmentation or Xenium cell mapping), assign discrete integer cell counts per type instead of continuous proportions.

**Method**:
- **E-step**: Integer Quadratic Programming (IQP) solver assigns discrete cell identities
- **M-step**: Ordinary Least Squares (OLS) optimizes per-marker beta coefficients
- Iterates until beta convergence or max iterations reached

**Key Functions**:
- `solve_discrete_cell_counts()`: Single-spot IQP optimization with Gurobi
- `optimize_discrete_cell_assignment_em()`: Full EM algorithm over all spots
- `run_discrete_cell_assignment()`: CitegeistModel entry point (replaces `run_cell_proportion_model()`)
- `preprocess_antibody_discrete()`: CLR normalization preserving cellularity signal

**Benchmark Results (Xenium pseudo-Visium, achievable-7 cell types)**:
| Metric | Discrete | Continuous Baseline | Improvement |
|--------|----------|---------------------|-------------|
| Proportion Pearson | 0.658 | 0.600 | +10% |
| Proportion RMSE | 0.128 | 0.167 | -23% |
| GEX Pearson | 0.407 | 0.35-0.37 | +10-16% |

**Usage**:
```python
# Requires nuclei counts (from StarDist or Xenium n_cells)
model.preprocess_antibody_discrete()  # NOT preprocess_antibody()
cell_counts_df = model.run_discrete_cell_assignment(nuclei_counts=nuclei_series)
model.run_cell_expression_pass1(use_discrete_mode=True, cell_counts=cell_counts_df)
```

**Tests**: `tests/test_discrete_assignment.py`
**Example**: `examples/run_module3_discrete.py`

---

### Module 4: Protein-Anchored Spatial Transcriptomic Program Discovery

**File**: `CITEgeist/model/anchored_program_discovery.py`

**Purpose**: Discovers gene expression programs from deconvolved cell-type-specific expression layers (Module 3 output), validated against protein profiles (Module 2).

**Method**:
1. For each cell type anchor (from Module 2 profiles):
   - Extract deconvolved expression layer from Module 3
   - Weight by cell type proportions
   - Non-negative matrix factorization (NMF) to discover K programs
   - Compute Moran's I for spatial coherence
   - Identify spatial subpopulations via k-means clustering

**Output**: `AnchoredProgramDiscoveryResult` containing:
- `SpatialProgram`: Top genes, loadings, variance explained, Moran's I
- `SpatialSubpopulation`: Spatially distinct subpopulations within each anchor
- Gene loadings (W matrix) and spot loadings (H matrix)

#### Module 4b: Bivariate Program Relationships
- Function: `analyze_program_relationships()`
- Computes bivariate Moran's I between pairs of programs
- Identifies co-localized (positive I) vs exclusive (negative I) programs
- Output: `BivariateProgramResult`

**Test**: `tests/test_module4_deconvolved.py`

---

### Module 5: Cross-Sample Integration

**File**: `CITEgeist/model/cross_sample_integration.py`

**Purpose**: Integrates gene expression programs across multiple patient samples, producing a similarity network of conserved programs and relationships.

**Method**:
1. Per-sample analysis from Module 4 (W matrices, H matrices, programs)
2. Module 4b provides bivariate Moran's I within each sample
3. Cross-sample integration:
   - Aligns programs using Harmony-style batch correction
   - PCA embedding in integrated latent space
   - Hierarchical clustering to match equivalent programs across patients
4. Relationship conservation:
   - Compares program pairs for consistent relationships across samples
   - Builds similarity network of aligned programs

**Output**: `IntegrationResult` containing:
- `AlignedProgram`: Programs aligned across samples with consensus signatures
- `ConservedRelationship`: Bivariate relationships that persist across patients
- Conservation scores and relationship types ('co-localized', 'exclusive', 'variable')

**Test**: `tests/test_cross_sample_integration.py`

---

### Module Summary Table

| Module | Name | Input | Output | Key File |
|--------|------|-------|--------|----------|
| **1** | Marker Interest Detection | Raw antibody data | Ranked interesting markers | `marker_interest.py` |
| **2a** | Marker Colocalization | Interesting markers | Scored marker pairs | `spatial_colocalization.py` |
| **2b** | Profile Discovery | Marker pairs | Cell type profiles | `spatial_colocalization.py` |
| **2c** | Profile Selection | Multiple profiles | Optimal profiles | `spatial_colocalization.py` |
| **3** | Deconvolution (2-Pass) | GEX + antibody profiles | Proportions + deconvolved layers | `gurobi_impl.py`, `citegeist_model.py` |
| **3-alt** | Discrete Cell Assignment | GEX + antibody + nuclei counts | Integer cell counts + deconvolved layers | `gurobi_impl.py`, `segmentation.py` |
| **4** | Anchored Programs | Deconvolved layers | Spatial programs per cell type | `anchored_program_discovery.py` |
| **4b** | Bivariate Relationships | Program H matrices | Program-pair correlations | `anchored_program_discovery.py` |
| **5** | Cross-sample Integration | Multi-sample programs | Conserved relationships | `cross_sample_integration.py` |

---

## Key Model Components

### CitegeistModel Class
Primary interface in `CITEgeist/model/citegeist_model.py`

**Initialization Modes:**
1. **Standard mode**: Single AnnData object with both GEX and antibody data
2. **Simulation mode**: Separate gene_expression_adata and antibody_capture_adata

**Key Methods:**
- `preprocess_gene_expression()` - Normalizes and filters gene expression data
- `preprocess_antibody()` - Processes antibody capture data
- `load_cell_profile_dict()` - Loads cell type to marker mappings (or use Module 2 for auto-discovery)
- `run_cell_proportion_model()` - Module 3 Pass 1: Cell type proportion optimization
- `run_cell_expression_pass1()` - Module 3 Pass 2: Gene expression deconvolution

### Gurobi Implementation
Located in `gurobi_impl.py`:
- `optimize_cell_proportions()` - Quadratic programming for proportions
- `optimize_cell_proportions_per_marker()` - Per-marker beta optimization variant
- `optimize_gene_expression()` - Gene expression optimization
- `compute_global_prior()` - Computes spatial priors for regularization

## Testing & Validation

### Test Structure
Located in `CITEgeist/tests/`:
- `test_citegeist_simulated.py` - Simulated data validation
- `visualize_prop_gridsearch.py` - Hyperparameter visualization
- `visualize_simulated_benchmarks.py` - Benchmark result visualization

### Running Tests

**IMPORTANT:** Load Gurobi module before running tests:
```bash
module load gurobi/12.0.3  # Run `module avail gurobi` if version unavailable
eval "$(conda shell.bash hook)" && conda activate ~/alc376_bgfs/envs/CITEgeist_env
cd CITEgeist/tests
python test_citegeist_simulated.py
# Or use pytest:
pytest tests/ -v
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
11. **Use sbatch for big jobs** - For computationally intensive tasks (loading large datasets, running benchmarks, etc.), submit via SLURM sbatch rather than running interactively. This prevents timeouts and allows proper resource allocation.

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

## Current Progress Summary (January 2026)

### Module Implementation Status

| Module | Name | Status |
|--------|------|--------|
| **1** | Marker Interest Detection | COMPLETE |
| **2a** | Marker Colocalization | COMPLETE |
| **2b** | Profile Discovery | COMPLETE |
| **2c** | Profile Selection | COMPLETE |
| **3** | Deconvolution (Proportions + GEX) | COMPLETE |
| **4** | Anchored Program Discovery | COMPLETE |
| **4b** | Bivariate Relationships | COMPLETE |
| **5** | Cross-Sample Integration | COMPLETE |

**All 5 core modules are fully implemented with comprehensive test coverage.**

### Xenium Benchmarking Status

The multi-method benchmarking framework compares CITEgeist against Cell2Location, Tangram, RCTD, and Seurat on 14 Xenium pseudo-Visium regions with 10 granular RNA-based cell type ground truth.

**Directory**: `Benchmarking/xenium_benchmarking/{CITEgeist,Cell2Location,RCTD,Seurat,Tangram,evaluation}/`

### Recent Cleanup (January 2026)

1. Removed 60 obsolete top-level benchmark files (refactored to method-specific directories)
2. Added SLURM mail directives to all 11 benchmark scripts
3. Fixed evaluation script for Seurat/Tangram output normalization

For full details, see: `docs/PROGRESS_REPORT.md`

---

**Last Updated**: 2026-01-09
**Repository**: https://github.com/leeoesterreich/CITEgeist
**Status**: Active Development (Pre-Publication)
