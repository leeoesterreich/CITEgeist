# CITEgeist Repository: Comprehensive Codebase Overview

## 1. PROJECT DESCRIPTION

### What is CITEgeist?
**CITEgeist** (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends) is a comprehensive computational framework for analyzing spatial multi-omic data.

**Key Innovation:** 
- Deconvolves spatial transcriptomics data using spatially-resolved CITE-seq measurements
- Performs both **cell-type proportion estimation** and **gene expression deconvolution** in a two-pass approach
- Leverages both protein (antibody capture) and RNA measurements from the same spatial locations
- Reference-independent alternative to traditional scRNA-seq based deconvolution

**Biological Context:**
- Integrates CITE-seq (protein) and spatial transcriptomics (RNA) data
- Enables spatial pattern analysis in tissues
- Provides robust cell-type deconvolution
- Supports flexible regularization options
- Scales to large datasets

**Citation:**
- Paper: Chang et al., 2025 - bioRxiv (10.1101/2025.02.15.638331)
- Focus: Breast cancer tissue analysis, translational pre/post-treatment samples

---

## 2. PROJECT STRUCTURE

### Directory Organization
```
/home/user/CITEgeist/
├── CITEgeist/                          # Main package
│   ├── model/                          # Core computational modules (2,667 lines)
│   │   ├── __init__.py                # Exports public API
│   │   ├── citegeist_model.py          # Main CitegeistModel class (834 lines)
│   │   ├── gurobi_impl.py              # Optimization algorithms (1,219 lines) [CRITICAL]
│   │   ├── checkpoints.py              # Checkpoint management (173 lines)
│   │   ├── utils.py                    # Utility functions (378 lines)
│   │   └── analysis_functions.py       # Analysis helpers (42 lines)
│   ├── examples/                       # Tutorial notebooks and scripts
│   │   ├── vignette_1_biopsy_heterogeneity.ipynb
│   │   ├── vignette_2_surgical_d538g.ipynb
│   │   ├── vignette_3_responder_macrophages.ipynb
│   │   ├── compute_sample.py
│   │   └── sbatch_sample.sh            # HPC job submission template
│   ├── tests/                          # Test and benchmark code
│   │   ├── test_citegeist_simulated.py # Main test script (380+ lines)
│   │   ├── visualize_prop_gridsearch.py
│   │   └── visualize_simulated_benchmarks.py
│   └── README.md                       # Package documentation
├── Benchmarking/                       # Benchmarking framework
│   ├── README.md                       # Benchmarking methodology
│   ├── benchmarking/                   # Method comparisons
│   │   ├── src/                        # Core benchmarking modules
│   │   │   ├── benchmarking_gex.py           # Gene expression metrics
│   │   │   ├── benchmarking_spot_deconv.py  # Cell type proportion metrics
│   │   │   ├── citegeist_bench_wrapper.py   # CITEgeist benchmarking wrapper
│   │   │   ├── cell2location_bench_wrapper.py
│   │   │   ├── tangram_bench_wrapper.py
│   │   │   ├── RCTD_bench_wrapper.py
│   │   │   └── seurat_bench_wrapper.py
│   │   ├── Cell2Location/              # Benchmarking scripts for Cell2Location
│   │   ├── Tangram/                    # Benchmarking scripts for Tangram
│   │   └── Seurat/                     # Benchmarking scripts for Seurat
│   └── simulation/                     # Synthetic data generation
│       └── README.md                   # Simulation methodology
├── CITEgeist_env.yml                   # Conda environment specification
├── README.md                           # Root documentation
├── LICENSE                             # BSD 3-Clause License
└── .github/workflows/
    └── pylint.yml                      # GitHub Actions CI configuration
```

---

## 3. SYSTEM REQUIREMENTS & DEPENDENCIES

### Python & Core Dependencies
- **Python Version:** 3.10 (required, though CI tests 3.8, 3.9, 3.10)
- **Critical Dependency:** Gurobi optimizer (v11.0.2) - requires academic license
- **Build/Package Tools:** setuptools, wheel, pip

### Key Scientific Libraries
| Package | Version | Purpose |
|---------|---------|---------|
| scanpy | 1.10.4 | Single-cell analysis, AnnData manipulation |
| anndata | 0.11.3 | Single-cell data structures |
| numpy | 1.26.4 | Numerical computation |
| pandas | 2.2.3 | Data manipulation |
| scipy | 1.13.1 | Scientific computing (optimization, stats) |
| scikit-learn | 1.6.1 | Machine learning utilities |
| gurobipy | 11.0.2 | Mathematical optimization solver |
| matplotlib | 3.10.0 | Visualization |
| seaborn | 0.13.2 | Statistical visualization |
| squidpy | 1.6.2 | Spatial transcriptomics analysis |
| spatialdata | 0.2.5.post0 | Spatial data handling |
| h5py | 3.12.1 | HDF5 file I/O |
| pyarrow | 19.0.0 | Parquet file I/O |
| statsmodels | 0.14.4 | Statistical models |

### Hardware Requirements
- **RAM:** Minimum 16GB, Recommended 32GB+, Optimal 64GB+
- **CPU:** Multi-core processor (8+ cores recommended)
- **Storage:** 16GB+ for installation and data
- **OS:** Linux (recommended), macOS, or Windows 10+ with WSL2

### Runtime Expectations
- Small dataset (~5000 spots): ~2 hours
- Medium dataset: ~4 hours
- Large dataset: ~10 hours
- (on 16 threads, 32GB RAM system)

---

## 4. MAIN MODULES & CORE FUNCTIONALITY

### 4.1 Core Class: `CitegeistModel` (citegeist_model.py, 834 lines)

**Purpose:** Main orchestrator for the CITEgeist analysis pipeline

**Key Methods:**

#### Data Management
- `split_adata()` - Separate gene expression and antibody capture data
- `load_cell_profile_dict(dict)` - Load cell type marker definitions
- `copy_gex_to_protein_adata()` - Synchronize spot indices between modalities

#### Preprocessing
- `filter_gex()` - Filter genes by expression criteria
- `preprocess_gex(target_sum)` - Count-preserving normalization
- `preprocess_antibody()` - Winsorization, CLR normalization
- `winsorize()`, `row_normalize()`, `global_clr()` - Normalization utilities

#### Core Analysis
- `run_cell_proportion_model()` - Estimate cell type proportions
- `run_cell_expression_pass1()` - First pass gene expression deconvolution
- `compute_expression_prior()` - Compute prior for second pass
- `run_cell_expression_pass2()` - Second pass with prior regularization

#### Utilities
- `register_gurobi(license_path)` - Configure Gurobi optimizer
- `append_proportions_to_adata()` - Store results in AnnData
- `append_gex_to_adata()` - Store gene expression layers
- `_save_profiles_to_parquet()` - Export results

### 4.2 Optimization Engine: `gurobi_impl.py` (1,219 lines) [CRITICAL]

**Purpose:** Mathematical optimization algorithms using Gurobi

**Key Functions:**

#### Antibody Mapping & Cell Profiling
- `map_antibodies_to_profiles()` - Map antibody data to cell type profiles
- `compute_global_prior()` - Compute global expression patterns from pass 1

#### Primary Optimization Functions
1. **`optimize_cell_proportions()`**
   - Performs global cell type proportion estimation
   - Uses linear optimization with Gurobi
   - Outputs initial Y (proportions) and β values
   
2. **`finetune_cell_proportions()`**
   - Spatial-aware refinement of proportions
   - Uses fixed-radius neighborhood information
   - Parallel processing capability
   - EM-algorithm with convergence criteria
   - Parameters: radius, lambda_reg, alpha (L1/L2 tradeoff), max_y_change

3. **`optimize_gene_expression()`**
   - Deconvolves gene expression profiles per cell type per spot
   - Incorporates enrichment weighting (global + local)
   - Integer programming formulation
   - Processes spots in parallel
   - Checkpoint system for long runs

4. **`normalize_counts()`**
   - Count-preserving normalization while maintaining integer values
   - Critical for discrete count optimization

#### Supporting Functions
- `compute_global_prior()` - Generate Bayesian priors from pass 1
- Expression enrichment calculations
- Neighborhood-based regularization

### 4.3 Checkpoint Management: `checkpoints.py` (173 lines)

**Purpose:** Save/load progress for long-running optimizations

**Key Methods:**
- `check_complete_run()` - Check if analysis is already complete
- `load_latest_checkpoint()` - Recover from interruptions
- `save_checkpoint()` - Periodic progress saving
- `mark_complete()` - Flag complete runs

### 4.4 Utilities: `utils.py` (378 lines)

**Purpose:** Helper functions for validation, logging, and spatial analysis

**Key Functions:**
- `validate_cell_profile_dict()` - Input validation
- `setup_logging()` - Initialize logging system
- `cleanup_memory()` - Force garbage collection
- `get_neighbors_with_fixed_radius()` - Spatial neighbor detection
- `find_fixed_radius_neighbors()` - Core neighbor search
- `assert_neighborhood_size()` - Validate neighborhood data
- `export_anndata_layers()` - Save expression layers
- `benchmark_cell_proportions()` - Calculate deconvolution metrics
- `calculate_expression_metrics()` - Calculate gene expression metrics

### 4.5 Analysis Functions: `analysis_functions.py` (42 lines)

**Purpose:** Post-analysis utilities

**Functions:**
- `filter_spots_by_celltype_proportion()` - Filter by proportion threshold
- `expand_prop_gex_adata()` - Expand to cell-type-level observations
- `get_celltype_expression_data()` - Extract cell type data

---

## 5. BENCHMARKING & TESTING INFRASTRUCTURE

### 5.1 Test Entry Point: `test_citegeist_simulated.py` (380+ lines)

**Type:** Simulation-based integration test

**Purpose:** 
- Run full pipeline on simulated data with known ground truth
- Calculate performance metrics
- Compare against other methods

**Key Components:**
```
Usage: python test_citegeist_simulated.py \
    --radius <float> \
    --lambda_reg <float> \
    --alpha_elastic <float> \
    --max_y_change <float> \
    --input_folder <path> \
    --output_folder <path> \
    --sample_prefix <str> \
    [--profiling_only] \
    [--skip_pass2]
```

**Workflow:**
1. Load paired CITE and GEX H5AD files
2. Initialize CitegeistModel with simulation mode
3. Preprocess data
4. Run cell proportion model with finetuning
5. Calculate proportion metrics (JSD, RMSE, MAE, correlation)
6. Run pass 1 gene expression optimization
7. Calculate gene expression metrics (RMSE, NRMSE, MAE)
8. Compute prior and prepare for pass 2
9. Track runtimes

**Metrics Computed:**
- **Proportions:** Jensen-Shannon Divergence, RMSE, MAE, Pearson correlation
- **Gene Expression:** RMSE, Normalized RMSE (NRMSE), MAE (after log1p transform)

### 5.2 Benchmarking Modules

#### `benchmarking_spot_deconv.py`
**Metrics for Cell Type Deconvolution:**
- Jensen-Shannon Divergence (JSD) - 0 (identical) to 1 (different)
- Per-cell-type RMSE & MAE
- Overall (sum) RMSE & MAE
- Pearson correlation

```python
benchmark_performance(predictions, ground_truth, cell_type_names)
# Returns: JSD, RMSE dict, Sum_RMSE, MAE dict, Sum_MAE, correlation
```

#### `benchmarking_gex.py`
**Metrics for Gene Expression Prediction:**
- Per-cell-type RMSE (after log1p)
- Normalized RMSE (NRMSE) - 'range' or 'mean' normalization
- Per-cell-type MAE
- Aggregate statistics (mean, median)

```python
calculate_rmse(ground_truth_dir, predictions_dir, normalize='range')
# Returns: per-cell-type metrics, averages, medians
```

#### Method-Specific Wrappers
- `citegeist_bench_wrapper.py` - CITEgeist evaluation
- `cell2location_bench_wrapper.py` - Cell2location comparison
- `tangram_bench_wrapper.py` - Tangram comparison
- `RCTD_bench_wrapper.py` - RCTD comparison
- `seurat_bench_wrapper.py` - Seurat comparison

### 5.3 Simulated Data

**Location:** `/home/user/CITEgeist/Benchmarking/simulation/`

**Data Source:**
- Single-cell reference: Wu et al. 2021 breast cancer scRNA-seq
- Spatial layout: 10x Visium platform (hex-based using scCube)

**Simulation Types:**
1. **Highly Segmented:** Distinct spatial regions with homogeneous composition
2. **Mixed:** Complex tissue with heterogeneous cell mixing

**Parameters:**
- ~5000 spots per simulation
- 9 cell types: Cancer Epithelial, Normal Epithelial, T-cells, B-cells, Myeloid, PVL, Plasmablasts, Endothelial, CAFs
- Replicates for statistical robustness

**Generated Data Formats:**
- Ground truth cell type proportions (CSV)
- Ground truth gene expression per cell type (CSV layers)
- Spatial coordinates

---

## 6. CI/CD INFRASTRUCTURE

### GitHub Actions Workflow: `.github/workflows/pylint.yml`

**Current Status:** Pylint code quality checks only (no pytest)

**Configuration:**
```yaml
- Python versions tested: 3.8, 3.9, 3.10
- Trigger: On push
- Tool: pylint
- Scope: CITEgeist folder
```

**Current Gaps:**
- No pytest configuration
- No unit tests in CI
- No integration test automation
- No test coverage reporting

---

## 7. DATA MANAGEMENT

### Data Formats

#### Input Data
- **H5AD format** (AnnData): Gene expression and antibody capture data
  - Required: `.obs` (spot metadata), `.X` (count matrix), `.var` (features)
  - Essential column: `.var['feature_types']` (categorizes as "Gene Expression" or "Antibody Capture")
  - Coordinates: `.obsm['spatial']` (2D coordinates)

#### Output Data
- **CSV files:** Cell type proportions (spots × cell types)
- **Parquet files:** Gene expression profiles (spotwise)
- **Layer exports:** Per-cell-type gene expression (CSV in per-celltype directories)

#### Cell Profile Dictionary
```python
cell_type_profiles = {
    "CellType1": {
        "Major": ["marker_protein_1", "marker_protein_2"],
        # Optional:
        "Minor": ["secondary_marker_1"]
    },
    "CellType2": { ... }
}
```

### Example Data Sources
- Simulation data: Generated via scCube (provided in repo)
- Real data: Clinical breast cancer samples (vignettes in examples/)
- Benchmarking comparisons: Multiple replicates per condition

---

## 8. KEY ENTRY POINTS & WORKFLOWS

### Standard Workflow

```python
# 1. Initialize model
model = CitegeistModel(
    sample_name="sample_1",
    adata=adata,  # or simulation=True with gene_expression_adata, antibody_capture_adata
    output_folder="/path/to/output"
)

# 2. Configure Gurobi
model.register_gurobi("/path/to/gurobi.lic")

# 3. Prepare data
model.split_adata()
model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
model.copy_gex_to_protein_adata()
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

# 4. Load cell profiles
model.load_cell_profile_dict(cell_type_profiles)

# 5. Cell proportion analysis
cell_props_global, cell_props_finetuned = model.run_cell_proportion_model(
    radius=4,
    lambda_reg=0.001,
    alpha=0.7,
    max_y_change=0.2
)

# 6. Gene expression pass 1
pass1_results = model.run_cell_expression_pass1(radius=4)

# 7. Gene expression pass 2 (optional)
prior_info = model.compute_expression_prior(...)
pass2_results = model.run_cell_expression_pass2(...)
```

### Command-Line Interface

```bash
python test_citegeist_simulated.py \
    --radius 4 \
    --lambda_reg 0.001 \
    --alpha_elastic 0.7 \
    --max_y_change 0.2 \
    --input_folder ./data \
    --output_folder ./results \
    --sample_prefix Wu_rep
```

### HPC Submission

```bash
sbatch examples/sbatch_sample.sh  # SLURM array jobs
```

---

## 9. TESTING NEEDS & GAPS

### Current Testing Status
- **Unit Tests:** None (0 exist)
- **Integration Tests:** Manual simulation test script only
- **Test Framework:** Not set up (no pytest.ini, conftest.py)
- **Coverage:** Unknown/untracked

### Critical Areas Requiring Tests

#### 1. **Preprocessing Pipeline**
- `filter_gex()` - Gene filtering logic
- `preprocess_gex()` - Count normalization preservation
- `preprocess_antibody()` - Winsorization, CLR normalization
- `split_adata()` - Data partitioning correctness

#### 2. **Optimization Algorithms** [HIGH PRIORITY]
- `optimize_cell_proportions()` - Gurobi model correctness
- `finetune_cell_proportions()` - Spatial regularization
- `optimize_gene_expression()` - Discrete optimization
- `compute_global_prior()` - Prior computation

#### 3. **Data Validation**
- `validate_cell_profile_dict()` - Input structure validation
- `map_antibodies_to_profiles()` - Antibody-profile mapping
- Neighborhood detection correctness

#### 4. **Spatial Analysis**
- `get_neighbors_with_fixed_radius()` - Neighbor identification
- Radius parameter sensitivity
- Edge cases (boundary spots)

#### 5. **Metrics & Benchmarking**
- `benchmark_cell_proportions()` - JSD, RMSE, MAE calculation
- `calculate_expression_metrics()` - Expression metric computation
- Correctness vs reference implementations

#### 6. **Checkpoint System**
- Save/load cycle correctness
- Recovery from interruptions
- Corrupted file handling

#### 7. **Numerical Stability**
- NaN/Inf handling in preprocessing
- Large dataset performance
- Memory efficiency under load

---

## 10. KEY CONFIGURATION PARAMETERS

### Cell Proportion Model Parameters
| Parameter | Default | Range | Purpose |
|-----------|---------|-------|---------|
| `radius` | None (required) | 2-10 (pixels) | Neighborhood size for spatial regularization |
| `lambda_reg` | 0.001 | 1e-5 to 1 | L1/L2 regularization strength |
| `alpha` | 0.7 | 0-1 | Elastic net mixing (0=L2 ridge, 1=L1 lasso) |
| `max_y_change` | 0.2 | 0-1 | Max proportion change per iteration |
| `tolerance` | 1e-4 | 1e-5 to 1e-3 | Convergence criterion |
| `max_iterations` | 20 | 10-100 | Maximum EM iterations |

### Gene Expression Model Parameters
| Parameter | Default | Purpose |
|-----------|---------|---------|
| `global_enrichment_weight` | 0.5 | Weight for global expression prior |
| `local_enrichment_weight` | 0.5 | Weight for local spatial enrichment |
| `lambda_reg_gex` | 0.001 | Gene expression regularization |
| `lambda_prior` | 1.0 | Prior influence strength (pass 2) |

---

## 11. EXTERNAL DEPENDENCIES & LICENSES

### Critical External Tools
- **Gurobi Optimizer** (11.0.2)
  - License: Academic (free for .edu) or Commercial
  - Required for optimization
  - No free alternative currently supported

### Optional Integration
- Tangram, Cell2location, RCTD, Seurat (for comparison)
- scCube (for simulation generation)

---

## 12. DOCUMENTATION REFERENCES

### Available Documentation
1. **Root README:** Project overview, installation, features
2. **CITEgeist/README.md:** Quick start, system requirements, parameter guide
3. **Examples/README.md:** Vignette descriptions, usage examples
4. **Benchmarking/README.md:** Methodology, metrics definitions
5. **Benchmarking/simulation/README.md:** Data generation details
6. **Vignettes:** Three Jupyter notebooks with real-world workflows
7. **Inline docstrings:** Class/method documentation

### Missing Documentation
- API reference documentation
- Parameter sensitivity analysis guide
- Troubleshooting guide
- Performance profiling tips

---

## SUMMARY TABLE

| Aspect | Details |
|--------|---------|
| **Language** | Python 3.10 |
| **Lines of Code (model/)** | 2,667 |
| **Main Classes** | 1 (CitegeistModel) |
| **Core Functions** | 4+ major optimization functions |
| **Test Files** | 1 integration test (manual) |
| **Unit Tests** | 0 |
| **CI/CD** | GitHub Actions (pylint only) |
| **Dependencies** | 100+ (heavy: scanpy, Gurobi, scipy) |
| **Benchmarks** | Comparison with 5 methods |
| **Real Data Sources** | Breast cancer CITE-seq samples |
| **Simulated Data** | 9 cell types, 2 spatial patterns |
| **Output Formats** | CSV, Parquet, H5AD layers |
| **Runtime** | 2-10 hours (depends on size) |

