# CITEgeist Development Progress Report

**Last Updated**: January 9, 2026
**Branch**: `hierarchical_approach`
**Status**: Active Development (Pre-Publication)

---

## Executive Summary

CITEgeist is a computational framework for spatial multi-omic deconvolution. This report documents the implementation status of all modules, benchmarking infrastructure, and remaining work.

**Key Achievement**: All 5 core modules are fully implemented with comprehensive test coverage.

---

## Module Implementation Status

### Module 1: Marker Interest Detection
**Status**: COMPLETE
**File**: `CITEgeist/model/marker_interest.py`
**Test**: `tests/test_marker_interest.py`

Identifies spatially-variable protein markers worth analyzing:
- Kurtosis gating for signal detection
- GMM-based signal/noise separation
- Moran's I spatial autocorrelation

### Module 2a: Marker Colocalization Analysis
**Status**: COMPLETE
**File**: `CITEgeist/model/spatial_colocalization.py`
**Test**: `tests/test_spatial_colocalization.py`

Analyzes spatial relationships between markers:
- Same-spot co-occurrence
- Expression correlation
- Adjacent-spot enrichment
- Bivariate Moran's I

### Module 2b: Profile Discovery
**Status**: COMPLETE
**File**: `CITEgeist/model/spatial_colocalization.py`
**Test**: `tests/test_module2_profile_discovery.py`

Discovers cell type profiles from colocalization patterns:
- Significance-filtered graph construction
- Connected component detection
- Hierarchical clustering with dynamic tree cutting

### Module 2c: Profile Selection
**Status**: COMPLETE
**File**: `CITEgeist/model/spatial_colocalization.py`
**Test**: `tests/test_module2c_profile_selection.py`

Selects optimal profiles using reconstruction accuracy:
- Dual variance checkpoints
- Coverage and non-redundancy optimization
- Spatial variance targeting

### Module 3: Deconvolution (Two-Pass)
**Status**: COMPLETE
**Files**: `CITEgeist/model/gurobi_impl.py`, `CITEgeist/model/citegeist_model.py`
**Test**: `tests/test_citegeist_simulated.py`

**Pass 1**: Cell-type proportion estimation
- Quadratic programming optimization
- Per-marker beta optimization
- Laplacian spatial smoothing
- Local finetuning with neighbors

**Pass 2**: Gene expression deconvolution
- Proportion-guided decomposition
- Global + local enrichment weights
- Per-cell-type expression layers

### Module 4: Anchored Program Discovery
**Status**: COMPLETE
**File**: `CITEgeist/model/anchored_program_discovery.py`
**Test**: `tests/test_module4_deconvolved.py`

Discovers gene programs from deconvolved layers:
- NMF-based program extraction
- Moran's I spatial coherence
- K-means subpopulation identification

### Module 4b: Bivariate Program Relationships
**Status**: COMPLETE
**File**: `CITEgeist/model/anchored_program_discovery.py`
**Test**: `tests/test_module4_deconvolved.py`

Analyzes program-program relationships:
- Bivariate Moran's I between programs
- Co-localized vs exclusive detection

### Module 5: Cross-Sample Integration
**Status**: COMPLETE
**File**: `CITEgeist/model/cross_sample_integration.py`
**Test**: `tests/test_cross_sample_integration.py`

Integrates programs across samples:
- Harmony-style batch correction
- Hierarchical clustering alignment
- Relationship conservation scoring

---

## Benchmarking Infrastructure

### Simulated Data Benchmarking

**Location**: `tests/` and `replicates/`

**Data**:
- 5 replicates per test set (`Wu_rep_0` to `Wu_rep_4`)
- 9 cell types with 18 specific + 82 nonspecific markers
- Ground truth proportions and GEX layers

**Test Framework**:
- Main test: `test_citegeist_simulated.py` (924 lines)
- Array jobs: `sbatch_CITEgeistArray.sh` (24 parameter combinations)
- Aggregator: `benchmark_array.py`

**Metrics**:
- Cell proportions: JSD, RMSE, MAE, Pearson correlation
- Gene expression: Holistic RMSE, per-celltype metrics
- Profile discovery: Recovery rate, precision, recall, FDR

### Xenium Benchmarking

**Location**: `Benchmarking/xenium_benchmarking/`

**Methods Compared**:
| Method | Proportions | GEX Layers | GPU Required |
|--------|-------------|------------|--------------|
| CITEgeist | Yes | Yes | No |
| Cell2Location | Yes | Yes | Yes |
| Tangram | Yes | Yes | Yes |
| RCTD | Yes | No | No |
| Seurat | Yes | No | No |

**Ground Truth**:
- 10 granular RNA-based cell types
- 14 Xenium pseudo-Visium regions

**Directory Structure**:
```
Benchmarking/xenium_benchmarking/
├── CITEgeist/
│   ├── slurm/          # 8 benchmark scripts
│   ├── src/            # 5 Python scripts
│   ├── output_autodiscovery/
│   ├── output_granular/
│   └── output_rna_gt/
├── Cell2Location/
│   ├── slurm/
│   ├── src/
│   ├── reference_model_granular/
│   └── output_granular/
├── RCTD/
├── Seurat/
├── Tangram/
├── evaluation/
│   ├── slurm/
│   ├── src/evaluate_all_methods.py
│   ├── results_granular/
│   └── results_autodiscovery/
└── reference_data/GSE156632/
```

---

## Recent Changes (January 2026)

### Cleanup Completed

1. **Removed obsolete files** (60 files):
   - Old top-level `slurm/` scripts replaced by method-specific versions
   - Old top-level `src/` code replaced by method-specific implementations
   - Old SLURM log files

2. **Added SLURM mail directives** (11 scripts):
   - All benchmark scripts now have required `--mail-type=ALL` and `--mail-user` directives

3. **Fixed evaluation inconsistencies**:
   - Seurat cell type name normalization (R dots to standard names)
   - Tangram proportion normalization (row-sum to 1.0)
   - Both output naming patterns supported

### Git Commits

1. `ea17b44` - chore: Remove obsolete top-level xenium benchmark files
2. `65f30f5` - fix: Add required SLURM mail directives to all benchmark scripts
3. `51507a9` - fix: Handle cell type naming inconsistencies and output format variations

---

## Test Coverage Summary

| Component | Test File | Status |
|-----------|-----------|--------|
| Module 1 (Marker Interest) | `test_marker_interest.py` | Passing |
| Module 2a (Colocalization) | `test_spatial_colocalization.py` | Passing |
| Module 2b (Profile Discovery) | `test_module2_profile_discovery.py` | Passing |
| Module 2c (Profile Selection) | `test_module2c_profile_selection.py` | Passing |
| Module 3 (Deconvolution) | `test_citegeist_simulated.py` | Passing |
| Module 4 (Anchored Programs) | `test_module4_deconvolved.py` | Passing |
| Module 4b (Bivariate) | `test_module4_deconvolved.py` | Passing |
| Module 5 (Integration) | `test_cross_sample_integration.py` | Passing |
| Gurobi Implementation | `test_gurobi_impl.py` | Passing |
| Checkpoints | `test_checkpoints.py` | Passing |
| Model Core | `test_model.py` | Passing |

---

## Remaining Work

### High Priority
- [ ] Run full Xenium benchmark comparison across all 5 methods
- [ ] Generate publication figures from benchmark results
- [ ] Validate Module 5 on multi-sample patient data

### Medium Priority
- [ ] Create project-managed Tangram environment (currently uses external)
- [ ] Integrate simulated benchmarks into CI/CD pipeline
- [ ] Add performance regression tracking

### Low Priority
- [ ] Standardize output file naming across all methods
- [ ] Create configuration file for method paths
- [ ] Update README with benchmark instructions

---

## Running Benchmarks

### Simulated Data (Quick Test)
```bash
cd tests
sbatch --array=0 sbatch_CITEgeistArray.sh
```

### Xenium Auto-Discovery
```bash
cd Benchmarking/xenium_benchmarking/CITEgeist/slurm
sbatch run_autodiscovery_benchmark.sh
```

### Full Evaluation
```bash
cd Benchmarking/xenium_benchmarking/evaluation/slurm
sbatch evaluate_all_methods.sh
```

---

## Contact

- **Developer**: Alexander Chang (alc376@pitt.edu)
- **Lab**: Lee/Oesterreich Lab
- **Repository**: https://github.com/leeoesterreich/CITEgeist
