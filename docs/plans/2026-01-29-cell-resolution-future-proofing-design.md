# Design: Cell-Resolution Future-Proofing for CITEgeist Modules 1-3

**Date**: 2026-01-29
**Branch**: TBD (feature branch from `dev`)
**Status**: Design

## Motivation

CITEgeist currently operates on spot-level spatial transcriptomics data (Visium), where each observation is a mixture of cells. As single-cell spatial platforms mature (Xenium, Visium HD, CosMx) and begin offering co-measured protein + RNA per cell, CITEgeist needs to provide value at cell-level resolution while preserving its optimization-based, protein-anchored philosophy.

**Core assumption preserved**: Protein markers are more stable and accurate than transcriptomic measurements. Protein anchors identity; RNA carries biological detail; spatial context resolves ambiguity.

## Architecture: Explicit Resolution Mode

A `resolution` parameter on `CitegeistModel` gates Module 3 behavior while Modules 1, 2, 4, and 5 remain unchanged in logic (only default parameters shift).

```python
class CitegeistModel:
    def __init__(self, ..., resolution: str = "spot"):
        """
        Args:
            resolution: "spot" for Visium/aggregated data (mixture deconvolution),
                        "cell" for Xenium/Visium HD single-cell data
                        (soft classification + expression optimization).
        """
        assert resolution in ("spot", "cell")
        self.resolution = resolution
```

No auto-detection. User declares their data type.

## Parameter Presets

Resolution-aware defaults resolved once at init. User can override any parameter individually.

```python
RESOLUTION_DEFAULTS = {
    "spot": {
        "neighbor_k": 8,
        "morans_k": 8,
        "coloc_k_range": [4, 6, 8, 10, 12, 14],
        "laplacian_k": 8,
        "lambda_spatial": 0.1,
        "lambda_sparse": 0.0,       # no sparsity for mixtures
        "pass2_library_slack": 1.0,  # strict count conservation
    },
    "cell": {
        "neighbor_k": 50,
        "morans_k": 50,
        "coloc_k_range": [20, 40, 60, 80, 100],
        "laplacian_k": 50,
        "lambda_spatial": 0.01,      # weaker — cells are locally heterogeneous
        "lambda_sparse": 0.1,        # encourage near-one-hot assignments
        "pass2_library_slack": 1.5,  # allow dropout recovery headroom
    },
}
```

Existing spot-level defaults are preserved exactly. No behavior change for `resolution="spot"`.

## Modules 1 & 2: Parameter Routing Only

No logic changes. Functions that currently hardcode `k` values accept them from the preset instead:

- `marker_interest.py`: `detect_interesting_markers()` accepts `morans_k` parameter
- `spatial_colocalization.py`: `analyze_colocalization()` accepts `k_range` parameter
- `utils.py`: `find_fixed_radius_neighbors()` and `build_spatial_laplacian()` accept `k` from preset

These modules are already resolution-agnostic — k-NN spatial graphs, Moran's I, colocalization statistics, and profile discovery all work identically at cell or spot scale.

## Module 3 Pass 1: Cell-Level Soft Classification

### Spot-level (unchanged)
Solves: given observed protein intensities per spot, find mixture proportions minimizing reconstruction error with spatial smoothing.

### Cell-level (new)
Same QP structure with an added **sparsity penalty** encouraging near-one-hot solutions:

```
minimize:
    sum_i sum_m (S_im - beta_m * sum_j Y_ij * P_jm)^2    # protein reconstruction
  + lambda_sparse * sum_i sum_j Y_ij                       # L1 sparsity penalty
  + lambda_spatial * Y^T L Y                                # spatial Laplacian coherence

subject to:
    Y_ij >= 0
    sum_j Y_ij >= 0.9,  sum_j Y_ij <= 1.1
```

**Key differences from spot-level:**
- `lambda_sparse > 0` pushes solutions toward near-one-hot vectors. A clean T-cell gets `Y = [0.0, 0.95, 0.05, ...]`. A doublet naturally gets `Y = [0.0, 0.45, 0.50, ...]`.
- Weaker `lambda_spatial` (cells are more locally heterogeneous than spots).
- Larger `k` for Laplacian (50 vs 8) to capture tissue-scale patterns at cellular density.
- Output format is identical: `(n_observations, n_types)` matrix summing to ~1. Downstream code is unaware of resolution.

**Implementation**: Add `lambda_sparse` parameter to existing `optimize_cell_proportions()` / `optimize_cell_proportions_per_marker()`. When `lambda_sparse=0` (spot mode), behavior is identical to current code.

## Module 3 Pass 2: Cell-Level True Count Estimation

### Spot-level (unchanged)
Partitions observed gene counts across cell types to maximize protein-informed enrichment weights. Count conservation is strict.

### Cell-level (new)
Optimizes estimated true expression per cell against observed counts, using protein-anchored priors and spatial neighbors for dropout recovery:

```
minimize:
    sum_g L(X_true[i,g], X_obs[i,g])                          # data fidelity
  + lambda_enrich * sum_g -E[j*,g] * X_true[i,g]              # enrichment prior
  + lambda_spatial * sum_g sum_{k in N_j*(i)} (X_true[i,g] - X_obs[k,g])^2  # spatial coherence

subject to:
    X_true[i,g] >= 0                                           # non-negative counts
    sum_g X_true[i,g] <= C * library_size_i                    # bounded total counts
```

Where:
- `L()` is quadratic loss keeping `X_true` close to `X_obs` where counts are non-zero
- `E[j*,g]` is the enrichment weight for gene `g` in cell type `j*` (dominant type from Pass 1) — same global + local prior computation as spot level
- `N_j*(i)` is spatial neighbors sharing the same dominant type — smoothing respects type boundaries
- `C` is `pass2_library_slack` (1.5 for cell mode) allowing modest count inflation for dropout recovery

**Key design choices:**
- **Optimization-based, not heuristic** — stays in CITEgeist's Gurobi QP framework
- **Zero-count genes**: Enrichment prior pulls `X_true` above zero for genes expected in type `j*`, but library size constraint prevents runaway imputation
- **Non-zero genes**: Data fidelity keeps values close to observed; enrichment prior only shifts when spatial neighbors strongly agree
- **Doublet handling**: Cells with max `Y[i,j] < 0.6` route to spot-level-style count partitioning across the top 2 types — the existing Pass 2 formulation handles this naturally

**Implementation**: New function `optimize_gene_expression_cell()` in `gurobi_impl.py`, parallel to existing `optimize_gene_expression()`.

## File Changes

| File | Change |
|------|--------|
| `citegeist_model.py` | Add `resolution` param, preset loading, dispatch logic in `run_cell_proportion_model()` and `run_cell_expression_pass1()` |
| `gurobi_impl.py` | Add `lambda_sparse` to `optimize_cell_proportions()`, add new `optimize_gene_expression_cell()` |
| `utils.py` | Accept `k` parameter from preset instead of hardcoded defaults |
| `marker_interest.py` | Accept `morans_k` parameter instead of hardcoded default |
| `spatial_colocalization.py` | Accept `k_range` parameter instead of hardcoded default |

5 files modified, 0 new files.

## Not In Scope

- Auto-detection of resolution from data
- Visium HD-specific binning utilities
- Preprocessing pipelines for Xenium/CosMx raw formats
- Changes to Modules 4 or 5
- Standalone doublet detection module
- Cell-level benchmarking framework (separate effort)

## Validation Plan

### Data

**10x Xenium Renal Cell Carcinoma Dataset**
- 465,534 single cells with co-measured protein (27 markers) + RNA (405 genes)
- 5 spatially distinct regions (regions 0-4)
- Cell-level ground truth: 10 granular RNA-derived cell types (461,886 labeled cells)
- Existing pseudo-Visium aggregation (7,054 spots) for spot-level comparison

**Data locations:**
- Single-cell raw: `/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma/`
- Cell-level GT: `Benchmarking/xenium_pseudovisium/data_rna_gt/cell_types.csv`
- Pseudo-Visium: `Benchmarking/xenium_pseudovisium/data_rna_gt/h5ad_objects/Xenium_region_{0-4}_CITE.h5ad`
- Spot-level GT proportions: `Benchmarking/xenium_pseudovisium/data_rna_gt/ground_truth/Xenium_region_{0-4}_prop.csv`

### Existing Baseline

An existing single-cell pipeline (`Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_pipeline.py`) has already been run on this data with outputs in `output_singlecell/`. This serves as the comparison baseline for the new `resolution="cell"` mode.

### Test Strategy

#### Test 1: Module 1/2 Consistency
- Run Modules 1 and 2 with `resolution="cell"` on Xenium regions 0-4
- Compare discovered profiles against `resolution="spot"` (pseudo-Visium) profiles
- **Success**: Same marker groupings emerge at both resolutions (profiles should be resolution-invariant)

#### Test 2: Pass 1 Cell Type Classification Accuracy
- Run Pass 1 with `resolution="cell"` on all 5 regions
- Compare dominant `Y[i,j]` assignments against RNA-derived ground truth labels
- **Metrics**: Accuracy, F1 per cell type, confusion matrix
- **Success**: Classification accuracy comparable to or better than the existing `run_singlecell_pipeline.py` cell assignment
- **Comparison**: Also evaluate against spot-level proportions aggregated from cell assignments vs spot-level ground truth proportions

#### Test 3: Pass 2 Expression Optimization Value
- Run Pass 2 with `resolution="cell"` on all 5 regions
- Compare `X_true` vs `X_obs` for known cell-type marker genes
- **Metrics**:
  - Dropout recovery rate: fraction of zero-count marker genes recovered to non-zero
  - Expression coherence: Pearson correlation of `X_true` profiles against type-average expression
  - Spatial smoothness: Moran's I of `X_true` vs `X_obs` for marker genes within cell type clusters
- **Success**: `X_true` has higher marker gene expression coherence than `X_obs` while maintaining biological variation (not over-smoothed)

#### Test 4: Doublet Robustness
- Identify cells with max `Y[i,j] < 0.6` (candidate doublets)
- Check if these overlap with cells at tissue boundaries (where mixing is biologically plausible)
- **Success**: Low-confidence cells are spatially enriched at type boundaries, not randomly distributed

#### Test 5: End-to-End Module 4 Impact
- Run Module 4 program discovery on `X_true` (cell-mode Pass 2 output) vs `X_obs` (raw counts)
- Compare number of discovered programs, spatial coherence (Moran's I), and variance explained
- **Success**: Programs from `X_true` show higher spatial coherence without loss of biological programs
