# Stage 2: Morphology-Based Single-Cell Assignment

**Date:** 2026-02-28
**Status:** Approved
**Author:** Claude + User

## Overview

Stage 2 assigns individual nuclei to cell types using morphology features, constrained by cell counts from Stage 1 (hybrid protein-based proportions).

**Key insight:** Global morphology classification fails (silhouette = -0.03), but constrained assignment works (46% accuracy vs 22% random). Instead of asking "what type is this cell?" we ask "given counts [2,2,1], which assignment maximizes likelihood?"

## Architecture

```
Stage 1: Hybrid Protein → Proportions (existing)
├── CITEgeist continuous model + detection post-filter
└── Output: proportions_df, cell_counts_df (per spot)
                    │
                    ▼
Stage 2: Morphology → Single-Cell Assignment (NEW)
├── Extract patches around cell centroids
├── Compute extended morphology features (12 features)
├── Learn GMM per cell type from high-purity spots
├── Hungarian assignment with count constraints
└── Output: per-nucleus cell type assignments
```

## Data Flow

### Inputs (Benchmark Mode - Approach B)

| Source | Purpose |
|--------|---------|
| Xenium cells.parquet | Cell centroids (x, y) |
| cell_type_assignments.csv | Per-cell GT labels |
| Morphology image (PNG) | DAPI + boundary channels |
| Pseudo-Visium h5ad | GEX + protein data |
| cell_to_spot_mapping.csv | Which cells belong to which spot |

### Outputs

- Per-nucleus cell type assignments
- Accuracy metrics vs GT
- Comparison to random baseline

### Production Mode (Approach A)

Same flow but uses Cellpose centroids + spatial registration to match Xenium cells for GT evaluation only. Users won't have Xenium GT in production.

## Components

### New Modules in `CITEgeist/model/`

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `morphology_features.py` | Extract extended features | `extract_extended_features(patch) → (12,)` |
| `morphology_classifier.py` | Learn GMM per cell type | `fit()`, `log_likelihood()` |
| `constrained_assignment.py` | Hungarian with constraints | `hungarian_assign(log_likes, counts)` |
| `stage2_pipeline.py` | Orchestrates Stage 2 | `Stage2Pipeline.run()` |

### Extended Features (12 total)

```python
# Basic (6) - from prototype
dapi_mean, dapi_std, dapi_area, boundary_mean, boundary_std, channel_corr

# Shape (3) - new
circularity, eccentricity, solidity

# Texture (3) - new
dapi_entropy, dapi_contrast, boundary_gradient_mag
```

### Benchmark Script

```
Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py
```

## Error Handling

| Edge Case | Handling |
|-----------|----------|
| Spot with 0 nuclei | Skip spot, log warning |
| Spot with 1 nucleus | Assign to highest-proportion type |
| Cell type with 0 samples for GMM | Use global mean/cov as fallback |
| NaN in features | Replace with zeros, flag in metadata |
| Count mismatch | Adjust using largest-remainder method |
| No high-purity spots for a type | Weight all spots by proportion |

## Testing & Success Criteria

### Unit Tests

```
CITEgeist/tests/test_stage2_morphology.py
├── test_extract_extended_features()
├── test_gmm_classifier_fit()
├── test_hungarian_respects_counts()
└── test_pipeline_end_to_end()
```

### Success Criteria

| Metric | Target | Rationale |
|--------|--------|-----------|
| Overall accuracy | > 35% | Prototype achieved 46% |
| vs random baseline | > 1.5x | Prototype was 2x (46/22) |
| Per-type accuracy | > random for all types | No regression |

### Deliverables

1. 4 new modules in `CITEgeist/model/`
2. 1 benchmark script
3. 1 test file
4. Results added to `canonical_results.json`

## Design Decisions

1. **Approach B for benchmark:** Use original Xenium cells directly for clean evaluation (no registration error)
2. **Design for Approach A:** Code structure supports Cellpose patches for production
3. **GMM per class:** Handles within-class heterogeneity better than single Gaussian
4. **Extended features:** Shape + texture features may improve discrimination
5. **Re-run Stage 1 inline:** Combined pipeline for reproducibility

## References

- Working prototype: `Benchmarking/xenium_benchmarking/CITEgeist/src/test_constrained_assignment.py`
- Stage 1 hybrid: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py`
- Xenium data: `Benchmarking/xenium_pseudovisium/data_protein_gt/`
