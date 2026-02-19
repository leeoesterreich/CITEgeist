# Design: Multi-Variant CITEgeist Xenium Benchmarking

**Date**: 2026-02-18
**Status**: Approved
**Author**: Claude (brainstorming session)

## Overview

Extend the existing Xenium benchmarking framework to compare all three CITEgeist variants (Continuous, Hybrid, Discrete) alongside other deconvolution methods, using both proportion and GEX metrics against protein-gated ground truth.

## Ground Truth

**Source**: `Benchmarking/xenium_pseudovisium/data_protein_gt/`

| Data Type | Location |
|-----------|----------|
| Proportions | `ground_truth/Xenium_region_{0-4}_prop.csv` |
| GEX | `ground_truth_gex/Xenium_region_{0-4}/{CellType}_GT.csv` |

**Cell Types (Achievable-7)**:
- B cells
- CD4+ T cells
- CD8+ T cells
- Macrophages
- Endothelial
- Epithelial
- Fibroblasts

**Fibroblasts Marker Definition**: alphaSMA only (VIM excluded due to ubiquity)

## CITEgeist Variants

| Variant | Output Path | Antibody Preprocessing |
|---------|-------------|------------------------|
| CITEgeist_Continuous | `CITEgeist/output/manual/` | `preprocess_antibody()` (CLR) |
| CITEgeist_Hybrid | `CITEgeist/output/hybrid_cellpose/` | `preprocess_antibody()` (CLR) |
| CITEgeist_Discrete | `CITEgeist/output_discrete_cellpose_fixed/` | `preprocess_antibody_discrete()` (per-marker scaling, no CLR) |

### Common Parameters (All Variants)

```python
CITEGEIST_COMMON_PARAMS = {
    "data_source": "data_protein_gt (protein-gated achievable-7)",
    "min_counts": 25,
    "gex_preprocessing": "filter_gex + preprocess_gex(target_sum=10000)",
    "cell_profile_dict": "ACHIEVABLE_7_CELL_PROFILE_DICT",
    "n_regions": 5,
}
```

### Preprocessing Difference (Documented)

The discrete variant uses different antibody preprocessing by design:

- **Continuous/Hybrid**: CLR normalization removes cellularity signal, proportions optimized via QP
- **Discrete**: Per-marker scaling preserves cellularity signal for IQP with nuclei count constraints

This difference is intentional and will be documented in the comparison output.

## Other Methods

| Method | Output Path |
|--------|-------------|
| Cell2Location | `Cell2Location/output_protein_gt/` |
| Tangram | `Tangram/output_protein_gt/` |
| RCTD | `RCTD/output_protein_gt/` |
| Seurat | `Seurat/output_protein_gt/` |
| CARD | `CARD/output_protein_gt_novim/` |

**Note**: CARD uses VIM-excluded results because VIM is pan-mesenchymal in kidney tissue.

## Metrics

### Proportion Evaluation
- Pearson correlation (r)
- RMSE
- MAE
- Jensen-Shannon Divergence (JSD)

### GEX Evaluation
- Per-cell-type Pearson correlation
- Overall mean across cell types

## Files to Modify

### 1. `compare_all_methods.py`

**Location**: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods.py`

**Changes**:
- Update METHODS dict to include all three CITEgeist variants
- Add METHOD_NOTES dict documenting preprocessing differences
- Add CITEGEIST_COMMON_PARAMS for reproducibility
- Print method notes in summary output
- Save notes to JSON

```python
METHODS = {
    # CITEgeist variants
    "CITEgeist_Continuous": BASE_DIR / "CITEgeist" / "output" / "manual",
    "CITEgeist_Hybrid": BASE_DIR / "CITEgeist" / "output" / "hybrid_cellpose",
    "CITEgeist_Discrete": BASE_DIR / "CITEgeist" / "output_discrete_cellpose_fixed",
    # Other methods
    "Cell2Location": BASE_DIR / "Cell2Location" / "output_protein_gt",
    "Tangram": BASE_DIR / "Tangram" / "output_protein_gt",
    "RCTD": BASE_DIR / "RCTD" / "output_protein_gt",
    "Seurat": BASE_DIR / "Seurat" / "output_protein_gt",
    "CARD": BASE_DIR / "CARD" / "output_protein_gt_novim",
}

METHOD_NOTES = {
    "CITEgeist_Continuous": "CLR antibody normalization, QP optimization",
    "CITEgeist_Hybrid": "CLR antibody normalization, QP optimization -> discretize via nuclei counts",
    "CITEgeist_Discrete": "Per-marker scaling (no CLR), IQP optimization with nuclei constraints",
}
```

### 2. `compare_all_methods_gex.py`

**Location**: `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py`

**Changes**:
- Add CITEGEIST_VARIANTS dict mapping variant names to relative paths
- Refactor `load_citegeist_gex()` to accept variant parameter
- Remove redundant `load_citegeist_discrete_gex()` function
- Update main comparison loop to iterate all CITEgeist variants
- Include same METHOD_NOTES in GEX output

```python
CITEGEIST_VARIANTS = {
    "CITEgeist_Continuous": "output/manual",
    "CITEgeist_Hybrid": "output/hybrid_cellpose",
    "CITEgeist_Discrete": "output_discrete_cellpose_fixed",
}

def load_citegeist_gex(region_id: int, variant: str = "CITEgeist_Continuous") -> Dict[str, pd.DataFrame]:
    """Load CITEgeist GEX layers for any variant."""
    sample_name = f"Xenium_region_{region_id}"
    variant_path = CITEGEIST_VARIANTS.get(variant)
    if variant_path is None:
        raise ValueError(f"Unknown CITEgeist variant: {variant}")

    base_dir = BASE_DIR / "CITEgeist" / variant_path
    layers_dir = base_dir / sample_name / f"{sample_name}_pass1" / "layers" / "pass1"

    if not layers_dir.exists():
        layers_dir = base_dir / sample_name / f"{sample_name}_pass1" / "layers"

    # ... rest of loading logic unchanged
```

## Output Format

### Proportion Comparison (`full_comparison.json`)

```
================================================================================
FULL METHOD COMPARISON - Protein-Gated Ground Truth (Achievable-7)
================================================================================

--- Method Notes ---
CITEgeist_Continuous: CLR antibody normalization, QP optimization
CITEgeist_Hybrid: CLR antibody normalization, QP optimization -> discretize via nuclei counts
CITEgeist_Discrete: Per-marker scaling (no CLR), IQP optimization with nuclei constraints

--- Common CITEgeist Parameters ---
data_source: data_protein_gt (protein-gated achievable-7)
min_counts: 25
gex_preprocessing: filter_gex + preprocess_gex(target_sum=10000)
cell_profile_dict: ACHIEVABLE_7_CELL_PROFILE_DICT (alphaSMA-only fibroblasts)

--- Overall Metrics ---
Method                Pearson r      RMSE       MAE       JSD
--------------------------------------------------------------
CITEgeist_Continuous     0.XXXX    0.XXXX    0.XXXX    0.XXXX
CITEgeist_Hybrid         0.XXXX    0.XXXX    0.XXXX    0.XXXX
CITEgeist_Discrete       0.XXXX    0.XXXX    0.XXXX    0.XXXX
Cell2Location            0.XXXX    0.XXXX    0.XXXX    0.XXXX
Tangram                  0.XXXX    0.XXXX    0.XXXX    0.XXXX
RCTD                     0.XXXX    0.XXXX    0.XXXX    0.XXXX
Seurat                   0.XXXX    0.XXXX    0.XXXX    0.XXXX
CARD                     0.XXXX    0.XXXX    0.XXXX    0.XXXX
```

### GEX Comparison (`full_comparison_gex.json`)

```
--- Per-Cell-Type GEX Pearson r ---
Cell Type           CITEgeist_Continuous  CITEgeist_Hybrid  CITEgeist_Discrete  ...
B cells                          0.XXX             0.XXX               0.XXX
CD4+ T cells                     0.XXX             0.XXX               0.XXX
CD8+ T cells                     0.XXX             0.XXX               0.XXX
Macrophages                      0.XXX             0.XXX               0.XXX
Endothelial                      0.XXX             0.XXX               0.XXX
Epithelial                       0.XXX             0.XXX               0.XXX
Fibroblasts                      0.XXX             0.XXX               0.XXX
```

## Verification Checklist

Before running the comparison:

- [ ] All 5 regions exist for CITEgeist_Continuous (`output/manual/`)
- [ ] All 5 regions exist for CITEgeist_Hybrid (`output/hybrid_cellpose/`)
- [ ] All 5 regions exist for CITEgeist_Discrete (`output_discrete_cellpose_fixed/`)
- [ ] Each region has `*_deconv_predictions.csv` (proportions)
- [ ] Each region has `*_pass1/layers/` directory (GEX)
- [ ] ACHIEVABLE_7_CELL_PROFILE_DICT uses alphaSMA-only for Fibroblasts

## Implementation Steps

1. Modify `compare_all_methods.py` with new METHODS dict and documentation
2. Modify `compare_all_methods_gex.py` with unified CITEgeist loader
3. Run proportion comparison via SLURM
4. Run GEX comparison via SLURM
5. Verify results in `evaluation/results/method_comparison/`
