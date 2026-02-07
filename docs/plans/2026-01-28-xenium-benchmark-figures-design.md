# Design: Xenium Benchmark Visualization Figures

**Date**: 2026-01-28
**Status**: Approved

## Goal

Create two publication-quality figures demonstrating the Xenium benchmarking setup and CITEgeist's performance.

## Deliverables

### Figure 1: Benchmark Schematic (`benchmark_schematic.png`)

Horizontal flowchart showing the benchmark design:

1. **Xenium Tissue** — "10x Xenium In Situ" with note "single-cell resolution RNA"
2. **→ Pseudo-Visium Conversion** — "Aggregate cells into ~55μm spots", "Region 4: 1,414 spots"
3. **→ Ground Truth Definition** — "Achievable-7 cell types from single-cell annotations" listing: B cells, CD4+ T cells, CD8+ T cells, Macrophages, Endothelial, Epithelial, Fibroblasts
4. **→ 5 Deconvolution Methods** — stacked: CITEgeist, Cell2Location, RCTD, Tangram, Seurat
5. **→ Evaluation** — "Pearson r, RMSE, MAE, JSD"

Implementation: matplotlib patches, FancyArrowPatch, text annotations. No external diagramming libraries.

### Figure 2: Spatial Comparison (`benchmark_spatial_comparison.png`)

Two-panel side-by-side spatial pie chart plots for Region 4:

- **Left panel**: Ground Truth (achievable-7 proportions)
- **Right panel**: CITEgeist Global (Pass 1) predictions
- **Shared legend**: 7 colors for achievable-7 cell types
- **Spot visualization**: Mini pie charts (~3-5pt radius) at all 1,414 spot positions
- **Coordinates**: from `obsm['spatial']` in the GEX h5ad file

## Data Sources

| Data | Path |
|------|------|
| Granular GT (10 types) | `Benchmarking/xenium_pseudovisium/data_granular_gt/ground_truth/Xenium_region_4_prop.csv` |
| CITEgeist predictions | `Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7/Xenium_region_4/Xenium_region_4_cell_prop_global_results.csv` |
| Spatial coordinates | `Benchmarking/xenium_pseudovisium/h5ad_objects/Xenium_region_4_GEX.h5ad` → `obsm['spatial']` |

## Achievable-7 Collapse Mapping

From 10 granular types:

```python
GT_TO_ACHIEVABLE_7 = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",
}
```

## Script Location

`Benchmarking/xenium_benchmarking/figures/plot_benchmark_overview.py`

## Dependencies

matplotlib, numpy, pandas, scanpy — all available in CITEgeist_env.
