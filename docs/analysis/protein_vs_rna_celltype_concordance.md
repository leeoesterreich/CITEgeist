# Protein vs RNA Cell Type Concordance Analysis

**Date**: 2026-02-19

## Summary

Overall cell-level concordance: **48.1%**

## Per-Cell-Type Concordance

| Cell Type | Concordance (%) |
|-----------|-----------------|
| B cells | 65.2 |
| CD4+ T cells | 1.6 |
| CD8+ T cells | 64.3 |
| Macrophages | 83.1 |
| Endothelial | 67.4 |
| Epithelial | 65.2 |
| Fibroblasts | 82.1 |
| Unknown | 0.0 |

## Spot-Level Correlations

| Cell Type | Pearson r | p-value |
|-----------|-----------|---------|
| B cells | 0.950 | 0.00e+00 |
| CD4+ T cells | 0.092 | 7.18e-15 |
| CD8+ T cells | 0.610 | 0.00e+00 |
| Macrophages | 0.702 | 0.00e+00 |
| Endothelial | 0.813 | 0.00e+00 |
| Epithelial | 0.620 | 0.00e+00 |
| Fibroblasts | 0.611 | 0.00e+00 |

## Interpretation

**Dual benchmarks required** for fair comparison (<60% concordance).

## Figures

- `concordance_heatmap.png` - Cell-level confusion matrix
- `spot_correlations.png` - Spot-level proportion scatter plots
- `concordance_matrix.csv` - Raw confusion matrix
- `spot_level_correlations.csv` - Correlation values
