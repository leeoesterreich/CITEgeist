# Protein vs RNA Cell Type Concordance Analysis

**Date**: 2026-02-20

## Summary

Overall cell-level concordance: **56.7%**

## Per-Cell-Type Concordance

| Cell Type | Concordance (%) |
|-----------|-----------------|
| B cells | 65.2 |
| T cells | 68.7 |
| Macrophages | 83.1 |
| Endothelial | 67.4 |
| Epithelial | 67.8 |
| Fibroblasts | 82.1 |
| Unknown | 0.0 |

## Spot-Level Correlations

| Cell Type | Pearson r | p-value |
|-----------|-----------|---------|
| B cells | 0.950 | 0.00e+00 |
| T cells | 0.847 | 0.00e+00 |
| Macrophages | 0.702 | 0.00e+00 |
| Endothelial | 0.813 | 0.00e+00 |
| Epithelial | 0.624 | 0.00e+00 |
| Fibroblasts | 0.611 | 0.00e+00 |

## Interpretation

**Dual benchmarks required** for fair comparison (<60% concordance).

## Figures

- `concordance_heatmap.png` - Cell-level confusion matrix
- `spot_correlations.png` - Spot-level proportion scatter plots
- `concordance_matrix.csv` - Raw confusion matrix
- `spot_level_correlations.csv` - Correlation values
