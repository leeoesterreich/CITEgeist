# Protein vs RNA Cell Type Concordance Analysis

**Date**: 2026-02-19
**Status**: Design approved
**Goal**: Investigate whether protein-gated cell types match RNA-defined cell types in Xenium data

## Background

The current Xenium benchmarking uses **protein-gated ground truth**: single cells are classified by hierarchical protein marker gating (CD68+ → Macrophage, etc.), then aggregated to spot-level proportions.

**Concern**: This may unfairly bias toward CITEgeist (which deconvolves the same protein markers) vs RNA-based methods (Cell2Location, Tangram, RCTD, Seurat) that predict RNA-defined cell types.

**Core question**: Are protein-defined cell types the same cells as RNA-defined cell types?

## Data Sources

| Data | Source | Format |
|------|--------|--------|
| Protein expression | `cell_feature_matrix.h5` | 465K cells × 22 proteins |
| Protein-gated cell types | `create_protein_gt.py` logic | Hierarchical gating |
| RNA k-means clusters | `analysis.tar.gz` → `gene_expression_kmeans_10_clusters/clusters.csv` | 10 clusters |
| RNA cluster annotations | `rna_cell_types.py` → `XENIUM_KIDNEY_RNA_CLUSTER_MAP` | Manual annotation |

**Cell types to compare** (7 common categories):
- B cells
- CD4+ T cells
- CD8+ T cells
- Macrophages
- Endothelial
- Epithelial
- Fibroblasts

## Methodology

### Step 1: Cell-Level Concordance Matrix

1. Apply protein gating to all 465K cells → `protein_cell_type`
2. Load RNA k-means clusters → map via `XENIUM_KIDNEY_RNA_CLUSTER_MAP` → `rna_cell_type`
3. Build confusion matrix: rows = protein types, columns = RNA types
4. Calculate per-cell-type concordance rate: `(agree / total)`

### Step 2: Disagreement Characterization

For cells where protein ≠ RNA assignment:
- Which mismatches are most common (e.g., protein=Macrophage, RNA=Mixed Immune)?
- Examine protein marker expression in discordant cells (are they borderline?)

### Step 3: Spot-Level Impact

1. Aggregate both cell type assignments to pseudo-Visium spots
2. Correlate protein-GT proportions vs RNA-GT proportions per cell type
3. Quantifies: "How much would benchmark results change if we used RNA-GT?"

## Output Deliverables

1. **Analysis script**: `Benchmarking/xenium_pseudovisium/src/analyze_celltype_concordance.py`
2. **Summary report**: `docs/analysis/protein_vs_rna_celltype_concordance.md`
3. **Figures**:
   - Confusion matrix heatmap (cell counts)
   - Per-cell-type agreement rates table
   - Spot-level correlation scatter plots

## Interpretation Guidelines

| Overall Concordance | Interpretation |
|---------------------|----------------|
| >80% | Current protein-GT benchmark is defensible |
| 60-80% | Consider reporting both benchmarks in supplement |
| <60% | Dual benchmarks required for fair comparison |

## Out of Scope (Future Work)

- Re-running all methods against RNA-GT
- Creating harmonized consensus GT (protein AND RNA agreement)
- Modifying main benchmark pipeline

## References

- `Benchmarking/xenium_pseudovisium/src/create_protein_gt.py` - current protein gating logic
- `Benchmarking/xenium_pseudovisium/src/rna_cell_types.py` - RNA cluster mapping
- `Benchmarking/xenium_pseudovisium/src/define_cell_types_validated.py` - validated profiles
- Zhao et al. (2025). BMC Bioinformatics. "Benchmarking cell type annotation methods for 10x Xenium"
