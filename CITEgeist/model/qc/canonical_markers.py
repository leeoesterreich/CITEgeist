"""
Canonical RNA markers for cell type validation.

These are RNA-level markers chosen for known protein-RNA concordance,
separate from the antibody panel markers used in QP deconvolution.
Cancer is excluded — its markers are tissue-specific.
"""

__all__ = ["CANONICAL_MARKERS", "get_available_markers"]

CANONICAL_MARKERS: dict[str, list[str]] = {
    "Macrophages": ["CD68", "CD163", "CSF1R", "MSR1", "MRC1"],
    "CD8_T_Cells": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7"],
    "CD4_T_Cells": ["CD4", "IL7R", "TCF7", "LEF1", "CCR7"],
    "B_Cells": ["MS4A1", "CD79A", "CD79B", "PAX5", "BANK1"],
    "Epithelial": ["EPCAM", "KRT18", "KRT8", "CDH1", "MUC1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "FLT1", "ERG"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "Monocytes": ["CD14", "FCGR3A", "S100A8", "S100A9", "VCAN"],
    "Dendritic_Cells": ["ITGAX", "FLT3", "IRF8", "BATF3", "CLEC9A"],
}


def get_available_markers(cell_type: str, gene_names: list[str]) -> list[str]:
    """Return canonical markers for a cell type that exist in the gene list.

    Args:
        cell_type: Cell type name (must match CANONICAL_MARKERS keys).
        gene_names: List of available gene names in the dataset.

    Returns:
        List of marker genes present in both the canonical table and gene_names.
    """
    markers = CANONICAL_MARKERS.get(cell_type, [])
    gene_set = set(gene_names)
    return [m for m in markers if m in gene_set]
