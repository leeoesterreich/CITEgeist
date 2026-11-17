"""
Xenium Pseudo-Visium Data Generation.

This module contains scripts for creating pseudo-Visium spots
from Xenium single-cell resolution data for benchmarking purposes.

Modules:
    load_xenium: Load raw Xenium single-cell data
    create_pseudo_spots: Create hexagonal grid and aggregate cells
    split_regions: Divide into spatial regions for cross-validation
    define_cell_types: Protein-based cell type classification
    rna_cell_types: RNA-based cell type classification
    generate_ground_truth: Calculate cell proportion ground truth
    generate_expression_weighted_gt: Expression-weighted proportions
    generate_gex_ground_truth: Gene expression ground truth per cell type
"""
