"""
Xenium Benchmarking Framework for CITEgeist

This module provides tools to benchmark CITEgeist on real Xenium data
by creating pseudo-Visium spots from single-cell resolution data.
"""

from .load_xenium import load_xenium_data, split_gex_protein
from .create_pseudo_spots import (
    create_hexagonal_grid,
    assign_cells_to_spots,
    aggregate_counts_per_spot,
)
from .define_cell_types import (
    XENIUM_CELL_PROFILE_DICT,
    classify_cells_by_protein,
)
from .generate_ground_truth import calculate_spot_proportions
from .split_regions import split_tissue_regions, save_region_datasets

__all__ = [
    "load_xenium_data",
    "split_gex_protein",
    "create_hexagonal_grid",
    "assign_cells_to_spots",
    "aggregate_counts_per_spot",
    "XENIUM_CELL_PROFILE_DICT",
    "classify_cells_by_protein",
    "calculate_spot_proportions",
    "split_tissue_regions",
    "save_region_datasets",
]
