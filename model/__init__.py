# Expose key classes and functions for easy access

from .citegeist_model import CitegeistModel
from .gurobi_impl import optimize_cell_proportions, optimize_gene_expression, map_antibodies_to_profiles
from .utils import (
    save_results_to_output,
    validate_cell_profile_dict,
    cleanup_memory,
    setup_logging
)

__all__ = [
    "CitegeistModel",
    "optimize_cell_proportions",
    "optimize_gene_expression",
    "map_antibodies_to_profiles",
    "save_results_to_output",
    "validate_cell_profile_dict",
    "cleanup_memory",
    "setup_logging"
]
