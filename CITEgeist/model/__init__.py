"""
CITEgeist model package - spatial transcriptomics deconvolution using CITE-seq.
"""
# Expose key classes and functions for easy access

from .citegeist_model import CitegeistModel
from .gurobi_impl import (
    map_antibodies_to_profiles,
    optimize_cell_proportions,
    optimize_gene_expression,
)
from .utils import cleanup_memory, save_results_to_output, setup_logging, validate_cell_profile_dict

# New marker analysis modules
from .marker_interest import MarkerInterest, MarkerInterestResult, identify_interesting_markers
from .spatial_colocalization import (
    MarkerPairColocalization,
    ColocalizationResult,
    analyze_marker_colocalization,
    LineageDendrogram,
    ProfileDiscoveryResult,
    discover_profiles,
    ProfileSelectionResult,
    select_profiles_by_reconstruction,
    SpatialVarianceSelectionResult,
    select_profiles_by_spatial_variance,
)

# Default profile selection method (spatial variance-based)
select_profiles = select_profiles_by_spatial_variance

__all__ = [
    # Core model
    "CitegeistModel",
    # Gurobi optimization
    "map_antibodies_to_profiles",
    "optimize_cell_proportions",
    "optimize_gene_expression",
    # Utilities
    "cleanup_memory",
    "save_results_to_output",
    "setup_logging",
    "validate_cell_profile_dict",
    # Marker interest analysis
    "MarkerInterest",
    "MarkerInterestResult",
    "identify_interesting_markers",
    # Spatial colocalization analysis
    "MarkerPairColocalization",
    "ColocalizationResult",
    "analyze_marker_colocalization",
    # Profile discovery
    "LineageDendrogram",
    "ProfileDiscoveryResult",
    "discover_profiles",
    # Profile selection (Module 2c)
    "ProfileSelectionResult",
    "select_profiles_by_reconstruction",
    "SpatialVarianceSelectionResult",
    "select_profiles_by_spatial_variance",
    "select_profiles",  # Default alias for spatial variance method
]
