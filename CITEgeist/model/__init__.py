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
    # Module 2c: Profile selection
    ProfileSelectionResult,
    select_profiles,
    # Module 2b enhancement: Hierarchical profile discovery
    ProfileTreeNode,
    ProfileTree,
    HierarchicalProfileResult,
    discover_hierarchical_profiles,
)

# Module 4: Protein-anchored program discovery
from .anchored_program_discovery import (
    SpatialSubpopulation,
    SpatialProgram,
    AnchoredProgramResult,
    AnchoredProgramDiscoveryResult,
    detect_spatial_subpopulations,
    discover_anchored_programs,  # Legacy: uses raw expression + contrastive
    discover_programs_from_layers,  # Recommended: uses deconvolved layers from Module 3
    store_results_in_adata,
    # Helper functions for deconvolved layers
    stack_deconvolved_layers,
    unstack_program_results,
    extract_celltype_expression,
    # Module 4b: Bivariate program relationships
    ProgramPairRelationship,
    BivariateProgramResult,
    analyze_program_relationships,
    # Module 4c: Region-aware program analysis
    analyze_program_regions,
    compare_programs_by_region,
    extract_program_context_genes,
)

# Module 5: Cross-sample integration
from .cross_sample_integration import (
    AlignedProgram,
    ConservedRelationship,
    IntegrationResult,
    load_multi_sample_results,
    integrate_samples,
    align_gene_sets,
    integrate_programs_harmony,
    match_programs_across_samples,
    compare_bivariate_relationships,
    build_similarity_network,
    save_integration_results,
)

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
    # Module 2c: Profile selection
    "ProfileSelectionResult",
    "select_profiles",
    # Module 2b enhancement: Hierarchical profile discovery
    "ProfileTreeNode",
    "ProfileTree",
    "HierarchicalProfileResult",
    "discover_hierarchical_profiles",
    # Module 4: Protein-anchored program discovery
    "SpatialSubpopulation",
    "SpatialProgram",
    "AnchoredProgramResult",
    "AnchoredProgramDiscoveryResult",
    "detect_spatial_subpopulations",
    "discover_anchored_programs",
    "discover_programs_from_layers",
    "store_results_in_adata",
    "stack_deconvolved_layers",
    "unstack_program_results",
    "extract_celltype_expression",
    # Module 4b: Bivariate program relationships
    "ProgramPairRelationship",
    "BivariateProgramResult",
    "analyze_program_relationships",
    # Module 4c: Region-aware program analysis
    "analyze_program_regions",
    "compare_programs_by_region",
    "extract_program_context_genes",
    # Module 5: Cross-sample integration
    "AlignedProgram",
    "ConservedRelationship",
    "IntegrationResult",
    "load_multi_sample_results",
    "integrate_samples",
    "align_gene_sets",
    "integrate_programs_harmony",
    "match_programs_across_samples",
    "compare_bivariate_relationships",
    "build_similarity_network",
    "save_integration_results",
]
