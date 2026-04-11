"""
CITEgeist model package - spatial transcriptomics deconvolution using CITE-seq.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS = {
    # Core model
    "CitegeistModel": (".citegeist_model", "CitegeistModel"),
    "RESOLUTION_DEFAULTS": (".citegeist_model", "RESOLUTION_DEFAULTS"),
    # Optimization (cuOPT backend)
    "map_antibodies_to_profiles": (".deconvolution.cuopt_impl", "map_antibodies_to_profiles"),
    "optimize_cell_proportions_per_marker": (".deconvolution.cuopt_impl", "optimize_cell_proportions_per_marker"),
    "optimize_gene_expression": (".deconvolution.cuopt_impl", "optimize_gene_expression"),
    "estimate_true_expression_cell": (".deconvolution.cuopt_impl", "estimate_true_expression_cell"),
    # GEX module-aware enrichment
    "discover_anchor_genes": (".gex.gex_modules", "discover_anchor_genes"),
    "compute_module_aware_enrichment": (".gex.gex_modules", "compute_module_aware_enrichment"),
    "compute_softmax_target": (".gex.gex_modules", "compute_softmax_target"),
    "compute_kl_penalty_coefficients": (".gex.gex_modules", "compute_kl_penalty_coefficients"),
    # Cell type detection
    "detect_cell_types": (".deconvolution.detection", "detect_cell_types"),
    # Utilities
    "cleanup_memory": (".utils", "cleanup_memory"),
    "save_results_to_output": (".utils", "save_results_to_output"),
    "setup_logging": (".utils", "setup_logging"),
    "validate_cell_profile_dict": (".utils", "validate_cell_profile_dict"),
    # Marker interest analysis
    "MarkerInterest": (".discovery.marker_interest", "MarkerInterest"),
    "MarkerInterestResult": (".discovery.marker_interest", "MarkerInterestResult"),
    "identify_interesting_markers": (".discovery.marker_interest", "identify_interesting_markers"),
    # Spatial colocalization analysis
    "MarkerPairColocalization": (".discovery.spatial_colocalization", "MarkerPairColocalization"),
    "ColocalizationResult": (".discovery.spatial_colocalization", "ColocalizationResult"),
    "analyze_marker_colocalization": (".discovery.spatial_colocalization", "analyze_marker_colocalization"),
    # Profile discovery
    "LineageDendrogram": (".discovery.spatial_colocalization", "LineageDendrogram"),
    "ProfileDiscoveryResult": (".discovery.spatial_colocalization", "ProfileDiscoveryResult"),
    "discover_profiles": (".discovery.spatial_colocalization", "discover_profiles_continuous"),
    "discover_profiles_continuous": (".discovery.spatial_colocalization", "discover_profiles_continuous"),
    "rescue_singletons": (".discovery.spatial_colocalization", "rescue_singletons"),
    "ProfileSelectionResult": (".discovery.spatial_colocalization", "ProfileSelectionResult"),
    "select_profiles": (".discovery.spatial_colocalization", "select_profiles"),
    "select_profiles_by_reconstruction": (".discovery.spatial_colocalization", "select_profiles"),
    "ProfileTreeNode": (".discovery.spatial_colocalization", "ProfileTreeNode"),
    "ProfileTree": (".discovery.spatial_colocalization", "ProfileTree"),
    "HierarchicalProfileResult": (".discovery.spatial_colocalization", "HierarchicalProfileResult"),
    "discover_hierarchical_profiles_continuous": (
        ".discovery.spatial_colocalization",
        "discover_hierarchical_profiles_continuous",
    ),
    # Module 4: Protein-anchored program discovery
    "SpatialSubpopulation": (".programs.anchored_program_discovery", "SpatialSubpopulation"),
    "SpatialProgram": (".programs.anchored_program_discovery", "SpatialProgram"),
    "AnchoredProgramResult": (".programs.anchored_program_discovery", "AnchoredProgramResult"),
    "AnchoredProgramDiscoveryResult": (".programs.anchored_program_discovery", "AnchoredProgramDiscoveryResult"),
    "detect_spatial_subpopulations": (".programs.anchored_program_discovery", "detect_spatial_subpopulations"),
    "discover_anchored_programs": (".programs.anchored_program_discovery", "discover_anchored_programs"),
    "discover_programs_from_layers": (".programs.anchored_program_discovery", "discover_programs_from_layers"),
    "store_results_in_adata": (".programs.anchored_program_discovery", "store_results_in_adata"),
    "stack_deconvolved_layers": (".programs.anchored_program_discovery", "stack_deconvolved_layers"),
    "unstack_program_results": (".programs.anchored_program_discovery", "unstack_program_results"),
    "extract_celltype_expression": (".programs.anchored_program_discovery", "extract_celltype_expression"),
    "ProgramPairRelationship": (".programs.anchored_program_discovery", "ProgramPairRelationship"),
    "BivariateProgramResult": (".programs.anchored_program_discovery", "BivariateProgramResult"),
    "analyze_program_relationships": (".programs.anchored_program_discovery", "analyze_program_relationships"),
    "analyze_program_regions": (".programs.anchored_program_discovery", "analyze_program_regions"),
    "compare_programs_by_region": (".programs.anchored_program_discovery", "compare_programs_by_region"),
    "extract_program_context_genes": (".programs.anchored_program_discovery", "extract_program_context_genes"),
    "JointProgram": (".programs.anchored_program_discovery", "JointProgram"),
    "JointDiscoveryResult": (".programs.anchored_program_discovery", "JointDiscoveryResult"),
    "discover_joint_programs": (".programs.anchored_program_discovery", "discover_joint_programs"),
    # Module 3b: Single-cell resolution
    "extract_nucleus_features": (".morphology.morphology_features", "extract_nucleus_features"),
    "extract_cell_features": (".morphology.morphology_features", "extract_cell_features"),
    "largest_remainder_discretize": (".morphology.morphology_features", "largest_remainder_discretize"),
    "SoftLabelClassifier": (".morphology.soft_label_classifier", "SoftLabelClassifier"),
    "assign_nuclei_to_types": (".assignment.hungarian_assignment", "assign_nuclei_to_types"),
    "run_nucleus_assignment": (".assignment.module3b_nucleus_assignment", "run_nucleus_assignment"),
    "run_nucleus_assignment_mil": (".assignment.module3b_nucleus_assignment", "run_nucleus_assignment_mil"),
    "run_nucleus_assignment_mil_em": (".assignment.module3b_nucleus_assignment", "run_nucleus_assignment_mil_em"),
    "NucleusAssignmentResult": (".assignment.module3b_nucleus_assignment", "NucleusAssignmentResult"),
    "distribute_gex_to_cells": (".gex.cell_level_gex", "distribute_gex_to_cells"),
    "allocate_gex_type_reference": (".gex.cell_level_gex", "allocate_gex_type_reference"),
    "create_single_cell_adata": (".assignment.single_cell_output", "create_single_cell_adata"),
    # Module 3.5 benchmark/projection helpers
    "aggregate_module3_5_results": (".annotation.module3_5_benchmark", "aggregate_module3_5_results"),
    "CoverageCheckResult": (".annotation.coverage_check", "CoverageCheckResult"),
    "check_module_coverage": (".annotation.coverage_check", "check_module_coverage"),
    "should_enrich_module3_5_output": (".annotation.module3_5_projection", "should_enrich_module3_5_output"),
    "should_enrich_single_cell_output": (".annotation.module3_5_projection", "should_enrich_single_cell_output"),
    "normalized_validated_pairs": (".annotation.module3_5_projection", "normalized_validated_pairs"),
    # Per-cell GEX pipeline components
    "MarkerGenes": (".gex.gex_modules", "MarkerGenes"),
    "GeneModule": (".gex.gex_modules", "GeneModule"),
    "GeneModules": (".gex.gex_modules", "GeneModules"),
    "discover_marker_genes": (".gex.gex_modules", "discover_marker_genes"),
    "build_gene_modules": (".gex.gex_modules", "build_gene_modules"),
    "IdentifiabilityReport": (".annotation.subtype_splitting", "IdentifiabilityReport"),
    "audit_gate_identifiability": (".annotation.subtype_splitting", "audit_gate_identifiability"),
    "build_subtype_proportions": (".annotation.subtype_splitting", "build_subtype_proportions"),
    # Module 5: Cross-sample integration
    "AlignedProgram": (".programs.cross_sample_integration", "AlignedProgram"),
    "ConservedRelationship": (".programs.cross_sample_integration", "ConservedRelationship"),
    "IntegrationResult": (".programs.cross_sample_integration", "IntegrationResult"),
    "load_multi_sample_results": (".programs.cross_sample_integration", "load_multi_sample_results"),
    "integrate_samples": (".programs.cross_sample_integration", "integrate_samples"),
    "align_gene_sets": (".programs.cross_sample_integration", "align_gene_sets"),
    "integrate_programs_harmony": (".programs.cross_sample_integration", "integrate_programs_harmony"),
    "match_programs_across_samples": (".programs.cross_sample_integration", "match_programs_across_samples"),
    "compare_bivariate_relationships": (".programs.cross_sample_integration", "compare_bivariate_relationships"),
    "build_similarity_network": (".programs.cross_sample_integration", "build_similarity_network"),
    "save_integration_results": (".programs.cross_sample_integration", "save_integration_results"),
    # QC framework
    "QCResult": (".qc", "QCResult"),
    "run_qc": (".qc", "run_qc"),
    # Per-type beta optimization
    "optimize_cell_proportions_per_type_beta": (".deconvolution.cuopt_impl", "optimize_cell_proportions_per_type_beta"),
    "MARKER_TYPE_TABLE": (".deconvolution.emission_init", "MARKER_TYPE_TABLE"),
    "CELL_TYPES": (".deconvolution.emission_init", "CELL_TYPES"),
    "build_marker_config": (".deconvolution.emission_init", "build_marker_config"),
    "initialize_beta_matrix": (".deconvolution.emission_init", "initialize_beta_matrix"),
    "build_beta_prior_sigma": (".deconvolution.emission_init", "build_beta_prior_sigma"),
    # Vision/SSL backbone
    "MorphologyBackbone": (".morphology.morphology_backbone", "MorphologyBackbone"),
    "DAPIBackbone": (".morphology.morphology_backbone", "DAPIBackbone"),
    "HEBackbone": (".morphology.morphology_backbone", "HEBackbone"),
    "ViTFeatureExtractor": (".morphology.vit_extractor", "ViTFeatureExtractor"),
    "load_uni_extractor": (".morphology.vit_extractor", "load_uni_extractor"),
    "ViTEncoder": (".morphology.vit_encoder", "ViTEncoder"),
    "create_vit_small": (".morphology.vit_encoder", "create_vit_small"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    try:
        module = import_module(module_name, __name__)
    except ImportError as exc:
        raise ImportError(
            f"Failed to import '{name}' from '{module_name}'. Install the optional "
            "dependencies needed for that feature, or use `pip install -e .[dev]` "
            "for the standard development environment."
        ) from exc

    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
