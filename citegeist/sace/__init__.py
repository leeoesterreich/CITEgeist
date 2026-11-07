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
    "map_antibodies_to_profiles": (".deconvolution.qp_solver", "map_antibodies_to_profiles"),
    "optimize_cell_proportions_per_marker": (".deconvolution.qp_solver", "optimize_cell_proportions_per_marker"),
    "optimize_gene_expression": (".deconvolution.qp_solver", "optimize_gene_expression"),
    "estimate_true_expression_cell": (".deconvolution.qp_solver", "estimate_true_expression_cell"),
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
    # Single-cell resolution: morphology features, GEX projection, AnnData output
    "extract_nucleus_features": (".morphology.morphology_features", "extract_nucleus_features"),
    "extract_cell_features": (".morphology.morphology_features", "extract_cell_features"),
    "largest_remainder_discretize": (".morphology.morphology_features", "largest_remainder_discretize"),
    "SoftLabelClassifier": (".morphology.soft_label_classifier", "SoftLabelClassifier"),
    "project_sace_to_cells": (".gex.sace_gex", "project_sace_to_cells"),
    "SaceInternals": (".gex.sace_gex", "SaceInternals"),
    "run_sace_layers": (".gex.sace_gex", "run_sace_layers"),
    # Module 3.5 benchmark/projection helpers
    "CoverageCheckResult": (".annotation.coverage_check", "CoverageCheckResult"),
    "check_module_coverage": (".annotation.coverage_check", "check_module_coverage"),
    # Per-cell GEX pipeline components
    "IdentifiabilityReport": (".annotation.subtype_splitting", "IdentifiabilityReport"),
    "audit_gate_identifiability": (".annotation.subtype_splitting", "audit_gate_identifiability"),
    "build_subtype_proportions": (".annotation.subtype_splitting", "build_subtype_proportions"),
    # Module 5: Cross-sample integration
    # QC framework
    # Per-type beta optimization
    "optimize_cell_proportions_per_type_beta": (".deconvolution.qp_solver", "optimize_cell_proportions_per_type_beta"),
    "MARKER_TYPE_TABLE": (".deconvolution.emission_init", "MARKER_TYPE_TABLE"),
    "CELL_TYPES": (".deconvolution.emission_init", "CELL_TYPES"),
    "build_marker_config": (".deconvolution.emission_init", "build_marker_config"),
    "initialize_beta_matrix": (".deconvolution.emission_init", "initialize_beta_matrix"),
    "build_beta_prior_sigma": (".deconvolution.emission_init", "build_beta_prior_sigma"),
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
            "dependencies needed for that feature — e.g. for imaging/segmentation support, "
            'run `pip install -e ".[imaging]"` from a cloned repo (CITEgeist is not yet on '
            'PyPI; once published, `pip install "citegeist[imaging]"` will also work).'
        ) from exc

    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
