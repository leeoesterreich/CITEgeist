"""
RNA-only pipeline for CITEgeist.

This module provides a unified interface for running the CITEgeist discovery
pipeline on RNA-only spatial transcriptomics data (no protein required).

The pipeline uses marker genes (curated or auto-discovered) to:
1. Identify spatially-interesting markers (Module 1)
2. Discover cell population profiles via colocalization (Module 2)
3. Assign cells/spots to profiles
4. Discover transcriptional programs anchored to marker-defined populations (Module 4)

This demonstrates CITEgeist as a general spatial analysis framework,
with protein being optimal but not required.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc

from .rna_marker_selection import (
    MarkerMode,
    RNAMarkerSelectionResult,
    select_rna_markers,
    get_curated_markers,
)
from .marker_interest import identify_interesting_markers, MarkerInterestResult
from .spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    ColocalizationResult,
    ProfileDiscoveryResult,
)

logger = logging.getLogger(__name__)


@dataclass
class RNAPipelineConfig:
    """Configuration for RNA-only pipeline."""

    # Marker selection
    marker_mode: MarkerMode = MarkerMode.HYBRID
    curated_markers: Optional[List[str]] = None
    n_hvgs: int = 2000
    max_discovered_markers: int = 50
    max_total_markers: int = 100
    redundancy_threshold: float = 0.7

    # Spatial parameters
    morans_k: int = 8
    smooth_k: int = 6
    neighbor_k: int = 15  # For colocalization

    # Profile discovery
    fdr_alpha: float = 0.05
    min_profile_size: int = 2
    max_profiles: int = 15

    # General
    n_permutations: int = 99
    seed: int = 42


@dataclass
class RNAPipelineResult:
    """Results from RNA-only pipeline."""

    # Marker selection
    marker_selection: RNAMarkerSelectionResult

    # Module 1: Spatial interest
    marker_interest: MarkerInterestResult

    # Module 2: Profiles
    colocalization: ColocalizationResult
    profiles: ProfileDiscoveryResult

    # Summary
    n_markers_selected: int
    n_interesting_markers: int
    n_profiles_discovered: int

    def summary(self) -> str:
        """Human-readable summary."""
        lines = [
            "RNA-only Pipeline Results",
            "=" * 40,
            f"Markers selected: {self.n_markers_selected}",
            f"  - Curated: {len(self.marker_selection.curated_markers)}",
            f"  - Discovered: {len(self.marker_selection.discovered_markers)}",
            f"Spatially interesting: {self.n_interesting_markers}",
            f"Profiles discovered: {self.n_profiles_discovered}",
            "",
            "Profiles:",
        ]
        for i, profile in enumerate(self.profiles.profiles):
            lines.append(f"  {i}: {profile}")

        return "\n".join(lines)


def run_rna_pipeline(
    adata: sc.AnnData,
    config: Optional[RNAPipelineConfig] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> RNAPipelineResult:
    """
    Run the full CITEgeist discovery pipeline on RNA-only data.

    This is the main entry point for RNA-only analysis. It:
    1. Selects marker genes (curated and/or auto-discovered)
    2. Runs Module 1 (spatial interest detection) on markers
    3. Runs Module 2 (colocalization-based profile discovery)

    The discovered profiles can then be used for:
    - Spot/cell assignment (Module 3 adaptation)
    - Transcriptional program discovery (Module 4)

    Args:
        adata: AnnData with:
            - .X: Gene expression (normalized, log-transformed recommended)
            - .var_names: Gene names
            - .obsm['spatial']: Spatial coordinates
        config: Pipeline configuration. If None, uses defaults.
        output_dir: Optional directory to save results.

    Returns:
        RNAPipelineResult with all intermediate and final results.

    Example:
        >>> # Basic usage with defaults
        >>> result = run_rna_pipeline(adata)
        >>> print(result.summary())
        >>>
        >>> # Custom configuration
        >>> config = RNAPipelineConfig(
        ...     marker_mode=MarkerMode.HYBRID,
        ...     curated_markers=["CD3E", "CD68", "EPCAM"],
        ...     max_discovered_markers=30,
        ... )
        >>> result = run_rna_pipeline(adata, config=config)
    """
    if config is None:
        config = RNAPipelineConfig()

    logger.info("=" * 60)
    logger.info("CITEgeist RNA-only Pipeline")
    logger.info("=" * 60)

    # Validate input
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] required for spatial analysis")

    # Step 1: Marker Selection
    logger.info("\n[Step 1] Marker Selection")
    logger.info("-" * 40)

    marker_result = select_rna_markers(
        adata=adata,
        mode=config.marker_mode,
        curated_markers=config.curated_markers,
        n_hvgs=config.n_hvgs,
        max_discovered=config.max_discovered_markers,
        max_total=config.max_total_markers,
        redundancy_threshold=config.redundancy_threshold,
        morans_k=config.morans_k,
        smooth_k=config.smooth_k,
        n_permutations=config.n_permutations,
        seed=config.seed,
    )

    selected_markers = marker_result.selected_markers
    logger.info(f"Selected {len(selected_markers)} markers for analysis")

    if len(selected_markers) < 3:
        raise ValueError(
            f"Only {len(selected_markers)} markers selected. "
            "Need at least 3 for profile discovery. "
            "Try adding more curated markers or using autodiscovery mode."
        )

    # Step 2: Module 1 - Spatial Interest Detection
    logger.info("\n[Step 2] Module 1: Spatial Interest Detection")
    logger.info("-" * 40)

    # Subset to selected markers
    adata_markers = adata[:, selected_markers].copy()

    X_markers = adata_markers.X
    if hasattr(X_markers, "toarray"):
        X_markers = X_markers.toarray()
    X_markers = np.asarray(X_markers, dtype=np.float64)

    coords = adata_markers.obsm["spatial"]

    interest_result = identify_interesting_markers(
        X=X_markers,
        coords=coords,
        marker_names=selected_markers,
        morans_k=config.morans_k,
        smooth_k=config.smooth_k,
        morans_n_perm=config.n_permutations,
        seed=config.seed,
        verbose=True,
    )

    interesting_markers = interest_result.interesting_markers
    logger.info(f"Found {len(interesting_markers)} spatially interesting markers")

    if len(interesting_markers) < 3:
        logger.warning(
            f"Only {len(interesting_markers)} interesting markers found. "
            "Profile discovery may be limited. Consider relaxing thresholds."
        )
        # Fall back to top markers by score
        interesting_markers = [m.name for m in interest_result.markers[:10]]

    # Step 3: Module 2 - Colocalization and Profile Discovery
    logger.info("\n[Step 3] Module 2: Colocalization & Profile Discovery")
    logger.info("-" * 40)

    # Module 2a: Colocalization analysis
    logger.info("Running colocalization analysis...")
    coloc_result = analyze_marker_colocalization(
        X=X_markers,
        coords=coords,
        marker_names=selected_markers,
        markers_to_analyze=interesting_markers,
        neighbor_k=config.neighbor_k,
        multi_scale_k=[config.neighbor_k // 2, config.neighbor_k, config.neighbor_k * 2],
        n_permutations=config.n_permutations,
        seed=config.seed,
    )
    logger.info(f"Analyzed {len(coloc_result.pairs)} marker pairs")

    # Module 2b: Profile discovery
    logger.info("Discovering profiles via hierarchical clustering...")
    profile_result = discover_profiles(
        colocalization_result=coloc_result,
        fdr_alpha=config.fdr_alpha,
        min_cluster_size=config.min_profile_size,
        top_k=config.max_profiles,
        seed=config.seed,
    )

    n_profiles = len(profile_result.profiles)
    logger.info(f"Discovered {n_profiles} marker profiles")

    for i, profile in enumerate(profile_result.profiles):
        logger.info(f"  Profile {i}: {profile}")

    # Save results if output_dir specified
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Marker selection
        marker_result.to_dataframe().to_csv(
            output_dir / "marker_selection.csv", index=False
        )

        # Module 1 results
        interest_result.to_dataframe().to_csv(
            output_dir / "marker_interest.csv", index=False
        )

        # Profile discovery
        profile_df = pd.DataFrame({
            "profile_id": range(n_profiles),
            "markers": [",".join(p) for p in profile_result.profiles],
            "n_markers": [len(p) for p in profile_result.profiles],
        })
        profile_df.to_csv(output_dir / "profiles.csv", index=False)

        logger.info(f"Results saved to {output_dir}")

    # Build result
    result = RNAPipelineResult(
        marker_selection=marker_result,
        marker_interest=interest_result,
        colocalization=coloc_result,
        profiles=profile_result,
        n_markers_selected=len(selected_markers),
        n_interesting_markers=len(interesting_markers),
        n_profiles_discovered=n_profiles,
    )

    logger.info("\n" + "=" * 60)
    logger.info("Pipeline complete!")
    logger.info("=" * 60)
    print(result.summary())

    return result


def compare_protein_vs_rna_profiles(
    adata_protein: sc.AnnData,
    adata_rna: sc.AnnData,
    protein_markers: List[str],
    rna_markers: List[str],
    **kwargs,
) -> pd.DataFrame:
    """
    Compare profiles discovered from protein vs RNA markers.

    Useful for validating that RNA-based profile discovery
    recovers similar structure to protein-based discovery.

    Args:
        adata_protein: AnnData with protein expression.
        adata_rna: AnnData with RNA expression.
        protein_markers: Protein marker names.
        rna_markers: RNA marker gene names.
        **kwargs: Additional arguments for pipeline.

    Returns:
        DataFrame comparing protein and RNA profiles.
    """
    # This is a stub for future implementation
    raise NotImplementedError(
        "Protein vs RNA comparison not yet implemented. "
        "See benchmarking scripts for manual comparison."
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    print("=" * 60)
    print("RNA Pipeline Module")
    print("=" * 60)
    print("\nThis module provides run_rna_pipeline() for RNA-only analysis.")
    print("\nExample usage:")
    print("""
    from CITEgeist.model.rna_pipeline import run_rna_pipeline, RNAPipelineConfig

    # Load your spatial RNA data
    adata = sc.read_h5ad("your_data.h5ad")

    # Run with defaults (hybrid mode)
    result = run_rna_pipeline(adata)

    # Or with custom config
    config = RNAPipelineConfig(
        marker_mode=MarkerMode.CURATED_ONLY,
        curated_markers=["CD3E", "CD68", "EPCAM", "COL1A1"],
    )
    result = run_rna_pipeline(adata, config=config)

    # Access results
    print(result.summary())
    profiles = result.profiles.profiles  # List of marker lists
    """)
