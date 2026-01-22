"""
Run CITEgeist deconvolution on Xenium pseudo-Visium data.

This module provides wrapper functions to run CITEgeist on the
pseudo-Visium datasets generated from Xenium single-cell data.
"""

import os
import sys
import argparse
import logging
import time
import json
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Add CITEgeist to path
# Path: Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py -> 5 levels to repo root
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

from model.citegeist_model import CitegeistModel
from model.marker_interest import identify_interesting_markers
from model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
)

logger = logging.getLogger(__name__)

# Cell profile dictionary for Xenium proteins (matching define_cell_types.py)
# Uses the same markers for deconvolution as for ground truth generation
XENIUM_CELL_PROFILE_DICT = {
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD45", "GranzymeB"],
    },
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45"],
    },
    "Plasma cells": {
        "Major": ["CD138"],
        "Minor": [],
    },
    "Macrophages": {
        "Major": ["CD68"],
        "Minor": ["CD163", "CD45", "HLA-DR"],
    },
    "Dendritic cells": {
        "Major": ["CD11c", "HLA-DR"],
        "Minor": ["CD45"],
    },
    "NK cells": {
        "Major": ["CD16", "GranzymeB"],
        "Minor": ["CD45"],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": ["E-Cadherin"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": ["Vimentin"],
    },
}

# RNA-based cell profile dictionary (matching RNA clustering ground truth)
# These 6 cell types match the simplified RNA-cluster-based ground truth.
# Using RNA-based clustering for ground truth avoids circular logic.
#
# Reference:
#   Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
#   Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
#   https://doi.org/10.1186/s12859-025-06044-0
RNA_CELL_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": [],
    },
    "T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD4"],
    },
    "Macrophages": {
        "Major": ["CD68", "HLA-DR"],
        "Minor": ["CD163"],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK", "E-Cadherin"],
        "Minor": [],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
}

# Achievable 8-cell type profile dictionary (realistic distinguishable types)
# These 8 types collapse the 10 granular types to only those distinguishable
# by the 27-antibody panel. Based on marker availability analysis 2026-01-20:
#
# Limitations of 27-antibody panel:
# - Stromal/Vascular Stromal: Both VIM+, missing DCN/LUM/PDGFRB/RGS5 to distinguish
# - Mixed Immune: CD3E+HLA-DR interface, maps to CD4+ T cells
# - Proliferating T: CD3E+PCNA, difficult to distinguish from CD8+ T cells
#
# GT Type Mapping for evaluation (10→8):
# - B cells → B cells
# - CD4+ T cells → Mixed Immune (activated CD4+ T cells show HLA-DR)
# - CD8+ T cells → CD8+ T cells + Proliferating T (both CD3E+, CD8A distinguishes)
# - Macrophages → Macrophages
# - Endothelial → Endothelial
# - Epithelial → Epithelial
# - Myofibroblasts → Myofibroblasts
# - Stromal → Stromal + Vascular Stromal (indistinguishable with panel)
ACHIEVABLE_CELL_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK", "E-Cadherin"],
        "Minor": ["Beta-catenin"],
    },
    "Myofibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": [],
    },
    "Stromal": {
        "Major": ["Vimentin"],
        "Minor": [],
    },
}

# Mapping from 10 GT types to 8 achievable types for evaluation
GT_TO_ACHIEVABLE_MAPPING = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",  # T/myeloid interface with HLA-DR
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",  # Merge with CD8+ (both CD3E+)
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Myofibroblasts",
    "Stromal": "Stromal",
    "Vascular Stromal": "Stromal",  # Merge with Stromal (both VIM+)
}

# Achievable 7-cell type profile dictionary (most conservative definition)
# Merges Myofibroblasts + Stromal → Fibroblasts because:
# - Both express Vimentin (99% of cells in both clusters)
# - alphaSMA is higher in Myofibroblasts (mean=108) but also present in Stromal (mean=20)
# - The overlap makes them indistinguishable at the spot level
#
# This is the most conservative achievable benchmark based on 2026-01-21 analysis.
ACHIEVABLE_7_CELL_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],  # E-Cadherin removed - does not co-localize (score=0.42)
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
}

# Mapping from 10 GT types to 7 achievable types for evaluation
GT_TO_ACHIEVABLE_7_MAPPING = {
    "B cells": "B cells",
    "Mixed Immune": "CD4+ T cells",  # T/myeloid interface with HLA-DR
    "CD8+ T cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells",  # Merge with CD8+ (both CD3E+)
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",  # CD31+ perivascular
    "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts",  # Merge - both VIM+
    "Stromal": "Fibroblasts",  # Merge - both VIM+
}

# Spatially-Validated 8-cell type profile dictionary (based on actual colocalization)
# These definitions are derived from bivariate Moran's I analysis across 5 Xenium regions
# (2026-01-22 analysis). Key findings:
#
# - PanCK + E-Cadherin: score = 0.42 (weak, do NOT co-localize in this data)
# - Vimentin + E-Cadherin: score = 0.85 (strong EMT signal)
# - alphaSMA + Vimentin: score = 0.83 (weaker than Vim+ECad)
#
# This captures cell STATES (EMT, exhausted T, M2 macs) not just cell TYPES.
# See docs/plans/2026-01-22-revised-gt-definitions-design.md for full analysis.
SPATIALLY_VALIDATED_CELL_PROFILE_DICT = {
    # Core immune populations (validated by colocalization)
    "B cells": {
        "Major": ["CD20", "CD45RA"],     # Score: 0.86
        "Minor": [],
    },
    "Macrophages": {
        "Major": ["CD68", "CD16"],        # Score: 0.90
        "Minor": [],
    },
    "T cells": {
        "Major": ["CD3E", "CD45RO"],      # Score: 0.88 (pan-T marker)
        "Minor": [],
    },
    # Functional immune states (validated by colocalization)
    "Exhausted T cells": {
        "Major": ["CD8A", "PD-1"],        # Score: 0.87
        "Minor": [],
    },
    "M2 Macrophages": {
        "Major": ["CD163"],               # Singleton, polarization marker
        "Minor": [],
    },
    # Structural populations
    "Endothelial": {
        "Major": ["CD31"],                # Singleton, consistent
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],               # True epithelial (NOT E-Cadherin)
        "Minor": [],
    },
    # EMT/Hybrid population
    "EMT cells": {
        "Major": ["Vimentin", "E-Cadherin"],  # Score: 0.85 (hybrid state)
        "Minor": [],
    },
}

# Mapping from 10 GT types to 8 spatially-validated types for evaluation
GT_TO_SPATIALLY_VALIDATED_MAPPING = {
    "B cells": "B cells",
    "Mixed Immune": "T cells",            # T cell interface
    "CD8+ T cells": "Exhausted T cells",  # Functional state
    "Proliferating T": "T cells",         # Pan-T
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial",    # CD31+ perivascular
    "Epithelial": "Epithelial",
    "Myofibroblasts": "EMT cells",        # EMT transitional state
    "Stromal": "EMT cells",               # Vimentin+ maps to EMT
}

# Negative marker competition rules for post-hoc redistribution
# Format: (loser_type, winner_type, discriminating_marker)
# When the discriminating marker is high, proportions are redistributed from loser to winner
NEGATIVE_MARKER_COMPETITIONS = [
    # Stromal (Vimentin+) is distinguishable from other types by negative markers:
    ('Stromal', 'Endothelial', 'CD31'),       # Endothelial: CD31+, Stromal: CD31-
    ('Stromal', 'Macrophages', 'CD68'),       # Macrophages: CD68+, Stromal: CD68-
    ('Stromal', 'CD8+ T cells', 'CD3E'),      # T cells: CD3E+, Stromal: CD3E-
    ('Stromal', 'CD4+ T cells', 'CD3E'),      # CD4+ T cells also CD3E+
    ('Stromal', 'Myofibroblasts', 'alphaSMA'),# Myofibroblasts: alphaSMA+, Stromal: alphaSMA-
    ('Stromal', 'Epithelial', 'PanCK'),       # Epithelial: PanCK+, Stromal: PanCK-
    # Additional competitions for marker disambiguation
    ('Myofibroblasts', 'Macrophages', 'CD68'),  # Macrophages: CD68+
    ('Myofibroblasts', 'CD8+ T cells', 'CD3E'), # T cells: CD3E+
]


def apply_negative_marker_redistribution(
    proportions: np.ndarray,
    protein_data: np.ndarray,
    protein_names: List[str],
    cell_type_names: List[str],
    competitions: List[Tuple[str, str, str]] = None,
    transfer_fraction: float = 0.6,
) -> np.ndarray:
    """
    Apply competition-based redistribution using negative marker logic.

    When a discriminating marker is high, proportion is transferred from the
    loser cell type to the winner cell type. This implements the negative
    marker constraint in a post-hoc manner.

    Args:
        proportions: Cell type proportions array (n_spots, n_types) or DataFrame
        protein_data: Protein expression matrix (n_spots, n_proteins)
        protein_names: List of protein names
        cell_type_names: List of cell type names (matching proportions columns)
        competitions: List of (loser, winner, marker) tuples
        transfer_fraction: Fraction of proportion to transfer when marker is high

    Returns:
        Adjusted proportions array (normalized to sum to 1)
    """
    if competitions is None:
        competitions = NEGATIVE_MARKER_COMPETITIONS

    # Convert DataFrame to numpy array if needed
    if isinstance(proportions, pd.DataFrame):
        adjusted = proportions.values.copy()
    else:
        adjusted = proportions.copy()

    # Ensure protein_data is numpy array
    if hasattr(protein_data, 'toarray'):
        protein_data = protein_data.toarray()
    elif isinstance(protein_data, pd.DataFrame):
        protein_data = protein_data.values

    protein_name_to_idx = {name: i for i, name in enumerate(protein_names)}
    cell_type_to_idx = {name: i for i, name in enumerate(cell_type_names)}

    for loser, winner, marker in competitions:
        # Skip if cell types or marker not present
        if loser not in cell_type_to_idx:
            continue
        if winner not in cell_type_to_idx:
            continue
        if marker not in protein_name_to_idx:
            continue

        loser_idx = cell_type_to_idx[loser]
        winner_idx = cell_type_to_idx[winner]
        marker_idx = protein_name_to_idx[marker]

        # Get marker expression
        marker_expr = protein_data[:, marker_idx]

        # Calculate threshold (median of nonzero values)
        nonzero = marker_expr[marker_expr > 0]
        if len(nonzero) == 0:
            continue
        threshold = np.percentile(nonzero, 50)

        # Transfer proportion where marker is high
        high_marker = marker_expr > threshold
        transfer = adjusted[high_marker, loser_idx] * transfer_fraction

        # Apply transfer
        adjusted[high_marker, loser_idx] -= transfer
        adjusted[high_marker, winner_idx] += transfer

    # Normalize to sum to 1
    row_sums = adjusted.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    adjusted = adjusted / row_sums

    return adjusted


# Granular 10-cell type profile dictionary (unsimplified RNA clustering)
# These 10 cell types match the unsimplified RNA k-means clusters.
# This provides maximum granularity to highlight CITEgeist's proteomic advantage.
#
# Cluster analysis performed 2026-01-02 (see analyze_cluster_profiles.py):
# Cluster 1: CD8+ T cells     - CD3E=378, CD8A=210
# Cluster 2: Macrophages      - CD68=430, CD163=88
# Cluster 3: Mixed Immune     - CD3E=142, CD8A=118, HLA-DR=142
# Cluster 4: Epithelial       - PanCK=39, Vimentin=311
# Cluster 5: Myofibroblasts   - alphaSMA=108, Vimentin=374
# Cluster 6: Stromal          - Mixed low markers
# Cluster 7: Endothelial      - CD31=168
# Cluster 8: B cells          - CD20=293, CD45RA=398
# Cluster 9: Proliferating T  - CD3E=679, PCNA=83
# Cluster 10: Vascular Stromal - CD31=53, Vimentin=209
GRANULAR_CELL_PROFILE_DICT: Dict[str, Dict[str, Any]] = {
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["CD45", "GranzymeB"],
    },
    "Macrophages": {
        "Major": ["CD68", "HLA-DR"],
        "Minor": ["CD163", "CD16"],
    },
    "Mixed Immune": {
        "Major": ["CD3E", "HLA-DR"],
        "Minor": ["CD8A", "CD45"],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": ["E-Cadherin", "Vimentin"],
    },
    "Myofibroblasts": {
        "Major": ["alphaSMA", "Vimentin"],
        "Minor": [],
    },
    "Stromal": {
        "Major": ["Vimentin"],
        "Minor": ["CD3E"],  # Low but detectable
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": ["Vimentin"],
    },
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45", "CD45RA"],
    },
    "Proliferating T": {
        "Major": ["CD3E"],
        "Minor": ["PCNA", "Ki-67"],
    },
    "Vascular Stromal": {
        "Major": ["Vimentin"],
        "Minor": ["CD31"],
    },
}


def run_citegeist_on_region(
    region_id: int,
    input_dir: str,
    output_dir: str,
    cell_profile_dict: Optional[Dict] = None,
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    lambda_laplacian: float = 0.1,
    laplacian_k: int = 8,
    min_counts: int = 25,
    nonzero_percentage: float = 0.01,
    mean_expression_threshold: float = 1.1,
    prefix: str = "Xenium",
    run_gex: bool = False,
    use_autodiscovery: bool = False,
    n_permutations: int = 199,
    fdr_threshold: float = 0.05,
    top_k: int = 5,
    variance_target: float = 0.90,
    use_negative_markers: bool = False,
    transfer_fraction: float = 0.6,
) -> Dict[str, Any]:
    """
    Run CITEgeist deconvolution on a single region.

    Args:
        region_id: Region identifier (0-4)
        input_dir: Directory containing h5ad_objects/
        output_dir: Directory to save outputs
        cell_profile_dict: Cell type to marker mapping (default: XENIUM_CELL_PROFILE_DICT)
        radius: Neighbor detection radius for spatial prior
        lambda_reg: Regularization strength
        alpha_elastic: Elastic net mixing parameter
        max_y_change: Maximum change in Y during optimization
        lambda_laplacian: Laplacian smoothing strength
        laplacian_k: Number of neighbors for Laplacian
        min_counts: Minimum counts per spot
        nonzero_percentage: Min percentage of spots with expression
        mean_expression_threshold: Min mean expression in nonzero spots
        prefix: Filename prefix
        run_gex: Whether to run gene expression deconvolution (Pass 1)
        use_autodiscovery: Use Module 1-2 to auto-discover cell profiles from protein data
        n_permutations: Number of permutations for colocalization analysis
        fdr_threshold: FDR threshold for profile discovery
        top_k: Top-k neighbors for mutual sparsification
        variance_target: Target variance for profile selection (dual checkpoint: spatial and protein)
        use_negative_markers: Apply negative marker redistribution post-hoc
        transfer_fraction: Fraction of proportion to transfer (0-1)

    Returns:
        Dict with run statistics and output paths
    """
    if cell_profile_dict is None and not use_autodiscovery:
        cell_profile_dict = XENIUM_CELL_PROFILE_DICT

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # File paths
    gex_path = input_dir / "h5ad_objects" / f"{prefix}_region_{region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"{prefix}_region_{region_id}_CITE.h5ad"

    logger.info(f"Loading region {region_id} data...")
    logger.info(f"  GEX: {gex_path}")
    logger.info(f"  Protein: {protein_path}")

    # Load data
    adata_gex = sc.read_h5ad(gex_path)
    adata_protein = sc.read_h5ad(protein_path)

    logger.info(f"  GEX shape: {adata_gex.shape}")
    logger.info(f"  Protein shape: {adata_protein.shape}")

    # Initialize model in simulation mode
    sample_name = f"{prefix}_region_{region_id}"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_protein,
    )

    # Auto-discover cell profiles from protein data (Module 1-2)
    discovery_stats = {}
    if use_autodiscovery:
        logger.info("=" * 60)
        logger.info("Running Auto Profile Discovery (Module 1-2)")
        logger.info("=" * 60)

        # Get raw protein data and coordinates
        X_protein = adata_protein.X
        if hasattr(X_protein, 'toarray'):
            X_protein = X_protein.toarray()
        coords = adata_protein.obsm.get('spatial', None)
        if coords is None:
            # Try to get from obs
            if 'x' in adata_protein.obs.columns and 'y' in adata_protein.obs.columns:
                coords = adata_protein.obs[['x', 'y']].values
            else:
                raise ValueError("No spatial coordinates found in protein data")
        marker_names = list(adata_protein.var_names)

        # Module 1: Identify interesting markers
        logger.info("Module 1: Identifying interesting markers...")
        interest_result = identify_interesting_markers(
            X=X_protein,
            coords=coords,
            marker_names=marker_names,
            morans_k=8,
            verbose=True,
        )
        interesting_markers = interest_result.interesting_markers
        logger.info(f"Module 1: Found {len(interesting_markers)} interesting markers")

        if len(interesting_markers) < 2:
            logger.warning("Not enough interesting markers. Using all markers.")
            interesting_markers = marker_names

        # Module 2a: Analyze marker colocalization
        logger.info("Module 2a: Analyzing marker colocalization...")
        coloc_result = analyze_marker_colocalization(
            X=X_protein,
            coords=coords,
            marker_names=marker_names,
            markers_to_analyze=interesting_markers,
            neighbor_k=8,
            n_permutations=n_permutations,
            seed=42,
            verbose=True,
        )
        logger.info(f"Module 2a: Found {len(coloc_result.pairs)} significant marker pairs")

        # Module 2b: Discover profiles
        logger.info("Module 2b: Discovering profiles...")
        discovery_result = discover_profiles(
            colocalization_result=coloc_result,
            fdr_alpha=fdr_threshold,
            top_k=top_k,
            seed=42,
            verbose=True,
        )
        logger.info(f"Module 2b: Discovered {len(discovery_result.profiles)} candidate profiles")
        for i, profile in enumerate(discovery_result.profiles):
            logger.info(f"  {i+1}. {profile}")

        # Module 2c: Select profiles by spatial variance
        logger.info("Module 2c: Selecting profiles by spatial variance...")
        selection_result = select_profiles(
            X=X_protein,
            coords=coords,
            marker_names=marker_names,
            profiles=discovery_result.profiles,
            interesting_markers=interesting_markers,
            colocalization_result=coloc_result,
            min_spatial_explained=variance_target,
            min_protein_explained=variance_target,  # Dual variance checkpoint
            verbose=True,
        )

        selected_profiles = selection_result.selected_profiles
        n_selected = selection_result.optimal_n
        total_ve = float(selection_result.variance_explained[n_selected - 1]) if n_selected > 0 else 0.0

        logger.info(f"Module 2c: Selected {len(selected_profiles)} profiles")
        logger.info(f"  Spatial variance explained: {total_ve:.1%}")
        logger.info(f"  Stopping reason: {selection_result.stopping_reason}")

        for i, profile in enumerate(selected_profiles):
            logger.info(f"  {i+1}. {profile}")

        # Convert to cell_profile_dict format
        cell_profile_dict = {}
        for i, profile in enumerate(selected_profiles):
            profile_name = f"Profile_{i+1}"
            markers_list = list(profile) if not isinstance(profile, list) else profile
            cell_profile_dict[profile_name] = {"Major": markers_list}

        logger.info(f"Created cell_profile_dict with {len(cell_profile_dict)} profiles")
        logger.info("=" * 60)

        # Store discovery stats
        discovery_stats = {
            "n_interesting_markers": len(interesting_markers),
            "n_colocalization_pairs": len(coloc_result.pairs),
            "n_candidate_profiles": len(discovery_result.profiles),
            "n_selected_profiles": len(selected_profiles),
            "variance_explained": total_ve,
            "stopping_reason": selection_result.stopping_reason,
            "profiles": {k: v["Major"] for k, v in cell_profile_dict.items()},
        }

    # Load cell profile dictionary
    model.load_cell_profile_dict(cell_profile_dict)

    # Preprocessing
    logger.info("Preprocessing gene expression...")
    model.filter_gex(
        nonzero_percentage=nonzero_percentage,
        mean_expression_threshold=mean_expression_threshold,
        min_counts=min_counts,
    )
    model.preprocess_gex(target_sum=10000)

    logger.info("Preprocessing antibody/protein...")
    model.preprocess_antibody()

    # Run cell proportion model
    logger.info("Running cell proportion optimization...")
    start_time = time.time()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        lambda_laplacian=lambda_laplacian,
        laplacian_k=laplacian_k,
        per_marker_beta=False,  # Disable per-marker beta (made results worse)
        validation_warn_only=True,  # Don't fail on high Unknown proportion
    )

    prop_time = time.time() - start_time
    logger.info(f"Cell proportion optimization completed in {prop_time:.1f}s")

    # Apply negative marker redistribution if requested
    if use_negative_markers and finetuned_props is not None:
        logger.info("Applying negative marker redistribution...")
        # Get protein data
        X_protein = adata_protein.X
        if hasattr(X_protein, 'toarray'):
            X_protein = X_protein.toarray()
        protein_names = list(adata_protein.var_names)
        # Note: cell_profile_dict.keys() already includes "Unknown" (added by model automatically)
        cell_type_names = list(cell_profile_dict.keys())

        finetuned_props = apply_negative_marker_redistribution(
            proportions=finetuned_props,
            protein_data=X_protein,
            protein_names=protein_names,
            cell_type_names=cell_type_names,
            competitions=NEGATIVE_MARKER_COMPETITIONS,
            transfer_fraction=transfer_fraction,
        )
        logger.info("Negative marker redistribution applied")

    # Run gene expression deconvolution if requested
    gex_time = 0.0
    if run_gex:
        logger.info("Running gene expression deconvolution (Pass 1)...")
        gex_start = time.time()

        try:
            pass1_results = model.run_cell_expression_pass1(
                radius=radius,
                alpha=0.5,
                lambda_reg_gex=0.001,
                global_enrichment_weight=0.5,
                local_enrichment_weight=0.5,
                checkpoint_interval=100,
                output_dir=str(output_dir / "checkpoints"),
                rerun=True,
            )
            gex_time = time.time() - gex_start
            logger.info(f"Gene expression deconvolution completed in {gex_time:.1f}s")
        except Exception as e:
            logger.error(f"Gene expression deconvolution failed: {e}")
            gex_time = -1.0  # Mark as failed

    # Get output paths - results are saved by the model during run
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)
    prop_csv = result_dir / "cell_proportions.csv"

    # Save proportions in benchmarking-compatible format
    # Note: cell_profile_dict.keys() already includes "Unknown" (added by model automatically)
    if finetuned_props is not None:
        props_df = pd.DataFrame(
            finetuned_props,
            index=model.antibody_capture_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )
        props_df.to_csv(result_dir / f"{sample_name}_deconv_predictions.csv")

    # Run statistics
    stats = {
        "region_id": region_id,
        "sample_name": sample_name,
        "n_spots": adata_gex.shape[0],
        "n_genes_filtered": model.gene_expression_adata.shape[1],
        "n_proteins": adata_protein.shape[1],
        "n_cell_types": len(cell_profile_dict) + 1,  # +1 for Unknown
        "runtime_proportions_sec": prop_time,
        "runtime_gex_sec": gex_time if run_gex else None,
        "output_dir": str(result_dir),
        "prop_csv": str(prop_csv),
        "gex_layers_dir": str(output_dir / f"{sample_name}_pass1" / "layers") if run_gex else None,
        "use_autodiscovery": use_autodiscovery,
        "discovery_stats": discovery_stats if use_autodiscovery else None,
    }

    # Save stats
    with open(result_dir / "run_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"Results saved to {result_dir}")

    return stats


def run_all_regions(
    input_dir: str,
    output_dir: str,
    n_regions: int = 5,
    **kwargs,
) -> Dict[str, Any]:
    """
    Run CITEgeist on all regions.

    Args:
        input_dir: Directory containing h5ad_objects/
        output_dir: Directory to save outputs
        n_regions: Number of regions
        **kwargs: Additional arguments for run_citegeist_on_region

    Returns:
        Dict with summary statistics
    """
    all_stats = []

    for region_id in range(n_regions):
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing region {region_id}/{n_regions-1}")
        logger.info(f"{'='*60}")

        stats = run_citegeist_on_region(
            region_id=region_id,
            input_dir=input_dir,
            output_dir=output_dir,
            **kwargs,
        )
        all_stats.append(stats)

    # Summary
    total_runtime = sum(s["runtime_proportions_sec"] for s in all_stats)
    total_spots = sum(s["n_spots"] for s in all_stats)

    summary = {
        "n_regions": n_regions,
        "total_spots": total_spots,
        "total_runtime_sec": total_runtime,
        "mean_runtime_per_region": total_runtime / n_regions,
        "regions": all_stats,
    }

    # Save summary
    output_dir = Path(output_dir)
    with open(output_dir / "benchmark_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\n{'='*60}")
    logger.info("All regions completed!")
    logger.info(f"Total runtime: {total_runtime:.1f}s ({total_runtime/60:.1f} min)")
    logger.info(f"Total spots: {total_spots}")
    logger.info(f"{'='*60}")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Run CITEgeist on Xenium pseudo-Visium data"
    )
    parser.add_argument(
        "--region-id",
        type=int,
        required=True,
        help="Region ID to process (0-4)",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="Benchmarking/xenium_pseudovisium/data",
        help="Input directory containing h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="Benchmarking/xenium_benchmarking/CITEgeist/output",
        help="Output directory for results",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=4.0,
        help="Neighbor detection radius",
    )
    parser.add_argument(
        "--lambda-reg",
        type=float,
        default=1.0,
        help="Regularization strength",
    )
    parser.add_argument(
        "--alpha-elastic",
        type=float,
        default=0.7,
        help="Elastic net mixing parameter",
    )
    parser.add_argument(
        "--max-y-change",
        type=float,
        default=0.4,
        help="Maximum change in Y during optimization",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=25,
        help="Minimum counts per spot",
    )
    parser.add_argument(
        "--use-rna-profiles",
        action="store_true",
        help="Use RNA-based cell profile dict (6 cell types matching simplified RNA clustering GT)",
    )
    parser.add_argument(
        "--use-granular-profiles",
        action="store_true",
        help="Use granular cell profile dict (10 cell types matching unsimplified RNA clustering GT)",
    )
    parser.add_argument(
        "--use-achievable-profiles",
        action="store_true",
        help="Use achievable cell profile dict (8 distinguishable types based on antibody panel)",
    )
    parser.add_argument(
        "--use-achievable-7-profiles",
        action="store_true",
        help="Use achievable-7 cell profile dict (7 types, merges Myofibroblasts+Stromal→Fibroblasts)",
    )
    parser.add_argument(
        "--run-gex",
        action="store_true",
        help="Run gene expression deconvolution (Pass 1) after cell proportions",
    )
    parser.add_argument(
        "--use-autodiscovery",
        action="store_true",
        help="Use Module 1-2 to auto-discover cell profiles from protein data (recommended)",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=199,
        help="Number of permutations for colocalization analysis",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for profile discovery",
    )
    parser.add_argument(
        "--variance-target",
        type=float,
        default=0.90,
        help="Target variance explained for profile selection (dual checkpoint: spatial and protein)",
    )
    parser.add_argument(
        "--use-negative-markers",
        action="store_true",
        help="Apply negative marker redistribution post-hoc (reduces Stromal over-prediction)",
    )
    parser.add_argument(
        "--transfer-fraction",
        type=float,
        default=0.6,
        help="Fraction of proportion to transfer during negative marker redistribution (0-1)",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Select cell profile dict based on flags
    if args.use_autodiscovery:
        cell_profile_dict = None  # Will be discovered automatically
        logger.info("Using AUTO PROFILE DISCOVERY (Module 1-2) - profiles will be discovered from protein data")
    elif args.use_achievable_7_profiles:
        cell_profile_dict = ACHIEVABLE_7_CELL_PROFILE_DICT
        logger.info("Using ACHIEVABLE-7 cell profile dict (7 types, most conservative)")
        logger.info("GT type mapping: Myofibroblasts+Stromal→Fibroblasts, Vascular Stromal→Endothelial")
    elif args.use_achievable_profiles:
        cell_profile_dict = ACHIEVABLE_CELL_PROFILE_DICT
        logger.info("Using ACHIEVABLE cell profile dict (8 distinguishable types based on antibody panel)")
        logger.info("GT type mapping: Stromal+Vascular Stromal→Stromal, CD8+ T+Proliferating T→CD8+ T, Mixed Immune→CD4+ T")
    elif args.use_granular_profiles:
        cell_profile_dict = GRANULAR_CELL_PROFILE_DICT
        logger.info("Using GRANULAR cell profile dict (10 cell types, unsimplified RNA clustering)")
    elif args.use_rna_profiles:
        cell_profile_dict = RNA_CELL_PROFILE_DICT
        logger.info("Using RNA-based cell profile dict (6 cell types, simplified)")
    else:
        cell_profile_dict = XENIUM_CELL_PROFILE_DICT
        logger.info("Using protein-based cell profile dict (10 cell types)")

    # Log negative marker status
    if args.use_negative_markers:
        logger.info(f"Negative marker redistribution ENABLED (transfer_fraction={args.transfer_fraction})")

    # Run
    stats = run_citegeist_on_region(
        region_id=args.region_id,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        cell_profile_dict=cell_profile_dict,
        radius=args.radius,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        max_y_change=args.max_y_change,
        min_counts=args.min_counts,
        run_gex=args.run_gex,
        use_autodiscovery=args.use_autodiscovery,
        n_permutations=args.n_permutations,
        fdr_threshold=args.fdr_threshold,
        variance_target=args.variance_target,
        use_negative_markers=args.use_negative_markers,
        transfer_fraction=args.transfer_fraction,
    )

    print(f"\nCompleted region {args.region_id}")
    print(f"  Spots: {stats['n_spots']}")
    print(f"  Runtime (proportions): {stats['runtime_proportions_sec']:.1f}s")
    if args.run_gex and stats.get('runtime_gex_sec'):
        print(f"  Runtime (GEX): {stats['runtime_gex_sec']:.1f}s")
        print(f"  GEX layers: {stats['gex_layers_dir']}")
    print(f"  Output: {stats['output_dir']}")


if __name__ == "__main__":
    main()
