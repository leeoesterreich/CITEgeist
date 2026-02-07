"""
Soft cell-to-profile assignment for single-cell spatial proteomics.

.. deprecated::
    This module is deprecated. Use ``cell_classification.py`` instead, which
    implements flow-cytometry-style Boolean gating -- the correct paradigm for
    single-cell classification. Soft scoring spreads probability across types,
    which is conceptually wrong when each cell IS one type.

    The gating-based approach is automatically dispatched when
    ``resolution="cell"`` in CitegeistModel.run_cell_proportion_model().

This module provides cell type assignment for single-cell data (like Xenium),
where each observation is a single cell rather than a spot containing multiple
cells. Unlike deconvolution (Module 3), which estimates mixtures within spots,
this assigns individual cells to protein-defined profiles.

Key differences from deconvolution:
- Input: Single-cell protein expression (one cell per row)
- Output: Soft assignment probabilities (not proportions)
- Each cell is primarily assigned to one profile (can have uncertainty)
- "Unknown" category for cells not matching any profile well

Methods:
1. Weighted sum: Score cells by sum of profile marker intensities
2. Correlation: Score by correlation with profile signatures
3. Spatial smoothing (optional): Regularize assignments using neighbors
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from scipy.special import softmax

logger = logging.getLogger(__name__)


@dataclass
class CellAssignmentResult:
    """Results container for cell-to-profile assignment."""

    soft_assignments: NDArray[np.floating]  # (n_cells, n_profiles) probabilities
    hard_assignments: NDArray[np.int_]  # (n_cells,) argmax profile index
    assignment_confidence: NDArray[np.floating]  # (n_cells,) max probability
    profile_names: List[str]  # Names of profiles (including "Unknown")
    cell_ids: List[str]  # Cell identifiers
    n_cells_per_profile: Dict[str, int]  # Counts per profile
    unknown_fraction: float  # Fraction assigned to "Unknown"

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert to DataFrame with cell IDs and assignments.

        Returns:
            DataFrame with columns for each profile probability and hard assignment.
        """
        df = pd.DataFrame(
            self.soft_assignments,
            columns=self.profile_names,
            index=self.cell_ids,
        )
        df["assigned_profile"] = [self.profile_names[i] for i in self.hard_assignments]
        df["assignment_confidence"] = self.assignment_confidence
        return df

    def get_cells_for_profile(
        self,
        profile_name: str,
        min_confidence: float = 0.5,
    ) -> List[str]:
        """
        Get cell IDs assigned to a specific profile.

        Args:
            profile_name: Name of the profile.
            min_confidence: Minimum assignment probability.

        Returns:
            List of cell IDs.
        """
        profile_idx = self.profile_names.index(profile_name)
        mask = (self.hard_assignments == profile_idx) & (
            self.assignment_confidence >= min_confidence
        )
        return [self.cell_ids[i] for i in np.where(mask)[0]]


def _build_spatial_weights(
    coords: NDArray[np.floating],
    k: int,
) -> Tuple[NDArray[np.int_], NDArray[np.floating]]:
    """
    Build k-NN spatial neighbor indices and weights.

    Args:
        coords: Spatial coordinates (n_cells, 2).
        k: Number of neighbors.

    Returns:
        Tuple of:
        - neighbor_indices: (n_cells, k) indices of neighbors
        - weights: (n_cells, k) distance-based weights (closer = higher)
    """
    n_cells = coords.shape[0]
    tree = cKDTree(coords)

    query_k = min(k + 1, n_cells)
    dists, indices = tree.query(coords, k=query_k)

    # Remove self (first index)
    if indices.ndim == 1:
        return np.array([]).reshape(n_cells, 0), np.array([]).reshape(n_cells, 0)

    neighbor_indices = indices[:, 1:]
    neighbor_dists = dists[:, 1:]

    # Distance-based weights (inverse distance, normalized)
    # Add small epsilon to avoid division by zero
    epsilon = 1e-6
    weights = 1.0 / (neighbor_dists + epsilon)
    weights = weights / weights.sum(axis=1, keepdims=True)

    return neighbor_indices, weights


def _compute_profile_scores_weighted_sum(
    X: NDArray[np.floating],
    profile_dict: Dict[str, List[str]],
    marker_names: List[str],
) -> NDArray[np.floating]:
    """
    Compute profile scores using weighted sum of marker intensities.

    For each profile, sum the expression of its markers for each cell.
    This is a simple and interpretable approach.

    Args:
        X: Expression matrix (n_cells, n_markers).
        profile_dict: {profile_name: [marker_list]}.
        marker_names: Names of markers in X columns.

    Returns:
        Scores matrix (n_cells, n_profiles).
    """
    n_cells = X.shape[0]
    profile_names = list(profile_dict.keys())
    n_profiles = len(profile_names)

    # Create marker name to index mapping
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    scores = np.zeros((n_cells, n_profiles))

    for profile_idx, profile_name in enumerate(profile_names):
        profile_markers = profile_dict[profile_name]

        # Find indices for profile markers
        marker_indices = []
        for marker in profile_markers:
            if marker in marker_to_idx:
                marker_indices.append(marker_to_idx[marker])
            else:
                logger.warning(f"Marker '{marker}' not found in data for profile '{profile_name}'")

        if not marker_indices:
            logger.warning(f"No valid markers for profile '{profile_name}'")
            continue

        # Sum marker intensities for this profile
        profile_expression = X[:, marker_indices]
        scores[:, profile_idx] = profile_expression.sum(axis=1)

    return scores


def _compute_profile_scores_correlation(
    X: NDArray[np.floating],
    profile_dict: Dict[str, List[str]],
    marker_names: List[str],
) -> NDArray[np.floating]:
    """
    Compute profile scores using correlation with profile signatures.

    For each profile, compute the Pearson correlation between the cell's
    marker expression and a binary profile signature (1 for profile markers,
    0 for others).

    Args:
        X: Expression matrix (n_cells, n_markers).
        profile_dict: {profile_name: [marker_list]}.
        marker_names: Names of markers in X columns.

    Returns:
        Scores matrix (n_cells, n_profiles).
    """
    n_cells, n_markers = X.shape
    profile_names = list(profile_dict.keys())
    n_profiles = len(profile_names)

    # Create marker name to index mapping
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    # Build profile signature matrix (n_markers, n_profiles)
    signatures = np.zeros((n_markers, n_profiles))
    for profile_idx, profile_name in enumerate(profile_names):
        for marker in profile_dict[profile_name]:
            if marker in marker_to_idx:
                signatures[marker_to_idx[marker], profile_idx] = 1.0

    # Compute correlation for each cell
    scores = np.zeros((n_cells, n_profiles))

    # Standardize X row-wise
    X_mean = X.mean(axis=1, keepdims=True)
    X_std = X.std(axis=1, keepdims=True)
    X_std = np.where(X_std < 1e-10, 1.0, X_std)  # Avoid division by zero
    X_norm = (X - X_mean) / X_std

    # Standardize signatures column-wise
    sig_mean = signatures.mean(axis=0, keepdims=True)
    sig_std = signatures.std(axis=0, keepdims=True)
    sig_std = np.where(sig_std < 1e-10, 1.0, sig_std)
    sig_norm = (signatures - sig_mean) / sig_std

    # Correlation is dot product of normalized vectors / n_markers
    scores = (X_norm @ sig_norm) / n_markers

    return scores


def _apply_spatial_smoothing(
    scores: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int,
    lambda_spatial: float,
) -> NDArray[np.floating]:
    """
    Apply spatial smoothing to profile scores.

    Regularize scores by averaging with neighbors, weighted by distance.

    Args:
        scores: Raw profile scores (n_cells, n_profiles).
        coords: Spatial coordinates (n_cells, 2).
        k: Number of neighbors for smoothing.
        lambda_spatial: Weight for spatial component (0-1).
            0 = no smoothing, 1 = only neighbors

    Returns:
        Smoothed scores (n_cells, n_profiles).
    """
    if lambda_spatial <= 0 or k <= 0:
        return scores

    n_cells = scores.shape[0]
    neighbor_indices, weights = _build_spatial_weights(coords, k)

    if neighbor_indices.shape[1] == 0:
        return scores

    # Compute neighbor-averaged scores
    neighbor_scores = np.zeros_like(scores)
    for i in range(n_cells):
        neighbors = neighbor_indices[i]
        w = weights[i]
        neighbor_scores[i] = np.sum(w[:, np.newaxis] * scores[neighbors], axis=0)

    # Blend original and neighbor scores
    smoothed = (1 - lambda_spatial) * scores + lambda_spatial * neighbor_scores

    return smoothed


def assign_cells_to_profiles(
    protein_data: NDArray[np.floating],
    coords: NDArray[np.floating],
    profile_dict: Dict[str, List[str]],
    marker_names: List[str],
    cell_ids: Optional[List[str]] = None,
    method: str = "weighted_sum",
    spatial_smoothing: float = 0.0,
    k_neighbors: int = 15,
    unknown_threshold: float = 0.3,
    temperature: float = 1.0,
) -> CellAssignmentResult:
    """
    Assign single cells to protein profiles using soft assignment.

    Unlike deconvolution (which estimates mixtures in spots), this assigns
    individual cells to profiles based on their protein expression pattern.
    The output is a soft assignment matrix giving the probability of each
    cell belonging to each profile.

    Args:
        protein_data: Protein expression matrix (n_cells, n_proteins).
        coords: Spatial coordinates (n_cells, 2).
        profile_dict: Dictionary mapping profile names to marker lists.
            Example: {"T cells": ["CD3E", "CD4"], "Macrophages": ["CD68", "CD163"]}
        marker_names: Names of proteins in protein_data columns.
        cell_ids: Optional list of cell identifiers.
        method: Scoring method:
            - "weighted_sum": Sum profile marker intensities (default)
            - "correlation": Correlation with profile signatures
        spatial_smoothing: Weight for spatial regularization (0-1).
            0 = no smoothing (default), 0.1-0.3 = mild smoothing
        k_neighbors: Number of neighbors for spatial smoothing.
        unknown_threshold: Max probability threshold for "Unknown" assignment.
            If best profile has probability < threshold, assign to "Unknown".
        temperature: Softmax temperature for probability conversion.
            Lower = more confident (sharper), higher = more uncertain (softer).

    Returns:
        CellAssignmentResult with soft and hard assignments.

    Example:
        >>> result = assign_cells_to_profiles(
        ...     protein_data=adata_protein.X.toarray(),
        ...     coords=adata_protein.obsm["spatial"],
        ...     profile_dict=discovered_profiles,
        ...     marker_names=list(adata_protein.var_names),
        ... )
        >>> # Get cells assigned to T cells with high confidence
        >>> t_cell_ids = result.get_cells_for_profile("T cells", min_confidence=0.7)
    """
    protein_data = np.asarray(protein_data, dtype=np.float64)
    coords = np.asarray(coords, dtype=np.float64)

    n_cells, n_markers = protein_data.shape
    profile_names = list(profile_dict.keys())
    n_profiles = len(profile_names)

    if cell_ids is None:
        cell_ids = [f"cell_{i}" for i in range(n_cells)]

    logger.info(f"Assigning {n_cells:,} cells to {n_profiles} profiles using '{method}' method")

    # Compute profile scores
    if method == "weighted_sum":
        scores = _compute_profile_scores_weighted_sum(
            protein_data, profile_dict, marker_names
        )
    elif method == "correlation":
        scores = _compute_profile_scores_correlation(
            protein_data, profile_dict, marker_names
        )
    else:
        raise ValueError(f"Unknown method: {method}. Use 'weighted_sum' or 'correlation'.")

    # Apply spatial smoothing if requested
    if spatial_smoothing > 0:
        logger.info(f"Applying spatial smoothing (lambda={spatial_smoothing}, k={k_neighbors})")
        scores = _apply_spatial_smoothing(scores, coords, k_neighbors, spatial_smoothing)

    # Convert scores to probabilities using softmax
    # Add "Unknown" profile with score 0 initially
    scores_with_unknown = np.zeros((n_cells, n_profiles + 1))
    scores_with_unknown[:, :n_profiles] = scores

    # Compute softmax probabilities
    probabilities = softmax(scores_with_unknown / temperature, axis=1)

    # The "Unknown" probability is implicitly from the zero score
    # But we want to explicitly assign to Unknown if max known profile prob < threshold
    known_probs = probabilities[:, :n_profiles]
    max_known_prob = known_probs.max(axis=1)

    # Create final soft assignments
    profile_names_with_unknown = profile_names + ["Unknown"]
    soft_assignments = np.zeros((n_cells, n_profiles + 1))
    soft_assignments[:, :n_profiles] = known_probs

    # Cells with low confidence get assigned to Unknown
    unknown_mask = max_known_prob < unknown_threshold
    soft_assignments[unknown_mask, :n_profiles] = 0.0
    soft_assignments[unknown_mask, n_profiles] = 1.0

    # Renormalize
    row_sums = soft_assignments.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums < 1e-10, 1.0, row_sums)
    soft_assignments = soft_assignments / row_sums

    # Hard assignments (argmax)
    hard_assignments = soft_assignments.argmax(axis=1)
    assignment_confidence = soft_assignments.max(axis=1)

    # Count cells per profile
    n_cells_per_profile = {}
    for profile_idx, profile_name in enumerate(profile_names_with_unknown):
        count = (hard_assignments == profile_idx).sum()
        n_cells_per_profile[profile_name] = int(count)

    unknown_count = n_cells_per_profile.get("Unknown", 0)
    unknown_fraction = unknown_count / n_cells

    # Log summary
    logger.info("Assignment summary:")
    for name, count in n_cells_per_profile.items():
        pct = 100 * count / n_cells
        logger.info(f"  {name}: {count:,} cells ({pct:.1f}%)")
    logger.info(f"  Unknown fraction: {unknown_fraction:.1%}")

    return CellAssignmentResult(
        soft_assignments=soft_assignments,
        hard_assignments=hard_assignments,
        assignment_confidence=assignment_confidence,
        profile_names=profile_names_with_unknown,
        cell_ids=cell_ids,
        n_cells_per_profile=n_cells_per_profile,
        unknown_fraction=unknown_fraction,
    )


def validate_assignments_against_clusters(
    assignment_result: CellAssignmentResult,
    cluster_labels: NDArray[np.int_],
    cluster_names: Optional[Dict[int, str]] = None,
) -> pd.DataFrame:
    """
    Validate cell assignments against external clustering (e.g., RNA-based).

    Computes contingency table showing overlap between protein-based
    profile assignments and RNA-based cluster labels.

    Args:
        assignment_result: Result from assign_cells_to_profiles.
        cluster_labels: External cluster labels for each cell.
        cluster_names: Optional mapping from cluster index to name.

    Returns:
        DataFrame contingency table (profiles × clusters).
    """
    hard_assignments = assignment_result.hard_assignments
    profile_names = assignment_result.profile_names

    # Get unique clusters
    unique_clusters = sorted(np.unique(cluster_labels))
    if cluster_names is None:
        cluster_names = {c: f"Cluster_{c}" for c in unique_clusters}

    # Build contingency table
    contingency = np.zeros((len(profile_names), len(unique_clusters)), dtype=int)
    for i, profile_idx in enumerate(hard_assignments):
        cluster_idx = cluster_labels[i]
        col = unique_clusters.index(cluster_idx)
        contingency[profile_idx, col] += 1

    df = pd.DataFrame(
        contingency,
        index=profile_names,
        columns=[cluster_names.get(c, f"Cluster_{c}") for c in unique_clusters],
    )

    return df


if __name__ == "__main__":
    # Simple test
    logging.basicConfig(level=logging.INFO)

    np.random.seed(42)
    n_cells = 1000
    n_markers = 10

    # Simulate data with 3 cell types
    X = np.random.exponential(1.0, (n_cells, n_markers))
    # Cell type 1: markers 0-2 high
    X[:300, :3] += 5.0
    # Cell type 2: markers 3-5 high
    X[300:600, 3:6] += 5.0
    # Cell type 3: markers 6-8 high
    X[600:900, 6:9] += 5.0
    # Rest is "unknown" (low everything)

    coords = np.random.rand(n_cells, 2) * 100

    profile_dict = {
        "Type1": ["M0", "M1", "M2"],
        "Type2": ["M3", "M4", "M5"],
        "Type3": ["M6", "M7", "M8"],
    }
    marker_names = [f"M{i}" for i in range(n_markers)]

    result = assign_cells_to_profiles(
        protein_data=X,
        coords=coords,
        profile_dict=profile_dict,
        marker_names=marker_names,
    )

    print("\nResult:")
    print(f"  Cells per profile: {result.n_cells_per_profile}")
    print(f"  Unknown fraction: {result.unknown_fraction:.1%}")

    # Expected: ~300 Type1, ~300 Type2, ~300 Type3, ~100 Unknown
    df = result.to_dataframe()
    print("\nFirst 10 rows:")
    print(df.head(10))
