"""
Gating-based cell classification for single-cell spatial proteomics.

Flow-cytometry-style Boolean gating for classifying individual cells into
protein-defined cell types. Each cell IS one type (or Unassigned), so this
uses hard Boolean gates rather than soft scoring or QP deconvolution.

Pipeline:
1. determine_thresholds() - BIC-based GMM model selection per marker
2. infer_negative_gates() - Auto-infer negative markers from profile exclusivity
3. classify_cells_gating() - Vectorized Boolean gate classification
4. compute_confidence() - Expression-based confidence scoring (margin + GMM posterior)
5. generate_threshold_report() - Per-marker diagnostic plots and summary
6. module2_to_gating_profiles() - Bridge from Module 2 discovery results

Replaces cell_assignment.py (soft scoring) which is conceptually wrong for
single-cell classification.
"""

import json
import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class AdaptiveThresholdInfo:
    """Per-marker adaptive thresholding metadata."""

    strategy: str  # "global" or "adaptive"
    n_clusters_used: int
    cluster_thresholds: Dict[int, float]  # cluster_id -> threshold
    cluster_methods: Dict[int, str]  # cluster_id -> "gmm_2comp" etc.
    cv_of_thresholds: float
    pct_positive_global: float
    pct_positive_adaptive: float


@dataclass
class SpatialComponentInfo:
    """Per-component spatial characterization for a marker's GMM."""

    n_components_tested: int
    component_morans_i: Dict[int, float]  # sorted_idx -> I value
    component_pvalues: Dict[int, float]  # sorted_idx -> p-value
    component_is_signal: Dict[int, bool]  # sorted_idx -> True=signal
    component_means: Dict[int, float]  # sorted_idx -> component mean
    lowest_signal_component: Optional[int] = None
    noise_components: Optional[List[int]] = None
    method: str = "spatial_gmm"  # or "bic_fallback"


@dataclass
class SpatialClusters:
    """Result of spatial expression-weighted Leiden clustering."""

    labels: NDArray[np.int_]  # (n_cells,) cluster assignment
    n_clusters: int
    cluster_sizes: Dict[int, int]
    resolution: float
    spatial_k: int


@dataclass
class MarkerThreshold:
    """Threshold info for a single marker."""

    marker_name: str
    threshold: float
    method: str  # "module1", "gmm_2comp", "gmm_3comp", "gmm_1comp", "percentile"
    n_components: int = 2
    gmm_means: Optional[Tuple] = None  # (bg,) or (bg, signal) or (bg, dim, bright)
    gmm_weights: Optional[Tuple] = None
    snr: float = 0.0
    bic_scores: Optional[Dict[int, float]] = None  # {1: bic, 2: bic, 3: bic}
    quality: str = "good"  # "good", "unimodal", "trimodal", "fallback"
    secondary_threshold: Optional[float] = None  # For 3-component (dim/bright boundary)
    adaptive_info: Optional[AdaptiveThresholdInfo] = None
    spatial_component_info: Optional[SpatialComponentInfo] = None


@dataclass
class ThresholdSet:
    """Container for all marker thresholds and pre-computed Boolean masks."""

    thresholds: Dict[str, MarkerThreshold]
    is_positive: Dict[str, NDArray[np.bool_]]  # marker -> (n_cells,) bool
    cluster_labels: Optional[NDArray[np.int_]] = None  # (n_cells,) spatial cluster IDs

    def save(self, path: str) -> None:
        """Export thresholds to JSON for reproducibility."""
        out = {}
        for name, mt in self.thresholds.items():
            entry = {
                "threshold": float(mt.threshold),
                "method": mt.method,
                "n_components": mt.n_components,
                "gmm_means": list(mt.gmm_means) if mt.gmm_means else None,
                "gmm_weights": list(mt.gmm_weights) if mt.gmm_weights else None,
                "snr": float(mt.snr),
                "bic_scores": {str(k): float(v) for k, v in mt.bic_scores.items()} if mt.bic_scores else None,
                "quality": mt.quality,
                "secondary_threshold": float(mt.secondary_threshold) if mt.secondary_threshold is not None else None,
            }
            # Serialize adaptive info if present
            if mt.adaptive_info is not None:
                ai = mt.adaptive_info
                ai_dict = {
                    "strategy": ai.strategy,
                    "n_clusters_used": ai.n_clusters_used,
                    "cluster_thresholds": {str(k): float(v) for k, v in ai.cluster_thresholds.items()},
                    "cluster_methods": {str(k): v for k, v in ai.cluster_methods.items()},
                    "cv_of_thresholds": float(ai.cv_of_thresholds),
                    "pct_positive_global": float(ai.pct_positive_global),
                    "pct_positive_adaptive": float(ai.pct_positive_adaptive),
                }
                entry["adaptive_info"] = ai_dict
            # Serialize spatial component info if present
            if mt.spatial_component_info is not None:
                sci = mt.spatial_component_info
                entry["spatial_component_info"] = {
                    "n_components_tested": sci.n_components_tested,
                    "component_morans_i": {str(k): float(v) for k, v in sci.component_morans_i.items()},
                    "component_pvalues": {str(k): float(v) for k, v in sci.component_pvalues.items()},
                    "component_is_signal": {str(k): v for k, v in sci.component_is_signal.items()},
                    "component_means": {str(k): float(v) for k, v in sci.component_means.items()},
                    "lowest_signal_component": sci.lowest_signal_component,
                    "noise_components": sci.noise_components,
                    "method": sci.method,
                }
            out[name] = entry

        # Save cluster labels alongside thresholds if present
        wrapper = {"thresholds": out}
        if self.cluster_labels is not None:
            wrapper["cluster_labels"] = self.cluster_labels.tolist()

        with open(path, "w") as f:
            json.dump(wrapper, f, indent=2)
        logger.info(f"Saved thresholds to {path}")

    @classmethod
    def load(
        cls,
        path: str,
        protein_data: NDArray[np.floating],
        marker_names: List[str],
    ) -> "ThresholdSet":
        """Load thresholds from JSON and recompute Boolean masks."""
        with open(path) as f:
            raw = json.load(f)

        # Handle both old format (flat dict) and new format (wrapper with cluster_labels)
        cluster_labels = None
        if "thresholds" in raw and isinstance(raw["thresholds"], dict):
            # New format
            if "cluster_labels" in raw and raw["cluster_labels"] is not None:
                cluster_labels = np.array(raw["cluster_labels"], dtype=np.int_)
            raw = raw["thresholds"]

        marker_to_idx = {name: i for i, name in enumerate(marker_names)}
        thresholds = {}
        is_positive = {}

        for name, info in raw.items():
            gmm_means = tuple(info["gmm_means"]) if info.get("gmm_means") else None
            gmm_weights = tuple(info["gmm_weights"]) if info.get("gmm_weights") else None
            bic_scores = {int(k): v for k, v in info["bic_scores"].items()} if info.get("bic_scores") else None

            # Reconstruct adaptive info if present
            adaptive_info = None
            if info.get("adaptive_info"):
                ai = info["adaptive_info"]
                adaptive_info = AdaptiveThresholdInfo(
                    strategy=ai.get("strategy", "global"),
                    n_clusters_used=ai["n_clusters_used"],
                    cluster_thresholds={int(k): v for k, v in ai["cluster_thresholds"].items()},
                    cluster_methods={int(k): v for k, v in ai["cluster_methods"].items()},
                    cv_of_thresholds=ai["cv_of_thresholds"],
                    pct_positive_global=ai["pct_positive_global"],
                    pct_positive_adaptive=ai["pct_positive_adaptive"],
                )

            # Reconstruct spatial component info if present
            spatial_component_info = None
            if info.get("spatial_component_info"):
                sci_raw = info["spatial_component_info"]
                spatial_component_info = SpatialComponentInfo(
                    n_components_tested=sci_raw["n_components_tested"],
                    component_morans_i={int(k): v for k, v in sci_raw["component_morans_i"].items()},
                    component_pvalues={int(k): v for k, v in sci_raw["component_pvalues"].items()},
                    component_is_signal={int(k): v for k, v in sci_raw["component_is_signal"].items()},
                    component_means={int(k): v for k, v in sci_raw["component_means"].items()},
                    lowest_signal_component=sci_raw.get("lowest_signal_component"),
                    noise_components=sci_raw.get("noise_components"),
                    method=sci_raw.get("method", "spatial_gmm"),
                )

            mt = MarkerThreshold(
                marker_name=name,
                threshold=info["threshold"],
                method=info.get("method", "gmm_2comp"),
                n_components=info.get("n_components", 2),
                gmm_means=gmm_means,
                gmm_weights=gmm_weights,
                snr=info.get("snr", 0.0),
                bic_scores=bic_scores,
                quality=info.get("quality", "good"),
                secondary_threshold=info.get("secondary_threshold"),
                adaptive_info=adaptive_info,
                spatial_component_info=spatial_component_info,
            )
            thresholds[name] = mt

            if name in marker_to_idx:
                idx = marker_to_idx[name]
                # For adaptive markers with cluster labels, reconstruct per-cluster masks
                if (adaptive_info is not None
                        and adaptive_info.strategy == "adaptive"
                        and cluster_labels is not None):
                    mask = np.zeros(protein_data.shape[0], dtype=bool)
                    expr = protein_data[:, idx]
                    for cl_id, cl_thresh in adaptive_info.cluster_thresholds.items():
                        cl_mask = cluster_labels == cl_id
                        mask[cl_mask] = expr[cl_mask] > cl_thresh
                    is_positive[name] = mask
                else:
                    is_positive[name] = protein_data[:, idx] > mt.threshold
            else:
                logger.warning(f"Marker '{name}' from threshold file not in data")

        return cls(thresholds=thresholds, is_positive=is_positive, cluster_labels=cluster_labels)


@dataclass
class GatingProfile:
    """Gate definition for a single cell type."""

    name: str
    major_markers: List[str]  # ALL must be positive (AND)
    minor_markers: List[str] = field(default_factory=list)  # Bonus, not required
    negative_markers: List[str] = field(default_factory=list)  # User-specified
    inferred_negatives: List[str] = field(default_factory=list)  # Auto-inferred
    priority: int = 0  # Higher = wins ties
    parent: Optional[str] = None  # Optional hierarchy


@dataclass
class GatingProfileSet:
    """Ordered set of gating profiles."""

    profiles: Dict[str, GatingProfile]
    gating_order: List[str]  # Sorted by priority descending

    @classmethod
    def from_profile_dict(
        cls,
        profile_dict: Dict[str, List[str]],
        priority_dict: Optional[Dict[str, int]] = None,
    ) -> "GatingProfileSet":
        """
        Build from a flat {name: [marker_list]} dictionary.

        All markers become Major. Priority is based on marker count
        (more markers = higher priority) unless overridden.
        """
        profiles = {}
        for name, markers in profile_dict.items():
            priority = len(markers)
            if priority_dict and name in priority_dict:
                priority = priority_dict[name]
            profiles[name] = GatingProfile(
                name=name,
                major_markers=list(markers),
                priority=priority,
            )

        gating_order = sorted(
            profiles.keys(), key=lambda n: profiles[n].priority, reverse=True
        )
        return cls(profiles=profiles, gating_order=gating_order)

    @classmethod
    def from_flat_dict(
        cls,
        flat_dict: Dict[str, Dict[str, Any]],
        priority_dict: Optional[Dict[str, int]] = None,
    ) -> "GatingProfileSet":
        """
        Build from structured dict with Major/Minor/Negative keys.

        Example:
            {"T cells": {"Major": ["CD3E"], "Minor": ["CD4"], "Negative": ["CD68"]}}
        """
        profiles = {}
        for name, spec in flat_dict.items():
            major = spec.get("Major", spec.get("major", []))
            minor = spec.get("Minor", spec.get("minor", []))
            negative = spec.get("Negative", spec.get("negative", []))
            parent = spec.get("parent", None)
            priority = len(major)
            if priority_dict and name in priority_dict:
                priority = priority_dict[name]
            profiles[name] = GatingProfile(
                name=name,
                major_markers=list(major),
                minor_markers=list(minor),
                negative_markers=list(negative),
                parent=parent,
                priority=priority,
            )

        gating_order = sorted(
            profiles.keys(), key=lambda n: profiles[n].priority, reverse=True
        )
        return cls(profiles=profiles, gating_order=gating_order)


@dataclass
class CellClassificationResult:
    """
    Results from gating-based cell classification.

    Compatible with downstream Pass 2 expectations: Y_assignments is a
    (n_cells, n_types) matrix usable as cell_props_values.
    """

    assignments: NDArray[np.int_]  # (n_cells,) int indices into cell_type_names
    cell_type_names: List[str]  # Including "Unassigned"
    confidence: NDArray[np.floating]  # (n_cells,) composite [0, 1]
    Y_assignments: NDArray[np.floating]  # (n_cells, n_types) one-hot, no Unassigned col
    thresholds: ThresholdSet
    doublet_flags: NDArray[np.bool_]  # (n_cells,) True = passed 2+ Major gates

    def to_dataframe(self, cell_ids: Optional[List[str]] = None) -> pd.DataFrame:
        """Convert to DataFrame with per-cell type assignments and confidence."""
        n_cells = len(self.assignments)
        if cell_ids is None:
            cell_ids = [f"cell_{i}" for i in range(n_cells)]

        # Type names excluding Unassigned for proportion columns
        type_names = [n for n in self.cell_type_names if n != "Unassigned"]

        df = pd.DataFrame(
            self.Y_assignments,
            columns=type_names,
            index=cell_ids,
        )
        df["assigned_type"] = [self.cell_type_names[i] for i in self.assignments]
        df["confidence"] = self.confidence
        df["doublet_flag"] = self.doublet_flags
        return df


# ---------------------------------------------------------------------------
# 1a. Spatial expression-weighted clustering (for adaptive thresholds)
# ---------------------------------------------------------------------------

def _build_spatial_expression_clusters(
    protein_data: NDArray[np.floating],
    spatial_coords: NDArray[np.floating],
    spatial_k: int = 30,
    min_cluster_size: int = 100,
    leiden_resolution: float = 1.0,
    seed: int = 42,
) -> SpatialClusters:
    """
    Build spatially coherent, expression-homogeneous cell neighborhoods via
    Leiden clustering on a spatial KNN graph weighted by protein expression
    similarity.

    Args:
        protein_data: Expression matrix (n_cells, n_markers).
        spatial_coords: Spatial coordinates (n_cells, 2).
        spatial_k: Number of spatial nearest neighbors.
        min_cluster_size: Merge clusters smaller than this into nearest neighbor.
        leiden_resolution: Leiden resolution parameter (higher = more clusters).
        seed: Random seed.

    Returns:
        SpatialClusters with cell cluster assignments.
    """
    from scipy.spatial import cKDTree
    from scipy import sparse
    import igraph as ig

    n_cells = protein_data.shape[0]
    protein_data = np.asarray(protein_data, dtype=np.float64)
    spatial_coords = np.asarray(spatial_coords, dtype=np.float64)

    # 1. Spatial KNN graph
    tree = cKDTree(spatial_coords)
    distances, indices = tree.query(spatial_coords, k=spatial_k + 1)
    # indices[:, 0] is self, skip it
    neighbor_indices = indices[:, 1:]

    # 2. Compute protein expression cosine similarity for edges
    # Normalize protein vectors for cosine similarity
    norms = np.linalg.norm(protein_data, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    protein_normed = protein_data / norms

    # Build sparse weighted adjacency (vectorized for performance)
    i_indices = np.repeat(np.arange(n_cells), spatial_k)
    j_indices = neighbor_indices.ravel()
    cos_sims = np.sum(protein_normed[i_indices] * protein_normed[j_indices], axis=1)
    cos_sims = np.maximum(cos_sims, 0.0)  # Clip negative similarities

    # Filter out zero-weight edges
    nonzero_mask = cos_sims > 0
    rows = i_indices[nonzero_mask]
    cols = j_indices[nonzero_mask]
    edge_weights_raw = cos_sims[nonzero_mask]

    adj = sparse.csr_matrix((edge_weights_raw, (rows, cols)), shape=(n_cells, n_cells))
    # Symmetrize: take max of (i,j) and (j,i) to avoid double-counting mutual KNN edges
    adj_t = adj.T
    adj = adj.maximum(adj_t)

    # 3. Convert to igraph and run Leiden
    sources, targets = adj.nonzero()
    edge_weights = np.array(adj[sources, targets]).flatten()

    g = ig.Graph(n=n_cells, edges=list(zip(sources.tolist(), targets.tolist())),
                 directed=False)
    g.es["weight"] = edge_weights.tolist()
    # Remove multi-edges (igraph may create them from symmetrization)
    g.simplify(combine_edges="max")

    resolution = leiden_resolution
    best_labels = None
    for attempt in range(4):
        partition = g.community_leiden(
            weights="weight",
            resolution=resolution,
            objective_function="modularity",
            n_iterations=3,
        )
        labels = np.array(partition.membership, dtype=np.int_)
        n_clusters = len(set(labels))

        if 15 <= n_clusters <= 60:
            best_labels = labels
            break
        elif n_clusters < 15:
            resolution += 0.3
        else:
            resolution = max(resolution - 0.3, 0.1)
        best_labels = labels  # Keep last attempt

    labels = best_labels
    n_clusters = len(set(labels))
    logger.info(f"Leiden clustering: {n_clusters} clusters (resolution={resolution:.2f})")

    # 4. Merge small clusters into nearest spatial neighbor cluster
    unique_labels = np.unique(labels)
    cluster_centroids = {}
    for cl in unique_labels:
        mask = labels == cl
        cluster_centroids[cl] = spatial_coords[mask].mean(axis=0)

    for cl in unique_labels:
        count = int((labels == cl).sum())
        if count < min_cluster_size:
            centroid = cluster_centroids[cl]
            best_dist = np.inf
            best_target = cl
            for other_cl, other_centroid in cluster_centroids.items():
                if other_cl == cl:
                    continue
                if int((labels == other_cl).sum()) < min_cluster_size:
                    continue  # Don't merge into another small cluster
                dist = np.linalg.norm(centroid - other_centroid)
                if dist < best_dist:
                    best_dist = dist
                    best_target = other_cl
            if best_target != cl:
                labels[labels == cl] = best_target
                logger.debug(f"Merged cluster {cl} ({count} cells) into {best_target}")

    # Relabel to consecutive integers
    unique_final = np.unique(labels)
    remap = {old: new for new, old in enumerate(unique_final)}
    labels = np.array([remap[l] for l in labels], dtype=np.int_)
    n_clusters = len(unique_final)

    cluster_sizes = {}
    for cl in range(n_clusters):
        cluster_sizes[cl] = int((labels == cl).sum())

    logger.info(
        f"After merging: {n_clusters} clusters, sizes: "
        f"min={min(cluster_sizes.values())}, max={max(cluster_sizes.values())}, "
        f"median={int(np.median(list(cluster_sizes.values())))}"
    )

    return SpatialClusters(
        labels=labels,
        n_clusters=n_clusters,
        cluster_sizes=cluster_sizes,
        resolution=resolution,
        spatial_k=spatial_k,
    )


# ---------------------------------------------------------------------------
# 1a2. Spatial GMM component characterization
# ---------------------------------------------------------------------------

def _compute_morans_i_permutation(
    x: NDArray[np.floating],
    neighbor_indices: NDArray[np.int_],
    n_perm: int = 199,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Compute Moran's I with binary spatial weights and permutation p-value.

    General-purpose function that takes pre-computed KNN neighbor indices
    to avoid rebuilding KDTree per component.

    Args:
        x: (n_cells,) continuous values to test for spatial autocorrelation.
        neighbor_indices: (n_cells, k) pre-built KNN indices (excluding self).
        n_perm: Number of permutations for p-value.
        seed: Random seed for permutations.

    Returns:
        (I_obs, p_value) where I_obs is observed Moran's I and p_value is
        from permutation test.
    """
    n_cells = len(x)
    spatial_k = neighbor_indices.shape[1]

    x_mean = x.mean()
    x_dev = x - x_mean
    ss = np.sum(x_dev ** 2)
    if ss < 1e-15:
        return 0.0, 1.0

    # Numerator: sum of w_ij * (x_i - mean) * (x_j - mean)
    neighbor_devs = x_dev[neighbor_indices]  # (n_cells, spatial_k)
    numerator = np.sum(x_dev[:, np.newaxis] * neighbor_devs)
    W = n_cells * spatial_k
    I_obs = float((n_cells / W) * (numerator / ss))

    # Permutation test
    rng = np.random.RandomState(seed)
    n_greater = 0
    for _ in range(n_perm):
        x_perm = rng.permutation(x)
        x_perm_dev = x_perm - x_perm.mean()
        ss_perm = np.sum(x_perm_dev ** 2)
        if ss_perm < 1e-15:
            continue
        neighbor_devs_perm = x_perm_dev[neighbor_indices]
        num_perm = np.sum(x_perm_dev[:, np.newaxis] * neighbor_devs_perm)
        I_perm = (n_cells / W) * (num_perm / ss_perm)
        if I_perm >= I_obs:
            n_greater += 1

    p_value = float((n_greater + 1) / (n_perm + 1))
    return I_obs, p_value


def _classify_components_by_spatial_structure(
    expr: NDArray[np.floating],
    gmm: GaussianMixture,
    neighbor_indices: NDArray[np.int_],
    n_perm: int = 199,
    morans_i_min: float = 0.05,
    pvalue_max: float = 0.05,
    seed: int = 42,
) -> Tuple[List[bool], List[float], List[float]]:
    """
    Classify each GMM component as noise or signal based on spatial structure.

    For each component (sorted by mean ascending):
    1. Compute posterior-weighted expression: z_k(i) = P(k|x_i) * x_i
    2. Compute Moran's I on z_k to test spatial autocorrelation
    3. Signal if I > morans_i_min and p < pvalue_max

    Args:
        expr: (n_cells,) full expression values.
        gmm: Fitted GaussianMixture model.
        neighbor_indices: (n_cells, k) pre-built KNN indices.
        n_perm: Number of permutations for Moran's I p-value.
        morans_i_min: Minimum Moran's I to classify as signal.
        pvalue_max: Maximum p-value for significance.
        seed: Random seed.

    Returns:
        (is_signal, morans_i_values, pvalues) per component in sorted order.
    """
    n_components = gmm.n_components
    means = gmm.means_.flatten()
    sort_idx = np.argsort(means)

    # Get posteriors for all cells
    posteriors = gmm.predict_proba(expr.reshape(-1, 1))

    is_signal = []
    morans_i_values = []
    pvalues = []

    for sorted_pos, comp_idx in enumerate(sort_idx):
        # Posterior-weighted expression
        w_k = posteriors[:, comp_idx]
        z_k = w_k * expr

        I_obs, p_val = _compute_morans_i_permutation(
            z_k, neighbor_indices, n_perm=n_perm, seed=seed + sorted_pos,
        )

        is_sig = (I_obs > morans_i_min) and (p_val < pvalue_max)
        is_signal.append(is_sig)
        morans_i_values.append(I_obs)
        pvalues.append(p_val)

    return is_signal, morans_i_values, pvalues


def _find_spatial_gmm_threshold(
    expr: NDArray[np.floating],
    spatial_coords: NDArray[np.floating],
    max_components: int = 5,
    spatial_k: int = 15,
    n_perm: int = 199,
    morans_i_min: float = 0.05,
    pvalue_max: float = 0.05,
    posterior_cutoff: float = 0.1,
    min_nonzero: int = 30,
    seed: int = 42,
) -> Optional[Tuple[MarkerThreshold, SpatialComponentInfo]]:
    """
    Find threshold using per-component spatial characterization.

    Algorithm:
    1. Build spatial KNN once
    2. Fit GMMs for k=1..max_components, select best by BIC
    3. Sort components by mean, classify each as noise/signal via Moran's I
    4. Walk from lowest-mean component upward to find noise/signal boundary
    5. If lowest component is already signal, try more components
    6. Threshold at P(lowest signal component | x) = posterior_cutoff

    Args:
        expr: (n_cells,) expression for one marker.
        spatial_coords: (n_cells, 2) spatial positions.
        max_components: Maximum GMM components to try.
        spatial_k: Number of spatial nearest neighbors.
        n_perm: Permutations for Moran's I test.
        morans_i_min: Minimum Moran's I for signal classification.
        pvalue_max: Maximum p-value for signal classification.
        posterior_cutoff: Posterior probability cutoff for threshold placement.
        min_nonzero: Minimum nonzero values required.
        seed: Random seed.

    Returns:
        (MarkerThreshold, SpatialComponentInfo) or None if no clean boundary found.
    """
    from scipy.spatial import cKDTree

    nonzero = expr[expr > 0]
    if len(nonzero) < min_nonzero:
        return None

    n_cells = len(expr)
    if n_cells < spatial_k + 1:
        return None

    # 1. Build spatial KNN once
    tree = cKDTree(spatial_coords)
    _, indices = tree.query(spatial_coords, k=spatial_k + 1)
    neighbor_indices = indices[:, 1:]  # exclude self

    # 2. Fit GMMs for k=1..max_components
    data = nonzero.reshape(-1, 1)
    bic_scores: Dict[int, float] = {}
    fitted_models: Dict[int, GaussianMixture] = {}

    for n_comp in range(1, max_components + 1):
        try:
            gmm = GaussianMixture(
                n_components=n_comp, random_state=seed, max_iter=200,
                covariance_type="full",
            )
            gmm.fit(data)
            bic_scores[n_comp] = gmm.bic(data)
            fitted_models[n_comp] = gmm
        except Exception as e:
            logger.debug(f"Spatial GMM with {n_comp} components failed: {e}")

    if not bic_scores:
        return None

    # BIC-selected k
    bic_best_k = min(bic_scores, key=bic_scores.get)

    # k=1 by BIC means unimodal — no noise/signal split possible
    if bic_best_k == 1:
        return None

    # 3. Try to find noise/signal boundary starting from BIC-selected k
    # Also try incrementing k if lowest component is already signal
    k_to_try = bic_best_k
    best_result = None

    while k_to_try <= max_components:
        if k_to_try not in fitted_models:
            break

        gmm = fitted_models[k_to_try]
        means = gmm.means_.flatten()
        sort_idx = np.argsort(means)
        sorted_means = means[sort_idx]

        # Classify components by spatial structure
        is_signal, morans_i_vals, pvals = _classify_components_by_spatial_structure(
            expr, gmm, neighbor_indices,
            n_perm=n_perm, morans_i_min=morans_i_min,
            pvalue_max=pvalue_max, seed=seed,
        )

        # Build component info dicts (sorted order)
        comp_morans = {i: morans_i_vals[i] for i in range(k_to_try)}
        comp_pvals = {i: pvals[i] for i in range(k_to_try)}
        comp_is_signal = {i: is_signal[i] for i in range(k_to_try)}
        comp_means = {i: float(sorted_means[i]) for i in range(k_to_try)}

        # Walk from lowest-mean component upward to find boundary
        lowest_signal = None
        noise_comps = []

        for sorted_pos in range(k_to_try):
            if is_signal[sorted_pos]:
                lowest_signal = sorted_pos
                break
            else:
                noise_comps.append(sorted_pos)

        # All components noise → no usable threshold
        if lowest_signal is None:
            return None

        # Lowest component is already signal → try more components
        if lowest_signal == 0:
            k_to_try += 1
            continue

        # Found a clean noise/signal boundary
        # Compute threshold at P(lowest signal comp | x) = posterior_cutoff
        lowest_signal_gmm_idx = sort_idx[lowest_signal]
        lowest_signal_mean = sorted_means[lowest_signal]

        # Search grid from 0 to lowest signal component mean
        x_min = max(0.0, sorted_means[0] - 1.0)
        x_max = lowest_signal_mean
        x_grid = np.linspace(x_min, x_max, 1000).reshape(-1, 1)
        posteriors = gmm.predict_proba(x_grid)
        signal_posterior = posteriors[:, lowest_signal_gmm_idx]

        # Find where posterior drops to posterior_cutoff (searching from high to low)
        threshold_val = None
        for i in range(len(x_grid) - 1, -1, -1):
            if signal_posterior[i] <= posterior_cutoff:
                if i < len(x_grid) - 1:
                    # Linear interpolation
                    x1, x2 = x_grid[i, 0], x_grid[i + 1, 0]
                    p1, p2 = signal_posterior[i], signal_posterior[i + 1]
                    if abs(p2 - p1) > 1e-10:
                        threshold_val = float(
                            x1 + (posterior_cutoff - p1) * (x2 - x1) / (p2 - p1)
                        )
                    else:
                        threshold_val = float((x1 + x2) / 2.0)
                else:
                    threshold_val = float(x_grid[i, 0])
                break

        if threshold_val is None:
            # Posterior never drops below cutoff — use lowest grid point
            threshold_val = float(x_grid[0, 0])

        # Clamp to data range
        threshold_val = float(np.clip(threshold_val, nonzero.min(), nonzero.max()))

        # Determine quality based on number of components needed
        if k_to_try <= 3:
            quality = "good"
        else:
            quality = "trimodal"

        # Extract GMM info for the threshold
        gmm_means_sorted = tuple(float(m) for m in sorted_means)
        gmm_weights_sorted = tuple(float(gmm.weights_[idx]) for idx in sort_idx)
        std_lo = max(
            np.sqrt(gmm.covariances_[sort_idx[0], 0, 0]
                    if gmm.covariances_.ndim == 3
                    else gmm.covariances_[sort_idx[0]]),
            1e-10,
        )
        snr = float((sorted_means[-1] - sorted_means[0]) / std_lo)

        mt = MarkerThreshold(
            marker_name="",  # Set by caller
            threshold=threshold_val,
            method=f"spatial_gmm_{k_to_try}comp",
            n_components=k_to_try,
            gmm_means=gmm_means_sorted,
            gmm_weights=gmm_weights_sorted,
            snr=snr,
            bic_scores=bic_scores,
            quality=quality,
        )

        sci = SpatialComponentInfo(
            n_components_tested=k_to_try,
            component_morans_i=comp_morans,
            component_pvalues=comp_pvals,
            component_is_signal=comp_is_signal,
            component_means=comp_means,
            lowest_signal_component=lowest_signal,
            noise_components=noise_comps,
            method="spatial_gmm",
        )

        return mt, sci

    # Exhausted max_components without finding clean boundary
    return None


# ---------------------------------------------------------------------------
# 1b. Per-cluster adaptive threshold determination
# ---------------------------------------------------------------------------

def _determine_thresholds_per_cluster(
    protein_data: NDArray[np.floating],
    marker_names: List[str],
    cluster_labels: NDArray[np.int_],
    global_thresholds: Dict[str, MarkerThreshold],
    global_is_positive: Dict[str, NDArray[np.bool_]],
    max_components: int = 3,
    min_nonzero: int = 50,
    seed: int = 42,
    cv_threshold: float = 0.15,
    gain_threshold: float = 0.10,
) -> Tuple[Dict[str, MarkerThreshold], Dict[str, NDArray[np.bool_]]]:
    """
    Fit per-cluster GMMs and decide global vs adaptive strategy per marker.

    For each marker:
    - Fit GMM independently within each spatial cluster
    - Compare cluster threshold variation (CV) and positive cell gain
    - Use adaptive only when it reveals meaningful differences

    Args:
        protein_data: Expression matrix (n_cells, n_markers).
        marker_names: Names of markers.
        cluster_labels: (n_cells,) cluster IDs.
        global_thresholds: Global MarkerThreshold dict from initial fit.
        global_is_positive: Global is_positive masks.
        max_components: Max GMM components per cluster fit.
        min_nonzero: Min nonzero values per cluster for GMM (stricter than global).
        seed: Random seed.
        cv_threshold: Coefficient of variation threshold for adaptive decision.
        gain_threshold: Minimum positive cell gain fraction to justify adaptive.

    Returns:
        Tuple of (updated thresholds dict, updated is_positive dict).
    """
    n_cells = protein_data.shape[0]
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}
    unique_clusters = np.unique(cluster_labels)

    updated_thresholds = dict(global_thresholds)
    updated_is_positive = dict(global_is_positive)

    n_adaptive = 0
    n_global = 0

    for marker_name in marker_names:
        mt_global = global_thresholds.get(marker_name)
        if mt_global is None:
            continue

        # Skip module1-derived thresholds (already have signal masks)
        if mt_global.method == "module1":
            continue

        m_idx = marker_to_idx[marker_name]
        expr = protein_data[:, m_idx]

        # Fit GMM per cluster
        cluster_thresholds: Dict[int, float] = {}
        cluster_methods: Dict[int, str] = {}

        for cl in unique_clusters:
            cl_mask = cluster_labels == cl
            cl_expr = expr[cl_mask]
            cl_nonzero = cl_expr[cl_expr > 0]

            if len(cl_nonzero) >= min_nonzero:
                result = _fit_gmm_bic(cl_nonzero, max_components, seed)
                if result is not None:
                    cluster_thresholds[int(cl)] = result.threshold
                    cluster_methods[int(cl)] = result.method
                    continue

            # Fall back to global threshold for this cluster
            cluster_thresholds[int(cl)] = mt_global.threshold
            cluster_methods[int(cl)] = "global_fallback"

        # Decision: adaptive vs global
        thresh_values = np.array(list(cluster_thresholds.values()))
        mean_thresh = thresh_values.mean()
        if mean_thresh > 0:
            cv = thresh_values.std() / mean_thresh
        else:
            cv = 0.0

        # Compute positive counts under each strategy
        n_pos_global = int(global_is_positive[marker_name].sum())
        pct_pos_global = n_pos_global / n_cells if n_cells > 0 else 0.0

        # Adaptive mask: each cell uses its cluster's threshold
        adaptive_mask = np.zeros(n_cells, dtype=bool)
        for cl in unique_clusters:
            cl_mask = cluster_labels == cl
            cl_thresh = cluster_thresholds[int(cl)]
            adaptive_mask[cl_mask] = expr[cl_mask] > cl_thresh

        n_pos_adaptive = int(adaptive_mask.sum())
        pct_pos_adaptive = n_pos_adaptive / n_cells if n_cells > 0 else 0.0

        # Gain: fractional increase in positive cells
        if n_pos_global > 0:
            gain = (n_pos_adaptive - n_pos_global) / n_pos_global
        elif n_pos_adaptive > 0:
            gain = 1.0  # Went from 0 to some positives
        else:
            gain = 0.0

        use_adaptive = (cv > cv_threshold) and (gain > gain_threshold)
        strategy = "adaptive" if use_adaptive else "global"

        adaptive_info = AdaptiveThresholdInfo(
            strategy=strategy,
            n_clusters_used=len(unique_clusters),
            cluster_thresholds=cluster_thresholds,
            cluster_methods=cluster_methods,
            cv_of_thresholds=float(cv),
            pct_positive_global=float(pct_pos_global),
            pct_positive_adaptive=float(pct_pos_adaptive),
        )

        if use_adaptive:
            n_adaptive += 1
            # Standard adaptive: median of cluster thresholds
            report_thresh = float(np.median(thresh_values))
            updated_mt = MarkerThreshold(
                marker_name=mt_global.marker_name,
                threshold=report_thresh,
                method=mt_global.method,
                n_components=mt_global.n_components,
                gmm_means=mt_global.gmm_means,
                gmm_weights=mt_global.gmm_weights,
                snr=mt_global.snr,
                bic_scores=mt_global.bic_scores,
                quality=mt_global.quality,
                secondary_threshold=mt_global.secondary_threshold,
                adaptive_info=adaptive_info,
                spatial_component_info=mt_global.spatial_component_info,
            )
            updated_thresholds[marker_name] = updated_mt
            updated_is_positive[marker_name] = adaptive_mask
            logger.info(
                f"  {marker_name}: ADAPTIVE - CV={cv:.3f}, gain={gain:+.1%}, "
                f"global_pos={pct_pos_global:.1%} -> adaptive_pos={pct_pos_adaptive:.1%}"
            )
        else:
            n_global += 1
            # Keep global threshold, just annotate with adaptive info
            updated_mt = MarkerThreshold(
                marker_name=mt_global.marker_name,
                threshold=mt_global.threshold,
                method=mt_global.method,
                n_components=mt_global.n_components,
                gmm_means=mt_global.gmm_means,
                gmm_weights=mt_global.gmm_weights,
                snr=mt_global.snr,
                bic_scores=mt_global.bic_scores,
                quality=mt_global.quality,
                secondary_threshold=mt_global.secondary_threshold,
                adaptive_info=adaptive_info,
                spatial_component_info=mt_global.spatial_component_info,
            )
            updated_thresholds[marker_name] = updated_mt
            logger.debug(
                f"  {marker_name}: GLOBAL (CV={cv:.3f}, gain={gain:+.1%})"
            )

    logger.info(
        f"Adaptive thresholding: {n_adaptive} adaptive, {n_global} global "
        f"(out of {n_adaptive + n_global} markers evaluated)"
    )

    return updated_thresholds, updated_is_positive


# ---------------------------------------------------------------------------
# 1c. Threshold determination (BIC-based GMM model selection)
# ---------------------------------------------------------------------------

def _find_posterior_crossing(gmm: GaussianMixture, lo_idx: int, hi_idx: int,
                             x_min: float, x_max: float) -> float:
    """
    Find the posterior crossing point where P(signal|x) = 0.5.

    This is the point where the posterior probability of belonging to the
    higher-mean component equals 0.5. More robust than midpoint of means
    for skewed distributions.

    Args:
        gmm: Fitted GaussianMixture model.
        lo_idx: Index of the lower-mean component.
        hi_idx: Index of the higher-mean component.
        x_min: Lower bound for search.
        x_max: Upper bound for search.

    Returns:
        Crossing point x where P(hi_component | x) = 0.5.
    """
    # Evaluate posteriors on a fine grid between the two means
    x_grid = np.linspace(x_min, x_max, 1000).reshape(-1, 1)
    posteriors = gmm.predict_proba(x_grid)

    # Find where posterior of high component crosses 0.5
    hi_posterior = posteriors[:, hi_idx]
    # Find first crossing from below 0.5 to above 0.5
    below_half = hi_posterior < 0.5
    if below_half.any() and (~below_half).any():
        # Find the transition point
        transitions = np.where(np.diff(below_half.astype(int)))[0]
        if len(transitions) > 0:
            cross_idx = transitions[0]
            # Linear interpolation for more precise crossing
            x1, x2 = x_grid[cross_idx, 0], x_grid[cross_idx + 1, 0]
            p1, p2 = hi_posterior[cross_idx], hi_posterior[cross_idx + 1]
            if abs(p2 - p1) > 1e-10:
                crossing = x1 + (0.5 - p1) * (x2 - x1) / (p2 - p1)
                return float(crossing)
            return float((x1 + x2) / 2.0)

    # Fallback: midpoint of means
    means = gmm.means_.flatten()
    return float((means[lo_idx] + means[hi_idx]) / 2.0)


def determine_thresholds(
    protein_data: NDArray[np.floating],
    marker_names: List[str],
    method: str = "auto",
    signal_masks: Optional[NDArray[np.bool_]] = None,
    signal_mask_marker_names: Optional[List[str]] = None,
    percentile_fallback: float = 75.0,
    max_components: int = 3,
    min_nonzero: int = 30,
    seed: int = 42,
    spatial_coords: Optional[NDArray[np.floating]] = None,
    adaptive: bool = True,
    spatial_k: int = 30,
    min_cluster_size: int = 100,
    leiden_resolution: float = 1.0,
    adaptive_cv_threshold: float = 0.15,
    adaptive_gain_threshold: float = 0.10,
) -> ThresholdSet:
    """
    Determine positive/negative thresholds for each marker using BIC-based
    GMM model selection, with optional adaptive spatial thresholding.

    For each marker:
    1. Fit GMM with 1, 2, and 3 components on nonzero values
    2. Select best model by BIC (lower = better)
    3. Set threshold based on selected model:
       - 1 component: unimodal, threshold at 90th percentile of fitted Gaussian
       - 2 components: posterior crossing point (P(signal|x) = 0.5)
       - 3 components: two thresholds (neg/dim boundary = primary)
    4. Percentile fallback only if GMM fails to converge

    If spatial_coords is provided and adaptive=True, additionally:
    5. Cluster cells into spatially coherent neighborhoods
    6. Fit per-cluster GMMs per marker
    7. Use adaptive thresholds where cluster variation is high and gains are meaningful

    No SNR cutoff is applied. SNR is stored as a quality metric for reporting.

    Args:
        protein_data: Expression matrix (n_cells, n_markers).
        marker_names: Names of markers in columns.
        method: "auto" (cascade), "gmm", or "percentile".
        signal_masks: Boolean (n_cells, n_markers) from Module 1.
        signal_mask_marker_names: Marker names for signal_masks columns.
        percentile_fallback: Percentile for fallback threshold.
        max_components: Maximum GMM components to test (default 3).
        min_nonzero: Minimum nonzero values required for GMM (default 30).
        seed: Random seed for GMM.
        spatial_coords: Spatial coordinates (n_cells, 2). If provided, enables
            adaptive spatial thresholding.
        adaptive: Whether to enable adaptive thresholding (requires spatial_coords).
        spatial_k: Number of spatial nearest neighbors for clustering.
        min_cluster_size: Minimum cluster size before merging.
        leiden_resolution: Leiden clustering resolution parameter.
        adaptive_cv_threshold: CV threshold for adaptive decision per marker.
        adaptive_gain_threshold: Minimum positive cell gain for adaptive decision.

    Returns:
        ThresholdSet with thresholds and pre-computed Boolean masks.
    """
    protein_data = np.asarray(protein_data, dtype=np.float64)
    n_cells, n_markers = protein_data.shape
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    thresholds: Dict[str, MarkerThreshold] = {}
    is_positive: Dict[str, NDArray[np.bool_]] = {}

    # Build signal_masks lookup if provided
    signal_mask_lookup: Dict[str, NDArray[np.bool_]] = {}
    if signal_masks is not None and signal_mask_marker_names is not None:
        for sm_idx, sm_name in enumerate(signal_mask_marker_names):
            if sm_name in marker_to_idx:
                signal_mask_lookup[sm_name] = signal_masks[:, sm_idx]

    for marker_name in marker_names:
        m_idx = marker_to_idx[marker_name]
        expr = protein_data[:, m_idx]

        # Path 1: Module 1 signal masks
        if method in ("auto", "module1") and marker_name in signal_mask_lookup:
            mask = signal_mask_lookup[marker_name]
            # Derive a threshold value from the mask boundary for reporting
            pos_vals = expr[mask]
            neg_vals = expr[~mask]
            if len(pos_vals) > 0 and len(neg_vals) > 0:
                thresh_val = (neg_vals.max() + pos_vals.min()) / 2.0
            elif len(pos_vals) > 0:
                thresh_val = pos_vals.min() * 0.9
            else:
                thresh_val = np.percentile(expr[expr > 0], percentile_fallback) if (expr > 0).any() else 0.0

            thresholds[marker_name] = MarkerThreshold(
                marker_name=marker_name,
                threshold=thresh_val,
                method="module1",
                n_components=0,
                quality="good",
            )
            is_positive[marker_name] = mask
            continue

        # Path 2a: Spatial GMM (primary when spatial_coords available)
        if method in ("auto", "gmm") and spatial_coords is not None:
            nonzero = expr[expr > 0]
            if len(nonzero) >= min_nonzero:
                spatial_result = _find_spatial_gmm_threshold(
                    expr, spatial_coords, max_components=5,
                    spatial_k=15, n_perm=199, seed=seed,
                )
                if spatial_result is not None:
                    mt, sci = spatial_result
                    mt.marker_name = marker_name
                    mt.spatial_component_info = sci
                    thresholds[marker_name] = mt
                    is_positive[marker_name] = expr > mt.threshold
                    logger.debug(
                        f"  {marker_name}: spatial_gmm threshold={mt.threshold:.3f}"
                    )
                    continue

        # Path 2b: BIC-based GMM model selection (fallback)
        if method in ("auto", "gmm"):
            nonzero = expr[expr > 0]
            if len(nonzero) >= min_nonzero:
                result = _fit_gmm_bic(nonzero, max_components, seed)
                if result is not None:
                    thresholds[marker_name] = result
                    is_positive[marker_name] = expr > result.threshold
                    logger.debug(
                        f"  {marker_name}: {result.method} threshold={result.threshold:.3f}, "
                        f"quality={result.quality}, SNR={result.snr:.2f}"
                    )
                    continue

        # Path 3: Percentile fallback
        nonzero = expr[expr > 0]
        if len(nonzero) > 0:
            thresh_val = float(np.percentile(nonzero, percentile_fallback))
        else:
            thresh_val = 0.0

        thresholds[marker_name] = MarkerThreshold(
            marker_name=marker_name,
            threshold=thresh_val,
            method="percentile",
            n_components=0,
            quality="fallback",
        )
        is_positive[marker_name] = expr > thresh_val

    # Log summary
    method_counts = {}
    for t in thresholds.values():
        method_counts[t.method] = method_counts.get(t.method, 0) + 1
    method_str = ", ".join(f"{count} {m}" for m, count in sorted(method_counts.items()))
    logger.info(f"Thresholds determined for {len(thresholds)} markers: {method_str}")

    quality_counts = {}
    for t in thresholds.values():
        quality_counts[t.quality] = quality_counts.get(t.quality, 0) + 1
    quality_str = ", ".join(f"{count} {q}" for q, count in sorted(quality_counts.items()))
    logger.info(f"Quality distribution: {quality_str}")

    # Adaptive spatial thresholding (optional)
    cluster_labels = None
    if spatial_coords is not None and adaptive:
        logger.info("Running adaptive spatial thresholding...")
        clusters = _build_spatial_expression_clusters(
            protein_data=protein_data,
            spatial_coords=spatial_coords,
            spatial_k=spatial_k,
            min_cluster_size=min_cluster_size,
            leiden_resolution=leiden_resolution,
            seed=seed,
        )
        cluster_labels = clusters.labels

        thresholds, is_positive = _determine_thresholds_per_cluster(
            protein_data=protein_data,
            marker_names=marker_names,
            cluster_labels=cluster_labels,
            global_thresholds=thresholds,
            global_is_positive=is_positive,
            max_components=max_components,
            min_nonzero=50,  # Stricter for per-cluster stability
            seed=seed,
            cv_threshold=adaptive_cv_threshold,
            gain_threshold=adaptive_gain_threshold,
        )

    return ThresholdSet(
        thresholds=thresholds,
        is_positive=is_positive,
        cluster_labels=cluster_labels,
    )


def _fit_gmm_bic(
    nonzero: NDArray[np.floating],
    max_components: int,
    seed: int,
) -> Optional[MarkerThreshold]:
    """
    Fit GMM with 1..max_components and select best by BIC.

    Returns:
        MarkerThreshold with BIC-selected threshold, or None if all fits fail.
    """
    data = nonzero.reshape(-1, 1)
    bic_scores: Dict[int, float] = {}
    fitted_models: Dict[int, GaussianMixture] = {}

    for n_comp in range(1, max_components + 1):
        try:
            gmm = GaussianMixture(
                n_components=n_comp, random_state=seed, max_iter=200,
                covariance_type="full",
            )
            gmm.fit(data)
            bic_scores[n_comp] = gmm.bic(data)
            fitted_models[n_comp] = gmm
        except Exception as e:
            logger.debug(f"GMM with {n_comp} components failed: {e}")

    if not bic_scores:
        return None

    # Select best model by BIC (lower = better)
    best_k = min(bic_scores, key=bic_scores.get)
    gmm = fitted_models[best_k]
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten() if gmm.covariances_.ndim == 1
                   else np.array([gmm.covariances_[i, 0, 0] for i in range(best_k)]))
    weights = gmm.weights_.flatten()

    # Sort components by mean
    sort_idx = np.argsort(means)
    means = means[sort_idx]
    stds = stds[sort_idx]
    weights = weights[sort_idx]

    if best_k == 1:
        # Unimodal: marker is not discriminative in this tissue
        # Threshold at 90th percentile of the fitted Gaussian
        mu, sigma = means[0], max(stds[0], 1e-10)
        from scipy.stats import norm
        thresh_val = float(norm.ppf(0.90, loc=mu, scale=sigma))
        # Clamp to data range
        thresh_val = float(np.clip(thresh_val, nonzero.min(), nonzero.max()))

        return MarkerThreshold(
            marker_name="",  # Will be set by caller
            threshold=thresh_val,
            method="gmm_1comp",
            n_components=1,
            gmm_means=(float(means[0]),),
            gmm_weights=(float(weights[0]),),
            snr=0.0,
            bic_scores=bic_scores,
            quality="unimodal",
        )

    elif best_k == 2:
        # Classic pos/neg. Threshold at posterior crossing point.
        lo_idx_sorted = 0
        hi_idx_sorted = 1

        # Map sorted indices back to GMM component indices
        lo_gmm_idx = sort_idx[lo_idx_sorted]
        hi_gmm_idx = sort_idx[hi_idx_sorted]

        mu_lo, mu_hi = means[0], means[1]
        std_lo = max(stds[0], 1e-10)
        snr = (mu_hi - mu_lo) / std_lo

        # Find posterior crossing point
        thresh_val = _find_posterior_crossing(
            gmm, lo_gmm_idx, hi_gmm_idx,
            x_min=mu_lo, x_max=mu_hi,
        )

        return MarkerThreshold(
            marker_name="",
            threshold=thresh_val,
            method="gmm_2comp",
            n_components=2,
            gmm_means=(float(mu_lo), float(mu_hi)),
            gmm_weights=(float(weights[0]), float(weights[1])),
            snr=float(snr),
            bic_scores=bic_scores,
            quality="good",
        )

    else:  # best_k == 3
        # neg/dim/bright. Two thresholds.
        # Primary threshold: neg vs dim+bright (lower boundary)
        # Secondary threshold: dim vs bright (upper boundary)
        lo_gmm_idx = sort_idx[0]
        mid_gmm_idx = sort_idx[1]
        hi_gmm_idx = sort_idx[2]

        mu_lo, mu_mid, mu_hi = means[0], means[1], means[2]
        std_lo = max(stds[0], 1e-10)
        snr = (mu_hi - mu_lo) / std_lo

        # Primary threshold: neg vs dim crossing
        primary_thresh = _find_posterior_crossing(
            gmm, lo_gmm_idx, mid_gmm_idx,
            x_min=mu_lo, x_max=mu_mid,
        )
        # Secondary threshold: dim vs bright crossing
        secondary_thresh = _find_posterior_crossing(
            gmm, mid_gmm_idx, hi_gmm_idx,
            x_min=mu_mid, x_max=mu_hi,
        )

        return MarkerThreshold(
            marker_name="",
            threshold=primary_thresh,
            method="gmm_3comp",
            n_components=3,
            gmm_means=(float(mu_lo), float(mu_mid), float(mu_hi)),
            gmm_weights=(float(weights[0]), float(weights[1]), float(weights[2])),
            snr=float(snr),
            bic_scores=bic_scores,
            quality="trimodal",
            secondary_threshold=secondary_thresh,
        )


# ---------------------------------------------------------------------------
# 2. Automatic negative gate inference
# ---------------------------------------------------------------------------

def infer_negative_gates(profiles: Dict[str, GatingProfile]) -> None:
    """
    Auto-infer negative markers from profile exclusivity. Mutates profiles in-place.

    Rules:
    - If marker X is Major for type A exclusively (not Major for any other type),
      then X is inferred-negative for all other types.
    - Skip markers that are Major for multiple types (shared markers).
    - Skip markers that are Minor for the target type.
    - User-specified negative_markers always take precedence.

    Args:
        profiles: Dict of GatingProfile to mutate.
    """
    # Count how many profiles claim each marker as Major
    major_counts: Dict[str, List[str]] = {}  # marker -> [profile_names]
    for name, profile in profiles.items():
        for marker in profile.major_markers:
            major_counts.setdefault(marker, []).append(name)

    # Identify exclusive markers (Major for exactly one type)
    exclusive_markers: Dict[str, str] = {}  # marker -> owning profile name
    for marker, owners in major_counts.items():
        if len(owners) == 1:
            exclusive_markers[marker] = owners[0]

    # Infer negatives
    for name, profile in profiles.items():
        inferred = []
        for marker, owner in exclusive_markers.items():
            if owner == name:
                continue  # This is our own marker
            if marker in profile.major_markers:
                continue  # Already Major for us (shouldn't happen, but guard)
            if marker in profile.minor_markers:
                continue  # Minor = bonus if present, don't negate
            if marker in profile.negative_markers:
                continue  # User already specified
            inferred.append(marker)
        profile.inferred_negatives = sorted(inferred)

    # Log summary
    total_inferred = sum(len(p.inferred_negatives) for p in profiles.values())
    logger.info(
        f"Inferred {total_inferred} negative gates across {len(profiles)} profiles "
        f"({len(exclusive_markers)} exclusive markers detected)"
    )


# ---------------------------------------------------------------------------
# 3. Vectorized Boolean gating classification
# ---------------------------------------------------------------------------

def classify_cells_gating(
    protein_data: NDArray[np.floating],
    marker_names: List[str],
    profile_set: GatingProfileSet,
    thresholds: ThresholdSet,
    use_negative_gates: bool = False,
) -> CellClassificationResult:
    """
    Classify cells using Boolean gating (flow-cytometry style).

    For each profile:
    - ALL Major markers must be positive (AND gate)
    - If use_negative_gates: ALL Negative + inferred-negative markers must be
      negative (AND NOT gate)
    - Minor markers add bonus score but are not required

    Cells not passing any Major gate are Unassigned.

    Args:
        protein_data: Expression matrix (n_cells, n_markers).
        marker_names: Names of markers in columns.
        profile_set: Ordered set of gating profiles.
        thresholds: ThresholdSet with is_positive masks.
        use_negative_gates: Whether to apply negative/inferred-negative gates.
            Default False. Set True to enable exclusion gating.

    Returns:
        CellClassificationResult with assignments, confidence, and doublet flags.
    """
    protein_data = np.asarray(protein_data, dtype=np.float64)
    n_cells = protein_data.shape[0]
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    profile_names = profile_set.gating_order
    n_profiles = len(profile_names)
    cell_scores = np.zeros((n_cells, n_profiles), dtype=np.float64)

    for p_idx, pname in enumerate(profile_names):
        profile = profile_set.profiles[pname]

        # Major gate: AND of all Major markers (no singleton boost)
        major_pass = np.ones(n_cells, dtype=bool)
        for marker in profile.major_markers:
            if marker not in marker_to_idx:
                logger.warning(f"Major marker '{marker}' for '{pname}' not found in data")
                major_pass[:] = False
                break

            if marker in thresholds.is_positive:
                major_pass &= thresholds.is_positive[marker]
            else:
                major_pass[:] = False

        # Negative gate: AND NOT of all negative + inferred-negative markers
        if use_negative_gates:
            all_negatives = profile.negative_markers + profile.inferred_negatives
            for marker in all_negatives:
                if marker in thresholds.is_positive:
                    major_pass &= ~thresholds.is_positive[marker]

        # Score: 1.0 for passing Major+Negative, plus Minor bonus
        score = major_pass.astype(np.float64)

        n_minor = max(len(profile.minor_markers), 1)
        for marker in profile.minor_markers:
            if marker in thresholds.is_positive:
                score += thresholds.is_positive[marker].astype(np.float64) * (0.1 / n_minor)

        # Priority tiebreaker (small epsilon)
        score += score * (profile.priority * 0.001)

        cell_scores[:, p_idx] = score

    # Assign best type per cell
    best_type_idx = np.argmax(cell_scores, axis=1)
    max_score = cell_scores.max(axis=1)

    # Build cell_type_names: profile names + "Unassigned"
    cell_type_names = list(profile_names) + ["Unassigned"]
    unassigned_idx = len(profile_names)

    # Cells with max_score == 0 -> Unassigned
    assignments = best_type_idx.copy()
    assignments[max_score == 0] = unassigned_idx

    # Doublet detection: cells passing Major gate for 2+ profiles
    n_passed = (cell_scores > 0).sum(axis=1)
    doublet_flags = n_passed >= 2

    # Build one-hot Y_assignments (n_cells, n_profiles) -- no Unassigned column
    Y_assignments = np.zeros((n_cells, n_profiles), dtype=np.float64)
    for i in range(n_profiles):
        Y_assignments[assignments == i, i] = 1.0
    # Unassigned cells have all-zero rows

    # Log summary
    for i, pname in enumerate(cell_type_names):
        count = int((assignments == i).sum())
        pct = 100.0 * count / n_cells
        logger.info(f"  {pname}: {count:,} cells ({pct:.1f}%)")

    n_doublet = int(doublet_flags.sum())
    logger.info(f"  Doublet flags: {n_doublet:,} cells ({100.0 * n_doublet / n_cells:.1f}%)")

    # Placeholder confidence (zeros); use compute_confidence() for real values
    confidence = np.zeros(n_cells, dtype=np.float64)
    # Set non-zero placeholder for assigned cells
    confidence[assignments != unassigned_idx] = 0.5

    return CellClassificationResult(
        assignments=assignments,
        cell_type_names=cell_type_names,
        confidence=confidence,
        Y_assignments=Y_assignments,
        thresholds=thresholds,
        doublet_flags=doublet_flags,
    )


# ---------------------------------------------------------------------------
# 4. Confidence scoring (expression-based, no spatial component)
# ---------------------------------------------------------------------------

def compute_confidence(
    protein_data: NDArray[np.floating],
    marker_names: List[str],
    result: CellClassificationResult,
    profile_set: GatingProfileSet,
    w_margin: float = 0.5,
    w_negative: float = 0.2,
    w_posterior: float = 0.3,
) -> NDArray[np.floating]:
    """
    Compute expression-based confidence score for classified cells.

    Components:
    - Threshold margin (w_margin): How far above threshold for Major markers.
    - Negative margin (w_negative): How far below threshold for negative markers.
    - GMM posterior (w_posterior): P(signal component | observed expression) for
      each Major marker. Directly quantifies how confident the GMM is about the
      positive call.

    Unassigned cells get confidence = 0.0.

    Args:
        protein_data: Expression matrix (n_cells, n_markers).
        marker_names: Names of markers in columns.
        result: CellClassificationResult from classify_cells_gating().
        profile_set: The GatingProfileSet used for classification.
        w_margin: Weight for threshold margin component.
        w_negative: Weight for negative margin component.
        w_posterior: Weight for GMM posterior component.

    Returns:
        confidence: (n_cells,) array of composite confidence [0, 1].
        Also updates result.confidence in place.
    """
    protein_data = np.asarray(protein_data, dtype=np.float64)
    n_cells = protein_data.shape[0]
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    unassigned_idx = result.cell_type_names.index("Unassigned")
    profile_names = [n for n in result.cell_type_names if n != "Unassigned"]
    thresholds = result.thresholds

    # Normalize weights to sum to 1
    total_w = w_margin + w_negative + w_posterior
    if total_w > 0:
        w_margin /= total_w
        w_negative /= total_w
        w_posterior /= total_w

    # Pre-compute max expression per marker for margin normalization
    max_per_marker = {}
    for marker in marker_names:
        m_idx = marker_to_idx[marker]
        max_per_marker[marker] = protein_data[:, m_idx].max()

    # Component 1: Threshold margin for Major markers
    margin_scores = np.zeros(n_cells, dtype=np.float64)
    for i in range(n_cells):
        a = result.assignments[i]
        if a == unassigned_idx:
            continue
        pname = profile_names[a]
        profile = profile_set.profiles[pname]
        margins = []
        for marker in profile.major_markers:
            if marker not in marker_to_idx or marker not in thresholds.thresholds:
                continue
            m_idx = marker_to_idx[marker]
            expr_val = protein_data[i, m_idx]
            thresh = thresholds.thresholds[marker].threshold
            denom = max(max_per_marker[marker] - thresh, 1e-10)
            m_val = (expr_val - thresh) / denom
            margins.append(np.clip(m_val, 0.0, 1.0))
        if margins:
            margin_scores[i] = float(np.mean(margins))

    # Component 2: Negative margin
    neg_scores = np.zeros(n_cells, dtype=np.float64)
    for i in range(n_cells):
        a = result.assignments[i]
        if a == unassigned_idx:
            continue
        pname = profile_names[a]
        profile = profile_set.profiles[pname]
        all_neg = profile.negative_markers + profile.inferred_negatives
        if not all_neg:
            neg_scores[i] = 1.0  # No negatives to violate
            continue
        margins = []
        for marker in all_neg:
            if marker not in marker_to_idx or marker not in thresholds.thresholds:
                continue
            m_idx = marker_to_idx[marker]
            expr_val = protein_data[i, m_idx]
            thresh = thresholds.thresholds[marker].threshold
            denom = max(thresh, 1e-10)
            m_val = (thresh - expr_val) / denom
            margins.append(np.clip(m_val, 0.0, 1.0))
        if margins:
            neg_scores[i] = float(np.mean(margins))
        else:
            neg_scores[i] = 1.0

    # Component 3: GMM posterior probability
    posterior_scores = np.zeros(n_cells, dtype=np.float64)
    # Pre-fit GMMs for markers that have them (reuse threshold info)
    gmm_cache: Dict[str, GaussianMixture] = {}
    for marker_name, mt in thresholds.thresholds.items():
        if mt.method in ("gmm_2comp", "gmm_3comp") and mt.gmm_means is not None:
            n_comp = len(mt.gmm_means)
            if marker_name in marker_to_idx:
                m_idx = marker_to_idx[marker_name]
                nonzero = protein_data[:, m_idx]
                nonzero = nonzero[nonzero > 0]
                if len(nonzero) >= 30:
                    try:
                        gmm = GaussianMixture(n_components=n_comp, random_state=42, max_iter=200)
                        gmm.fit(nonzero.reshape(-1, 1))
                        gmm_cache[marker_name] = gmm
                    except Exception:
                        pass

    for i in range(n_cells):
        a = result.assignments[i]
        if a == unassigned_idx:
            continue
        pname = profile_names[a]
        profile = profile_set.profiles[pname]
        posteriors = []
        for marker in profile.major_markers:
            if marker not in marker_to_idx:
                continue
            m_idx = marker_to_idx[marker]
            expr_val = protein_data[i, m_idx]
            if marker in gmm_cache and expr_val > 0:
                gmm = gmm_cache[marker]
                probs = gmm.predict_proba(np.array([[expr_val]]))[0]
                # Probability of highest-mean component
                hi_comp = np.argmax(gmm.means_.flatten())
                posteriors.append(float(probs[hi_comp]))
            else:
                # No GMM available: use simple sigmoid-like mapping from margin
                mt = thresholds.thresholds.get(marker)
                if mt:
                    margin = (expr_val - mt.threshold) / max(mt.threshold, 1e-10)
                    posteriors.append(float(np.clip(0.5 + 0.5 * np.tanh(margin), 0.0, 1.0)))
        if posteriors:
            posterior_scores[i] = float(np.mean(posteriors))

    # Composite
    confidence = (
        w_margin * margin_scores
        + w_negative * neg_scores
        + w_posterior * posterior_scores
    )
    # Unassigned stays at 0
    confidence[result.assignments == unassigned_idx] = 0.0

    result.confidence = confidence

    assigned_mask = result.assignments != unassigned_idx
    if assigned_mask.any():
        logger.info(
            f"Confidence computed: mean={confidence[assigned_mask].mean():.3f} "
            f"(assigned cells only)"
        )
    return confidence


# ---------------------------------------------------------------------------
# 5. Threshold report generation
# ---------------------------------------------------------------------------

def generate_threshold_report(
    thresholds: ThresholdSet,
    protein_data: NDArray[np.floating],
    marker_names: List[str],
    output_dir: str,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generate per-marker diagnostic plots and summary table.

    For each marker, produces:
    - Histogram of nonzero expression with GMM components overlaid
    - Threshold line(s) marked
    - BIC scores for 1/2/3 components
    - Summary table (CSV) of all markers with thresholds and quality metrics

    Args:
        thresholds: ThresholdSet from determine_thresholds().
        protein_data: Expression matrix (n_cells, n_markers).
        marker_names: Names of markers in columns.
        output_dir: Directory to save plots and CSV.
        seed: Random seed for refitting GMMs for visualization.

    Returns:
        Summary DataFrame with per-marker threshold info.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(output_dir, exist_ok=True)
    figures_dir = os.path.join(output_dir, "marker_plots")
    os.makedirs(figures_dir, exist_ok=True)

    protein_data = np.asarray(protein_data, dtype=np.float64)
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    summary_rows = []

    for marker_name, mt in thresholds.thresholds.items():
        if marker_name not in marker_to_idx:
            continue

        m_idx = marker_to_idx[marker_name]
        expr = protein_data[:, m_idx]
        nonzero = expr[expr > 0]
        n_positive = int((expr > mt.threshold).sum())
        n_total = len(expr)
        pct_positive = 100.0 * n_positive / n_total if n_total > 0 else 0.0

        row = {
            "marker": marker_name,
            "threshold": mt.threshold,
            "method": mt.method,
            "n_components": mt.n_components,
            "quality": mt.quality,
            "snr": mt.snr,
            "n_positive": n_positive,
            "n_total": n_total,
            "pct_positive": pct_positive,
            "n_nonzero": len(nonzero),
            "secondary_threshold": mt.secondary_threshold,
        }
        if mt.bic_scores:
            for k, v in mt.bic_scores.items():
                row[f"bic_{k}comp"] = v
        if mt.gmm_means:
            row["gmm_means"] = str(mt.gmm_means)
        if mt.gmm_weights:
            row["gmm_weights"] = str(mt.gmm_weights)
        # Adaptive thresholding info
        if mt.adaptive_info is not None:
            ai = mt.adaptive_info
            row["adaptive_strategy"] = ai.strategy
            row["adaptive_cv"] = ai.cv_of_thresholds
            row["pct_positive_global"] = ai.pct_positive_global
            row["pct_positive_adaptive"] = ai.pct_positive_adaptive
        # Spatial GMM component info
        if mt.spatial_component_info is not None:
            sci = mt.spatial_component_info
            row["spatial_method"] = sci.method
            row["spatial_n_components"] = sci.n_components_tested
            row["spatial_lowest_signal"] = sci.lowest_signal_component
            if sci.lowest_signal_component is not None:
                row["spatial_lowest_signal_mean"] = sci.component_means.get(
                    sci.lowest_signal_component
                )
            row["spatial_component_morans_i"] = str(sci.component_morans_i)

        summary_rows.append(row)

        # Generate diagnostic plot
        if len(nonzero) < 10:
            continue

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left panel: Histogram with GMM overlay
        ax = axes[0]
        ax.hist(nonzero, bins=80, density=True, alpha=0.6, color="steelblue",
                edgecolor="none", label="Nonzero values")

        # Overlay GMM components if available
        if mt.gmm_means is not None and len(nonzero) >= 30:
            try:
                n_comp = len(mt.gmm_means)
                gmm = GaussianMixture(n_components=n_comp, random_state=seed, max_iter=200)
                gmm.fit(nonzero.reshape(-1, 1))

                x_plot = np.linspace(nonzero.min(), nonzero.max(), 500).reshape(-1, 1)
                log_dens = gmm.score_samples(x_plot)
                total_density = np.exp(log_dens)
                ax.plot(x_plot.flatten(), total_density, 'k-', lw=2, label="GMM total")

                # Individual components
                colors = ["#e74c3c", "#2ecc71", "#3498db"]
                from scipy.stats import norm as scipy_norm
                means_sorted = np.sort(gmm.means_.flatten())
                sort_idx = np.argsort(gmm.means_.flatten())
                for j, comp_idx in enumerate(sort_idx):
                    mu = gmm.means_[comp_idx, 0]
                    sigma = np.sqrt(gmm.covariances_[comp_idx, 0, 0]) if gmm.covariances_.ndim == 3 else np.sqrt(gmm.covariances_[comp_idx])
                    w = gmm.weights_[comp_idx]
                    comp_y = w * scipy_norm.pdf(x_plot.flatten(), mu, sigma)
                    label_str = f"Comp {j+1}: mu={mu:.1f}, w={w:.2f}"
                    ax.plot(x_plot.flatten(), comp_y, '--', color=colors[j % len(colors)],
                            lw=1.5, label=label_str)
            except Exception:
                pass

        # Threshold line(s)
        ax.axvline(mt.threshold, color="red", lw=2, linestyle="-",
                   label=f"Threshold: {mt.threshold:.2f}")
        if mt.secondary_threshold is not None:
            ax.axvline(mt.secondary_threshold, color="orange", lw=1.5, linestyle="--",
                       label=f"Secondary: {mt.secondary_threshold:.2f}")

        ax.set_xlabel("Expression (nonzero)")
        ax.set_ylabel("Density")
        ax.set_title(f"{marker_name}\n{mt.method} ({mt.quality}), SNR={mt.snr:.2f}")
        ax.legend(fontsize=8, loc="upper right")

        # Right panel: BIC scores
        ax2 = axes[1]
        if mt.bic_scores:
            ks = sorted(mt.bic_scores.keys())
            bics = [mt.bic_scores[k] for k in ks]
            colors_bar = ["#95a5a6"] * len(ks)
            best_idx = ks.index(mt.n_components) if mt.n_components in ks else None
            if best_idx is not None:
                colors_bar[best_idx] = "#2ecc71"
            ax2.bar([str(k) for k in ks], bics, color=colors_bar, edgecolor="black")
            ax2.set_xlabel("Number of Components")
            ax2.set_ylabel("BIC (lower = better)")
            ax2.set_title(f"BIC Model Selection\nBest: {mt.n_components} components")
        else:
            ax2.text(0.5, 0.5, "No BIC data\n(Module 1 or percentile)",
                     ha="center", va="center", fontsize=12, transform=ax2.transAxes)
            ax2.set_title("BIC Model Selection")

        # Add summary text
        summary_text = (
            f"Positive: {n_positive:,}/{n_total:,} ({pct_positive:.1f}%) | "
            f"Nonzero: {len(nonzero):,}"
        )
        if mt.adaptive_info is not None:
            ai = mt.adaptive_info
            summary_text += (
                f" | Strategy: {ai.strategy.upper()}"
                f" (CV={ai.cv_of_thresholds:.3f},"
                f" global={ai.pct_positive_global:.1%},"
                f" adaptive={ai.pct_positive_adaptive:.1%})"
            )
        fig.text(0.5, 0.01, summary_text, ha="center", fontsize=10, style="italic")

        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        safe_name = marker_name.replace("/", "_").replace(" ", "_")
        plt.savefig(os.path.join(figures_dir, f"{safe_name}.png"), dpi=150, bbox_inches="tight")
        plt.close()

        # Extra plot: per-cluster threshold box plot for adaptive markers
        if (mt.adaptive_info is not None
                and mt.adaptive_info.strategy == "adaptive"
                and mt.adaptive_info.cluster_thresholds):
            fig_ad, ax_ad = plt.subplots(1, 1, figsize=(8, 4))
            cl_thresh_vals = list(mt.adaptive_info.cluster_thresholds.values())
            ax_ad.boxplot([cl_thresh_vals], vert=True, widths=0.5)
            ax_ad.axhline(mt.threshold, color="blue", ls="--", lw=1.5,
                          label=f"Median adaptive: {mt.threshold:.1f}")
            ax_ad.set_ylabel("Threshold value")
            ax_ad.set_title(
                f"{marker_name} — Per-cluster thresholds\n"
                f"CV={mt.adaptive_info.cv_of_thresholds:.3f}, "
                f"n_clusters={mt.adaptive_info.n_clusters_used}"
            )
            ax_ad.legend(fontsize=9)
            plt.tight_layout()
            plt.savefig(
                os.path.join(figures_dir, f"{safe_name}_adaptive.png"),
                dpi=150, bbox_inches="tight",
            )
            plt.close()

    # Save summary CSV
    summary_df = pd.DataFrame(summary_rows)
    csv_path = os.path.join(output_dir, "threshold_summary.csv")
    summary_df.to_csv(csv_path, index=False)
    logger.info(f"Threshold report saved: {csv_path} ({len(summary_rows)} markers)")
    logger.info(f"Diagnostic plots saved: {figures_dir}/")

    return summary_df


# ---------------------------------------------------------------------------
# 6. Bridge from Module 2
# ---------------------------------------------------------------------------

def module2_to_gating_profiles(
    discovery_result: Any,
    user_additions: Optional[Dict[str, Dict[str, List[str]]]] = None,
    user_renames: Optional[Dict[str, str]] = None,
    reference_signatures: Optional[Dict[str, List[str]]] = None,
    priority_dict: Optional[Dict[str, int]] = None,
) -> GatingProfileSet:
    """
    Convert Module 2 ProfileDiscoveryResult into GatingProfileSet.

    Shared markers (claimed by multiple profiles) are demoted to Minor.
    Exclusive markers become Major. User can override with additions/renames.

    Args:
        discovery_result: ProfileDiscoveryResult from Module 2 (or HierarchicalProfileResult).
        user_additions: Extra profiles to merge.
            Example: {"Epithelial": {"Major": ["PanCK"], "Negative": ["CD45"]}}
        user_renames: Rename discovered profiles.
            Example: {"CD3E+CD8A": "CD8+ T cells"}
        reference_signatures: Known cell type -> marker lists for auto-annotation
            via Jaccard matching.
        priority_dict: Override priorities. Default: more Major = higher priority.

    Returns:
        GatingProfileSet ready for classify_cells_gating().
    """
    # Get flat dict from discovery result
    if hasattr(discovery_result, "to_profile_dict"):
        flat_dict = discovery_result.to_profile_dict()
    elif isinstance(discovery_result, dict):
        flat_dict = discovery_result
    else:
        raise TypeError(
            f"Expected ProfileDiscoveryResult or dict, got {type(discovery_result)}"
        )

    # Count marker occurrences across profiles to find shared markers
    marker_owners: Dict[str, List[str]] = {}
    for pname, markers in flat_dict.items():
        for m in markers:
            marker_owners.setdefault(m, []).append(pname)

    shared_markers = {m for m, owners in marker_owners.items() if len(owners) > 1}

    # Build structured profiles: exclusive -> Major, shared -> Minor
    structured: Dict[str, Dict[str, List[str]]] = {}
    for pname, markers in flat_dict.items():
        major = [m for m in markers if m not in shared_markers]
        minor = [m for m in markers if m in shared_markers]
        structured[pname] = {"Major": major, "Minor": minor}

    if shared_markers:
        logger.info(f"Shared markers demoted to Minor: {sorted(shared_markers)}")

    # Apply renames
    if user_renames:
        for old_name, new_name in user_renames.items():
            if old_name in structured:
                structured[new_name] = structured.pop(old_name)

    # Merge user additions
    if user_additions:
        for name, spec in user_additions.items():
            if name in structured:
                # Merge into existing
                structured[name]["Major"] = list(
                    set(structured[name].get("Major", []) + spec.get("Major", []))
                )
                structured[name]["Minor"] = list(
                    set(structured[name].get("Minor", []) + spec.get("Minor", []))
                )
                structured[name].setdefault("Negative", []).extend(
                    spec.get("Negative", [])
                )
            else:
                structured[name] = spec

    # Auto-annotate using reference signatures (Jaccard matching)
    if reference_signatures:
        rename_map = {}
        for pname, spec in structured.items():
            all_markers = set(spec.get("Major", []) + spec.get("Minor", []))
            best_name = None
            best_jaccard = 0.0
            for ref_name, ref_markers in reference_signatures.items():
                ref_set = set(ref_markers)
                intersection = len(all_markers & ref_set)
                union = len(all_markers | ref_set)
                if union > 0:
                    jaccard = intersection / union
                    if jaccard > best_jaccard:
                        best_jaccard = jaccard
                        best_name = ref_name
            if best_name and best_jaccard > 0.3:
                rename_map[pname] = best_name

        for old_name, new_name in rename_map.items():
            if old_name in structured and new_name not in structured:
                structured[new_name] = structured.pop(old_name)
                logger.info(f"Auto-annotated '{old_name}' -> '{new_name}' (Jaccard)")

    return GatingProfileSet.from_flat_dict(structured, priority_dict=priority_dict)
