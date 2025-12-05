"""
Automatic Multi-Order Antibody Profile Discovery for CITEgeist.

This module provides unsupervised discovery of cell type profiles from antibody
capture data, eliminating the need for manual marker specification. It uses an
iterative greedy search with EM refinement and BIC-based model selection.

The algorithm discovers profiles of mixed sizes (singles, doubles, triplets)
and learns marker-specific scaling factors (β) that handle shared markers
appropriately (e.g., CD3D shared between CD4+ and CD8+ T cells).

Author: A. Chang
License: BSD 3-Clause
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from itertools import combinations
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from scipy.spatial import cKDTree

if TYPE_CHECKING:
    from numpy.random import Generator

logger = logging.getLogger(__name__)


@dataclass
class ProfileDiscoveryResult:
    """Container for profile discovery results."""

    profiles: Dict[str, Dict[str, List[str]]]
    """Discovered profiles in CITEgeist-compatible format."""

    beta: Dict[str, float]
    """Learned marker specificity weights."""

    proportions: NDArray[np.floating]
    """Estimated profile proportions per spot (n_spots × n_profiles)."""

    bic_trace: List[float]
    """BIC values at each iteration for diagnostics."""

    n_iterations: int
    """Number of discovery iterations performed."""

    metadata: Dict = field(default_factory=dict)
    """Additional metadata (seed, convergence info, etc.)."""


def discover_profiles(
    X: NDArray[np.floating] | sparse.spmatrix,
    marker_names: List[str],
    max_k: int = 3,
    seed: int = 1234,
    *,
    robust_zscore: bool = False,
    verbose: bool = True,
    n_perm: int = 500,
    alpha: float = 0.05,
    max_profiles: int = 20,
    coords: Optional[NDArray[np.floating]] = None,
    morans_k: int = 8,
) -> ProfileDiscoveryResult:
    """
    Discover antibody profiles from spatial transcriptomics data.

    Uses iterative greedy search with EM refinement and BIC model selection
    to find the optimal set of marker combinations that explain the data,
    automatically filtering out markers with no significant signal and
    prioritizing spatially structured markers via Moran's I.

    Parameters
    ----------
    X : array-like of shape (n_spots, n_markers)
        Antibody capture matrix, typically after winsorization and CLR normalization.
        Can be dense numpy array or scipy sparse matrix.
    marker_names : list of str
        Names of antibody markers corresponding to columns of X.
    max_k : int, default=3
        Maximum number of markers per profile. Biological constraint.
    seed : int, default=1234
        Random seed for reproducibility of permutation tests.
    robust_zscore : bool, default=False
        If True, use median/MAD instead of mean/std for z-score normalization.
        More robust to outliers.
    verbose : bool, default=True
        If True, log progress during discovery.
    n_perm : int, default=500
        Number of permutations for null distribution when scoring candidates.
    alpha : float, default=0.05
        Significance threshold for accepting candidates.
    max_profiles : int, default=20
        Safety cap on number of discovered profiles.
    coords : array-like of shape (n_spots, 2), optional
        Spatial coordinates (e.g., spot_x/spot_y). When provided, Moran's I
        is used to identify spatially informative markers before combination
        search.
    morans_k : int, default=8
        Number of nearest neighbors for Moran's I computation.

    Returns
    -------
    ProfileDiscoveryResult
        Contains discovered profiles, beta weights, proportions, and diagnostics.

    Examples
    --------
    >>> from CITEgeist.model import discover_profiles
    >>> result = discover_profiles(adata.X, list(adata.var_names), max_k=3)
    >>> print(result.profiles)
    {'EPCAM': {'Major': ['EPCAM']},
     'CD3D_CD4': {'Major': ['CD3D', 'CD4']},
     'CD3D_CD8': {'Major': ['CD3D', 'CD8']}}
    """
    rng = np.random.default_rng(seed)

    # Convert sparse to dense (antibody panels are small, ~10-50 markers)
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    n_spots, n_markers = X.shape

    # Spatial coordinates are REQUIRED for profile discovery
    # The algorithm uses Moran's I to identify spatially-structured markers
    if coords is None:
        raise ValueError(
            "coords is required for profile discovery. "
            "Spatial coordinates are used for Moran's I filtering to identify "
            "markers with significant spatial structure."
        )
    coords = np.asarray(coords)
    if coords.shape != (n_spots, 2):
        raise ValueError("coords must have shape (n_spots, 2)")

    if len(marker_names) != n_markers:
        raise ValueError(
            f"marker_names length ({len(marker_names)}) != X columns ({n_markers})"
        )

    # Standardize markers
    Z, valid_mask = _standardize_markers(X, robust=robust_zscore)
    valid_indices = np.where(valid_mask)[0]
    valid_names = [marker_names[i] for i in valid_indices]

    if len(valid_indices) < n_markers:
        dropped = [marker_names[i] for i in range(n_markers) if not valid_mask[i]]
        logger.warning(f"Dropped {len(dropped)} zero-variance markers: {dropped}")

    Z_valid = Z[:, valid_mask]

    # UPFRONT Moran's I filtering: identify spatially-structured markers ONCE
    # This ensures only markers with real spatial signal can participate in profiles
    eligible_markers = _identify_significant_single_markers(
        Z_valid,
        np.ones(len(valid_indices)),  # uniform beta for initial screening
        rng,
        n_perm,
        alpha,
        coords=coords,
        morans_k=morans_k,
    )

    if verbose:
        logger.info(
            f"Moran's I filter: {len(eligible_markers)}/{len(valid_indices)} "
            f"markers passed spatial autocorrelation test"
        )
        if eligible_markers:
            eligible_names = [valid_names[i] for i in eligible_markers]
            logger.debug(f"Eligible markers: {eligible_names}")

    if not eligible_markers:
        logger.warning("No markers passed Moran's I filter. Returning empty result.")
        return ProfileDiscoveryResult(
            profiles={},
            beta={name: 0.0 for name in marker_names},
            proportions=np.zeros((n_spots, 0)),
            bic_trace=[],
            n_iterations=0,
            metadata={
                "seed": seed,
                "max_k": max_k,
                "robust_zscore": robust_zscore,
                "alpha": alpha,
                "n_perm": n_perm,
                "morans_k": morans_k,
                "coords_provided": coords is not None,
                "eligible_markers": 0,
            },
        )

    # Score ALL candidates jointly on the original data
    beta = np.ones(len(valid_indices))

    if verbose:
        logger.info("Scoring all candidate profiles on original data...")

    all_candidates = _score_all_candidates(
        Z_valid,
        eligible_markers,
        beta,
        max_k,
        rng,
        n_perm=n_perm,
        alpha=alpha,
    )

    if verbose:
        logger.info(f"Found {len(all_candidates)} significant candidates")

    # Greedy selection of non-overlapping profiles
    profiles = _greedy_select_profiles(all_candidates, max_profiles)

    if verbose:
        logger.info(f"Selected {len(profiles)} non-overlapping profiles")
        for p in profiles:
            profile_str = "_".join(valid_names[i] for i in sorted(p))
            logger.info(f"  {profile_str}")

    # EM refinement on selected profiles
    bic_trace: List[float] = []
    if profiles:
        Y, beta, log_lik = _em_refine(Z_valid, profiles, beta)

        # Compute final BIC
        n_params = sum(len(p) for p in profiles) + len(profiles)
        bic = -2 * log_lik + n_params * np.log(n_spots)
        bic_trace.append(bic)
    else:
        Y = np.zeros((n_spots, 0))
        if verbose:
            logger.warning("No profiles discovered. Returning empty result.")

    # Format output
    profiles_dict, beta_dict = _format_output(
        profiles, beta, valid_names, marker_names, valid_mask
    )

    return ProfileDiscoveryResult(
        profiles=profiles_dict,
        beta=beta_dict,
        proportions=Y,
        bic_trace=bic_trace,
        n_iterations=len(bic_trace),
        metadata={
            "seed": seed,
            "max_k": max_k,
            "robust_zscore": robust_zscore,
            "alpha": alpha,
            "n_perm": n_perm,
            "morans_k": morans_k,
            "coords_provided": coords is not None,
            "eligible_markers": len(eligible_markers),
            "total_markers": len(valid_indices),
        },
    )


def _standardize_markers(
    X: NDArray[np.floating], robust: bool = False
) -> Tuple[NDArray[np.floating], NDArray[np.bool_]]:
    """
    Z-score standardize markers, handling zero-variance cases.

    Parameters
    ----------
    X : ndarray of shape (n_spots, n_markers)
    robust : bool
        Use median/MAD instead of mean/std.

    Returns
    -------
    Z : ndarray
        Standardized matrix.
    valid_mask : ndarray of bool
        Mask indicating which markers have non-zero variance.
    """
    n_spots, n_markers = X.shape
    Z = np.zeros_like(X)
    valid_mask = np.ones(n_markers, dtype=bool)
    eps = 1e-9

    for i in range(n_markers):
        col = X[:, i]
        if robust:
            center = np.median(col)
            scale = np.median(np.abs(col - center)) * 1.4826  # MAD to std
        else:
            center = col.mean()
            scale = col.std()

        if scale < eps:
            valid_mask[i] = False
            Z[:, i] = 0.0
        else:
            Z[:, i] = (col - center) / (scale + eps)

    return Z, valid_mask


def _score_all_candidates(
    Z: NDArray[np.floating],
    eligible_markers: List[int],
    beta: NDArray[np.floating],
    max_k: int,
    rng: Generator,
    n_perm: int = 500,
    alpha: float = 0.05,
) -> List[Tuple[Set[int], float, float]]:
    """
    Score ALL candidate profiles on the original data jointly.

    Evaluates all k-sets (k=1..max_k) from eligible_markers.
    Returns list of (candidate, score, p_value) for significant candidates.

    Parameters
    ----------
    Z : ndarray
        Standardized expression matrix (original, not residual).
    eligible_markers : list of int
        Marker indices that passed upfront Moran's I filter.
    beta : ndarray
        Marker weights.
    max_k : int
        Maximum markers per profile.
    rng : Generator
        Random number generator.
    n_perm : int
        Number of permutations for significance testing.
    alpha : float
        Significance threshold.

    Returns
    -------
    List of (candidate_set, score, p_value) tuples, sorted by score descending.
    Only includes candidates with p_value < alpha.
    """
    if not eligible_markers:
        return []

    significant_candidates: List[Tuple[Set[int], float, float]] = []
    candidate_markers = sorted(eligible_markers)

    # Score ALL candidates on original data
    for k in range(1, max_k + 1):
        for combo in combinations(candidate_markers, k):
            candidate = set(combo)

            score = _score_candidate(Z, candidate, beta)
            null_scores = _compute_null_distribution(Z, candidate, beta, rng, n_perm)
            p_value = (1 + np.sum(null_scores >= score)) / (1 + n_perm)

            if p_value < alpha:
                significant_candidates.append((candidate, score, p_value))

    # Sort by score descending
    significant_candidates.sort(key=lambda x: -x[1])

    return significant_candidates


def _greedy_select_profiles(
    candidates: List[Tuple[Set[int], float, float]],
    max_profiles: int,
) -> List[Set[int]]:
    """
    Greedily select non-overlapping profiles from scored candidates.

    Takes the highest-scoring candidate, excludes it, then takes next highest
    that doesn't share markers with already selected profiles.

    Parameters
    ----------
    candidates : list of (candidate_set, score, p_value)
        Significant candidates sorted by score descending.
    max_profiles : int
        Maximum number of profiles to select.

    Returns
    -------
    List of selected profile sets.
    """
    selected: List[Set[int]] = []
    used_markers: Set[int] = set()

    for candidate, score, p_value in candidates:
        if len(selected) >= max_profiles:
            break

        # Check if this candidate overlaps with already selected
        if candidate & used_markers:
            continue

        selected.append(candidate)
        used_markers.update(candidate)
        logger.debug(f"Selected profile {candidate}, score={score:.4f}, p={p_value:.4f}")

    return selected


def _score_candidate(
    Z: NDArray[np.floating],
    marker_set: Set[int],
    beta: NDArray[np.floating],
) -> float:
    """
    Compute β-weighted co-expression score.

    Score = mean over spots of product of β_i * max(z_i, 0) for markers in set.
    High when all markers are jointly elevated in the same spots.
    """
    markers = list(marker_set)
    joint = np.ones(Z.shape[0])
    for m in markers:
        joint *= beta[m] * np.maximum(Z[:, m], 0)
    return float(joint.mean())


def _compute_null_distribution(
    Z: NDArray[np.floating],
    marker_set: Set[int],
    beta: NDArray[np.floating],
    rng: Generator,
    n_perm: int,
) -> NDArray[np.floating]:
    """Generate null distribution via permutation."""
    markers = list(marker_set)
    n_spots = Z.shape[0]
    nulls = np.zeros(n_perm)

    if len(markers) == 1:
        m0 = markers[0]
        col = Z[:, m0]
        scale = col.std() if col.std() > 0 else 1.0
        for b in range(n_perm):
            noise = rng.normal(loc=0.0, scale=scale, size=n_spots)
            joint = beta[m0] * np.maximum(noise, 0)
            nulls[b] = joint.mean()
    else:
        for b in range(n_perm):
            perm_idx = rng.permutation(n_spots)
            joint = beta[markers[0]] * np.maximum(Z[perm_idx, markers[0]], 0)
            for m in markers[1:]:
                perm_idx = rng.permutation(n_spots)
                joint *= beta[m] * np.maximum(Z[perm_idx, m], 0)
            nulls[b] = joint.mean()

    return nulls


def _compute_morans_i_batch(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int,
    rng: Generator,
    n_perm: int,
) -> List[Tuple[float, float]]:
    """
    Compute Moran's I and permutation p-value for each marker.

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k : int
        Number of nearest neighbors.

    Returns
    -------
    List of tuples (I_obs, p_value) for each marker.
    """
    n_spots, n_markers = Z.shape

    # Edge case: need at least 2 spots for spatial autocorrelation
    if n_spots < 2:
        return [(np.nan, 1.0)] * n_markers

    tree = cKDTree(coords)
    # Query k+1 to include self then drop
    query_k = min(k + 1, n_spots)
    dists, idx = tree.query(coords, k=query_k)

    # Handle case where query returns 1D array (when query_k == 1)
    if idx.ndim == 1:
        return [(np.nan, 1.0)] * n_markers

    weights = np.ones_like(idx, dtype=float)
    # Remove self (distance 0)
    idx = idx[:, 1:]
    weights = weights[:, 1:]

    # Check if we have any neighbors after removing self
    if idx.shape[1] == 0:
        return [(np.nan, 1.0)] * n_markers

    S0 = weights.sum()
    if S0 == 0:
        return [(np.nan, 1.0)] * n_markers

    results: List[Tuple[float, float]] = []
    for m in range(n_markers):
        values = Z[:, m]
        if np.allclose(values, 0):
            results.append((np.nan, 1.0))
            continue

        z = values - values.mean()
        denom = np.sum(z**2)
        if denom == 0:
            results.append((np.nan, 1.0))
            continue

        neighbor_vals = z[idx]
        num = np.sum(weights * (z[:, None] * neighbor_vals))
        I_obs = (n_spots / S0) * (num / denom)

        # Permutation null
        null_I = np.zeros(n_perm)
        for b in range(n_perm):
            perm = rng.permutation(z)
            neighbor_perm = perm[idx]
            num_perm = np.sum(weights * (perm[:, None] * neighbor_perm))
            null_I[b] = (n_spots / S0) * (num_perm / denom)
        p_value = (1 + np.sum(null_I >= I_obs)) / (1 + len(null_I))
        results.append((I_obs, p_value))

    return results


def _identify_significant_single_markers(
    Z: NDArray[np.floating],
    beta: NDArray[np.floating],
    rng: Generator,
    n_perm: int,
    alpha: float,
    coords: Optional[NDArray[np.floating]] = None,
    morans_k: int = 8,
) -> List[int]:
    """
    Return markers with significant single-marker signal.

    If spatial coordinates are provided, Moran's I is used first to select
    spatially structured markers; if none pass, fall back to expression
    permutation testing.
    """
    n_markers = Z.shape[1]
    significant: List[int] = []

    if coords is not None:
        moran_stats = _compute_morans_i_batch(Z, coords, morans_k, rng, n_perm)
        for m, (I_obs, p_value) in enumerate(moran_stats):
            if np.isnan(I_obs):
                continue
            if I_obs > 0 and p_value < alpha:
                significant.append(m)
        if significant:
            return significant
        logger.debug("No markers passed Moran's I; falling back to expression permutation test.")

    for m in range(n_markers):
        score = _score_candidate(Z, {m}, beta)
        null_scores = _compute_null_distribution(Z, {m}, beta, rng, n_perm)
        p_value = (1 + np.sum(null_scores >= score)) / (1 + n_perm)
        if p_value < alpha:
            significant.append(m)

    return significant


def _em_refine(
    Z: NDArray[np.floating],
    profiles: List[Set[int]],
    beta_init: NDArray[np.floating],
    max_iter: int = 20,
    tol: float = 1e-4,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], float]:
    """
    EM algorithm to refine Y (proportions) and β (marker weights).

    Mirrors CITEgeist's existing EM for deconvolution, adapted for
    profile discovery phase.
    """
    n_spots, n_markers = Z.shape
    K = len(profiles)

    if K == 0:
        return np.zeros((n_spots, 0)), beta_init.copy(), 0.0

    A = _build_profile_matrix(profiles, n_markers)
    beta = beta_init.copy()
    Y = np.full((n_spots, K), 1.0 / K)

    prev_ll = -np.inf

    for _ in range(max_iter):
        # E-step: Update Y given β
        for k, profile in enumerate(profiles):
            markers = list(profile)
            # Gaussian likelihood for profile markers
            diff = Z[:, markers] - beta[markers]
            profile_fit = np.exp(-0.5 * np.sum(diff**2, axis=1))
            Y[:, k] = profile_fit

        # Normalize to proportions
        row_sums = Y.sum(axis=1, keepdims=True)
        Y = Y / np.maximum(row_sums, 1e-9)

        # M-step: Update β given Y
        for i in range(n_markers):
            containing = [k for k in range(K) if A[k, i] > 0]

            if len(containing) == 0:
                beta[i] = 0.1
            elif len(containing) == 1:
                beta[i] = 1.0
            else:
                # Shared marker: weight by concentration
                high_mask = Z[:, i] > 0
                n_high = high_mask.sum()
                if n_high > 10:
                    Y_high = Y[high_mask][:, containing]
                    concentration = Y_high.max(axis=1).mean()
                    beta[i] = 0.3 + 0.7 * concentration
                else:
                    beta[i] = 0.5

        # Log-likelihood for convergence
        expected = Y @ A @ np.diag(beta)
        ll = -0.5 * np.sum((Z - expected) ** 2)

        if abs(ll - prev_ll) < tol:
            break
        prev_ll = ll

    return Y, beta, ll


def _build_profile_matrix(
    profiles: List[Set[int]], n_markers: int
) -> NDArray[np.floating]:
    """Build binary (K × n_markers) profile definition matrix."""
    K = len(profiles)
    A = np.zeros((K, n_markers))
    for k, profile in enumerate(profiles):
        for i in profile:
            A[k, i] = 1.0
    return A


def _format_output(
    profiles: List[Set[int]],
    beta: NDArray[np.floating],
    valid_names: List[str],
    all_names: List[str],
    valid_mask: NDArray[np.bool_],
) -> Tuple[Dict[str, Dict[str, List[str]]], Dict[str, float]]:
    """Format profiles and beta into CITEgeist-compatible dictionaries."""
    profiles_dict: Dict[str, Dict[str, List[str]]] = {}

    for profile in profiles:
        markers = sorted([valid_names[i] for i in profile])
        name = "_".join(markers)
        profiles_dict[name] = {"Major": markers}

    beta_dict: Dict[str, float] = {}
    valid_idx = 0
    for i, name in enumerate(all_names):
        if valid_mask[i]:
            beta_dict[name] = float(beta[valid_idx])
            valid_idx += 1
        else:
            beta_dict[name] = 0.0

    return profiles_dict, beta_dict


# Convenience function for integration with CitegeistModel
def integrate_with_model(
    model: "CitegeistModel",  # type: ignore[name-defined]
    max_k: int = 3,
    seed: int = 1234,
    **kwargs,
) -> ProfileDiscoveryResult:
    """
    Convenience function to run profile discovery on a CitegeistModel.

    Parameters
    ----------
    model : CitegeistModel
        Initialized model with antibody_capture_adata loaded.
    max_k : int
        Maximum profile size.
    seed : int
        Random seed.
    **kwargs
        Additional arguments passed to discover_profiles.

    Returns
    -------
    ProfileDiscoveryResult
        Discovery results. Also sets model.cell_type_profiles and model.marker_beta.
    """
    if model.antibody_capture_adata is None:
        raise ValueError("Model must have antibody_capture_adata loaded.")

    X = model.antibody_capture_adata.X
    marker_names = list(model.antibody_capture_adata.var_names)

    result = discover_profiles(X, marker_names, max_k=max_k, seed=seed, **kwargs)

    # Set on model for downstream use
    model.cell_type_profiles = result.profiles
    model.marker_beta = result.beta

    logger.info(
        f"Discovered {len(result.profiles)} profiles: {list(result.profiles.keys())}"
    )

    return result
