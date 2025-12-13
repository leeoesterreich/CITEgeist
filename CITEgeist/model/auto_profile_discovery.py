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
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import kurtosis as scipy_kurtosis

# Local spatial statistics from esda/libpysal
try:
    from esda.moran import Moran_Local
    from esda.getisord import G_Local
    from libpysal.weights import KNN as LibPySAL_KNN
    HAS_LOCAL_SPATIAL = True
except ImportError:
    HAS_LOCAL_SPATIAL = False
    Moran_Local = None
    G_Local = None
    LibPySAL_KNN = None

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

    # Spatial Mixture Model results (if SMM was applied)
    smm_applied: bool = False
    """Whether SMM background correction was applied."""

    smm_filtered_markers: List[str] = field(default_factory=list)
    """Markers filtered out by SMM due to low SNR (legacy, now empty)."""

    smm_snr_values: Optional[Dict[str, float]] = None
    """Per-marker SNR values from SMM for ALL markers (None if SMM not applied)."""

    smm_signal_fractions: Optional[Dict[str, float]] = None
    """Per-marker signal fractions from SMM for ALL markers (None if SMM not applied)."""

    smm_beta_learned: Optional[float] = None
    """Learned spatial regularization parameter from SMM (None if SMM not applied)."""

    # NEW: Data for diagnostic plotting
    smm_raw_matrix: Optional[NDArray[np.floating]] = None
    """Raw expression matrix before correction (n_spots, n_markers)."""

    smm_corrected_matrix: Optional[NDArray[np.floating]] = None
    """Background-corrected matrix (X * P(signal)) (n_spots, n_markers)."""

    smm_smoothed_matrix: Optional[NDArray[np.floating]] = None
    """Smoothed corrected matrix (n_spots, n_markers)."""

    smm_signal_posteriors: Optional[NDArray[np.floating]] = None
    """Per-spot signal probability P(signal | x) (n_spots, n_markers)."""

    smm_gmm_params: Optional[Dict[str, NDArray[np.floating]]] = None
    """GMM parameters {'mu': (n_markers, 2), 'sigma': (n_markers, 2)}."""


def discover_profiles(
    X: NDArray[np.floating] | sparse.spmatrix,
    marker_names: List[str],
    max_k: int = 3,
    seed: int = 1234,
    *,
    robust_zscore: bool = False,
    verbose: bool = True,
    n_perm: int = 500,
    alpha: float = 0.1,
    max_profiles: int = 20,
    coords: Optional[NDArray[np.floating]] = None,
    morans_k: Union[int, List[int]] = 8,
    morans_i_threshold: float = 0.1,
    per_marker_fallback: bool = True,
    model_selection: str = "cv",
    cv_folds: int = 5,
    min_profiles: int = 2,
    hierarchical: bool = False,
    # NEW: Local spatial statistics parameters
    use_local_spatial: bool = True,
    local_method: str = "both",
    hotspot_threshold: float = 0.10,
    local_permutations: int = 99,
    # NEW: Intensity rescue parameters
    use_intensity_rescue: bool = True,
    intensity_percentile: float = 0.75,
    min_hotspot_fraction: float = 0.03,
    # NEW: Abundance-adaptive parameters
    use_abundance_adaptive: bool = True,
    ubiquitous_cv_threshold: float = 0.5,
    ubiquitous_presence_threshold: float = 0.7,
    rare_presence_threshold: float = 0.10,
    rare_intensity_factor: float = 2.0,
    # NEW: Profile selection method
    selection_method: str = "reconstruction",
    min_reconstruction_improvement: float = 0.05,
    redundancy_threshold: float = 0.9,
    allow_overlap: Union[bool, str] = "auto",
    # NEW: MIQP-specific parameters
    miqp_lambda_spatial: float = 0.1,
    miqp_lambda_complexity: float = 0.01,
    miqp_time_limit: float = 300.0,
    miqp_gap: float = 0.01,
    # NEW: Hierarchical MIQP parameters
    miqp_lambda_overlap: float = 0.5,
    miqp_lambda_orphan: float = 0.2,
    miqp_lambda_sparsity: float = 0.3,
    miqp_enforce_hierarchy: bool = False,
    miqp_sparsity_aggregation: str = "mean",
    # NEW: Spatial Mixture Model (SMM) background correction parameters
    use_smm: bool = False,  # Default OFF for backward compatibility; enable for background filtering
    smm_k_neighbors: int = 6,
    smm_snr_threshold: float = 1.5,
    smm_min_signal_fraction: float = 0.05,
    smm_max_signal_fraction: float = 0.95,
    smm_beta_init: float = 1.0,
    smm_max_iter: int = 50,
    # NEW: Spatial smoothing parameters (applied before SMM)
    smm_apply_smoothing: bool = True,  # Apply spatial smoothing before SMM
    smm_smoothing_sigma: float = 1.5,  # Gaussian kernel bandwidth
    smm_smoothing_k: int = 6,  # Number of neighbors for smoothing
    smm_smoothing_self_weight: float = 0.5,  # Weight for spot's own value
) -> ProfileDiscoveryResult:
    """
    Discover antibody profiles from spatial transcriptomics data.

    Uses a three-tier detection hierarchy to identify markers with significant
    spatial signal, then combines them into cell type profiles via greedy
    search with EM refinement and model selection.

    Three-Tier Marker Detection:
    1. **Tier 1 (Global Spatial)**: Global Moran's I >= threshold at ANY scale
    2. **Tier 2 (Local Spatial)**: >hotspot_threshold of spots in significant
       LISA HH clusters or Getis-Ord hotspots (finds dispersed but strong signals)
    3. **Tier 3 (Intensity Rescue)**: High variance/kurtosis + minimal local structure

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
    alpha : float, default=0.1
        Significance threshold for accepting candidates.
    max_profiles : int, default=20
        Safety cap on number of discovered profiles.
    coords : array-like of shape (n_spots, 2), optional
        Spatial coordinates (e.g., spot_x/spot_y). Required for spatial testing.
    morans_k : int or list of int, default=8
        Number of nearest neighbors for Moran's I computation. Can be a single
        int (backward compatible) or a list for multi-scale testing.
        For mixed/dysplastic tumor environments, use [3, 5, 8, 12].
    morans_i_threshold : float, default=0.1
        Minimum Moran's I value for Tier 1. Higher threshold rejects more noise.
        For mixed tissues with weak spatial signals, set to 0.05 or lower.
    per_marker_fallback : bool, default=True
        If True, markers failing all spatial tiers get expression-based testing.
    model_selection : str, default='cv'
        Model selection strategy:
        - 'cv': Use spatial cross-validation (recommended for mixed tissues)
        - 'bic': Use BIC (original behavior, may stop too early)
        - 'greedy': No model selection, use all significant candidates
    cv_folds : int, default=5
        Number of folds for cross-validation (only used when model_selection='cv').
    min_profiles : int, default=2
        Minimum number of profiles before model selection stopping applies.
        Prevents premature termination after discovering one dominant profile.
    hierarchical : bool, default=False
        If True, use hierarchical profile discovery that allows shared markers.
    use_local_spatial : bool, default=True
        If True, enable Tier 2 (LISA/Getis-Ord local spatial statistics).
        This is the primary mechanism for detecting dispersed but strong signals
        like T-cells scattered in tumor stroma.
    local_method : str, default="both"
        Which local statistic to use: "lisa", "getis", or "both".
    hotspot_threshold : float, default=0.10
        Minimum fraction of spots in significant hotspots for Tier 2.
        Lower values (0.05) for rare cell types, higher (0.15) for specificity.
    local_permutations : int, default=99
        Number of permutations for local statistics p-values.
    use_intensity_rescue : bool, default=True
        If True, enable Tier 3 (intensity-based rescue) for high-variance markers.
    intensity_percentile : float, default=0.75
        Minimum variance/kurtosis percentile for Tier 3.
    min_hotspot_fraction : float, default=0.03
        Minimum hotspot fraction for Tier 3 (relaxed threshold).
    miqp_lambda_spatial : float, default=0.1
        Weight for spatial penalty in MIQP optimization. Higher values prefer
        profiles with stronger spatial signals. Only used when selection_method="miqp".
    miqp_lambda_complexity : float, default=0.01
        Weight for profile count penalty in MIQP. Higher values prefer fewer profiles.
    miqp_time_limit : float, default=300.0
        Time limit in seconds for MIQP solver.
    miqp_gap : float, default=0.01
        Acceptable MIP optimality gap (1% default).

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

    # Initialize SMM tracking variables
    smm_result = None
    snr_weights = None
    smm_filtered_markers = []
    original_marker_names = list(marker_names)  # Keep original for tracking

    # ========================================================================
    # TIER -1: Spatial Mixture Model (SMM) Background Correction
    # ========================================================================
    # Pipeline: GMM on RAW data -> Per-spot soft-scale correction -> Smoothing
    # All markers proceed through the pipeline (no pre-filtering by SNR).
    # Filtering decisions are made downstream based on Moran's I.
    if use_smm and coords is not None:
        from .background_correction import fit_spatial_mixture_model

        if verbose:
            logger.info("Running Spatial Mixture Model (SMM) background correction...")
            logger.info("Pipeline: GMM(raw) -> Soft-scale(X * P(signal)) -> Smooth")

        smm_result = fit_spatial_mixture_model(
            X, coords, marker_names,
            k_neighbors=smm_k_neighbors,
            snr_threshold=smm_snr_threshold,
            min_signal_fraction=smm_min_signal_fraction,
            max_signal_fraction=smm_max_signal_fraction,
            beta_init=smm_beta_init,
            max_iter=smm_max_iter,
            apply_smoothing=smm_apply_smoothing,
            smoothing_sigma=smm_smoothing_sigma,
            smoothing_k_neighbors=smm_smoothing_k,
            smoothing_self_weight=smm_smoothing_self_weight,
            verbose=verbose,
        )

        # Use the SMOOTHED corrected matrix for downstream analysis
        # All markers proceed - no SNR-based filtering
        X = smm_result.smoothed_matrix
        snr_weights = smm_result.snr_values

        if verbose:
            logger.info(
                f"SMM complete: All {n_markers} markers proceed to Moran's I filtering, "
                f"β={smm_result.beta_learned:.3f}"
            )
            # Log SNR info for diagnostics
            signal_fracs = smm_result.metadata.get("signal_fractions", {})
            for m_name in marker_names[:5]:  # First 5 markers
                snr = smm_result.snr_values.get(m_name, 0)
                sig_frac = signal_fracs.get(m_name, 0)
                logger.debug(f"  {m_name}: SNR={snr:.2f}, sig_frac={sig_frac:.2f}")

    # Standardize markers
    Z, valid_mask = _standardize_markers(X, robust=robust_zscore)
    valid_indices = np.where(valid_mask)[0]
    valid_names = [marker_names[i] for i in valid_indices]

    if len(valid_indices) < n_markers:
        dropped = [marker_names[i] for i in range(n_markers) if not valid_mask[i]]
        logger.warning(f"Dropped {len(dropped)} zero-variance markers: {dropped}")

    Z_valid = Z[:, valid_mask]

    # UPFRONT marker filtering: identify significant markers using abundance-adaptive detection
    # This ensures only markers with real signal can participate in profiles
    # Tier 0 (abundance-adaptive) catches ubiquitous and rare markers
    # Tiers 1-3 catch spatially-structured markers at different scales
    eligible_markers = _identify_significant_single_markers(
        Z_valid,
        np.ones(len(valid_indices)),  # uniform beta for initial screening
        rng,
        n_perm,
        alpha,
        coords=coords,
        morans_k=morans_k,
        morans_i_threshold=morans_i_threshold,
        per_marker_fallback=per_marker_fallback,
        # Local spatial statistics parameters
        use_local_spatial=use_local_spatial,
        local_method=local_method,
        hotspot_threshold=hotspot_threshold,
        local_permutations=local_permutations,
        # Intensity rescue parameters
        use_intensity_rescue=use_intensity_rescue,
        intensity_percentile=intensity_percentile,
        min_hotspot_fraction=min_hotspot_fraction,
        # Abundance-adaptive parameters
        use_abundance_adaptive=use_abundance_adaptive,
        ubiquitous_cv_threshold=ubiquitous_cv_threshold,
        ubiquitous_presence_threshold=ubiquitous_presence_threshold,
        rare_presence_threshold=rare_presence_threshold,
        rare_intensity_factor=rare_intensity_factor,
    )

    if verbose:
        logger.info(
            f"Abundance-adaptive filter: {len(eligible_markers)}/{len(valid_indices)} "
            f"markers passed significance testing"
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
                "morans_i_threshold": morans_i_threshold,
                "per_marker_fallback": per_marker_fallback,
                "model_selection": model_selection,
                "cv_folds": cv_folds,
                "min_profiles": min_profiles,
                "hierarchical": hierarchical,
                "cv_score": None,
                "coords_provided": coords is not None,
                "eligible_markers": 0,
            },
            smm_applied=smm_result is not None,
            smm_filtered_markers=smm_filtered_markers,
            smm_snr_values=smm_result.snr_values if smm_result else None,
            smm_signal_fractions=smm_result.metadata.get("signal_fractions") if smm_result else None,
            smm_beta_learned=smm_result.beta_learned if smm_result else None,
            # NEW: Data for diagnostic plotting
            smm_raw_matrix=smm_result.raw_matrix if smm_result else None,
            smm_corrected_matrix=smm_result.corrected_matrix if smm_result else None,
            smm_smoothed_matrix=smm_result.smoothed_matrix if smm_result else None,
            smm_signal_posteriors=smm_result.signal_posteriors if smm_result else None,
            smm_gmm_params=smm_result.gmm_params if smm_result else None,
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

    # Select profiles based on selection_method and model_selection strategy
    # selection_method determines HOW profiles are scored/selected
    # model_selection determines stopping criterion (for permutation-based methods)

    if selection_method == "reconstruction":
        # Reconstruction-based selection: include profiles that reduce reconstruction error
        # This naturally handles both flat and hierarchical structures
        #
        # IMPORTANT: Only use candidates that passed the permutation significance test.
        # This prevents spurious noise correlations from being selected.
        # The reconstruction score is used to RANK among significant candidates,
        # not to determine significance.
        if verbose:
            logger.info(
                f"Using reconstruction-based selection "
                f"(min_improvement={min_reconstruction_improvement})"
            )

        # Only use candidates that passed permutation significance test
        # Do NOT add extra profiles that didn't pass the test
        candidate_sets = [cand for cand, score, pval in all_candidates]

        if not candidate_sets:
            logger.warning(
                "No candidates passed permutation significance test. "
                "Reconstruction selection has no candidates to work with."
            )
            profiles = []
            detected_hierarchical = False
            shared_markers_detected = []
        else:
            # Automatic detection of hierarchical vs flat structure
            if allow_overlap == "auto":
                detected_hierarchical, shared_markers_detected = _detect_hierarchical_structure(
                    candidate_sets
                )
                allow_overlap_bool = detected_hierarchical
                if verbose:
                    if detected_hierarchical:
                        shared_names = [valid_names[i] for i in shared_markers_detected]
                        logger.info(
                            f"Auto-detected HIERARCHICAL structure: "
                            f"{len(shared_markers_detected)} shared markers: {shared_names}"
                        )
                    else:
                        logger.info(
                            "Auto-detected FLAT structure: no shared markers found"
                        )
            else:
                allow_overlap_bool = bool(allow_overlap)
                detected_hierarchical = allow_overlap_bool
                shared_markers_detected = []

            if verbose:
                logger.info(
                    f"Reconstruction selection working with {len(candidate_sets)} "
                    f"permutation-significant candidates (allow_overlap={allow_overlap_bool})"
                )
            profiles = _select_profiles_reconstruction_based(
                Z_valid,
                candidate_sets,
                min_improvement=min_reconstruction_improvement,
                max_profiles=max_profiles,
                redundancy_threshold=redundancy_threshold,
                allow_overlap=allow_overlap_bool,
            )

    elif selection_method == "miqp":
        # MIQP-based selection using Gurobi joint optimization
        from .gurobi_impl import optimize_profiles_miqp

        if verbose:
            logger.info("Using MIQP-based joint optimization (Gurobi)")

        # Only use permutation-significant candidates
        candidate_sets = [cand for cand, score, pval in all_candidates]
        if not candidate_sets:
            logger.warning(
                "No candidates passed permutation significance test. "
                "MIQP has no candidates to work with."
            )
            profiles = []
            detected_hierarchical = False
            shared_markers_detected = []
        else:
            # Compute composite spatial scores (permissive MAX approach)
            k_scales = morans_k if isinstance(morans_k, list) else [morans_k]
            if verbose:
                logger.info(f"Computing composite spatial scores at scales {k_scales}")
            marker_spatial_scores = compute_composite_spatial_score(
                Z_valid, coords, k_scales=k_scales, n_perm=n_perm, alpha=alpha, seed=seed
            )

            # Auto-detect hierarchical structure
            if allow_overlap == "auto":
                detected_hierarchical, shared_markers_detected = _detect_hierarchical_structure(
                    candidate_sets
                )
                allow_overlap_bool = detected_hierarchical
                if verbose:
                    if detected_hierarchical:
                        shared_names = [valid_names[i] for i in shared_markers_detected]
                        logger.info(
                            f"Auto-detected HIERARCHICAL structure: "
                            f"{len(shared_markers_detected)} shared markers: {shared_names}"
                        )
                    else:
                        logger.info("Auto-detected FLAT structure: no shared markers found")
            else:
                allow_overlap_bool = bool(allow_overlap)
                detected_hierarchical = allow_overlap_bool
                shared_markers_detected = []

            # Run MIQP optimization
            try:
                profiles, miqp_proportions, miqp_metadata = optimize_profiles_miqp(
                    Z_valid,
                    candidate_sets,
                    marker_spatial_scores,
                    lambda_spatial=miqp_lambda_spatial,
                    lambda_complexity=miqp_lambda_complexity,
                    max_profiles=max_profiles,
                    min_profiles=min_profiles,
                    allow_overlap=allow_overlap_bool,
                    time_limit=miqp_time_limit,
                    mip_gap=miqp_gap,
                    seed=seed,
                    verbose=verbose,
                )

                if miqp_metadata.get("status") in ("failed", "gurobi_error", "error", "no_candidates"):
                    raise RuntimeError(f"MIQP failed: {miqp_metadata.get('status')}")

                if verbose:
                    logger.info(f"MIQP optimization status: {miqp_metadata.get('status')}")
                    logger.info(f"MIQP objective value: {miqp_metadata.get('objective', 'N/A')}")

            except Exception as e:
                logger.warning(f"MIQP optimization failed ({e}), falling back to reconstruction")
                profiles = _select_profiles_reconstruction_based(
                    Z_valid,
                    candidate_sets,
                    min_improvement=min_reconstruction_improvement,
                    max_profiles=max_profiles,
                    redundancy_threshold=redundancy_threshold,
                    allow_overlap=allow_overlap_bool,
                )

    elif selection_method == "miqp_hierarchical":
        # Hierarchical MIQP-based selection with spot-level sparsity
        from .gurobi_impl import optimize_profiles_miqp_hierarchical

        if verbose:
            logger.info("Using Hierarchical MIQP optimization (Gurobi)")
            logger.info(
                f"  Hierarchy params: lambda_overlap={miqp_lambda_overlap}, "
                f"lambda_orphan={miqp_lambda_orphan}, lambda_sparsity={miqp_lambda_sparsity}"
            )

        # Only use permutation-significant candidates
        candidate_sets = [cand for cand, score, pval in all_candidates]
        if not candidate_sets:
            logger.warning(
                "No candidates passed permutation significance test. "
                "Hierarchical MIQP has no candidates to work with."
            )
            profiles = []
            detected_hierarchical = True  # Hierarchical method was requested
            shared_markers_detected = []
        else:
            # Compute composite spatial scores
            k_scales = morans_k if isinstance(morans_k, list) else [morans_k]
            if verbose:
                logger.info(f"Computing composite spatial scores at scales {k_scales}")
            marker_spatial_scores = compute_composite_spatial_score(
                Z_valid, coords, k_scales=k_scales, n_perm=n_perm, alpha=alpha, seed=seed
            )

            # Hierarchical MIQP always allows overlap (it uses soft hierarchy penalties)
            detected_hierarchical = True
            shared_markers_detected = _find_shared_markers(candidate_sets)

            if verbose and shared_markers_detected:
                shared_names = [valid_names[i] for i in shared_markers_detected]
                logger.info(f"Shared markers detected: {shared_names}")

            # Run hierarchical MIQP optimization
            try:
                profiles, miqp_proportions, miqp_metadata = optimize_profiles_miqp_hierarchical(
                    Z_valid,
                    candidate_sets,
                    marker_spatial_scores,
                    lambda_spatial=miqp_lambda_spatial,
                    lambda_complexity=miqp_lambda_complexity,
                    lambda_overlap=miqp_lambda_overlap,
                    lambda_orphan=miqp_lambda_orphan,
                    lambda_sparsity=miqp_lambda_sparsity,
                    max_profiles=max_profiles,
                    min_profiles=min_profiles,
                    enforce_hierarchy=miqp_enforce_hierarchy,
                    sparsity_aggregation=miqp_sparsity_aggregation,
                    time_limit=miqp_time_limit,
                    mip_gap=miqp_gap,
                    seed=seed,
                    verbose=verbose,
                )

                if miqp_metadata.get("status") in ("failed", "gurobi_error", "error", "no_candidates"):
                    raise RuntimeError(f"Hierarchical MIQP failed: {miqp_metadata.get('status')}")

                if verbose:
                    logger.info(f"Hierarchical MIQP status: {miqp_metadata.get('status')}")
                    logger.info(f"Hierarchical MIQP objective: {miqp_metadata.get('objective', 'N/A')}")
                    if "hierarchy_levels" in miqp_metadata:
                        level_counts = {
                            lv: len(indices) for lv, indices in miqp_metadata["hierarchy_levels"].items()
                        }
                        logger.info(f"Candidate hierarchy levels: {level_counts}")
                    if "selected_levels" in miqp_metadata:
                        from collections import Counter
                        selected_level_counts = Counter(miqp_metadata["selected_levels"].values())
                        logger.info(f"Selected profiles by level: {dict(selected_level_counts)}")

            except Exception as e:
                logger.warning(f"Hierarchical MIQP failed ({e}), falling back to reconstruction")
                profiles = _select_profiles_reconstruction_based(
                    Z_valid,
                    candidate_sets,
                    min_improvement=min_reconstruction_improvement,
                    max_profiles=max_profiles,
                    redundancy_threshold=redundancy_threshold,
                    allow_overlap=True,  # Hierarchical allows overlap
                )

    else:  # selection_method == "permutation" or default
        # Original permutation-based selection
        # For permutation method, hierarchical structure is controlled by the
        # `hierarchical` parameter, not auto-detected
        detected_hierarchical = hierarchical
        shared_markers_detected = []
        if hierarchical:
            if verbose:
                logger.info("Using hierarchical profile discovery")
            profiles = _discover_hierarchical_profiles(
                Z_valid,
                coords,
                eligible_markers,
                beta,
                rng,
                max_k=max_k,
                n_perm=n_perm,
                alpha=alpha,
                max_profiles=max_profiles,
            )
        elif model_selection == "cv":
            if verbose:
                logger.info(f"Using CV-based model selection ({cv_folds}-fold)")
            profiles = _discover_with_cv_selection(
                Z_valid,
                coords,
                all_candidates,
                beta,
                max_profiles=max_profiles,
                min_profiles=min_profiles,
                cv_folds=cv_folds,
                seed=seed,
            )
        elif model_selection == "greedy":
            if verbose:
                logger.info("Using greedy selection (no model selection)")
            profiles = _greedy_select_profiles(all_candidates, max_profiles)
        else:  # "bic" or default
            if verbose:
                logger.info("Using greedy selection with BIC")
            profiles = _greedy_select_profiles(all_candidates, max_profiles)

    selection_type = selection_method if selection_method in ("reconstruction", "miqp", "miqp_hierarchical") else (
        "hierarchical" if hierarchical else "permutation"
    )
    if verbose:
        logger.info(f"Selected {len(profiles)} {selection_type} profiles")
        for p in profiles:
            profile_str = "_".join(valid_names[i] for i in sorted(p))
            logger.info(f"  {profile_str}")

    # EM refinement on selected profiles
    bic_trace: List[float] = []
    cv_score: Optional[float] = None
    if profiles:
        Y, beta, log_lik = _em_refine(Z_valid, profiles, beta)

        # Compute final BIC for diagnostics
        n_params = sum(len(p) for p in profiles) + len(profiles)
        bic = -2 * log_lik + n_params * np.log(n_spots)
        bic_trace.append(bic)

        # Also compute CV score for diagnostics
        if model_selection == "cv":
            cv_score = _spatial_cv_likelihood(
                Z_valid, coords, profiles, beta, n_folds=cv_folds, seed=seed
            )
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
            "morans_i_threshold": morans_i_threshold,
            "per_marker_fallback": per_marker_fallback,
            "model_selection": model_selection,
            "cv_folds": cv_folds,
            "min_profiles": min_profiles,
            "hierarchical": hierarchical,
            "cv_score": cv_score,
            "coords_provided": coords is not None,
            "eligible_markers": len(eligible_markers),
            "total_markers": len(valid_indices),
            # Abundance-adaptive parameters
            "use_abundance_adaptive": use_abundance_adaptive,
            "ubiquitous_cv_threshold": ubiquitous_cv_threshold,
            "rare_presence_threshold": rare_presence_threshold,
            # Selection method parameters
            "selection_method": selection_method,
            "min_reconstruction_improvement": min_reconstruction_improvement,
            "redundancy_threshold": redundancy_threshold,
            "allow_overlap": allow_overlap,
            # Auto-detection results
            "detected_hierarchical": detected_hierarchical,
            "shared_markers_detected": [valid_names[i] for i in shared_markers_detected] if shared_markers_detected else [],
            # SNR weights for downstream MIQP
            "snr_weights": snr_weights,
        },
        smm_applied=smm_result is not None,
        smm_filtered_markers=smm_filtered_markers,
        smm_snr_values=smm_result.snr_values if smm_result else None,
        smm_signal_fractions=smm_result.metadata.get("signal_fractions") if smm_result else None,
        smm_beta_learned=smm_result.beta_learned if smm_result else None,
        # NEW: Data for diagnostic plotting
        smm_raw_matrix=smm_result.raw_matrix if smm_result else None,
        smm_corrected_matrix=smm_result.corrected_matrix if smm_result else None,
        smm_smoothed_matrix=smm_result.smoothed_matrix if smm_result else None,
        smm_signal_posteriors=smm_result.signal_posteriors if smm_result else None,
        smm_gmm_params=smm_result.gmm_params if smm_result else None,
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
    alpha: float = 0.1,
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


def _compute_morans_i_multiscale(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k_scales: List[int],
    rng: "Generator",
    n_perm: int,
) -> Dict[int, List[Tuple[float, float]]]:
    """
    Compute Moran's I at multiple spatial scales for all markers.

    This allows detection of:
    - Fine-grained clustering (small k): TILs, scattered fibroblasts
    - Medium clustering (medium k): T-cell aggregates
    - Broad domains (large k): Cancer nests, stromal regions

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k_scales : list of int
        List of k values (neighbor counts) to test, e.g., [3, 5, 8, 12].
    rng : Generator
        Random number generator for permutation tests.
    n_perm : int
        Number of permutations for null distribution.

    Returns
    -------
    Dict mapping k -> list of (I_obs, p_value) tuples for each marker.
    """
    n_spots, n_markers = Z.shape
    results: Dict[int, List[Tuple[float, float]]] = {}

    # Edge case: need at least 2 spots
    if n_spots < 2:
        for k in k_scales:
            results[k] = [(np.nan, 1.0)] * n_markers
        return results

    # Build KDTree once
    tree = cKDTree(coords)

    # Sort scales to process largest first (reuse neighbor queries)
    sorted_scales = sorted(k_scales, reverse=True)
    max_k = sorted_scales[0]

    # Query for maximum k+1 neighbors (includes self)
    query_k = min(max_k + 1, n_spots)
    _, all_idx = tree.query(coords, k=query_k)

    # Handle 1D array case (single neighbor)
    if all_idx.ndim == 1:
        for k in k_scales:
            results[k] = [(np.nan, 1.0)] * n_markers
        return results

    # Remove self (first column is always self with distance 0)
    all_idx = all_idx[:, 1:]  # Shape: (n_spots, max_k)

    for k in k_scales:
        # Subset to k neighbors
        effective_k = min(k, all_idx.shape[1])
        if effective_k == 0:
            results[k] = [(np.nan, 1.0)] * n_markers
            continue

        idx = all_idx[:, :effective_k]
        weights = np.ones_like(idx, dtype=float)

        S0 = weights.sum()
        if S0 == 0:
            results[k] = [(np.nan, 1.0)] * n_markers
            continue

        scale_results: List[Tuple[float, float]] = []
        for m in range(n_markers):
            values = Z[:, m]
            if np.allclose(values, 0):
                scale_results.append((np.nan, 1.0))
                continue

            z = values - values.mean()
            denom = np.sum(z**2)
            if denom == 0:
                scale_results.append((np.nan, 1.0))
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
            scale_results.append((I_obs, p_value))

        results[k] = scale_results

    return results


# =============================================================================
# LOCAL SPATIAL STATISTICS (LISA, Getis-Ord Gi*)
# =============================================================================


def _compute_lisa_batch(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int,
    alpha: float = 0.05,
    permutations: int = 99,
) -> List[Dict[str, float]]:
    """
    Compute Local Moran's I (LISA) for each marker.

    LISA identifies local clusters (High-High, Low-Low) and outliers
    (High-Low, Low-High) at each spatial location. This enables detection
    of markers with significant local spatial structure even when the
    global Moran's I is weak or negative.

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k : int
        Number of nearest neighbors for spatial weights.
    alpha : float, default=0.05
        Significance threshold for local clusters.
    permutations : int, default=99
        Number of permutations for pseudo p-values. Use 99 for screening,
        999 for final analysis.

    Returns
    -------
    List of dicts, one per marker, containing:
        - 'significant_HH_fraction': Proportion of spots in significant HH clusters
        - 'significant_LL_fraction': Proportion of spots in significant LL clusters
        - 'significant_any_fraction': Proportion in any significant cluster
    """
    if not HAS_LOCAL_SPATIAL:
        logger.warning(
            "esda/libpysal not available; LISA computation skipped. "
            "Install with: pip install esda libpysal"
        )
        return [{'significant_HH_fraction': 0.0, 'significant_LL_fraction': 0.0,
                 'significant_any_fraction': 0.0}] * Z.shape[1]

    n_spots, n_markers = Z.shape

    # Edge case: need sufficient spots for local statistics
    if n_spots < k + 1:
        return [{'significant_HH_fraction': 0.0, 'significant_LL_fraction': 0.0,
                 'significant_any_fraction': 0.0}] * n_markers

    # Build spatial weights using libpysal
    try:
        w = LibPySAL_KNN.from_array(coords, k=k)
        w.transform = 'r'  # Row-standardize
    except Exception as e:
        logger.warning(f"Failed to build spatial weights: {e}")
        return [{'significant_HH_fraction': 0.0, 'significant_LL_fraction': 0.0,
                 'significant_any_fraction': 0.0}] * n_markers

    results: List[Dict[str, float]] = []
    for m in range(n_markers):
        values = Z[:, m]

        # Skip constant markers
        if np.allclose(values, values[0]):
            results.append({'significant_HH_fraction': 0.0,
                           'significant_LL_fraction': 0.0,
                           'significant_any_fraction': 0.0})
            continue

        try:
            lisa = Moran_Local(values, w, permutations=permutations)
            # lisa.q: quadrant (1=HH, 2=LH, 3=LL, 4=HL)
            # lisa.p_sim: pseudo p-values from permutation
            significant_mask = lisa.p_sim < alpha

            # Count fractions by quadrant
            hh_fraction = np.mean((lisa.q == 1) & significant_mask)
            ll_fraction = np.mean((lisa.q == 3) & significant_mask)
            any_fraction = np.mean(significant_mask)

            results.append({
                'significant_HH_fraction': float(hh_fraction),
                'significant_LL_fraction': float(ll_fraction),
                'significant_any_fraction': float(any_fraction),
            })
        except Exception as e:
            logger.debug(f"LISA failed for marker {m}: {e}")
            results.append({'significant_HH_fraction': 0.0,
                           'significant_LL_fraction': 0.0,
                           'significant_any_fraction': 0.0})

    return results


def _compute_getis_ord_batch(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int,
    alpha: float = 0.05,
    permutations: int = 99,
) -> List[Dict[str, float]]:
    """
    Compute Getis-Ord Gi* for each marker.

    Gi* identifies local hotspots (high values clustered together) and
    coldspots (low values clustered together). Unlike LISA, Gi* includes
    the focal point in the computation (star variant).

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k : int
        Number of nearest neighbors for spatial weights.
    alpha : float, default=0.05
        Significance threshold for hotspots.
    permutations : int, default=99
        Number of permutations for pseudo p-values.

    Returns
    -------
    List of dicts, one per marker, containing:
        - 'significant_hotspot_fraction': Proportion of spots that are significant hotspots
        - 'significant_coldspot_fraction': Proportion of significant coldspots
        - 'significant_any_fraction': Proportion in any significant cluster
    """
    if not HAS_LOCAL_SPATIAL:
        logger.warning(
            "esda/libpysal not available; Getis-Ord Gi* computation skipped."
        )
        return [{'significant_hotspot_fraction': 0.0, 'significant_coldspot_fraction': 0.0,
                 'significant_any_fraction': 0.0}] * Z.shape[1]

    n_spots, n_markers = Z.shape

    if n_spots < k + 1:
        return [{'significant_hotspot_fraction': 0.0, 'significant_coldspot_fraction': 0.0,
                 'significant_any_fraction': 0.0}] * n_markers

    try:
        w = LibPySAL_KNN.from_array(coords, k=k)
        w.transform = 'B'  # Binary weights for Gi*
    except Exception as e:
        logger.warning(f"Failed to build spatial weights: {e}")
        return [{'significant_hotspot_fraction': 0.0, 'significant_coldspot_fraction': 0.0,
                 'significant_any_fraction': 0.0}] * n_markers

    results: List[Dict[str, float]] = []
    for m in range(n_markers):
        values = Z[:, m]

        if np.allclose(values, values[0]):
            results.append({'significant_hotspot_fraction': 0.0,
                           'significant_coldspot_fraction': 0.0,
                           'significant_any_fraction': 0.0})
            continue

        try:
            # Use Gi* (star=True includes focal point)
            gi = G_Local(values, w, star=True, permutations=permutations)
            # gi.Zs: standardized Gi* statistics
            # gi.p_sim: pseudo p-values (two-sided)

            significant_mask = gi.p_sim < alpha

            # Hotspots: significant AND positive z-score
            hotspot_fraction = np.mean(significant_mask & (gi.Zs > 0))
            # Coldspots: significant AND negative z-score
            coldspot_fraction = np.mean(significant_mask & (gi.Zs < 0))
            any_fraction = np.mean(significant_mask)

            results.append({
                'significant_hotspot_fraction': float(hotspot_fraction),
                'significant_coldspot_fraction': float(coldspot_fraction),
                'significant_any_fraction': float(any_fraction),
            })
        except Exception as e:
            logger.debug(f"Getis-Ord Gi* failed for marker {m}: {e}")
            results.append({'significant_hotspot_fraction': 0.0,
                           'significant_coldspot_fraction': 0.0,
                           'significant_any_fraction': 0.0})

    return results


def _compute_local_spatial_multiscale(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k_scales: List[int],
    alpha: float = 0.05,
    permutations: int = 99,
    method: str = "both",
) -> Dict[int, List[Dict[str, float]]]:
    """
    Compute local spatial statistics at multiple scales.

    Enables detection of patterns at different spatial granularities:
    - Small k (3-5): Fine-grained clusters (TILs, scattered cells)
    - Medium k (8): Standard aggregates (T-cell patches)
    - Large k (12+): Broad domains (cancer nests, stroma)

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k_scales : list of int
        List of k values to test, e.g., [3, 5, 8, 12].
    alpha : float, default=0.05
        Significance threshold for local clusters.
    permutations : int, default=99
        Number of permutations for pseudo p-values.
    method : str, default="both"
        Which statistics to compute: "lisa", "getis", or "both".

    Returns
    -------
    Dict mapping k -> list of result dicts per marker.
    Each dict contains keys from LISA and/or Getis-Ord depending on method.
    """
    n_markers = Z.shape[1]
    results: Dict[int, List[Dict[str, float]]] = {}

    for k in k_scales:
        scale_results: List[Dict[str, float]] = []

        # Compute LISA
        if method in ("lisa", "both"):
            lisa_results = _compute_lisa_batch(Z, coords, k, alpha, permutations)
        else:
            lisa_results = [{}] * n_markers

        # Compute Getis-Ord
        if method in ("getis", "both"):
            getis_results = _compute_getis_ord_batch(Z, coords, k, alpha, permutations)
        else:
            getis_results = [{}] * n_markers

        # Merge results for each marker
        for m in range(n_markers):
            merged = {}
            if method in ("lisa", "both"):
                merged.update(lisa_results[m])
            if method in ("getis", "both"):
                merged.update(getis_results[m])
            scale_results.append(merged)

        results[k] = scale_results

    return results


# =============================================================================
# ABUNDANCE-ADAPTIVE CLASSIFICATION
# =============================================================================


def _classify_marker_abundance(
    marker_values: NDArray[np.floating],
    ubiquitous_cv_threshold: float = 0.5,
    ubiquitous_presence_threshold: float = 0.7,
    rare_presence_threshold: float = 0.10,
    rare_intensity_factor: float = 2.0,
) -> str:
    """
    Classify a marker by its abundance pattern.

    This enables detection of cell types that fail spatial autocorrelation tests:
    - Ubiquitous types (Epithelial, CAFs): present everywhere, low Moran's I
    - Rare types (Plasmablasts, PVL): too few spots for spatial power

    Parameters
    ----------
    marker_values : ndarray (n_spots,)
        Expression values for a single marker (z-scored).
    ubiquitous_cv_threshold : float, default=0.5
        Maximum coefficient of variation to be classified as ubiquitous.
        Lower CV means more uniform expression.
    ubiquitous_presence_threshold : float, default=0.7
        Minimum fraction of spots with expression above median to be ubiquitous.
    rare_presence_threshold : float, default=0.10
        Maximum fraction of spots with high expression to be classified as rare.
    rare_intensity_factor : float, default=2.0
        Factor by which max expression must exceed mean for rare classification.

    Returns
    -------
    str
        One of: "ubiquitous", "rare", or "standard"
    """
    n_spots = len(marker_values)

    # Basic statistics
    mean_expr = marker_values.mean()
    std_expr = marker_values.std()
    max_expr = marker_values.max()

    # Coefficient of variation (handle near-zero mean)
    cv = std_expr / (abs(mean_expr) + 1e-6)

    # Presence fraction: spots above median (for z-scored data, median ≈ 0)
    median_expr = np.median(marker_values)
    presence_frac = (marker_values > median_expr).sum() / n_spots

    # High expression fraction: spots > 1 std above mean
    high_threshold = mean_expr + std_expr
    high_frac = (marker_values > high_threshold).sum() / n_spots

    # Ubiquitous classification:
    # - Low coefficient of variation (uniform expression)
    # - High presence fraction (expressed almost everywhere)
    if cv < ubiquitous_cv_threshold and presence_frac > ubiquitous_presence_threshold:
        return "ubiquitous"

    # Rare classification:
    # - Very low fraction with high expression (truly sparse)
    # - When present, expression is significantly above background
    # - High kurtosis (peaked distribution with heavy tails, not noise)
    if high_frac < rare_presence_threshold:
        # For rare markers, require TRUE intensity signal:
        # 1. Max expression well above the distribution (not just 2x mean)
        # 2. High kurtosis indicating peaked distribution (real cells, not noise)
        # 3. Some minimum fraction expressing (not just one outlier spot)
        min_expressing_frac = (marker_values > median_expr + 0.5 * std_expr).sum() / n_spots

        # Kurtosis: normal distribution has kurtosis=3, excess kurtosis=0
        # Peaked distributions (rare cell types) have excess kurtosis > 2
        excess_kurtosis = scipy_kurtosis(marker_values, fisher=True)

        # Require: peaked distribution AND meaningful expression fraction AND intensity contrast
        # Upper bound 0.30 allows rare cell types like Plasmablasts that may have ~25% expression
        if (excess_kurtosis > 2.0 and
            min_expressing_frac > 0.02 and min_expressing_frac < 0.30 and
            max_expr > 3.0 * std_expr):
            return "rare"

    return "standard"


def _compute_intensity_metrics(
    Z: NDArray[np.floating],
) -> Dict[str, NDArray[np.floating]]:
    """
    Compute intensity-based metrics for each marker.

    These metrics capture signal strength independent of spatial pattern,
    allowing rescue of markers with high expression but dispersed spatial
    distribution (e.g., scattered plasmablasts).

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.

    Returns
    -------
    Dict containing:
        - 'variance': Variance per marker
        - 'kurtosis': Excess kurtosis per marker (peakedness)
        - 'iqr': Interquartile range per marker (robust spread)
        - 'variance_percentile': Percentile rank of variance (0-1)
        - 'kurtosis_percentile': Percentile rank of kurtosis (0-1)
    """
    n_markers = Z.shape[1]

    # Compute raw metrics
    variance = np.var(Z, axis=0)
    kurt = np.array([scipy_kurtosis(Z[:, m], fisher=True) for m in range(n_markers)])
    iqr = np.array([
        np.percentile(Z[:, m], 75) - np.percentile(Z[:, m], 25)
        for m in range(n_markers)
    ])

    # Compute percentile ranks (useful for adaptive thresholding)
    def percentile_rank(arr: NDArray) -> NDArray:
        """Convert values to their percentile rank (0-1)."""
        # Handle NaN and constant arrays
        valid = ~np.isnan(arr)
        if not np.any(valid) or np.allclose(arr[valid], arr[valid][0]):
            return np.zeros_like(arr)
        ranks = np.zeros_like(arr)
        ranks[valid] = (np.argsort(np.argsort(arr[valid])) + 1) / np.sum(valid)
        return ranks

    variance_percentile = percentile_rank(variance)
    kurtosis_percentile = percentile_rank(kurt)

    return {
        'variance': variance,
        'kurtosis': kurt,
        'iqr': iqr,
        'variance_percentile': variance_percentile,
        'kurtosis_percentile': kurtosis_percentile,
    }


def compute_composite_spatial_score(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k_scales: Optional[List[int]] = None,
    n_perm: int = 99,
    alpha: float = 0.1,
    seed: int = 1234,
) -> NDArray[np.floating]:
    """
    Compute permissive composite spatial score for each marker.

    Uses MAX aggregation: a marker qualifies if ANY spatial metric shows signal.
    This is more permissive than requiring all metrics to pass.

    The composite score combines:
    - Global Moran's I (scaled to 0-1 if > 0 and significant)
    - Local LISA HH fraction × 2 (scaled since typically 0-0.5)
    - Getis-Ord hotspot fraction × 2 (scaled similarly)
    - Intensity score × 0.5 (variance/kurtosis percentile, weaker signal)

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k_scales : list of int, optional
        Neighbor scales to test. Default: [3, 5, 8, 12].
    n_perm : int, default=99
        Number of permutations for significance testing.
    alpha : float, default=0.1
        Significance level.
    seed : int, default=1234
        Random seed for reproducibility.

    Returns
    -------
    ndarray (n_markers,)
        Composite spatial score for each marker in [0, 1].
        Higher = stronger spatial signal.
    """
    if k_scales is None:
        k_scales = [3, 5, 8, 12]

    n_spots, n_markers = Z.shape
    scores = np.zeros(n_markers)

    rng = np.random.default_rng(seed)

    # Compute Moran's I at multiple scales
    morans_results = _compute_morans_i_multiscale(
        Z, coords, k_scales=k_scales, rng=rng, n_perm=n_perm
    )

    # Compute local spatial statistics at multiple scales
    local_results = _compute_local_spatial_multiscale(
        Z, coords, k_scales=k_scales, alpha=alpha, method="both", permutations=n_perm
    )

    # Compute intensity metrics
    intensity_metrics = _compute_intensity_metrics(Z)

    for m in range(n_markers):
        # Component 1: Global Moran's I (best across scales)
        morans_scores = []
        for k in k_scales:
            if k in morans_results and m < len(morans_results[k]):
                I_obs, p_val = morans_results[k][m]
                if I_obs > 0 and p_val < alpha:
                    # Scale Moran's I to [0, 1] (cap at 1.0)
                    morans_scores.append(min(I_obs, 1.0))
        morans_score = max(morans_scores) if morans_scores else 0.0

        # Component 2: Local LISA HH fraction (best across scales)
        lisa_scores = []
        for k in k_scales:
            if k in local_results and m < len(local_results[k]):
                local_dict = local_results[k][m]
                hh_frac = local_dict.get('significant_HH_fraction', 0.0)
                lisa_scores.append(hh_frac)
        lisa_score = max(lisa_scores) if lisa_scores else 0.0

        # Component 3: Getis-Ord hotspot fraction (best across scales)
        getis_scores = []
        for k in k_scales:
            if k in local_results and m < len(local_results[k]):
                local_dict = local_results[k][m]
                hotspot_frac = local_dict.get('significant_hotspot_fraction', 0.0)
                getis_scores.append(hotspot_frac)
        getis_score = max(getis_scores) if getis_scores else 0.0

        # Component 4: Intensity score (variance/kurtosis percentile)
        var_pct = intensity_metrics['variance_percentile'][m]
        kurt_pct = intensity_metrics['kurtosis_percentile'][m]
        intensity_score = max(var_pct, kurt_pct)

        # PERMISSIVE COMBINATION: take MAXIMUM across all components
        # Each component is scaled appropriately:
        # - Moran's I already in [0, 1]
        # - LISA HH fraction typically 0-0.5, so scale by 2
        # - Getis-Ord fraction typically 0-0.5, so scale by 2
        # - Intensity is a weaker signal, so scale by 0.5
        scores[m] = max(
            morans_score,
            lisa_score * 2.0,
            getis_score * 2.0,
            intensity_score * 0.5,
        )

        # Clip to [0, 1]
        scores[m] = min(1.0, scores[m])

    return scores


# =============================================================================
# RECONSTRUCTION-BASED SCORING AND SELECTION
# =============================================================================


def _score_candidate_by_reconstruction(
    X: NDArray[np.floating],
    candidate_markers: Set[int],
    residual_X: Optional[NDArray[np.floating]] = None,
) -> float:
    """
    Score a candidate profile by how much it reduces reconstruction error.

    This is the core of the reconstruction-based approach: a profile is valuable
    if including it reduces the error in explaining the observed antibody data.

    Parameters
    ----------
    X : ndarray (n_spots, n_markers)
        Original expression matrix (z-scored).
    candidate_markers : set of int
        Marker indices in this candidate profile.
    residual_X : ndarray (n_spots, n_markers), optional
        Residual expression after accounting for previously selected profiles.
        If None, uses the original X.

    Returns
    -------
    float
        Relative reconstruction improvement (0-1 scale).
        Higher is better.
    """
    if residual_X is None:
        residual_X = X

    markers = list(candidate_markers)
    n_markers_in_profile = len(markers)

    if n_markers_in_profile == 0:
        return 0.0

    # Extract the relevant columns
    marker_data = residual_X[:, markers]

    # Before error: sum of squared residuals for these markers
    before_error = np.sum(marker_data ** 2)

    if before_error < 1e-9:
        return 0.0  # No error to reduce

    # Profile expression: average across markers in the profile
    # This represents the "presence" of this cell type at each spot
    profile_expr = marker_data.mean(axis=1)

    # Predicted values: profile expression broadcasted to all markers
    # Each marker in the profile gets the same predicted value
    predicted = np.outer(profile_expr, np.ones(n_markers_in_profile))

    # After error: residual after explaining with this profile
    after_error = np.sum((marker_data - predicted) ** 2)

    # Relative improvement
    improvement = (before_error - after_error) / before_error

    return float(max(0.0, improvement))


def _update_residual(
    residual_X: NDArray[np.floating],
    candidate_markers: Set[int],
) -> NDArray[np.floating]:
    """
    Update residual matrix after accounting for a selected profile.

    Parameters
    ----------
    residual_X : ndarray (n_spots, n_markers)
        Current residual expression matrix.
    candidate_markers : set of int
        Marker indices in the selected profile.

    Returns
    -------
    ndarray (n_spots, n_markers)
        Updated residual matrix.
    """
    markers = list(candidate_markers)
    updated = residual_X.copy()

    # Profile expression: average across markers in the profile
    profile_expr = residual_X[:, markers].mean(axis=1)

    # Subtract the explained portion from each marker
    for m in markers:
        updated[:, m] = residual_X[:, m] - profile_expr

    return updated


def _is_profile_redundant(
    candidate: Set[int],
    selected_profiles: List[Set[int]],
    X: NDArray[np.floating],
    correlation_threshold: float = 0.9,
) -> bool:
    """
    Check if a candidate profile is redundant with already selected profiles.

    Redundancy is measured by expression correlation, NOT marker overlap.
    This allows shared markers (hierarchy) while preventing duplicate profiles.

    Parameters
    ----------
    candidate : set of int
        Marker indices in candidate profile.
    selected_profiles : list of sets
        Already selected profiles.
    X : ndarray (n_spots, n_markers)
        Expression matrix.
    correlation_threshold : float, default=0.9
        Correlation above which profiles are considered redundant.

    Returns
    -------
    bool
        True if candidate is redundant with any selected profile.
    """
    if not selected_profiles:
        return False

    # Compute candidate's expression pattern
    cand_markers = list(candidate)
    cand_expr = X[:, cand_markers].mean(axis=1)

    for sel in selected_profiles:
        sel_markers = list(sel)
        sel_expr = X[:, sel_markers].mean(axis=1)

        # Pearson correlation
        corr_matrix = np.corrcoef(cand_expr, sel_expr)
        corr = corr_matrix[0, 1]

        if np.isfinite(corr) and abs(corr) > correlation_threshold:
            return True

    return False


def _detect_hierarchical_structure(
    candidates: List[Set[int]],
    min_shared_markers: int = 1,
) -> Tuple[bool, List[int]]:
    """
    Automatically detect if candidate profiles suggest hierarchical structure.

    Hierarchical structure is indicated when the same marker appears in multiple
    significant candidates with DIFFERENT partners. This is characteristic of
    real biological data where markers like CD3D (pan-T-cell) appear in both
    {CD3D, CD4} and {CD3D, CD8} profiles.

    For flat/simulated data, each marker should appear in only one profile
    (e.g., {B-cells_Protein_1, B-cells_Protein_2}).

    Parameters
    ----------
    candidates : list of sets
        Significant candidate profiles from permutation testing.
    min_shared_markers : int, default=1
        Minimum number of markers that must show sharing pattern to declare
        hierarchical structure.

    Returns
    -------
    is_hierarchical : bool
        True if hierarchical structure detected.
    shared_markers : list of int
        Marker indices that appear in multiple profiles with different partners.

    Examples
    --------
    Flat structure (simulated):
        candidates = [{0, 1}, {2, 3}, {4, 5}]
        -> is_hierarchical=False (no marker appears with multiple partners)

    Hierarchical structure (real T-cells):
        candidates = [{CD3D}, {CD3D, CD4}, {CD3D, CD8}, {CD4}, {CD8}]
        -> is_hierarchical=True (CD3D appears with CD4 and CD8 as partners)
    """
    from collections import defaultdict

    if not candidates:
        return False, []

    # Map each marker to the profiles it appears in
    marker_to_profiles: Dict[int, List[Set[int]]] = defaultdict(list)
    for cand in candidates:
        for marker in cand:
            marker_to_profiles[marker].append(cand)

    shared_markers = []

    for marker, profiles in marker_to_profiles.items():
        if len(profiles) <= 1:
            continue

        # Check if this marker has different partners across profiles
        # (not just the same profile appearing multiple times or subsets)
        unique_partner_sets = set()
        for profile in profiles:
            partners = frozenset(profile - {marker})
            # Only count profiles where marker has partners (not single-marker profiles)
            # OR count single-marker profiles as a distinct "no partner" case
            unique_partner_sets.add(partners)

        # If marker appears with multiple different partner combinations -> hierarchy
        if len(unique_partner_sets) > 1:
            shared_markers.append(marker)

    is_hierarchical = len(shared_markers) >= min_shared_markers

    return is_hierarchical, shared_markers


def _find_shared_markers(candidates: List[Set[int]]) -> List[int]:
    """
    Find markers that appear in multiple candidate profiles.

    Simple helper that extracts just the shared markers list.

    Parameters
    ----------
    candidates : list of sets
        Candidate profiles (each is a set of marker indices).

    Returns
    -------
    shared_markers : list of int
        Marker indices that appear in more than one profile.
    """
    from collections import Counter

    if not candidates:
        return []

    # Count occurrences of each marker
    marker_counts: Counter = Counter()
    for cand in candidates:
        for marker in cand:
            marker_counts[marker] += 1

    # Return markers that appear in multiple profiles
    return [marker for marker, count in marker_counts.items() if count > 1]


def _select_profiles_reconstruction_based(
    X: NDArray[np.floating],
    candidates: List[Set[int]],
    min_improvement: float = 0.05,
    max_profiles: int = 20,
    redundancy_threshold: float = 0.9,
    allow_overlap: bool = False,
) -> List[Set[int]]:
    """
    Select profiles by reconstruction improvement.

    This algorithm can handle both flat and hierarchical structures:
    - Flat (allow_overlap=False): Enforces non-overlapping markers like greedy selection
    - Hierarchical (allow_overlap=True): Uses correlation-based redundancy check

    Parameters
    ----------
    X : ndarray (n_spots, n_markers)
        Expression matrix (z-scored).
    candidates : list of sets
        Candidate profiles (each is a set of marker indices).
    min_improvement : float, default=0.05
        Minimum relative reconstruction improvement to accept a profile.
    max_profiles : int, default=20
        Maximum number of profiles to select.
    redundancy_threshold : float, default=0.9
        Correlation threshold for redundancy check (only used if allow_overlap=True).
    allow_overlap : bool, default=False
        If False, enforce non-overlapping markers (for flat/simulated data).
        If True, allow shared markers and use correlation-based redundancy check
        (for hierarchical/real data).

    Returns
    -------
    List of selected profile sets.
    """
    selected_profiles: List[Set[int]] = []
    used_markers: Set[int] = set()  # Track used markers for non-overlapping mode
    residual_X = X.copy()

    # Sort candidates by initial reconstruction score (descending)
    scored_candidates = [
        (cand, _score_candidate_by_reconstruction(X, cand))
        for cand in candidates
    ]
    scored_candidates.sort(key=lambda x: -x[1])

    for iteration in range(max_profiles):
        best_candidate = None
        best_score = -np.inf

        for cand, initial_score in scored_candidates:
            # Check overlap/redundancy based on mode
            if allow_overlap:
                # Hierarchical mode: use correlation-based redundancy
                if _is_profile_redundant(cand, selected_profiles, X, redundancy_threshold):
                    continue
            else:
                # Flat mode: enforce non-overlapping markers
                if cand & used_markers:
                    continue

            # Score on current residual (not original)
            current_score = _score_candidate_by_reconstruction(X, cand, residual_X)

            if current_score > best_score:
                best_score = current_score
                best_candidate = cand

        if best_candidate is None:
            break

        if best_score < min_improvement:
            logger.debug(
                f"Stopping reconstruction selection: best score {best_score:.4f} "
                f"< threshold {min_improvement}"
            )
            break

        selected_profiles.append(best_candidate)
        used_markers.update(best_candidate)  # Track used markers
        residual_X = _update_residual(residual_X, best_candidate)

        logger.debug(
            f"Selected profile {best_candidate} with reconstruction "
            f"improvement {best_score:.4f}"
        )

    return selected_profiles


# =============================================================================
# MARKER SIGNIFICANCE TESTING
# =============================================================================


def _identify_significant_single_markers(
    Z: NDArray[np.floating],
    beta: NDArray[np.floating],
    rng: "Generator",
    n_perm: int,
    alpha: float,
    coords: Optional[NDArray[np.floating]] = None,
    morans_k: Union[int, List[int]] = 8,
    morans_i_threshold: float = 0.1,
    per_marker_fallback: bool = True,
    # NEW: Local spatial statistics parameters
    use_local_spatial: bool = True,
    local_method: str = "both",
    hotspot_threshold: float = 0.10,
    local_permutations: int = 99,
    # NEW: Intensity rescue parameters
    use_intensity_rescue: bool = True,
    intensity_percentile: float = 0.75,
    min_hotspot_fraction: float = 0.03,
    # NEW: Abundance-adaptive parameters
    use_abundance_adaptive: bool = True,
    ubiquitous_cv_threshold: float = 0.5,
    ubiquitous_presence_threshold: float = 0.7,
    rare_presence_threshold: float = 0.10,
    rare_intensity_factor: float = 2.0,
) -> List[int]:
    """
    Return markers with significant signal using abundance-adaptive detection.

    Detection Hierarchy:
    0. **Tier 0 (Abundance-Adaptive)**: Classify by abundance pattern (ubiquitous/rare)
       - Ubiquitous markers pass without spatial test (e.g., Epithelial everywhere)
       - Rare markers pass with intensity test (e.g., scattered Plasmablasts)
    1. **Tier 1 (Global Spatial)**: Global Moran's I >= threshold at ANY scale
    2. **Tier 2 (Local Spatial)**: >hotspot_threshold of spots in significant
       local clusters (LISA HH or Getis-Ord hotspots)
    3. **Tier 3 (Intensity Rescue)**: High variance/kurtosis + minimal local structure

    A marker passes if it passes ANY tier.

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    beta : ndarray
        Marker weights for scoring.
    rng : Generator
        Random number generator.
    n_perm : int
        Number of permutations for global Moran's I testing.
    alpha : float
        Significance threshold (p-value cutoff).
    coords : ndarray (n_spots, 2), optional
        Spatial coordinates. Required for spatial testing.
    morans_k : int or list of int, default=8
        Single k value or list of scales for Moran's I.
        For mixed/dysplastic tissues, use [3, 5, 8, 12].
    morans_i_threshold : float, default=0.1
        Minimum Moran's I value for Tier 1. Higher threshold rejects more noise.
        For mixed tissues with weak spatial signals, set to 0.05 or lower.
    per_marker_fallback : bool, default=True
        If True, markers failing all spatial tiers get expression testing.
    use_local_spatial : bool, default=True
        If True, enable Tier 2 (LISA/Getis-Ord local statistics).
    local_method : str, default="both"
        Which local statistic: "lisa", "getis", or "both".
    hotspot_threshold : float, default=0.10
        Minimum fraction of spots in significant hotspots for Tier 2.
    local_permutations : int, default=99
        Number of permutations for local statistics p-values.
    use_intensity_rescue : bool, default=True
        If True, enable Tier 3 (intensity-based rescue).
    intensity_percentile : float, default=0.75
        Minimum variance/kurtosis percentile for Tier 3.
    min_hotspot_fraction : float, default=0.03
        Minimum hotspot fraction for Tier 3 (relaxed threshold).
    use_abundance_adaptive : bool, default=True
        If True, enable Tier 0 (abundance-adaptive classification).
    ubiquitous_cv_threshold : float, default=0.5
        Maximum CV for ubiquitous classification.
    ubiquitous_presence_threshold : float, default=0.7
        Minimum presence fraction for ubiquitous classification.
    rare_presence_threshold : float, default=0.10
        Maximum high-expression fraction for rare classification.
    rare_intensity_factor : float, default=2.0
        Max/mean ratio required for rare classification.

    Returns
    -------
    List of marker indices that passed significance testing.

    Notes
    -----
    Tier 0 (abundance-adaptive) enables detection of cell types that fail
    spatial autocorrelation tests:
    - Ubiquitous types (Epithelial, CAFs): Present everywhere, Moran's I ≈ 0
    - Rare types (Plasmablasts, PVL): Too few spots for spatial power
    """
    n_markers = Z.shape[1]
    significant: List[int] = []
    markers_passed_tier0_ubiq: Set[int] = set()
    markers_passed_tier0_rare: Set[int] = set()
    markers_passed_tier1: Set[int] = set()
    markers_passed_tier2: Set[int] = set()
    markers_passed_tier3: Set[int] = set()

    # Normalize k_scales to list
    if isinstance(morans_k, int):
        k_scales = [morans_k]
    else:
        k_scales = list(morans_k)

    # =========================================================================
    # TIER 0: Abundance-Adaptive Classification
    # =========================================================================
    if use_abundance_adaptive:
        for m in range(n_markers):
            abundance_class = _classify_marker_abundance(
                Z[:, m],
                ubiquitous_cv_threshold=ubiquitous_cv_threshold,
                ubiquitous_presence_threshold=ubiquitous_presence_threshold,
                rare_presence_threshold=rare_presence_threshold,
                rare_intensity_factor=rare_intensity_factor,
            )

            if abundance_class == "ubiquitous":
                significant.append(m)
                markers_passed_tier0_ubiq.add(m)
                logger.debug(
                    f"Marker {m} passed Tier 0 (Ubiquitous): present everywhere, "
                    f"no spatial structure required"
                )
            elif abundance_class == "rare":
                # For rare markers, check intensity (already done in classification)
                significant.append(m)
                markers_passed_tier0_rare.add(m)
                logger.debug(
                    f"Marker {m} passed Tier 0 (Rare): sparse but intense expression"
                )

        logger.info(
            f"Tier 0 (Abundance-Adaptive): {len(markers_passed_tier0_ubiq)} ubiquitous, "
            f"{len(markers_passed_tier0_rare)} rare markers passed"
        )

    # =========================================================================
    # TIER 1: Global Spatial (Moran's I)
    # =========================================================================
    # Skip markers that already passed Tier 0
    markers_passed_tier0 = markers_passed_tier0_ubiq | markers_passed_tier0_rare
    markers_for_tier1 = set(range(n_markers)) - markers_passed_tier0

    if coords is not None and markers_for_tier1:
        multiscale_results = _compute_morans_i_multiscale(
            Z, coords, k_scales, rng, n_perm
        )

        for m in markers_for_tier1:
            passed = False
            best_I = -np.inf
            best_k = None

            for k in k_scales:
                I_obs, p_value = multiscale_results[k][m]

                if np.isnan(I_obs):
                    continue

                if I_obs > best_I:
                    best_I = I_obs
                    best_k = k

                if I_obs >= morans_i_threshold and p_value < alpha:
                    passed = True
                    break

            if passed:
                significant.append(m)
                markers_passed_tier1.add(m)
                logger.debug(
                    f"Marker {m} passed Tier 1 (Global Moran's I) at k={best_k}: "
                    f"I={best_I:.3f}, p<{alpha}"
                )

        logger.info(
            f"Tier 1 (Global Moran's I): {len(markers_passed_tier1)}/{len(markers_for_tier1)} "
            f"remaining markers passed"
        )

    # =========================================================================
    # TIER 2: Local Spatial (LISA / Getis-Ord Gi*)
    # =========================================================================
    if use_local_spatial and coords is not None and HAS_LOCAL_SPATIAL:
        markers_for_tier2 = set(range(n_markers)) - markers_passed_tier0 - markers_passed_tier1

        if markers_for_tier2:
            # Compute local statistics at all scales
            local_results = _compute_local_spatial_multiscale(
                Z, coords, k_scales,
                alpha=alpha,
                permutations=local_permutations,
                method=local_method,
            )

            for m in markers_for_tier2:
                passed = False
                best_fraction = 0.0
                best_k = None

                for k in k_scales:
                    result = local_results[k][m]

                    # Check LISA HH fraction
                    if local_method in ("lisa", "both"):
                        hh_frac = result.get('significant_HH_fraction', 0.0)
                        if hh_frac > best_fraction:
                            best_fraction = hh_frac
                            best_k = k
                        if hh_frac >= hotspot_threshold:
                            passed = True

                    # Check Getis-Ord hotspot fraction
                    if local_method in ("getis", "both"):
                        hotspot_frac = result.get('significant_hotspot_fraction', 0.0)
                        if hotspot_frac > best_fraction:
                            best_fraction = hotspot_frac
                            best_k = k
                        if hotspot_frac >= hotspot_threshold:
                            passed = True

                    if passed:
                        break  # No need to check other scales

                if passed:
                    significant.append(m)
                    markers_passed_tier2.add(m)
                    logger.debug(
                        f"Marker {m} passed Tier 2 (Local Spatial) at k={best_k}: "
                        f"hotspot_fraction={best_fraction:.3f} >= {hotspot_threshold}"
                    )

            logger.info(
                f"Tier 2 (Local Spatial): {len(markers_passed_tier2)}/{len(markers_for_tier2)} "
                f"remaining markers passed"
            )

    elif use_local_spatial and not HAS_LOCAL_SPATIAL:
        logger.warning(
            "Tier 2 (Local Spatial) skipped: esda/libpysal not available. "
            "Install with: pip install esda libpysal"
        )

    # =========================================================================
    # TIER 3: Intensity Rescue
    # =========================================================================
    if use_intensity_rescue and coords is not None:
        markers_for_tier3 = (
            set(range(n_markers)) - markers_passed_tier0 - markers_passed_tier1 - markers_passed_tier2
        )

        if markers_for_tier3:
            # Compute intensity metrics
            intensity_metrics = _compute_intensity_metrics(Z)
            variance_pct = intensity_metrics['variance_percentile']
            kurtosis_pct = intensity_metrics['kurtosis_percentile']

            # Get best local hotspot fraction for each marker (from Tier 2 results)
            # If Tier 2 wasn't run, compute it now with relaxed threshold
            if use_local_spatial and HAS_LOCAL_SPATIAL:
                # Reuse local_results from Tier 2
                pass
            else:
                # Need to compute local stats for Tier 3
                local_results = _compute_local_spatial_multiscale(
                    Z, coords, k_scales,
                    alpha=alpha * 2,  # More lenient for Tier 3
                    permutations=local_permutations,
                    method="both",
                )

            for m in markers_for_tier3:
                # Check intensity requirement
                intensity_pass = (
                    variance_pct[m] >= intensity_percentile or
                    kurtosis_pct[m] >= intensity_percentile
                )

                if not intensity_pass:
                    continue

                # Check minimal local structure (relaxed threshold)
                best_fraction = 0.0
                for k in k_scales:
                    result = local_results[k][m]
                    hh_frac = result.get('significant_HH_fraction', 0.0)
                    hotspot_frac = result.get('significant_hotspot_fraction', 0.0)
                    best_fraction = max(best_fraction, hh_frac, hotspot_frac)

                if best_fraction >= min_hotspot_fraction:
                    significant.append(m)
                    markers_passed_tier3.add(m)
                    logger.debug(
                        f"Marker {m} passed Tier 3 (Intensity Rescue): "
                        f"variance_pct={variance_pct[m]:.2f}, "
                        f"kurtosis_pct={kurtosis_pct[m]:.2f}, "
                        f"hotspot_fraction={best_fraction:.3f}"
                    )

            logger.info(
                f"Tier 3 (Intensity Rescue): {len(markers_passed_tier3)}/{len(markers_for_tier3)} "
                f"remaining markers passed"
            )

    # =========================================================================
    # FALLBACK: Expression-based permutation test
    # =========================================================================
    all_passed = markers_passed_tier0 | markers_passed_tier1 | markers_passed_tier2 | markers_passed_tier3

    if per_marker_fallback:
        markers_to_test = set(range(n_markers)) - all_passed
    else:
        if all_passed:
            return significant
        markers_to_test = set(range(n_markers))
        logger.debug(
            "No markers passed any tier; "
            "falling back to expression permutation test for all markers."
        )

    fallback_passed = 0
    for m in markers_to_test:
        score = _score_candidate(Z, {m}, beta)
        null_scores = _compute_null_distribution(Z, {m}, beta, rng, n_perm)
        p_value = (1 + np.sum(null_scores >= score)) / (1 + n_perm)

        if p_value < alpha:
            significant.append(m)
            fallback_passed += 1
            logger.debug(
                f"Marker {m} passed fallback (Expression test): "
                f"score={score:.3f}, p={p_value:.4f}"
            )

    if fallback_passed > 0:
        logger.info(
            f"Fallback (Expression test): {fallback_passed}/{len(markers_to_test)} "
            f"remaining markers passed"
        )

    logger.info(
        f"Total eligible markers: {len(significant)}/{n_markers} "
        f"(Tier0_Ubiq={len(markers_passed_tier0_ubiq)}, Tier0_Rare={len(markers_passed_tier0_rare)}, "
        f"Tier1={len(markers_passed_tier1)}, Tier2={len(markers_passed_tier2)}, "
        f"Tier3={len(markers_passed_tier3)}, Fallback={fallback_passed})"
    )

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


def _spatial_cv_likelihood(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    profiles: List[Set[int]],
    beta: NDArray[np.floating],
    n_folds: int = 5,
    seed: int = 1234,
) -> float:
    """
    Compute cross-validated log-likelihood using spatial stratification.

    Spatial stratification ensures train/test spots are spatially separated
    to avoid information leakage from spatial autocorrelation.

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    profiles : list of sets
        List of profile marker index sets.
    beta : ndarray
        Marker specificity weights.
    n_folds : int, default=5
        Number of cross-validation folds.
    seed : int, default=1234
        Random seed for k-means stratification.

    Returns
    -------
    float
        Mean held-out log-likelihood across folds.
    """
    from sklearn.cluster import KMeans

    n_spots, n_markers = Z.shape
    K = len(profiles)

    if K == 0:
        return -np.inf

    # Handle edge case: too few spots for CV
    if n_spots < n_folds * 2:
        logger.debug(f"Too few spots ({n_spots}) for {n_folds}-fold CV, using 2-fold")
        n_folds = max(2, n_spots // 10)
        if n_folds < 2:
            # Can't do CV, return EM likelihood instead
            _, _, ll = _em_refine(Z, profiles, beta.copy())
            return ll

    # Spatial stratification: cluster spots into n_folds groups using k-means
    kmeans = KMeans(n_clusters=n_folds, random_state=seed, n_init=10)
    fold_assignment = kmeans.fit_predict(coords)

    cv_log_liks = []
    A = _build_profile_matrix(profiles, n_markers)

    for fold in range(n_folds):
        test_mask = fold_assignment == fold
        train_mask = ~test_mask

        # Need at least 10 spots in each set
        if train_mask.sum() < 10 or test_mask.sum() < 5:
            continue

        # Fit EM on training spots
        Z_train = Z[train_mask]
        Y_train, beta_cv, _ = _em_refine(Z_train, profiles, beta.copy())

        # Evaluate on test spots
        Z_test = Z[test_mask]
        n_test = test_mask.sum()

        # E-step on test data to get Y_test
        Y_test = np.full((n_test, K), 1.0 / K)
        for k, profile in enumerate(profiles):
            markers = list(profile)
            diff = Z_test[:, markers] - beta_cv[markers]
            Y_test[:, k] = np.exp(-0.5 * np.sum(diff**2, axis=1))

        # Normalize to proportions
        row_sums = Y_test.sum(axis=1, keepdims=True)
        Y_test = Y_test / np.maximum(row_sums, 1e-9)

        # Compute test log-likelihood
        expected = Y_test @ A @ np.diag(beta_cv)
        test_ll = -0.5 * np.sum((Z_test - expected) ** 2)

        # Normalize by number of test spots for comparability
        cv_log_liks.append(test_ll / n_test)

    if not cv_log_liks:
        # Couldn't do any CV, fall back to EM likelihood
        _, _, ll = _em_refine(Z, profiles, beta.copy())
        return ll

    return float(np.mean(cv_log_liks))


def _discover_with_cv_selection(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    candidates: List[Tuple[Set[int], float, float]],
    beta: NDArray[np.floating],
    max_profiles: int = 20,
    min_profiles: int = 2,
    cv_folds: int = 5,
    seed: int = 1234,
) -> List[Set[int]]:
    """
    Iteratively add profiles, keeping those that improve CV score.

    Uses spatial cross-validation instead of BIC for model selection.
    This prevents premature stopping that occurs with BIC when the
    penalty term overwhelms likelihood improvements.

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    candidates : list of (set, score, p_value) tuples
        Candidate profiles sorted by score descending.
    beta : ndarray
        Initial marker weights.
    max_profiles : int, default=20
        Maximum number of profiles to discover.
    min_profiles : int, default=2
        Minimum profiles before CV stopping applies.
    cv_folds : int, default=5
        Number of cross-validation folds.
    seed : int, default=1234
        Random seed for CV stratification.

    Returns
    -------
    List of selected profile sets.
    """
    profiles: List[Set[int]] = []
    used_markers: Set[int] = set()
    prev_cv = -np.inf

    for candidate, score, p_value in candidates:
        if len(profiles) >= max_profiles:
            break

        # Check if this candidate overlaps with already selected
        # Allow overlapping until min_profiles is reached (for shared markers)
        if len(profiles) >= min_profiles and (candidate & used_markers):
            continue

        # Tentatively add the candidate
        profiles.append(candidate)
        used_markers.update(candidate)

        # Compute CV score with the new profile
        cv_score = _spatial_cv_likelihood(
            Z, coords, profiles, beta.copy(), n_folds=cv_folds, seed=seed
        )

        # Apply CV stopping only after min_profiles
        if len(profiles) > min_profiles:
            if cv_score <= prev_cv:
                # This profile doesn't improve held-out likelihood
                profiles.pop()
                used_markers -= candidate
                logger.debug(
                    f"Profile {candidate} rejected by CV: "
                    f"score={cv_score:.2f} <= prev={prev_cv:.2f}"
                )
                # Continue checking other candidates (might find one that helps)
                continue

        prev_cv = cv_score
        logger.debug(
            f"Profile {candidate} accepted: CV score={cv_score:.2f}"
        )

    return profiles


def is_parent_marker(
    marker_idx: int,
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    eligible_markers: Set[int],
    morans_i_threshold: float = 0.15,
    min_coexpression_partners: int = 2,
    coexpression_threshold: float = 0.3,
) -> bool:
    """
    Detect markers that could be parent markers for hierarchical profiles.

    A parent marker (like CD3D for T-cells, CD68 for myeloid cells) should have:
    1. Strong spatial clustering (Moran's I > threshold)
    2. Multiple co-expression partners (potential subtypes)

    Parameters
    ----------
    marker_idx : int
        Index of the marker to test.
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    eligible_markers : set of int
        Indices of markers that passed spatial filtering.
    morans_i_threshold : float, default=0.15
        Minimum Moran's I for strong spatial clustering.
    min_coexpression_partners : int, default=2
        Minimum number of co-expression partners to be a parent.
    coexpression_threshold : float, default=0.3
        Minimum Pearson correlation to count as co-expression.

    Returns
    -------
    bool
        True if marker meets parent marker criteria.
    """
    if marker_idx not in eligible_markers:
        return False

    n_spots = Z.shape[0]

    # Compute Moran's I for this marker
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=9)  # k=8 neighbors + self
    k_neighbors = indices[:, 1:]

    marker_values = Z[:, marker_idx]
    marker_mean = marker_values.mean()
    marker_centered = marker_values - marker_mean
    variance = np.sum(marker_centered**2)

    if variance < 1e-9:
        return False

    lag_sum = 0.0
    for i in range(n_spots):
        for j in k_neighbors[i]:
            lag_sum += marker_centered[i] * marker_centered[j]

    n_edges = n_spots * 8
    morans_i = (n_spots / n_edges) * (lag_sum / variance)

    if morans_i < morans_i_threshold:
        return False

    # Count co-expression partners
    n_partners = 0
    for other_idx in eligible_markers:
        if other_idx == marker_idx:
            continue
        # Compute Pearson correlation in high-expression spots
        high_expr_mask = Z[:, marker_idx] > np.percentile(Z[:, marker_idx], 50)
        if high_expr_mask.sum() < 10:
            continue
        corr = np.corrcoef(Z[high_expr_mask, marker_idx],
                          Z[high_expr_mask, other_idx])[0, 1]
        if np.isfinite(corr) and abs(corr) > coexpression_threshold:
            n_partners += 1

    return n_partners >= min_coexpression_partners


def _discover_hierarchical_profiles(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    eligible_markers: Set[int],
    beta: NDArray[np.floating],
    rng: "Generator",
    max_k: int = 3,
    n_perm: int = 500,
    alpha: float = 0.1,
    max_profiles: int = 20,
    max_parent_profiles: int = 5,
    max_subtypes_per_parent: int = 4,
) -> List[Set[int]]:
    """
    Discover profiles hierarchically, allowing shared parent markers.

    Phase 1: Discover parent profiles (single markers with strong spatial signal
             and multiple co-expression partners).
    Phase 2: Within parent-expressing spots, discover subtypes (child markers
             with spatial structure and mutual exclusivity with siblings).
    Phase 3: For remaining markers, use non-hierarchical discovery.

    Parameters
    ----------
    Z : ndarray (n_spots, n_markers)
        Standardized expression matrix.
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    eligible_markers : set of int
        Indices of markers that passed spatial filtering.
    beta : ndarray
        Marker weights.
    rng : Generator
        Random number generator.
    max_k : int, default=3
        Maximum markers per profile.
    n_perm : int, default=500
        Permutations for significance testing.
    alpha : float, default=0.1
        Significance threshold.
    max_profiles : int, default=20
        Maximum total profiles.
    max_parent_profiles : int, default=5
        Maximum number of parent profiles.
    max_subtypes_per_parent : int, default=4
        Maximum subtypes per parent.

    Returns
    -------
    List of profile sets (may contain shared markers).
    """
    # Convert eligible_markers to set for set operations
    eligible_set = set(eligible_markers) if not isinstance(eligible_markers, set) else eligible_markers

    profiles: List[Set[int]] = []
    used_as_child: Set[int] = set()

    # Phase 1: Identify potential parent markers
    parent_markers: List[int] = []
    for marker_idx in eligible_set:
        if is_parent_marker(marker_idx, Z, coords, eligible_set):
            parent_markers.append(marker_idx)

    parent_markers = parent_markers[:max_parent_profiles]
    logger.debug(f"Identified {len(parent_markers)} potential parent markers")

    # Phase 2: For each parent, find subtypes
    for parent_idx in parent_markers:
        if len(profiles) >= max_profiles:
            break

        # Get spots where parent is highly expressed
        parent_values = Z[:, parent_idx]
        parent_threshold = np.percentile(parent_values, 50)
        parent_spots = parent_values > parent_threshold

        if parent_spots.sum() < 20:
            # Not enough spots, just add parent as single profile
            profiles.append({parent_idx})
            logger.debug(f"Parent marker {parent_idx} added as single profile")
            continue

        # Find child markers with spatial structure WITHIN parent spots
        Z_parent = Z[parent_spots, :]
        coords_parent = coords[parent_spots, :]

        child_candidates: List[Tuple[int, float]] = []
        for other_idx in eligible_set:
            if other_idx == parent_idx:
                continue
            if other_idx in used_as_child:
                continue

            # Check spatial structure within parent spots
            child_values = Z_parent[:, other_idx]
            if child_values.std() < 0.1:
                continue

            # Compute Moran's I within parent spots
            n_parent_spots = len(coords_parent)
            if n_parent_spots < 20:
                continue

            tree = cKDTree(coords_parent)
            k_neighbors = min(8, n_parent_spots - 1)
            if k_neighbors < 3:
                continue
            _, indices = tree.query(coords_parent, k=k_neighbors + 1)
            indices = indices[:, 1:]

            child_mean = child_values.mean()
            child_centered = child_values - child_mean
            variance = np.sum(child_centered**2)

            if variance < 1e-9:
                continue

            lag_sum = 0.0
            for i in range(n_parent_spots):
                for j in indices[i]:
                    if j < n_parent_spots:
                        lag_sum += child_centered[i] * child_centered[j]

            n_edges = n_parent_spots * k_neighbors
            child_morans_i = (n_parent_spots / n_edges) * (lag_sum / variance)

            if child_morans_i > 0.1:
                # Child has spatial structure within parent region
                signal_strength = child_values[child_values > 0].mean() if (child_values > 0).any() else 0
                child_candidates.append((other_idx, child_morans_i * signal_strength))

        # Sort by signal strength
        child_candidates.sort(key=lambda x: -x[1])

        # Select subtypes (checking for mutual exclusivity)
        subtypes_found: List[Set[int]] = []
        for child_idx, score in child_candidates[:max_subtypes_per_parent]:
            if len(profiles) >= max_profiles:
                break

            # Check mutual exclusivity with existing subtypes
            child_high = Z_parent[:, child_idx] > np.percentile(Z_parent[:, child_idx], 70)
            is_exclusive = True

            for existing_subtype in subtypes_found:
                for existing_child in existing_subtype:
                    if existing_child == parent_idx:
                        continue
                    existing_high = Z_parent[:, existing_child] > np.percentile(
                        Z_parent[:, existing_child], 70
                    )
                    overlap = (child_high & existing_high).sum()
                    min_count = min(child_high.sum(), existing_high.sum())
                    if min_count > 0 and overlap / min_count > 0.2:
                        is_exclusive = False
                        break
                if not is_exclusive:
                    break

            if is_exclusive:
                subtype_profile = {parent_idx, child_idx}
                subtypes_found.append(subtype_profile)
                profiles.append(subtype_profile)
                used_as_child.add(child_idx)
                logger.debug(
                    f"Subtype profile {subtype_profile} added "
                    f"(parent={parent_idx}, child={child_idx})"
                )

        # If no subtypes found, add parent as single profile
        if not subtypes_found:
            profiles.append({parent_idx})
            logger.debug(f"No subtypes found for parent {parent_idx}, added as single")

    # Phase 3: Non-hierarchical discovery for remaining markers
    remaining_markers = eligible_set - set(parent_markers) - used_as_child
    if remaining_markers and len(profiles) < max_profiles:
        # Score remaining candidates
        remaining_candidates = _score_all_candidates(
            Z, remaining_markers, beta, max_k, rng, n_perm=n_perm, alpha=alpha
        )

        # Greedy non-overlapping selection for remaining
        used_in_remaining: Set[int] = set()
        for candidate, score, p_value in remaining_candidates:
            if len(profiles) >= max_profiles:
                break
            if candidate & used_in_remaining:
                continue
            profiles.append(candidate)
            used_in_remaining.update(candidate)
            logger.debug(f"Non-hierarchical profile {candidate} added")

    return profiles


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

    # Set on model for downstream use via the proper method
    # This sets model.cell_profile_dict which is what run_cell_proportion_model() checks
    model.load_cell_profile_dict(result.profiles)
    model.marker_beta = result.beta

    logger.info(
        f"Discovered {len(result.profiles)} profiles: {list(result.profiles.keys())}"
    )

    return result
