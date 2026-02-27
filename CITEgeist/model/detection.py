# CITEgeist/model/detection.py
"""Cell type detection using multivariate Gaussian Mixture Models.

This module implements Stage 1 of the two-stage detection + estimation model.
For each cell type, a 2-component GMM is fit in the joint space of its markers
to classify spots as signal (cell type present) or background (absent).

The multivariate approach captures marker covariance - e.g., CD4+ T cells
should have BOTH CD3 and CD4 elevated, not just one.
"""
import logging
from typing import Dict, List, Optional

import numpy as np
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


def _compute_adaptive_threshold(
    gmm: GaussianMixture,
    signal_cluster: int,
    base_threshold: float = 0.5,
) -> float:
    """
    Compute adaptive threshold based on GMM cluster weights.

    Rare types (small signal cluster) get lower threshold for sensitivity.
    Common types (large signal cluster) get higher threshold to avoid over-detection.

    The intuition is:
    - If signal cluster weight is small (rare type) → GMM found true signal/background
      split → use lower threshold to catch all signal
    - If signal cluster weight is large (common type) → GMM may be splitting within
      signal (moderate vs high expression) → use higher threshold to avoid over-detection

    Args:
        gmm: Fitted GaussianMixture model.
        signal_cluster: Index of the signal cluster (higher mean).
        base_threshold: Default threshold when signal weight is moderate (0.5-0.7).

    Returns:
        Adaptive threshold value between 0.3 and 0.6.
    """
    w_signal = gmm.weights_[signal_cluster]

    if w_signal < 0.3:
        # Very rare type - GMM found true separation, be sensitive
        return 0.3
    elif w_signal < 0.5:
        # Rare to moderate - slightly lower threshold
        return 0.4
    elif w_signal < 0.7:
        # Common type - use base threshold
        return base_threshold
    else:
        # Very common type - GMM may be splitting within signal, be conservative
        return 0.6


def detect_cell_types(
    X: np.ndarray,
    marker_groups: Dict[str, List[int]],
    threshold: float = 0.5,
    random_state: int = 42,
    log_transform: bool = True,
    adaptive_threshold: bool = True,
) -> np.ndarray:
    """
    Binary detection per cell type using multivariate GMM.

    For each cell type, fits a 2-component GMM (background vs signal) in the
    joint space of its markers. Spots are classified as "detected" if the
    posterior probability of belonging to the signal cluster exceeds threshold.

    Args:
        X: (n_spots, n_markers) antibody signal matrix.
        marker_groups: Dict mapping cell_type_name -> list of marker indices.
            Example: {"CD4+ T cells": [0, 3], "B cells": [5]}
        threshold: Posterior probability threshold for detection (default 0.5).
            If adaptive_threshold=True, this serves as the base threshold.
        random_state: Random seed for GMM initialization.
        log_transform: If True, apply log1p transform before GMM fitting.
            This compresses heavy-tailed distributions and improves separation
            of true signal from background (recommended for antibody data).
        adaptive_threshold: If True, adjust threshold based on GMM cluster weights.
            Rare types (small signal cluster weight) get lower thresholds for
            sensitivity; common types (large signal cluster weight) get higher
            thresholds to avoid over-detection. Recommended for heterogeneous
            cell type prevalences.

    Returns:
        detected: (n_spots, n_types) boolean mask where detected[i,k]=True
            means cell type k is present in spot i.
    """
    n_spots = X.shape[0]
    n_types = len(marker_groups)
    cell_type_names = list(marker_groups.keys())

    detected = np.zeros((n_spots, n_types), dtype=bool)

    for k, cell_type in enumerate(cell_type_names):
        marker_indices = marker_groups[cell_type]

        # Extract markers for this cell type
        marker_data = X[:, marker_indices]  # (n_spots, n_markers_k)

        # Handle edge case: single marker
        if marker_data.ndim == 1:
            marker_data = marker_data.reshape(-1, 1)

        # Log-transform to compress dynamic range (heavy right tail)
        # This helps GMM find true signal/background boundary rather than
        # splitting moderate from high expression
        if log_transform:
            marker_data = np.log1p(marker_data)

        # Fit 2-component GMM
        gmm = GaussianMixture(
            n_components=2,
            covariance_type='full',
            random_state=random_state,
            n_init=3,  # multiple initializations for stability
        )

        try:
            gmm.fit(marker_data)
        except Exception as e:
            logger.warning(f"GMM fit failed for {cell_type}: {e}. Marking all as not detected.")
            continue

        # Identify signal cluster (higher mean across markers)
        cluster_means = gmm.means_.sum(axis=1)
        signal_cluster = int(np.argmax(cluster_means))

        # Compute threshold (adaptive or fixed)
        if adaptive_threshold:
            effective_threshold = _compute_adaptive_threshold(
                gmm, signal_cluster, base_threshold=threshold
            )
        else:
            effective_threshold = threshold

        # Get posterior probability of signal
        posteriors = gmm.predict_proba(marker_data)[:, signal_cluster]

        # Binary detection
        detected[:, k] = posteriors > effective_threshold

        n_detected = detected[:, k].sum()
        w_signal = gmm.weights_[signal_cluster]
        logger.debug(
            f"{cell_type}: {n_detected}/{n_spots} spots detected ({100*n_detected/n_spots:.1f}%), "
            f"signal_weight={w_signal:.2f}, threshold={effective_threshold:.2f}"
        )

    return detected
