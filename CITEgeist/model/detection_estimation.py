# CITEgeist/model/detection_estimation.py
"""Combined detection + estimation pipeline with learned noise variance.

This module combines Stage 1 (GMM detection) and Stage 2 (masked IQP) into
a single pipeline with EM-style iteration to learn per-marker noise variance.

The EM algorithm:
- E-step: Solve IQP for counts, alpha, beta given current weights (1/sigma^2)
- M-step: Update sigma^2 from residuals

This learns which markers are noisy and down-weights them appropriately.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np

from .detection import detect_cell_types
from .masked_iqp import solve_masked_iqp

logger = logging.getLogger(__name__)


def solve_detection_estimation(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    marker_groups: Dict[str, List[int]],
    max_iter: int = 10,
    detection_threshold: float = 0.5,
    convergence_rtol: float = 0.01,
    use_robust_variance: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Two-stage detection + estimation with learned noise variance.

    Args:
        X: (n_spots, n_markers) antibody signal matrix.
        nuclei_counts: (n_spots,) from Cellpose segmentation.
        profile: (n_types, n_markers) binary assignment matrix.
        marker_groups: Dict mapping cell_type -> marker indices.
        max_iter: Maximum EM iterations.
        detection_threshold: Posterior threshold for GMM detection.
        convergence_rtol: Relative tolerance for sigma_sq convergence.
        use_robust_variance: If True, use MAD-based variance (robust to outliers).

    Returns:
        detected: (n_spots, n_types) binary presence mask.
        counts: (n_spots, n_types) integer cell counts.
        alpha: (n_markers,) learned baselines.
        beta: (n_markers,) learned signal-per-cell.
        sigma_sq: (n_markers,) learned noise variances.
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    logger.info("Stage 1: Cell type detection via multivariate GMM")
    detected = detect_cell_types(X, marker_groups, threshold=detection_threshold)

    # Log detection summary
    for k, cell_type in enumerate(marker_groups.keys()):
        n_det = detected[:, k].sum()
        logger.info(f"  {cell_type}: {n_det}/{n_spots} spots detected ({100*n_det/n_spots:.1f}%)")

    # Check for edge case: no types detected anywhere
    if not detected.any():
        logger.warning("No cell types detected in any spot. Returning zeros.")
        return (
            detected,
            np.zeros((n_spots, n_types), dtype=int),
            np.zeros(n_markers),
            np.ones(n_markers),
            np.ones(n_markers),
        )

    logger.info("Stage 2: Global IQP estimation with learned variance")

    # Initialize with uniform weights
    sigma_sq = np.ones(n_markers)

    for iteration in range(max_iter):
        logger.debug(f"EM iteration {iteration + 1}/{max_iter}")

        # E-step: Solve IQP with current weights
        weights = 1.0 / sigma_sq
        counts, alpha, beta = solve_masked_iqp(
            X, nuclei_counts, profile, detected, weights
        )

        # M-step: Update sigma_sq from residuals
        # predicted = alpha + counts @ profile * beta
        predicted = alpha + (counts @ profile) * beta
        residuals = X - predicted

        if use_robust_variance:
            # MAD-based robust variance estimation
            # sigma = 1.4826 * median(|residuals|)
            mad = np.median(np.abs(residuals), axis=0)
            sigma_sq_new = (1.4826 * mad) ** 2
        else:
            # Standard variance
            sigma_sq_new = (residuals ** 2).mean(axis=0)

        # Floor for numerical stability
        sigma_sq_new = np.maximum(sigma_sq_new, 1e-6)

        # Check convergence
        if np.allclose(sigma_sq, sigma_sq_new, rtol=convergence_rtol):
            logger.info(f"Converged at iteration {iteration + 1}")
            break

        sigma_sq = sigma_sq_new

    logger.info(f"Learned sigma: {np.sqrt(sigma_sq)}")
    logger.info(f"Learned alpha (baseline): {alpha}")
    logger.info(f"Learned beta (signal/cell): {beta}")

    return detected, counts, alpha, beta, sigma_sq
