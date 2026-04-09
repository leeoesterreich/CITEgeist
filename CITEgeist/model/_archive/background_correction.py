"""
Spatial Mixture Model (SMM) for background correction.

This module implements a Hidden Markov Random Field (HMRF) with Expectation-Maximization
for spatial background correction of protein marker data. The algorithm combines:
- Gaussian Mixture Model (GMM) for intensity-based classification
- Hidden Markov Random Field (Potts model) for spatial regularization
- Cross-entropy minimization for automatic beta (spatial regularization) estimation

The primary purpose is to identify and filter out nonspecific/background markers
before profile discovery, improving downstream cell type deconvolution.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
import logging

logger = logging.getLogger(__name__)


@dataclass
class SMMResult:
    """Container for Spatial Mixture Model results.

    Attributes:
        corrected_matrix: Background-corrected expression matrix (n_spots, n_markers).
            Values are scaled by P(signal | x) per spot.
        raw_matrix: Original expression matrix before any correction (n_spots, n_markers).
        smoothed_matrix: Spatially smoothed corrected matrix (n_spots, n_markers).
        signal_posteriors: Per-spot signal probability P(signal | x) (n_spots, n_markers).
        background_estimates: Per-marker estimated background mean from GMM.
        snr_values: Per-marker signal-to-noise ratio (for diagnostics).
        passing_mask: Boolean mask of markers passing the SNR filter (legacy, kept for compat).
        passing_marker_names: Names of markers that passed filtering (legacy).
        filtered_marker_names: Names of markers that were filtered out (legacy).
        spot_labels: Per-spot binary labels (0=background, 1=signal) from HMRF.
        beta_learned: Learned spatial regularization parameter.
        gmm_params: Dictionary with 'mu' and 'sigma' arrays from GMM fitting.
        metadata: Additional information (convergence, iterations, etc.).
    """
    corrected_matrix: NDArray[np.floating]
    raw_matrix: NDArray[np.floating]
    smoothed_matrix: NDArray[np.floating]
    signal_posteriors: NDArray[np.floating]
    background_estimates: Dict[str, float]
    snr_values: Dict[str, float]
    passing_mask: NDArray[np.bool_]
    passing_marker_names: List[str]
    filtered_marker_names: List[str]
    spot_labels: NDArray[np.int_]
    beta_learned: float
    gmm_params: Dict[str, NDArray[np.floating]]
    metadata: Dict = field(default_factory=dict)


def spatial_smooth_markers(
    X: np.ndarray,
    coords: np.ndarray,
    sigma: float = 1.5,
    k_neighbors: int = 6,
    self_weight: float = 0.5,
) -> np.ndarray:
    """Apply spatial smoothing to marker expression matrix using Gaussian KNN.

    This function smooths marker expression values by averaging with spatially
    neighboring spots, weighted by distance using a Gaussian kernel. Smoothing
    helps reduce point-wise noise and enhances spatial signal patterns before
    downstream analysis.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        coords: Spatial coordinates of shape (n_spots, 2).
        sigma: Gaussian kernel bandwidth in spatial units (default 1.5).
            Larger values produce more smoothing.
        k_neighbors: Number of nearest neighbors to use (default 6).
        self_weight: Weight given to the spot's own value (default 0.5).
            The neighbor average gets weight (1 - self_weight).

    Returns:
        Smoothed expression matrix of shape (n_spots, n_markers).

    Example:
        >>> X_smoothed = spatial_smooth_markers(X, coords, sigma=1.5, k_neighbors=6)
    """
    n_spots, n_markers = X.shape

    # Build KNN tree
    tree = cKDTree(coords)
    # Query k_neighbors + 1 because result includes the point itself
    distances, indices = tree.query(coords, k=k_neighbors + 1)

    X_smoothed = np.zeros_like(X, dtype=np.float64)

    for i in range(n_spots):
        # Exclude self (index 0) from neighbors
        neighbor_dists = distances[i, 1:]
        neighbor_indices = indices[i, 1:]

        # Compute Gaussian weights based on distance
        # Avoid division by zero if sigma is very small
        if sigma > 1e-10:
            weights = np.exp(-0.5 * (neighbor_dists / sigma) ** 2)
        else:
            # Uniform weights if sigma is essentially zero
            weights = np.ones_like(neighbor_dists)

        # Normalize weights to sum to 1
        weight_sum = weights.sum()
        if weight_sum > 1e-10:
            weights = weights / weight_sum
        else:
            weights = np.ones_like(weights) / len(weights)

        # Compute smoothed values: weighted combination of self and neighbors
        for m in range(n_markers):
            neighbor_values = X[neighbor_indices, m]
            neighbor_avg = np.dot(weights, neighbor_values)
            X_smoothed[i, m] = self_weight * X[i, m] + (1 - self_weight) * neighbor_avg

    return X_smoothed


def _build_spatial_graph(
    coords: np.ndarray,
    k: int = 6,
) -> List[List[int]]:
    """Build k-nearest neighbor graph from spatial coordinates.

    Args:
        coords: Spatial coordinates array of shape (n_spots, 2).
        k: Number of nearest neighbors to consider.

    Returns:
        List of neighbor indices for each spot.
    """
    tree = cKDTree(coords)
    # Query k+1 because the point itself is included
    _, indices = tree.query(coords, k=k + 1)
    # Remove self-references (first column)
    neighbors = [list(idx[1:]) for idx in indices]
    return neighbors


def _initialize_gmm_per_marker(
    X: np.ndarray,
    n_components: int = 2,
    random_state: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Initialize GMM parameters for each marker using sklearn.

    Fits a 2-component GMM to each marker independently to get initial
    estimates of background (component 0) and signal (component 1) distributions.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        n_components: Number of GMM components (default 2 for bg/signal).
        random_state: Random seed for reproducibility.

    Returns:
        Tuple of (mu, sigma, initial_labels):
            - mu: Means array of shape (n_markers, n_components)
            - sigma: Standard deviations array of shape (n_markers, n_components)
            - initial_labels: Initial spot labels of shape (n_spots,)
    """
    n_spots, n_markers = X.shape
    mu = np.zeros((n_markers, n_components))
    sigma = np.zeros((n_markers, n_components))

    # Track per-marker GMM posteriors for consensus labeling
    all_posteriors = np.zeros((n_spots, n_markers))

    for m in range(n_markers):
        x_m = X[:, m].reshape(-1, 1)

        # Handle edge case of constant values
        if np.std(x_m) < 1e-10:
            mu[m, :] = [np.mean(x_m), np.mean(x_m)]
            sigma[m, :] = [1e-6, 1e-6]
            all_posteriors[:, m] = 0.5
            continue

        gmm = GaussianMixture(
            n_components=n_components,
            random_state=random_state,
            max_iter=100,
        )
        gmm.fit(x_m)

        # Get component means and order by mean (lower = background)
        means = gmm.means_.flatten()
        order = np.argsort(means)

        mu[m, :] = means[order]
        sigma[m, :] = np.sqrt(gmm.covariances_.flatten()[order])

        # Store posterior probability of signal class
        posteriors = gmm.predict_proba(x_m)
        signal_idx = order[1]  # Higher mean component
        all_posteriors[:, m] = posteriors[:, signal_idx]

    # Consensus initial labels: average posterior across markers
    mean_posterior = np.mean(all_posteriors, axis=1)
    initial_labels = (mean_posterior > 0.5).astype(int)

    return mu, sigma, initial_labels


def _compute_gmm_log_likelihood(
    x: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> np.ndarray:
    """Compute GMM log-likelihood for each spot and component.

    Args:
        x: Expression values for one marker, shape (n_spots,).
        mu: Component means, shape (n_components,).
        sigma: Component standard deviations, shape (n_components,).

    Returns:
        Log-likelihood array of shape (n_spots, n_components).
    """
    n_spots = len(x)
    n_components = len(mu)
    ll = np.zeros((n_spots, n_components))

    for k in range(n_components):
        # Gaussian log-likelihood
        if sigma[k] < 1e-10:
            sigma_k = 1e-10
        else:
            sigma_k = sigma[k]

        ll[:, k] = -0.5 * np.log(2 * np.pi * sigma_k**2) - \
                   0.5 * ((x - mu[k]) / sigma_k)**2

    return ll


def _e_step_icm(
    X: np.ndarray,
    labels: np.ndarray,
    neighbors: List[List[int]],
    mu: np.ndarray,
    sigma: np.ndarray,
    beta: float,
    max_icm_iter: int = 10,
) -> np.ndarray:
    """E-step using Iterated Conditional Modes (ICM).

    Updates spot labels using MAP estimation with spatial prior.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        labels: Current spot labels of shape (n_spots,).
        neighbors: List of neighbor indices for each spot.
        mu: Component means of shape (n_markers, n_components).
        sigma: Component std devs of shape (n_markers, n_components).
        beta: Spatial regularization parameter.
        max_icm_iter: Maximum ICM iterations.

    Returns:
        Updated spot labels of shape (n_spots,).
    """
    n_spots, n_markers = X.shape
    labels = labels.copy()

    for _ in range(max_icm_iter):
        labels_old = labels.copy()

        for i in range(n_spots):
            # Compute GMM log-likelihood for each component
            ll_components = np.zeros(2)
            for k in range(2):
                for m in range(n_markers):
                    ll_m = _compute_gmm_log_likelihood(
                        np.array([X[i, m]]),
                        mu[m, :],
                        sigma[m, :],
                    )
                    ll_components[k] += ll_m[0, k]

            # Compute spatial prior (Potts model)
            neighbor_labels = labels[neighbors[i]]
            spatial_prior = np.zeros(2)
            for k in range(2):
                # Count neighbors with same label
                n_same = np.sum(neighbor_labels == k)
                spatial_prior[k] = beta * n_same

            # MAP assignment
            posterior = ll_components + spatial_prior
            labels[i] = np.argmax(posterior)

        # Check convergence
        if np.all(labels == labels_old):
            break

    return labels


def _m_step_gmm(
    X: np.ndarray,
    labels: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """M-step: update GMM parameters from current labels.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        labels: Current spot labels of shape (n_spots,).

    Returns:
        Tuple of (mu, sigma) with updated parameters.
    """
    n_markers = X.shape[1]
    mu = np.zeros((n_markers, 2))
    sigma = np.zeros((n_markers, 2))

    for k in range(2):
        mask = labels == k
        if np.sum(mask) < 2:
            # Not enough samples for this class
            continue

        for m in range(n_markers):
            x_m_k = X[mask, m]
            mu[m, k] = np.mean(x_m_k)
            sigma[m, k] = np.std(x_m_k) + 1e-6  # Add small constant for stability

    return mu, sigma


def _compute_gmm_posterior(
    X: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> np.ndarray:
    """Compute GMM posterior probability for signal class.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        mu: Component means of shape (n_markers, n_components).
        sigma: Component std devs of shape (n_markers, n_components).

    Returns:
        Posterior probability of signal class for each spot, shape (n_spots,).
    """
    n_spots, n_markers = X.shape

    # Sum log-likelihoods across markers
    ll_bg = np.zeros(n_spots)
    ll_signal = np.zeros(n_spots)

    for m in range(n_markers):
        ll_m = _compute_gmm_log_likelihood(X[:, m], mu[m, :], sigma[m, :])
        ll_bg += ll_m[:, 0]
        ll_signal += ll_m[:, 1]

    # Convert to posterior using softmax
    max_ll = np.maximum(ll_bg, ll_signal)
    p_bg = np.exp(ll_bg - max_ll)
    p_signal = np.exp(ll_signal - max_ll)
    posterior = p_signal / (p_bg + p_signal + 1e-10)

    return posterior


def _compute_spatial_posterior(
    labels: np.ndarray,
    neighbors: List[List[int]],
    beta: float,
) -> np.ndarray:
    """Compute spatial posterior probability for signal class.

    Uses pseudo-likelihood approximation of the Potts model.

    Args:
        labels: Current spot labels of shape (n_spots,).
        neighbors: List of neighbor indices for each spot.
        beta: Spatial regularization parameter.

    Returns:
        Spatial posterior probability of signal class, shape (n_spots,).
    """
    n_spots = len(labels)
    spatial_posterior = np.zeros(n_spots)

    for i in range(n_spots):
        neighbor_labels = labels[neighbors[i]]
        n_signal = np.sum(neighbor_labels == 1)
        n_bg = len(neighbor_labels) - n_signal

        # Potts model pseudo-likelihood
        energy_bg = -beta * n_bg
        energy_signal = -beta * n_signal

        # Softmax to get probability
        max_energy = max(energy_bg, energy_signal)
        p_bg = np.exp(-(energy_bg - max_energy))
        p_signal = np.exp(-(energy_signal - max_energy))
        spatial_posterior[i] = p_signal / (p_bg + p_signal + 1e-10)

    return spatial_posterior


def _cross_entropy(p: np.ndarray, q: np.ndarray) -> float:
    """Compute cross-entropy H(p, q) = -sum(p * log(q)).

    Args:
        p: True distribution.
        q: Predicted distribution.

    Returns:
        Cross-entropy value.
    """
    # Clip to avoid log(0)
    q = np.clip(q, 1e-10, 1 - 1e-10)
    return -np.mean(p * np.log(q) + (1 - p) * np.log(1 - q))


def _estimate_beta_cross_entropy(
    X: np.ndarray,
    labels: np.ndarray,
    neighbors: List[List[int]],
    mu: np.ndarray,
    sigma: np.ndarray,
    beta_range: Tuple[float, float] = (0.1, 5.0),
) -> float:
    """Estimate optimal beta via pseudo-likelihood cross-entropy minimization.

    Finds beta that minimizes the cross-entropy between GMM posteriors
    and spatial posteriors.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        labels: Current spot labels of shape (n_spots,).
        neighbors: List of neighbor indices for each spot.
        mu: Component means of shape (n_markers, n_components).
        sigma: Component std devs of shape (n_markers, n_components).
        beta_range: Search range for beta.

    Returns:
        Optimal beta value.
    """
    # Compute GMM posterior (target distribution)
    p_gmm = _compute_gmm_posterior(X, mu, sigma)

    def objective(beta: float) -> float:
        # Compute spatial posterior for this beta
        q_spatial = _compute_spatial_posterior(labels, neighbors, beta)
        return _cross_entropy(p_gmm, q_spatial)

    # Minimize cross-entropy
    result = minimize_scalar(
        objective,
        bounds=beta_range,
        method='bounded',
    )

    return result.x


def _compute_marker_snr(
    X: np.ndarray,
    labels: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    marker_names: List[str],
) -> Dict[str, float]:
    """Compute signal-to-noise ratio for each marker.

    SNR = (mu_signal - mu_background) / sqrt(sigma_bg^2 + sigma_signal^2)

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        labels: Spot labels of shape (n_spots,).
        mu: Component means of shape (n_markers, n_components).
        sigma: Component std devs of shape (n_markers, n_components).
        marker_names: List of marker names.

    Returns:
        Dictionary mapping marker names to SNR values.
    """
    n_markers = X.shape[1]
    snr_values = {}

    for m in range(n_markers):
        mu_bg, mu_signal = mu[m, 0], mu[m, 1]
        sigma_bg, sigma_signal = sigma[m, 0], sigma[m, 1]

        # Compute SNR
        denominator = np.sqrt(sigma_bg**2 + sigma_signal**2)
        if denominator < 1e-10:
            snr = 0.0
        else:
            snr = (mu_signal - mu_bg) / denominator

        snr_values[marker_names[m]] = max(0.0, snr)  # SNR should be non-negative

    return snr_values


def _compute_signal_fraction(
    X: np.ndarray,
    labels: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> np.ndarray:
    """Compute per-marker fraction of spots assigned to signal class.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        labels: Spot labels of shape (n_spots,).
        mu: Component means of shape (n_markers, n_components).
        sigma: Component std devs of shape (n_markers, n_components).

    Returns:
        Signal fraction for each marker, shape (n_markers,).
    """
    n_spots, n_markers = X.shape
    signal_fractions = np.zeros(n_markers)

    for m in range(n_markers):
        # Use per-marker GMM to classify spots
        ll = _compute_gmm_log_likelihood(X[:, m], mu[m, :], sigma[m, :])
        marker_labels = np.argmax(ll, axis=1)
        signal_fractions[m] = np.mean(marker_labels == 1)

    return signal_fractions


def _apply_background_correction(
    X: np.ndarray,
    mu: np.ndarray,
    method: str = "subtract_mean",
) -> np.ndarray:
    """Apply background correction to expression matrix.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        mu: Component means of shape (n_markers, n_components).
        method: Correction method ("subtract_mean" or "subtract_quantile").

    Returns:
        Background-corrected expression matrix.
    """
    X_corrected = X.copy()

    if method == "subtract_mean":
        # Subtract background mean from each marker
        for m in range(X.shape[1]):
            X_corrected[:, m] = X[:, m] - mu[m, 0]
            # Clip negative values to zero
            X_corrected[:, m] = np.maximum(X_corrected[:, m], 0)

    return X_corrected


def _apply_per_spot_background_correction(
    X: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply soft-scale background correction per spot using GMM posteriors.

    For each spot i and marker m:
    1. Compute P(signal | X[i,m]) using GMM parameters
    2. X_corrected[i,m] = X[i,m] * P(signal | X[i,m])

    This keeps all spots but downweights background values proportionally
    to their probability of being background signal.

    Args:
        X: Raw expression matrix of shape (n_spots, n_markers).
        mu: GMM component means of shape (n_markers, 2) where [:, 0] is
            background and [:, 1] is signal.
        sigma: GMM component standard deviations of shape (n_markers, 2).

    Returns:
        Tuple of (corrected_matrix, signal_posteriors):
            - corrected_matrix: Background-scaled values (n_spots, n_markers)
            - signal_posteriors: P(signal | x) per spot/marker (n_spots, n_markers)
    """
    n_spots, n_markers = X.shape
    signal_posteriors = np.zeros((n_spots, n_markers), dtype=np.float64)

    for m in range(n_markers):
        # Get GMM parameters for this marker
        mu_bg, mu_signal = mu[m, 0], mu[m, 1]
        sigma_bg = max(sigma[m, 0], 1e-10)  # Avoid division by zero
        sigma_signal = max(sigma[m, 1], 1e-10)

        # Compute Gaussian likelihoods for background and signal components
        ll_bg = norm.pdf(X[:, m], mu_bg, sigma_bg)
        ll_signal = norm.pdf(X[:, m], mu_signal, sigma_signal)

        # Posterior: P(signal | x) = ll_signal / (ll_bg + ll_signal)
        # Using equal priors for background and signal
        total = ll_bg + ll_signal + 1e-10
        signal_posteriors[:, m] = ll_signal / total

    # Soft-scale: multiply original values by signal probability
    # Background spots (low P(signal)) get downweighted toward zero
    # Signal spots (high P(signal)) retain most of their original value
    corrected_matrix = X * signal_posteriors

    return corrected_matrix, signal_posteriors


def fit_spatial_mixture_model(
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    *,
    k_neighbors: int = 6,
    max_iter: int = 50,
    tol: float = 1e-4,
    beta_init: float = 1.0,
    beta_range: Tuple[float, float] = (0.1, 5.0),
    snr_threshold: float = 1.5,
    min_signal_fraction: float = 0.05,
    max_signal_fraction: float = 0.95,
    apply_smoothing: bool = True,
    smoothing_sigma: float = 1.5,
    smoothing_k_neighbors: int = 6,
    smoothing_self_weight: float = 0.5,
    verbose: bool = True,
) -> SMMResult:
    """Fit Spatial Mixture Model for background correction.

    This function implements HMRF-EM (Hidden Markov Random Field with EM)
    for spatial background correction with the following pipeline:

    1. GMM fitting on RAW data (per marker) to identify background vs signal
    2. HMRF-EM iterations to refine classifications with spatial context
    3. Per-spot soft-scale background correction: X * P(signal | x)
    4. Spatial smoothing on CORRECTED data (if enabled)

    All markers proceed through the pipeline (no pre-filtering).
    Filtering decisions should be made downstream based on Moran's I.

    Args:
        X: Expression matrix of shape (n_spots, n_markers).
        coords: Spatial coordinates of shape (n_spots, 2).
        marker_names: List of marker names.
        k_neighbors: Number of neighbors for spatial graph (default 6).
        max_iter: Maximum EM iterations (default 50).
        tol: Convergence tolerance (default 1e-4).
        beta_init: Initial spatial regularization (default 1.0).
        beta_range: Search range for beta optimization (default (0.1, 5.0)).
        snr_threshold: Minimum SNR for diagnostic logging (default 1.5).
            Note: Markers are NOT filtered by SNR; this is for diagnostics only.
        min_signal_fraction: Min signal fraction for diagnostic logging (default 0.05).
        max_signal_fraction: Max signal fraction for diagnostic logging (default 0.95).
        apply_smoothing: Whether to apply spatial smoothing AFTER correction (default True).
        smoothing_sigma: Gaussian kernel bandwidth for smoothing (default 1.5).
        smoothing_k_neighbors: Number of neighbors for smoothing (default 6).
        smoothing_self_weight: Weight for spot's own value in smoothing (default 0.5).
        verbose: Whether to print progress (default True).

    Returns:
        SMMResult containing:
            - raw_matrix: Original input data
            - corrected_matrix: Background-scaled values (X * P(signal))
            - smoothed_matrix: Spatially smoothed corrected data
            - signal_posteriors: P(signal | x) per spot/marker
            - gmm_params: Learned GMM parameters (mu, sigma)
            - snr_values: Per-marker SNR (for diagnostics)
            - And other metadata
    """
    n_spots, n_markers = X.shape

    if len(marker_names) != n_markers:
        raise ValueError(
            f"Number of marker names ({len(marker_names)}) must match "
            f"number of markers in X ({n_markers})"
        )

    # Store raw data before any modifications
    X_raw = X.copy()

    if verbose:
        logger.info(f"Fitting SMM on {n_spots} spots, {n_markers} markers")
        logger.info("Pipeline: GMM(raw) -> HMRF-EM -> Soft-scale correction -> Smooth")

    # Step 1: Build spatial graph
    if verbose:
        logger.info(f"Building spatial graph with k={k_neighbors} neighbors")
    neighbors = _build_spatial_graph(coords, k=k_neighbors)

    # Step 2: Initialize GMM parameters on RAW data (no pre-smoothing)
    if verbose:
        logger.info("Initializing GMM parameters per marker on RAW data")
    mu, sigma, labels = _initialize_gmm_per_marker(X_raw)
    beta = beta_init

    # Step 3: HMRF-EM iterations on RAW data
    convergence_history = []
    iteration = 0

    for iteration in range(max_iter):
        labels_old = labels.copy()

        # E-step: ICM update
        labels = _e_step_icm(X_raw, labels, neighbors, mu, sigma, beta)

        # M-step: Update GMM parameters
        mu, sigma = _m_step_gmm(X_raw, labels)

        # Beta update: Cross-entropy minimization
        beta = _estimate_beta_cross_entropy(
            X_raw, labels, neighbors, mu, sigma, beta_range
        )

        # Check convergence
        changed_fraction = np.mean(labels != labels_old)
        convergence_history.append(changed_fraction)

        if verbose and (iteration + 1) % 10 == 0:
            logger.info(
                f"Iteration {iteration + 1}: {changed_fraction:.4f} labels changed, "
                f"beta={beta:.3f}"
            )

        if changed_fraction < tol:
            if verbose:
                logger.info(f"Converged at iteration {iteration + 1}")
            break

    # Step 4: Per-spot background correction using soft-scale
    # X_corrected = X * P(signal | x)
    if verbose:
        logger.info("Applying per-spot soft-scale background correction")
    corrected_matrix, signal_posteriors = _apply_per_spot_background_correction(
        X_raw, mu, sigma
    )

    # Step 5: Spatial smoothing on CORRECTED data (if enabled)
    if apply_smoothing:
        if verbose:
            logger.info(
                f"Applying spatial smoothing on corrected data: sigma={smoothing_sigma}, "
                f"k={smoothing_k_neighbors}, self_weight={smoothing_self_weight}"
            )
        smoothed_matrix = spatial_smooth_markers(
            corrected_matrix, coords,
            sigma=smoothing_sigma,
            k_neighbors=smoothing_k_neighbors,
            self_weight=smoothing_self_weight,
        )
    else:
        smoothed_matrix = corrected_matrix.copy()

    # Step 6: Compute marker statistics (for diagnostics, NOT for filtering)
    snr_values = _compute_marker_snr(X_raw, labels, mu, sigma, marker_names)
    signal_fractions = _compute_signal_fraction(X_raw, labels, mu, sigma)

    # Log SNR info for diagnostics (all markers proceed, no filtering)
    passing_mask = np.ones(n_markers, dtype=bool)  # All markers pass
    passing_marker_names = list(marker_names)
    filtered_marker_names = []

    if verbose:
        # Log markers that would have been filtered under old system (for comparison)
        would_filter = []
        for m, name in enumerate(marker_names):
            snr = snr_values[name]
            sig_frac = signal_fractions[m]
            passes_snr = snr >= snr_threshold
            passes_min_frac = sig_frac >= min_signal_fraction
            passes_max_frac = sig_frac <= max_signal_fraction
            if not (passes_snr and passes_min_frac and passes_max_frac):
                would_filter.append(name)
                logger.debug(
                    f"Marker {name}: SNR={snr:.2f}, sig_frac={sig_frac:.2f} "
                    f"(would have been filtered under old system)"
                )
        if would_filter:
            logger.info(
                f"Note: {len(would_filter)} markers have low SNR but will proceed "
                f"to Moran's I filtering: {would_filter[:5]}{'...' if len(would_filter) > 5 else ''}"
            )

    # Step 7: Compute background estimates
    background_estimates = {
        marker_names[m]: mu[m, 0] for m in range(n_markers)
    }

    if verbose:
        logger.info(
            f"SMM complete: All {n_markers} markers proceed to downstream analysis, "
            f"beta={beta:.3f}"
        )

    return SMMResult(
        corrected_matrix=corrected_matrix,
        raw_matrix=X_raw,
        smoothed_matrix=smoothed_matrix,
        signal_posteriors=signal_posteriors,
        background_estimates=background_estimates,
        snr_values=snr_values,
        passing_mask=passing_mask,
        passing_marker_names=passing_marker_names,
        filtered_marker_names=filtered_marker_names,
        spot_labels=labels,
        beta_learned=beta,
        gmm_params={'mu': mu, 'sigma': sigma},
        metadata={
            "n_iterations": iteration + 1,
            "convergence_history": convergence_history,
            "signal_fractions": {
                marker_names[m]: signal_fractions[m] for m in range(n_markers)
            },
            "k_neighbors": k_neighbors,
            "snr_threshold": snr_threshold,
            "min_signal_fraction": min_signal_fraction,
            "max_signal_fraction": max_signal_fraction,
            "smoothing_applied": apply_smoothing,
            "smoothing_sigma": smoothing_sigma if apply_smoothing else None,
            "smoothing_k_neighbors": smoothing_k_neighbors if apply_smoothing else None,
            "smoothing_self_weight": smoothing_self_weight if apply_smoothing else None,
        },
    )
