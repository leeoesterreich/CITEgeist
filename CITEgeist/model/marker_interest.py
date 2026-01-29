"""
Marker interest scoring for spatial transcriptomics data.

Identifies markers with significant signal variation using:
1. Kurtosis gating: Peaked distributions indicate real signal vs flat noise
2. GMM signal/noise separation: SNR-based classification
3. Moran's I on RAW data: Spatial autocorrelation (must use raw, not smoothed)

This module provides informational ranking only - no automatic profile selection.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from scipy.stats import kurtosis as scipy_kurtosis
from sklearn.mixture import GaussianMixture


@dataclass
class MarkerInterest:
    """Interest metrics for a single marker."""

    name: str
    interest_score: float  # Combined score (higher = more interesting)
    kurtosis: float  # Excess kurtosis (>2.0 indicates peaked)
    gmm_snr: float  # Signal-to-noise ratio from GMM
    gmm_signal_fraction: float  # Fraction of spots classified as signal
    morans_i: float  # Global Moran's I on raw data
    morans_i_pvalue: float  # Permutation p-value for Moran's I
    passed_kurtosis: bool  # True if kurtosis > threshold
    passed_gmm: bool  # True if in high-SNR GMM component
    passed_morans: bool  # True if Moran's I significant


@dataclass
class MarkerInterestResult:
    """Results container for marker interest analysis."""

    markers: List[MarkerInterest]
    kurtosis_threshold: float
    morans_threshold: float
    morans_k: int
    morans_alpha: float
    signal_masks: Optional[NDArray[np.bool_]] = field(default=None, repr=False)  # (n_spots, n_markers) boolean, True = signal component
    signal_mask_marker_names: Optional[List[str]] = None  # marker names corresponding to columns

    def to_dataframe(self) -> pd.DataFrame:
        """Convert results to ranked DataFrame sorted by interest_score descending."""
        records = [
            {
                "marker": m.name,
                "interest_score": m.interest_score,
                "kurtosis": m.kurtosis,
                "gmm_snr": m.gmm_snr,
                "gmm_signal_fraction": m.gmm_signal_fraction,
                "morans_i": m.morans_i,
                "morans_i_pvalue": m.morans_i_pvalue,
                "passed_kurtosis": m.passed_kurtosis,
                "passed_gmm": m.passed_gmm,
                "passed_morans": m.passed_morans,
                "passed_either": m.passed_kurtosis or m.passed_morans,
            }
            for m in self.markers
        ]
        df = pd.DataFrame(records)
        return df.sort_values("interest_score", ascending=False).reset_index(drop=True)

    @property
    def interesting_markers(self) -> List[str]:
        """Return names of markers passing EITHER kurtosis OR Moran's I gate (plus GMM)."""
        return [
            m.name for m in self.markers
            if (m.passed_kurtosis or m.passed_morans) and m.passed_gmm
        ]

    @property
    def boring_markers(self) -> List[str]:
        """Return names of markers failing both kurtosis AND Moran's I gates (or GMM)."""
        return [
            m.name for m in self.markers
            if not ((m.passed_kurtosis or m.passed_morans) and m.passed_gmm)
        ]


def _compute_kurtosis(X: NDArray[np.floating]) -> NDArray[np.floating]:
    """
    Compute excess kurtosis for each marker.

    Excess kurtosis (Fisher=True) is 0 for normal distribution.
    Peaked distributions (bimodal signals) have kurtosis > 2.
    Flat distributions (noise) have kurtosis < 2.

    Args:
        X: Expression matrix (n_spots, n_markers).

    Returns:
        Array of kurtosis values (n_markers,).
    """
    n_markers = X.shape[1]
    kurtosis_values = np.zeros(n_markers)

    for m in range(n_markers):
        values = X[:, m]
        if np.std(values) < 1e-10:
            kurtosis_values[m] = 0.0
        else:
            kurtosis_values[m] = scipy_kurtosis(values, fisher=True)

    return kurtosis_values


def _fit_gmm_per_marker(
    X: NDArray[np.floating],
    seed: int,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.bool_]]:
    """
    Fit 2-component GMM to each marker to separate signal from noise.

    Args:
        X: Expression matrix (n_spots, n_markers).
        seed: Random seed for GMM initialization.

    Returns:
        Tuple of:
        - snr_values: SNR = (mu_signal - mu_background) / sigma_background (n_markers,)
        - signal_fractions: Fraction of spots in signal component (n_markers,)
        - signal_masks: Boolean array (n_spots, n_markers) where True means the spot
          is classified in the signal (higher-mean) GMM component.
    """
    n_spots, n_markers = X.shape
    snr_values = np.zeros(n_markers)
    signal_fractions = np.zeros(n_markers)
    signal_masks = np.zeros((n_spots, n_markers), dtype=bool)

    for m in range(n_markers):
        values = X[:, m].reshape(-1, 1)

        # Skip markers with no variance
        if np.std(values) < 1e-10:
            snr_values[m] = 0.0
            signal_fractions[m] = 0.0
            continue

        try:
            gmm = GaussianMixture(
                n_components=2,
                covariance_type="full",
                random_state=seed,
                n_init=3,
                max_iter=100,
            )
            gmm.fit(values)

            # Identify signal (higher mean) vs background (lower mean)
            means = gmm.means_.flatten()
            stds = np.sqrt(gmm.covariances_.flatten())
            weights = gmm.weights_

            if means[0] > means[1]:
                signal_idx, bg_idx = 0, 1
            else:
                signal_idx, bg_idx = 1, 0

            mu_signal = means[signal_idx]
            mu_bg = means[bg_idx]
            sigma_bg = max(stds[bg_idx], 1e-6)  # Avoid division by zero

            snr_values[m] = (mu_signal - mu_bg) / sigma_bg
            signal_fractions[m] = weights[signal_idx]

            # Per-spot signal classification
            posteriors = gmm.predict_proba(values)
            signal_masks[:, m] = posteriors[:, signal_idx] > 0.5

        except Exception as e:
            logging.debug(f"GMM fitting failed for marker {m}: {e}")
            snr_values[m] = 0.0
            signal_fractions[m] = 0.0

    return snr_values, signal_fractions, signal_masks


def _spatial_smooth(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int,
) -> NDArray[np.floating]:
    """
    Apply spatial smoothing by averaging each spot with its k nearest neighbors.

    This creates local averages that amplify true spatial patterns while
    reducing spot-level noise.

    Args:
        X: Expression matrix (n_spots, n_markers).
        coords: Spatial coordinates (n_spots, 2).
        k: Number of nearest neighbors for smoothing.

    Returns:
        Smoothed expression matrix (n_spots, n_markers).
    """
    n_spots = X.shape[0]
    tree = cKDTree(coords)

    # Query k+1 to include self
    query_k = min(k + 1, n_spots)
    _, idx = tree.query(coords, k=query_k)

    # Handle 1D case (very few spots)
    if idx.ndim == 1:
        return X.copy()

    # Average over neighbors (including self)
    smoothed = np.zeros_like(X)
    for i in range(n_spots):
        neighbor_indices = idx[i]
        smoothed[i] = np.mean(X[neighbor_indices], axis=0)

    return smoothed


def _compute_morans_i_batch(
    Z: NDArray[np.floating],
    coords: NDArray[np.floating],
    k: int,
    rng: np.random.Generator,
    n_perm: int,
) -> List[Tuple[float, float]]:
    """
    Compute Moran's I and permutation p-value for each marker.

    Args:
        Z: Standardized expression (n_spots, n_markers). Can be smoothed.
        coords: Spatial coordinates (n_spots, 2).
        k: Number of nearest neighbors.
        rng: Random number generator for permutations.
        n_perm: Number of permutations for null distribution.

    Returns:
        List of tuples (I_obs, p_value) for each marker.
    """
    n_spots, n_markers = Z.shape

    # Edge case: need at least 2 spots
    if n_spots < 2:
        return [(np.nan, 1.0)] * n_markers

    tree = cKDTree(coords)
    query_k = min(k + 1, n_spots)
    dists, idx = tree.query(coords, k=query_k)

    # Handle 1D array case
    if idx.ndim == 1:
        return [(np.nan, 1.0)] * n_markers

    weights = np.ones_like(idx, dtype=float)
    # Remove self (distance 0)
    idx = idx[:, 1:]
    weights = weights[:, 1:]

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

        # Permutation null distribution
        null_I = np.zeros(n_perm)
        for b in range(n_perm):
            perm = rng.permutation(z)
            neighbor_perm = perm[idx]
            num_perm = np.sum(weights * (perm[:, None] * neighbor_perm))
            null_I[b] = (n_spots / S0) * (num_perm / denom)

        p_value = (1 + np.sum(null_I >= I_obs)) / (1 + len(null_I))
        results.append((float(I_obs), float(p_value)))

    return results


def _fit_kurtosis_gmm(
    kurtosis_values: NDArray[np.floating],
    seed: int,
) -> Tuple[float, NDArray[np.bool_]]:
    """
    Fit 2-component GMM to kurtosis values to separate high/low kurtosis markers.

    Instead of a hard threshold, this learns the boundary between "interesting"
    (high kurtosis, peaked distributions) and "boring" (low kurtosis, flat
    distributions) markers from the data itself.

    Args:
        kurtosis_values: Array of kurtosis values (n_markers,).
        seed: Random seed for GMM initialization.

    Returns:
        Tuple of:
        - learned_threshold: The kurtosis value at the decision boundary
        - passed_kurtosis: Boolean array indicating high-kurtosis markers (n_markers,)
    """
    # Filter out NaN/infinite values for fitting
    valid_mask = np.isfinite(kurtosis_values)
    valid_kurtosis = kurtosis_values[valid_mask]

    if len(valid_kurtosis) < 4:
        # Not enough data, fall back to all pass
        logging.warning("Too few valid kurtosis values for GMM fitting, skipping kurtosis gate")
        return 0.0, np.ones(len(kurtosis_values), dtype=bool)

    try:
        gmm = GaussianMixture(
            n_components=2,
            covariance_type="full",
            random_state=seed,
            n_init=5,
            max_iter=200,
        )
        gmm.fit(valid_kurtosis.reshape(-1, 1))

        # Identify high vs low component
        means = gmm.means_.flatten()
        if means[0] > means[1]:
            high_idx, low_idx = 0, 1
        else:
            high_idx, low_idx = 1, 0

        # Decision boundary is where posterior probability = 0.5
        # For Gaussian: midpoint weighted by variances
        mu_high = means[high_idx]
        mu_low = means[low_idx]
        var_high = gmm.covariances_.flatten()[high_idx]
        var_low = gmm.covariances_.flatten()[low_idx]

        # Approximate threshold as weighted midpoint
        # More precise would solve the quadratic, but this is good enough
        w_high = 1.0 / max(np.sqrt(var_high), 0.01)
        w_low = 1.0 / max(np.sqrt(var_low), 0.01)
        learned_threshold = (mu_high * w_low + mu_low * w_high) / (w_high + w_low)

        # Assign markers to high component based on posterior probability
        posteriors = gmm.predict_proba(valid_kurtosis.reshape(-1, 1))
        p_high = posteriors[:, high_idx]

        # Marker passes if >50% probability of being in high-kurtosis component
        passed_valid = p_high > 0.5

        # Map back to full array
        passed_kurtosis = np.zeros(len(kurtosis_values), dtype=bool)
        passed_kurtosis[valid_mask] = passed_valid

        logging.info(
            f"Kurtosis GMM: low_mean={mu_low:.2f}, high_mean={mu_high:.2f}, "
            f"threshold={learned_threshold:.2f}, n_high={passed_kurtosis.sum()}"
        )

        return float(learned_threshold), passed_kurtosis

    except Exception as e:
        logging.warning(f"Kurtosis GMM fitting failed: {e}, falling back to threshold=2.0")
        return 2.0, kurtosis_values >= 2.0


def _fit_morans_gmm(
    morans_values: NDArray[np.floating],
    seed: int,
) -> Tuple[float, NDArray[np.bool_]]:
    """
    Fit 2-component GMM to Moran's I values to separate high/low spatial autocorrelation.

    Similar to kurtosis GMM, this learns the boundary between markers with
    strong spatial patterns (high Moran's I) and those without (low/zero Moran's I).

    Args:
        morans_values: Array of Moran's I values (n_markers,).
        seed: Random seed for GMM initialization.

    Returns:
        Tuple of:
        - learned_threshold: The Moran's I value at the decision boundary
        - passed_morans: Boolean array indicating high-Moran's I markers (n_markers,)
    """
    # Filter out NaN/infinite values for fitting
    valid_mask = np.isfinite(morans_values)
    valid_morans = morans_values[valid_mask]

    if len(valid_morans) < 4:
        # Not enough data, fall back to all pass
        logging.warning("Too few valid Moran's I values for GMM fitting, skipping Moran's gate")
        return 0.0, np.ones(len(morans_values), dtype=bool)

    try:
        gmm = GaussianMixture(
            n_components=2,
            covariance_type="full",
            random_state=seed,
            n_init=5,
            max_iter=200,
        )
        gmm.fit(valid_morans.reshape(-1, 1))

        # Identify high vs low component
        means = gmm.means_.flatten()
        if means[0] > means[1]:
            high_idx, low_idx = 0, 1
        else:
            high_idx, low_idx = 1, 0

        # Decision boundary
        mu_high = means[high_idx]
        mu_low = means[low_idx]
        var_high = gmm.covariances_.flatten()[high_idx]
        var_low = gmm.covariances_.flatten()[low_idx]

        # Approximate threshold as weighted midpoint
        w_high = 1.0 / max(np.sqrt(var_high), 0.01)
        w_low = 1.0 / max(np.sqrt(var_low), 0.01)
        learned_threshold = (mu_high * w_low + mu_low * w_high) / (w_high + w_low)

        # Assign markers to high component based on posterior probability
        posteriors = gmm.predict_proba(valid_morans.reshape(-1, 1))
        p_high = posteriors[:, high_idx]

        # Marker passes if >50% probability of being in high-Moran's component
        passed_valid = p_high > 0.5

        # Map back to full array
        passed_morans = np.zeros(len(morans_values), dtype=bool)
        passed_morans[valid_mask] = passed_valid

        logging.info(
            f"Moran's I GMM: low_mean={mu_low:.3f}, high_mean={mu_high:.3f}, "
            f"threshold={learned_threshold:.3f}, n_high={passed_morans.sum()}"
        )

        return float(learned_threshold), passed_morans

    except Exception as e:
        logging.warning(f"Moran's I GMM fitting failed: {e}, falling back to threshold=0.1")
        return 0.1, morans_values >= 0.1


def _compute_interest_score(
    kurtosis: float,
    gmm_snr: float,
    morans_i: float,
) -> float:
    """
    Combine metrics into single interest score.

    Score = kurtosis * gmm_snr * max(morans_i, 0)

    This multiplicative form ensures markers need signal in all three
    dimensions to score highly. Markers with zero in any dimension
    will have a low score.

    Args:
        kurtosis: Excess kurtosis value.
        gmm_snr: Signal-to-noise ratio from GMM.
        morans_i: Moran's I spatial autocorrelation.

    Returns:
        Combined interest score.
    """
    # Clamp negative values to small positive to avoid sign issues
    k = max(kurtosis, 0.01)
    s = max(gmm_snr, 0.01)
    i = max(morans_i, 0.01)

    return k * s * i


def identify_interesting_markers(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    *,
    kurtosis_threshold: Optional[float] = None,
    morans_threshold: Optional[float] = None,
    gmm_snr_threshold: float = 1.0,
    morans_k: int = 8,
    smooth_k: int = 6,
    morans_n_perm: int = 199,
    seed: int = 1234,
    verbose: bool = True,
) -> MarkerInterestResult:
    """
    Identify markers with interesting signal magnitude variation.

    Uses three complementary approaches with OR logic:
    1. Kurtosis: Detects peaked distributions (real signal has peaks where
       specific cell types are located)
    2. Moran's I: Detects spatial autocorrelation on SMOOTHED data
    3. GMM: Separates signal from background noise using 2-component mixture

    A marker is "interesting" if it passes EITHER:
    - High kurtosis population (peaked distribution), OR
    - High Moran's I population (spatial structure)
    AND passes the GMM signal/noise threshold.

    Args:
        X: Raw expression matrix (n_spots, n_markers).
        coords: Spatial coordinates for each spot (n_spots, 2).
        marker_names: Names for each marker column.
        kurtosis_threshold: Kurtosis threshold for passing.
            - None (default): Use adaptive GMM-based threshold learned from data
            - float: Use fixed threshold (e.g., 2.0)
        morans_threshold: Moran's I threshold for passing.
            - None (default): Use adaptive GMM-based threshold learned from data
            - float: Use fixed threshold (e.g., 0.1)
        gmm_snr_threshold: Minimum SNR from GMM to pass (default: 1.0).
        morans_k: Number of neighbors for Moran's I spatial weights (default: 8).
        smooth_k: Number of neighbors for spatial smoothing before Moran's I (default: 6).
        morans_n_perm: Number of permutations for Moran's I null (default: 199).
        seed: Random seed for reproducibility (default: 1234).
        verbose: Log progress information (default: True).

    Returns:
        MarkerInterestResult with ranked list of markers and interest scores.

    Example:
        >>> # Use adaptive thresholds (recommended)
        >>> result = identify_interesting_markers(X, coords, marker_names)
        >>> df = result.to_dataframe()
        >>> df.head(10)  # Top 10 most interesting markers
    """
    X = np.asarray(X, dtype=np.float64)
    coords = np.asarray(coords, dtype=np.float64)

    n_spots, n_markers = X.shape

    if len(marker_names) != n_markers:
        raise ValueError(
            f"Number of marker names ({len(marker_names)}) must match "
            f"number of columns in X ({n_markers})"
        )

    if verbose:
        logging.info(f"Analyzing {n_markers} markers across {n_spots} spots")

    rng = np.random.default_rng(seed)

    # Step 1: Compute kurtosis (on raw data)
    if verbose:
        logging.info("Computing kurtosis for each marker...")
    kurtosis_values = _compute_kurtosis(X)

    # Step 2: Determine kurtosis threshold (adaptive or fixed)
    use_adaptive_kurtosis = kurtosis_threshold is None
    if use_adaptive_kurtosis:
        if verbose:
            logging.info("Fitting GMM to kurtosis values for adaptive threshold...")
        kurtosis_thresh_learned, kurtosis_passed_array = _fit_kurtosis_gmm(kurtosis_values, seed)
    else:
        kurtosis_thresh_learned = kurtosis_threshold
        kurtosis_passed_array = kurtosis_values >= kurtosis_threshold

    # Step 3: Fit GMM per marker for signal/noise separation
    if verbose:
        logging.info("Fitting GMM to separate signal/noise per marker...")
    snr_values, signal_fractions, signal_masks = _fit_gmm_per_marker(X, seed)

    # Step 4: Spatial smoothing before Moran's I
    if verbose:
        logging.info(f"Applying spatial smoothing (k={smooth_k} neighbors)...")
    X_smoothed = _spatial_smooth(X, coords, smooth_k)

    # Step 5: Compute Moran's I on smoothed, z-scored data
    if verbose:
        logging.info(f"Computing Moran's I on smoothed data (k={morans_k}, {morans_n_perm} permutations)...")

    # Z-score the smoothed data
    Z = np.zeros_like(X_smoothed)
    for m in range(n_markers):
        col = X_smoothed[:, m]
        std = np.std(col)
        if std > 1e-10:
            Z[:, m] = (col - np.mean(col)) / std
        else:
            Z[:, m] = 0.0

    morans_results = _compute_morans_i_batch(Z, coords, morans_k, rng, morans_n_perm)

    # Extract Moran's I values for GMM fitting
    morans_values = np.array([r[0] for r in morans_results])

    # Step 6: Determine Moran's I threshold (adaptive or fixed)
    use_adaptive_morans = morans_threshold is None
    if use_adaptive_morans:
        if verbose:
            logging.info("Fitting GMM to Moran's I values for adaptive threshold...")
        morans_thresh_learned, morans_passed_array = _fit_morans_gmm(morans_values, seed)
    else:
        morans_thresh_learned = morans_threshold
        morans_passed_array = morans_values >= morans_threshold

    # Build results
    markers = []
    for m, name in enumerate(marker_names):
        kurt = kurtosis_values[m]
        snr = snr_values[m]
        sig_frac = signal_fractions[m]
        morans_i, morans_p = morans_results[m]

        # Apply gates (using adaptive GMM classifications)
        passed_kurt = bool(kurtosis_passed_array[m])
        passed_gmm = snr >= gmm_snr_threshold
        passed_morans = bool(morans_passed_array[m])

        # Compute combined score
        interest = _compute_interest_score(kurt, snr, morans_i if not np.isnan(morans_i) else 0.0)

        markers.append(MarkerInterest(
            name=name,
            interest_score=interest,
            kurtosis=kurt,
            gmm_snr=snr,
            gmm_signal_fraction=sig_frac,
            morans_i=morans_i if not np.isnan(morans_i) else 0.0,
            morans_i_pvalue=morans_p,
            passed_kurtosis=passed_kurt,
            passed_gmm=passed_gmm,
            passed_morans=passed_morans,
        ))

    # Sort by interest score
    markers.sort(key=lambda x: x.interest_score, reverse=True)

    result = MarkerInterestResult(
        markers=markers,
        kurtosis_threshold=kurtosis_thresh_learned,
        morans_threshold=morans_thresh_learned,
        morans_k=morans_k,
        morans_alpha=0.0,  # Not used with adaptive GMM
        signal_masks=signal_masks,
        signal_mask_marker_names=list(marker_names),
    )

    if verbose:
        n_interesting = len(result.interesting_markers)
        n_kurtosis_only = sum(1 for m in markers if m.passed_kurtosis and not m.passed_morans and m.passed_gmm)
        n_morans_only = sum(1 for m in markers if m.passed_morans and not m.passed_kurtosis and m.passed_gmm)
        n_both = sum(1 for m in markers if m.passed_kurtosis and m.passed_morans and m.passed_gmm)
        logging.info(
            f"Found {n_interesting}/{n_markers} interesting markers "
            f"(kurtosis_only={n_kurtosis_only}, morans_only={n_morans_only}, both={n_both})"
        )

    return result
