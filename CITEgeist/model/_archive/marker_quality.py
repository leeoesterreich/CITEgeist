"""
Marker quality classification via max-SNR approach.

Classifies markers using a simple rule:
  - Compute SNR at local (k1) and regional (k3) scales
  - Fit GMM on max(SNR_local, SNR_regional)
  - PASS if in high-SNR component, FAIL otherwise

This captures both:
  - Abundant markers: high local SNR
  - Rare markers: low local SNR but high regional SNR (spatial coherence)

No delta, no rescue paths, no thresholds - just max(SNR) and GMM.
"""

import logging
from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
from scipy.stats import entropy, kurtosis as scipy_kurtosis
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import kneighbors_graph

logger = logging.getLogger(__name__)


@dataclass
class MarkerQuality:
    """Quality metrics for a single marker."""
    name: str
    passed: bool           # True if PASS (high-SNR + kurtosis gate)
    snr_local: float       # SNR at k1 (local scale)
    snr_regional: float    # SNR at k3 (regional scale)
    max_snr: float         # max(snr_local, snr_regional)
    gmm_p_high: float      # GMM probability of high-SNR component
    kurtosis: float        # Excess kurtosis (peaked distribution indicator)
    passed_snr: bool       # True if passed SNR gate (high-SNR component)
    passed_kurtosis: bool  # True if passed kurtosis gate (>2.0)


@dataclass
class MarkerQualityResult:
    """Results for all markers."""
    markers: List[MarkerQuality]

    @property
    def passed_names(self) -> List[str]:
        """Names of markers that passed."""
        return [m.name for m in self.markers if m.passed]

    @property
    def passed_indices(self) -> List[int]:
        """Indices of markers that passed."""
        return [i for i, m in enumerate(self.markers) if m.passed]

    @property
    def pass_fail_counts(self) -> Tuple[int, int]:
        """Count of PASS and FAIL markers."""
        n_pass = sum(1 for m in self.markers if m.passed)
        n_fail = len(self.markers) - n_pass
        return n_pass, n_fail

    @property
    def gate_counts(self) -> Tuple[int, int, int]:
        """Count markers passing each gate: (snr_only, kurtosis_only, both)."""
        n_snr = sum(1 for m in self.markers if m.passed_snr)
        n_kurt = sum(1 for m in self.markers if m.passed_kurtosis)
        n_both = sum(1 for m in self.markers if m.passed)
        return n_snr, n_kurt, n_both

    def summary(self) -> str:
        """Summary string of classification results."""
        n_pass, n_fail = self.pass_fail_counts
        n_snr, n_kurt, n_both = self.gate_counts
        return (f"PASS: {n_pass}, FAIL: {n_fail} "
                f"(SNR gate: {n_snr}, Kurtosis gate: {n_kurt})")


def filter_markers(
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    k1: Optional[int] = None,
    k3: Optional[int] = None,
    n_smooth_iter: int = 3,
    verbose: bool = True,
) -> MarkerQualityResult:
    """
    Classify markers using max-SNR approach.

    Simple rule: Pass if max(SNR_local, SNR_regional) is in high-SNR GMM component.

    Parameters
    ----------
    X : ndarray (n_spots, n_markers)
        Expression matrix (raw counts or log-transformed).
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    marker_names : list of str
        Names for each marker.
    k1 : int, optional
        Local neighborhood size. If None, scales with n_spots (~1%).
    k3 : int, optional
        Regional neighborhood size. If None, scales with n_spots (~20%).
        This should capture a major region of the tissue.
    n_smooth_iter : int, default=3
        Number of mean-field smoothing iterations.
    verbose : bool, default=True
        Whether to log progress.

    Returns
    -------
    MarkerQualityResult
        Classification results for all markers.
    """
    n_spots, n_markers = X.shape

    # Scale k values with sample size if not provided
    # k1: local (~1% of spots, minimum 6)
    # k3: regional (~20% of spots, captures major tissue regions)
    if k1 is None:
        k1 = max(6, int(n_spots * 0.01))
    if k3 is None:
        k3 = max(30, int(n_spots * 0.20))

    # Ensure k values don't exceed n_spots - 1
    k1 = min(k1, n_spots - 1)
    k3 = min(k3, n_spots - 1)

    if verbose:
        logger.info(f"Classifying {n_markers} markers using max-SNR approach")
        logger.info(f"n_spots={n_spots}, k1={k1} (local), k3={k3} (regional)")

    # Log transform if not already
    x_log = np.log1p(X)

    # Build neighborhood graphs at 2 scales
    if verbose:
        logger.info("Building k-NN neighbor graphs...")
    neighbors_k1 = _build_neighbor_list(coords, k1)
    neighbors_k3 = _build_neighbor_list(coords, k3)

    # Compute SNR at both scales for all markers
    if verbose:
        logger.info("Computing SNR at local and regional scales...")
    all_snr_1 = []
    all_snr_3 = []
    for m in range(n_markers):
        x = x_log[:, m]
        snr_1 = _compute_entropy_snr(x, neighbors_k1, n_smooth_iter)
        snr_3 = _compute_entropy_snr(x, neighbors_k3, n_smooth_iter)
        all_snr_1.append(snr_1)
        all_snr_3.append(snr_3)

    all_snr_1 = np.array(all_snr_1)
    all_snr_3 = np.array(all_snr_3)

    # Compute max SNR for each marker
    all_max_snr = np.maximum(all_snr_1, all_snr_3)

    # Fit GMM on max SNR values
    if verbose:
        logger.info("Fitting GMM on max(SNR_local, SNR_regional)...")
    gmm, high_idx = _fit_snr_gmm(all_max_snr)
    if verbose:
        means = gmm.means_.flatten()
        logger.info(f"GMM components: low={means[1-high_idx]:.4f}, high={means[high_idx]:.4f}")

    # Compute kurtosis for each marker (diagnostic only, not used for gating)
    # Note: In simulated data, GT markers often have LOWER kurtosis than NS markers
    # because real signals are spread across larger spatial domains
    if verbose:
        logger.info("Computing kurtosis for each marker...")
    all_kurtosis = []
    for m in range(n_markers):
        x = x_log[:, m]
        # Excess kurtosis (Fisher=True): normal distribution = 0
        kurt = scipy_kurtosis(x, fisher=True)
        all_kurtosis.append(kurt)
    all_kurtosis = np.array(all_kurtosis)

    # Classify each marker using SNR only
    # Kurtosis is kept as metadata but not used for gating
    if verbose:
        logger.info("Classifying markers...")
    results = []
    for m in range(n_markers):
        snr_local = all_snr_1[m]
        snr_regional = all_snr_3[m]
        max_snr = all_max_snr[m]
        kurt = all_kurtosis[m]

        passed_snr, p_high = _classify_max_snr(max_snr, gmm, high_idx)

        # Classification based on SNR only
        # Kurtosis gate disabled - in simulated data, GT markers have lower kurtosis
        passed = passed_snr

        results.append(MarkerQuality(
            name=marker_names[m],
            passed=passed,
            snr_local=snr_local,
            snr_regional=snr_regional,
            max_snr=max_snr,
            gmm_p_high=p_high,
            kurtosis=kurt,
            passed_snr=passed_snr,
            passed_kurtosis=True,  # Always True since gate is disabled
        ))

    result = MarkerQualityResult(markers=results)

    if verbose:
        logger.info(f"Classification complete: {result.summary()}")

    return result


def _build_neighbor_list(coords: np.ndarray, k: int) -> List[List[int]]:
    """
    Build k-NN neighbor list from coordinates.

    Parameters
    ----------
    coords : ndarray (n_spots, 2)
        Spatial coordinates.
    k : int
        Number of neighbors.

    Returns
    -------
    list of list of int
        neighbors[i] contains indices of neighbors for spot i.
    """
    n_spots = coords.shape[0]

    # Ensure k doesn't exceed available neighbors
    k_actual = min(k, n_spots - 1)

    # Build k-NN graph
    A = kneighbors_graph(coords, k_actual, mode='connectivity', include_self=False)

    # Convert to neighbor list
    neighbors = []
    A_csr = A.tocsr()
    for i in range(n_spots):
        start, end = A_csr.indptr[i], A_csr.indptr[i + 1]
        neighbors.append(list(A_csr.indices[start:end]))

    return neighbors


def _compute_entropy_snr(
    x: np.ndarray,
    neighbors: List[List[int]],
    n_iter: int = 3,
) -> float:
    """
    Compute SNR via GMM posterior entropy + spatial smoothing.

    SNR = -mean(H(posterior)), where H is entropy.
    Higher SNR means more confident signal/background separation.

    Parameters
    ----------
    x : ndarray (n_spots,)
        Log-transformed expression values for one marker.
    neighbors : list of list of int
        Neighbor indices for each spot.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    float
        Spatial SNR (negative mean entropy).
    """
    # Handle constant or near-constant markers
    if np.std(x) < 1e-6:
        return 0.0  # No signal

    # Fit 2-component GMM
    try:
        gmm = GaussianMixture(n_components=2, random_state=0, max_iter=100)
        gmm.fit(x.reshape(-1, 1))
        post = gmm.predict_proba(x.reshape(-1, 1))  # (n_spots, 2)
    except Exception:
        # GMM failed - return low SNR
        return 0.0

    # Spatial smoothing of posteriors
    post_spatial = _spatially_smooth_posteriors(post, neighbors, n_iter)

    # SNR = negative mean entropy
    # entropy expects (n_classes, n_samples), so transpose
    H_spatial = entropy(post_spatial.T)  # (n_spots,)

    return -np.mean(H_spatial)


def _spatially_smooth_posteriors(
    post: np.ndarray,
    neighbors: List[List[int]],
    n_iter: int = 3,
) -> np.ndarray:
    """
    Mean-field spatial smoothing of posteriors.

    No beta parameter - implicitly determined by data.
    This is equivalent to one iteration of a Potts MRF.

    Parameters
    ----------
    post : ndarray (n_spots, n_components)
        GMM posterior probabilities.
    neighbors : list of list of int
        Neighbor indices for each spot.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    ndarray (n_spots, n_components)
        Spatially smoothed posteriors.
    """
    post_sm = post.copy()

    for _ in range(n_iter):
        new_post = post_sm.copy()
        for i, nbrs in enumerate(neighbors):
            if len(nbrs) == 0:
                continue
            # Average neighbor posteriors (Bayesian consensus)
            nbr_mean = post_sm[nbrs].mean(axis=0)
            # Multiplicative update
            new_post[i] *= nbr_mean
            # Renormalize
            new_post[i] /= new_post[i].sum() + 1e-10
        post_sm = new_post

    return post_sm


def _fit_snr_gmm(all_max_snr: np.ndarray) -> Tuple[GaussianMixture, int]:
    """
    Fit 2-component GMM to max-SNR values to find natural boundary.

    Parameters
    ----------
    all_max_snr : ndarray
        max(SNR_local, SNR_regional) for all markers.

    Returns
    -------
    gmm : GaussianMixture
        Fitted GMM model.
    high_idx : int
        Index of the high-SNR component (0 or 1).
    """
    gmm = GaussianMixture(n_components=2, random_state=0, max_iter=100)
    gmm.fit(all_max_snr.reshape(-1, 1))

    # Identify which component is "high SNR" (less negative = higher)
    high_idx = int(np.argmax(gmm.means_.flatten()))

    return gmm, high_idx


def _classify_max_snr(
    max_snr: float,
    gmm: GaussianMixture,
    high_idx: int,
) -> Tuple[bool, float]:
    """
    Classify marker by max SNR using GMM.

    Simple rule: If in high-SNR component → PASS, else FAIL.
    No delta, no rescue paths, no thresholds.

    Parameters
    ----------
    max_snr : float
        max(SNR_local, SNR_regional) for this marker.
    gmm : GaussianMixture
        Fitted GMM on all max-SNR values.
    high_idx : int
        Index of the high-SNR component.

    Returns
    -------
    passed : bool
        True if marker passes (in high-SNR component).
    p_high : float
        GMM probability of belonging to high-SNR component.
    """
    assignment = gmm.predict(np.array([[max_snr]]))[0]
    is_high = (assignment == high_idx)

    # Get probability for metadata
    probs = gmm.predict_proba(np.array([[max_snr]]))[0]
    p_high = probs[high_idx]

    return is_high, p_high


def plot_snr_distribution(
    result: MarkerQualityResult,
    output_path: Optional[str] = None,
) -> None:
    """
    Plot distribution of max-SNR values, colored by PASS/FAIL.

    Parameters
    ----------
    result : MarkerQualityResult
        Classification results.
    output_path : str, optional
        Path to save figure. If None, displays interactively.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))

    passed_snr = [m.max_snr for m in result.markers if m.passed]
    failed_snr = [m.max_snr for m in result.markers if not m.passed]

    if passed_snr:
        ax.hist(passed_snr, bins=20, alpha=0.7, label='PASS', color='green')
    if failed_snr:
        ax.hist(failed_snr, bins=20, alpha=0.7, label='FAIL', color='gray')

    ax.set_xlabel('max(SNR_local, SNR_regional)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Max-SNR Distribution by Classification', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_snr_scatter(
    result: MarkerQualityResult,
    output_path: Optional[str] = None,
    highlight_names: Optional[List[str]] = None,
) -> None:
    """
    Scatter plot of SNR_local vs SNR_regional, colored by PASS/FAIL.

    Parameters
    ----------
    result : MarkerQualityResult
        Classification results.
    output_path : str, optional
        Path to save figure.
    highlight_names : list of str, optional
        Marker names to highlight with labels.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))

    for m in result.markers:
        color = 'green' if m.passed else 'gray'
        alpha = 0.8 if m.passed else 0.4
        ax.scatter(m.snr_local, m.snr_regional, c=color, alpha=alpha, s=50)

        # Highlight specific markers
        if highlight_names and m.name in highlight_names:
            ax.annotate(m.name, (m.snr_local, m.snr_regional), fontsize=8)

    # Add diagonal line (max = local or regional)
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.3, label='SNR_local = SNR_regional')

    ax.set_xlabel('SNR_local (k1)', fontsize=12)
    ax.set_ylabel('SNR_regional (k3)', fontsize=12)
    ax.set_title('Local vs Regional SNR', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
