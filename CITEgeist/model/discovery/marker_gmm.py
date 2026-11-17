"""BIC-selected 1-D Gaussian mixture per marker.

This backs the opt-in brightness correction in marker share
(``brightness_correction=True``), which was evaluated on held-out data and
rejected as the default. Dictionary discovery does not use it: it uses the
negative-binomial mixture in ``nb_mixture.py``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.mixture import GaussianMixture


@dataclass
class MarkerGMMFit:
    """Result of fitting a 1-D GMM to one marker channel.

    means are sorted ascending. positive_mean is the mean of the highest
    component (nan if only one component); positive_mask marks spots assigned
    to that component (all False if only one component).
    """

    n_components: int
    means: np.ndarray
    weights: np.ndarray
    positive_mean: float
    positive_mask: np.ndarray


def fit_marker_gmm(values: np.ndarray, max_components: int = 3, seed: int = 1234) -> MarkerGMMFit:
    """Fit GMMs with 1..max_components components and keep the lowest BIC."""
    v = np.asarray(values, dtype=float).ravel()
    n = v.shape[0]
    if np.std(v) < 1e-10:
        return MarkerGMMFit(1, np.array([float(np.mean(v))]), np.array([1.0]), float("nan"), np.zeros(n, dtype=bool))
    col = v.reshape(-1, 1)
    best, best_bic = None, np.inf
    for k in range(1, max_components + 1):
        if k > n:
            break
        gmm = GaussianMixture(n_components=k, covariance_type="full", random_state=seed, n_init=3, max_iter=200).fit(
            col
        )
        bic = gmm.bic(col)
        if bic < best_bic:
            best, best_bic = gmm, bic
    order = np.argsort(best.means_.ravel())
    means = best.means_.ravel()[order]
    weights = best.weights_[order]
    k = len(means)
    if k == 1:
        return MarkerGMMFit(1, means, weights, float("nan"), np.zeros(n, dtype=bool))
    top = int(order[-1])
    mask = best.predict(col) == top
    return MarkerGMMFit(k, means, weights, float(means[-1]), mask)
