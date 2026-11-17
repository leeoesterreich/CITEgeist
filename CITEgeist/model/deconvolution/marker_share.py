"""Marker-share cell-type proportion estimator (SPOTS-style, CPU only).

Follows Ben-Chetrit et al. 2023 (Nat Biotechnol): cosine-normalize each spot's
antibody vector, take the mean of each cell type's Major marker channels, then
divide by the sum over types so each spot's fractions sum to 1. Spots whose
marker channels are all zero get a uniform 1/K composition and are counted.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


def cosine_normalize(X: np.ndarray) -> np.ndarray:
    """Scale each row (spot) to unit L2 norm; all-zero rows stay zero."""
    X = np.asarray(X, dtype=float)
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    return np.divide(X, norm, out=np.zeros_like(X), where=norm > 0)


def _major_markers(entry) -> list:
    if isinstance(entry, Mapping):
        return list(entry["Major"])
    return list(entry)


def marker_share_proportions(
    X: np.ndarray,
    var_names: Sequence[str],
    profile: Mapping[str, Mapping[str, list] | list],
    brightness: Mapping[str, float] | None = None,
) -> tuple[pd.DataFrame, int]:
    """Estimate per-spot cell-type fractions from marker-channel shares.

    Args:
        X: Raw antibody counts (spots x all channels).
        var_names: Channel names for the columns of X.
        profile: type -> {"Major": [markers], ...} or type -> [markers].
            Only Major markers are used.
        brightness: Optional marker -> b_m. The cosine-normalized channel is
            divided by b_m before the per-type mean.

    Returns:
        (spots x types DataFrame of fractions, number of all-zero-marker spots).

    Raises:
        ValueError: if a profile marker is not in var_names.
    """
    var_names = list(var_names)
    types = list(profile.keys())
    Xc = cosine_normalize(X)
    scores = np.zeros((Xc.shape[0], len(types)))
    for j, ct in enumerate(types):
        markers = _major_markers(profile[ct])
        idx = [var_names.index(m) for m in markers]  # ValueError on missing: loud
        cols = Xc[:, idx]
        if brightness is not None:
            b = np.array([brightness.get(m, 1.0) for m in markers], dtype=float)
            cols = cols / b
        scores[:, j] = cols.mean(axis=1)
    denom = scores.sum(axis=1)
    zero = denom == 0
    frac = np.divide(scores, denom[:, None], out=np.zeros_like(scores), where=denom[:, None] > 0)
    frac[zero] = 1.0 / len(types)
    return pd.DataFrame(frac, columns=types), int(zero.sum())


def estimate_marker_brightness_with_diagnostics(
    X: np.ndarray,
    var_names: Sequence[str],
    markers: Iterable[str],
    seed: int = 1234,
    max_components: int = 3,
) -> tuple[dict[str, float], dict[str, dict]]:
    """Per-marker brightness plus per-marker fallback diagnostics.

    Diagnostics per marker: ``used_brightness``, ``fallback`` (True when the GMM gave no
    usable positive mean; b is then the top-decile mean, or 1.0 if that is not positive),
    ``fallback_value`` ("top_decile_mean" or "one"; fallbacks only) and ``reason`` (None, "one_component", "non_finite_positive_mean" or
    "non_positive_positive_mean").
    """
    from CITEgeist.model.discovery.marker_gmm import fit_marker_gmm

    var_names = list(var_names)
    Xc = cosine_normalize(X)
    out: dict[str, float] = {}
    diag: dict[str, dict] = {}
    for m in markers:
        fit = fit_marker_gmm(Xc[:, var_names.index(m)], max_components=max_components, seed=seed)
        pm = fit.positive_mean
        reason = None
        if fit.n_components == 1:
            reason = "one_component"
        elif not np.isfinite(pm):
            reason = "non_finite_positive_mean"
        elif pm <= 0:
            reason = "non_positive_positive_mean"
        entry = {"fallback": reason is not None, "reason": reason}
        if reason is None:
            b = float(pm)
        else:
            # Fallback: the marker's own top-decile (>= 90th percentile) mean on the same
            # cosine scale. 1.0 would shrink a marker whose real brightness is 0.004-0.81.
            col = Xc[:, var_names.index(m)]
            top = col[col >= np.percentile(col, 90)]
            tm = float(top.mean()) if top.size else float("nan")
            usable = np.isfinite(tm) and tm > 0
            b = tm if usable else 1.0
            entry["fallback_value"] = "top_decile_mean" if usable else "one"
        out[m] = b
        diag[m] = {"used_brightness": b, **entry}
    return out, diag


def estimate_marker_brightness(
    X: np.ndarray,
    var_names: Sequence[str],
    markers: Iterable[str],
    seed: int = 1234,
) -> dict[str, float]:
    """Per-marker brightness b_m = mean of the positive GMM mode of the cosine-normalized channel.

    For markers with no positive population (one BIC component or a non-positive
    positive mean) returns the mean of the marker's own >=90th-percentile values
    (1.0 if that is not positive, e.g. an all-zero channel).
    """
    return estimate_marker_brightness_with_diagnostics(X, var_names, markers, seed)[0]
