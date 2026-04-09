"""Utilities for cellularity-scaled QP deconvolution.

Handles nuclei count preparation, winsorization, and discrete rounding.
"""
import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def round_counts_largest_remainder(
    c: np.ndarray,
    N: int,
) -> np.ndarray:
    """Round continuous cell counts to integers summing to N.

    Uses the largest-remainder (Hamilton) method from apportionment theory.
    Guaranteed closest integer solution to the continuous optimum in L1 norm.

    Args:
        c: Continuous cell counts, shape (T,). Must be non-negative.
        N: Target integer sum (nuclei count for this spot).

    Returns:
        Integer counts, shape (T,), summing to exactly N.
    """
    c = np.asarray(c, dtype=np.float64)
    T = len(c)
    N = int(N)

    if N <= 0:
        return np.zeros(T, dtype=np.int64)

    # Normalize c to sum to N (handle case where solver gives slightly off sum)
    c_sum = c.sum()
    if c_sum > 0:
        c_scaled = c * (N / c_sum)
    else:
        # All zeros: uniform allocation
        c_scaled = np.full(T, N / T)

    # Floor allocation
    floors = np.floor(c_scaled).astype(np.int64)
    remainders = c_scaled - floors

    # Distribute deficit to largest remainders
    deficit = N - floors.sum()
    if deficit > 0:
        # Break ties by original magnitude (stable sort)
        indices = np.argsort(-remainders, kind="stable")
        for k in range(int(deficit)):
            floors[indices[k]] += 1
    elif deficit < 0:
        # Over-allocated (rounding edge case): remove from smallest remainders
        indices = np.argsort(remainders, kind="stable")
        for k in range(int(-deficit)):
            if floors[indices[k]] > 0:
                floors[indices[k]] -= 1

    return floors


def winsorize_cellularity(
    N_raw: np.ndarray,
    percentile: float = 95,
) -> np.ndarray:
    """Winsorize nuclei counts to cap outliers and clamp zeros.

    Args:
        N_raw: Raw nuclei counts, shape (n_spots,).
        percentile: Cap values above this percentile (default: 95).

    Returns:
        Winsorized counts with zeros clamped to 1.
    """
    N = np.asarray(N_raw, dtype=np.float64).copy()
    # Clamp zeros to 1 (avoid division by zero in c/N)
    N = np.maximum(N, 1.0)
    # Cap at percentile
    cap = np.percentile(N, percentile)
    n_capped = int((N > cap).sum())
    N = np.minimum(N, cap)
    if n_capped > 0:
        logger.info(
            "Cellularity winsorized: %d spots capped at N=%.0f (%.0fth percentile)",
            n_capped, cap, percentile,
        )
    return N


def prepare_cellularity(
    nuclei_counts: Optional[pd.Series],
    spot_names: pd.Index,
    percentile: float = 95,
) -> np.ndarray:
    """Prepare cellularity array from optional nuclei counts.

    Aligns to spot names, fills missing with 1, winsorizes outliers.

    Args:
        nuclei_counts: Per-spot nuclei counts (optional). None → all ones.
        spot_names: Spot identifiers to align to.
        percentile: Winsorization percentile.

    Returns:
        Array of shape (n_spots,) with prepared cellularity values.
    """
    if nuclei_counts is None:
        return np.ones(len(spot_names), dtype=np.float64)

    aligned = nuclei_counts.reindex(spot_names, fill_value=0).to_numpy(dtype=np.float64)
    n_missing = int((aligned == 0).sum())
    if n_missing > 0:
        n_total = len(spot_names)
        frac = n_missing / n_total
        if frac > 0.05:
            logger.warning(
                "%.1f%% of spots (%d/%d) have no nuclei counts — falling back to N=1 for those",
                frac * 100, n_missing, n_total,
            )
    return winsorize_cellularity(aligned, percentile=percentile)
