"""Constrained assignment using Hungarian algorithm."""
import numpy as np
from scipy.optimize import linear_sum_assignment

from .morphology_features import largest_remainder_discretize


def hungarian_assign(log_likes: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """
    Assign samples to types using Hungarian algorithm with count constraints.

    Args:
        log_likes: (N, K) log-likelihoods for each sample and type
        counts: (K,) integer counts for each type

    Returns:
        (N,) assignments - type index for each sample
    """
    N = log_likes.shape[0]
    K = len(counts)

    # Adjust counts to sum to N
    counts = counts.copy().astype(int)
    while counts.sum() < N:
        counts[np.argmax(counts)] += 1
    while counts.sum() > N:
        idx = np.where(counts > 0)[0]
        counts[idx[np.argmin(counts[idx])]] -= 1

    # Build expanded cost matrix
    # Each column represents one "slot" for a cell type
    expanded_cols = []
    col_to_type = []

    for k in range(K):
        for _ in range(int(counts[k])):
            # Negative because we maximize likelihood but minimize cost
            expanded_cols.append(-log_likes[:, k])
            col_to_type.append(k)

    if len(expanded_cols) == 0:
        return np.zeros(N, dtype=int)

    cost_matrix = np.column_stack(expanded_cols)

    # Solve assignment problem
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Map column indices back to type indices
    assignments = np.array([col_to_type[c] for c in col_ind])

    return assignments


def random_assign(counts: np.ndarray, n_samples: int) -> np.ndarray:
    """
    Random assignment respecting count constraints.

    Args:
        counts: (K,) integer counts for each type
        n_samples: total number of samples

    Returns:
        (N,) random assignments
    """
    counts = counts.copy().astype(int)

    # Adjust counts to sum to n_samples
    while counts.sum() < n_samples:
        counts[np.argmax(counts)] += 1
    while counts.sum() > n_samples:
        idx = np.where(counts > 0)[0]
        counts[idx[np.argmin(counts[idx])]] -= 1

    assignments = []
    for k, c in enumerate(counts):
        assignments.extend([k] * int(c))

    np.random.shuffle(assignments)
    return np.array(assignments)


# Alias for API consistency - delegates to the canonical implementation
proportions_to_counts = largest_remainder_discretize
