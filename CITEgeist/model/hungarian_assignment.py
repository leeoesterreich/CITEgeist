"""Hungarian algorithm for optimal nucleus-to-celltype assignment."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from typing import Dict


def assign_nuclei_to_types(
    probs: np.ndarray,
    counts: np.ndarray,
    nucleus_ids: np.ndarray,
) -> Dict[int, int]:
    """
    Assign nuclei to cell types using Hungarian algorithm.

    Given probability predictions for each nucleus and target counts per type,
    find the optimal assignment that maximizes total probability.

    Args:
        probs: (n_nuclei, n_types) probability matrix from classifier
        counts: (n_types,) integer counts per cell type
        nucleus_ids: (n_nuclei,) unique identifier for each nucleus

    Returns:
        Dict mapping nucleus_id -> cell_type_index
    """
    n_nuclei, n_types = probs.shape
    n_slots = int(counts.sum())

    # Expand counts into slot list
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * int(count))

    # Handle edge case: more nuclei than slots
    if n_slots == 0:
        # No cells to assign - assign each nucleus to its most probable type
        assignments = {}
        for i, nid in enumerate(nucleus_ids):
            assignments[int(nid)] = int(np.argmax(probs[i]))
        return assignments

    # Build cost matrix: -log(prob) for maximum probability matching
    cost_matrix = np.zeros((n_nuclei, max(n_slots, n_nuclei)))

    for i in range(n_nuclei):
        for j in range(n_slots):
            type_idx = slots[j]
            cost_matrix[i, j] = -np.log(probs[i, type_idx] + 1e-10)
        # For padding columns (if n_nuclei > n_slots), use high cost
        for j in range(n_slots, cost_matrix.shape[1]):
            cost_matrix[i, j] = -np.log(probs[i].max() + 1e-10) + 0.01

    # Pad rows if n_slots > n_nuclei
    if n_slots > n_nuclei:
        padding = np.full((n_slots - n_nuclei, cost_matrix.shape[1]), 1e6)
        cost_matrix = np.vstack([cost_matrix, padding])

    # Solve assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Build result dictionary
    assignments = {}
    for i, j in zip(row_ind, col_ind):
        if i < n_nuclei:
            nid = int(nucleus_ids[i])
            if j < n_slots:
                assignments[nid] = slots[j]
            else:
                assignments[nid] = int(np.argmax(probs[i]))

    return assignments
