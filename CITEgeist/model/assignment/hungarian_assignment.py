"""Hungarian algorithm for optimal nucleus-to-celltype assignment."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from typing import Dict, Optional


def assign_nuclei_to_types(
    probs: np.ndarray,
    counts: np.ndarray,
    nucleus_ids: np.ndarray,
    lambda_prior: float = 0.0,
    proportions: Optional[np.ndarray] = None,
) -> Dict[int, int]:
    """
    Assign nuclei to cell types using Hungarian algorithm.

    Cost matrix combines morphology likelihood and optional proportion prior:
        cost[i,k] = -log(probs[i,k]) - lambda_prior * log(proportions[k])

    Args:
        probs: (n_nuclei, n_types) probability matrix (e.g., MIL attention)
        counts: (n_types,) integer counts per cell type
        nucleus_ids: (n_nuclei,) unique identifier for each nucleus
        lambda_prior: Weight for proportion prior (0 = morphology only)
        proportions: (n_types,) continuous proportions from Module 3.
                     If None, lambda_prior is ignored.

    Returns:
        Dict mapping nucleus_id -> cell_type_index
    """
    n_nuclei, n_types = probs.shape
    n_slots = int(counts.sum())

    # Expand counts into slot list
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * int(count))

    # Coerce nucleus IDs to native Python types (handles both int and str)
    def _coerce_id(nid):
        try:
            return int(nid)
        except (ValueError, TypeError):
            return nid

    # Handle edge case: no cells to assign
    if n_slots == 0:
        assignments = {}
        for i, nid in enumerate(nucleus_ids):
            assignments[_coerce_id(nid)] = int(np.argmax(probs[i]))
        return assignments

    # Compute proportion prior (broadcast across nuclei)
    prop_prior = np.zeros(n_types)
    if proportions is not None and lambda_prior > 0:
        safe_props = np.clip(proportions, 1e-10, None)
        prop_prior = -lambda_prior * np.log(safe_props)

    # Build cost matrix with expanded slots
    cost_matrix = np.zeros((max(n_nuclei, n_slots), max(n_nuclei, n_slots)))
    cost_matrix[:] = 1e6  # High cost for padding

    for i in range(n_nuclei):
        for j in range(n_slots):
            type_idx = slots[j]
            cost_matrix[i, j] = -np.log(probs[i, type_idx] + 1e-10) + prop_prior[type_idx]

    # Solve assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Build result dictionary
    assignments = {}
    for i, j in zip(row_ind, col_ind):
        if i < n_nuclei:
            nid = _coerce_id(nucleus_ids[i])
            if j < n_slots:
                assignments[nid] = slots[j]
            else:
                assignments[nid] = int(np.argmax(probs[i]))

    return assignments
