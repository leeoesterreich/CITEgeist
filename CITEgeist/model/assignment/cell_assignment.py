"""Discrete cell assignment with optional morphology nudge.

Post-processing pipeline that converts spot-level proportions into
per-cell type assignments using:
  Stage 1: Largest-remainder rounding (proportions -> integer counts)
  Stage 2: Optional morphology scoring from precomputed scores (bayesian only)
  Stage 3: Hungarian assignment (cells -> types within each spot)
"""

__all__ = [
    "assign_cells",
    "bayesian_assign_cells",
    "discretize_proportions",
    "hungarian_assign_spot",
    "assign_cells_to_types",
]

import logging
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from .cellularity_utils import round_counts_largest_remainder

logger = logging.getLogger(__name__)


def discretize_proportions(
    cell_prop: pd.DataFrame,
    nuclei_counts: pd.Series,
    existing_integer_counts: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Discretize continuous proportions to integer cell counts per spot.

    Args:
        cell_prop: (I, T) spot-level proportions. Columns are cell types.
        nuclei_counts: (I,) nuclei count per spot, aligned to cell_prop index.
        existing_integer_counts: If provided, reuse directly (skip rounding).

    Returns:
        (I, T) DataFrame of integer cell counts summing to N_i per spot.
    """
    if existing_integer_counts is not None:
        logger.info("Reusing existing integer counts (%d spots)", len(existing_integer_counts))
        return existing_integer_counts

    cell_types = cell_prop.columns.tolist()
    n_spots = len(cell_prop)
    int_counts = np.zeros((n_spots, len(cell_types)), dtype=np.int64)

    for i, spot_id in enumerate(cell_prop.index):
        N_i = int(nuclei_counts.loc[spot_id])
        if N_i == 0:
            continue
        p_i = cell_prop.loc[spot_id].values.astype(np.float64)
        continuous_counts = N_i * p_i
        int_counts[i] = round_counts_largest_remainder(continuous_counts, N_i)

    return pd.DataFrame(int_counts, index=cell_prop.index, columns=cell_types)


def hungarian_assign_spot(
    n_cells: int,
    integer_counts: np.ndarray,
    cell_types: List[str],
    morph_scores: Optional[np.ndarray],
    morphology_weight: float,
) -> List[str]:
    """Assign n_cells to types for a single spot using Hungarian matching.

    Builds a column list by replicating each type according to its integer
    count, then optionally uses morphology scores to find the minimum-cost
    assignment via the linear sum assignment algorithm.

    Args:
        n_cells: Number of cells in this spot.
        integer_counts: (T,) integer count per cell type (must sum to n_cells
            after reconciliation).
        cell_types: Ordered list of T cell type names.
        morph_scores: (n_cells, T) morphology score matrix, or None.
        morphology_weight: Weight in [0, 1] for morphology scores.  If 0 or
            morph_scores is None, assignment is purely deterministic.

    Returns:
        List of n_cells cell-type name strings.
    """
    if n_cells == 0:
        return []

    # Build ordered target list: [type_0] * k0 + [type_1] * k1 + ...
    col_to_type: List[str] = []
    for t_idx, ct in enumerate(cell_types):
        col_to_type.extend([ct] * int(integer_counts[t_idx]))

    # Reconcile: if rounding produced a different total, trim / pad with the
    # most-represented type.
    while len(col_to_type) < n_cells:
        dominant = cell_types[int(np.argmax(integer_counts))]
        col_to_type.append(dominant)
    col_to_type = col_to_type[:n_cells]

    if morph_scores is not None and morphology_weight > 0.0 and n_cells > 1:
        # Cost matrix: rows = cells, cols = slots in col_to_type.
        # Lower cost = better match → negate the score.
        n_slots = len(col_to_type)
        cost = np.zeros((n_cells, n_slots), dtype=np.float64)
        for slot_idx, ct in enumerate(col_to_type):
            t_idx = cell_types.index(ct)
            cost[:, slot_idx] = -morphology_weight * morph_scores[:, t_idx]

        row_ind, col_ind = linear_sum_assignment(cost)
        assignments = [""] * n_cells
        for r, c in zip(row_ind, col_ind):
            assignments[r] = col_to_type[c]
        return assignments

    # Deterministic (no morphology): assign cells in order.
    return list(col_to_type)


def assign_cells_to_types(
    integer_counts: pd.DataFrame,
    cell_to_spot: np.ndarray,
    cell_types: List[str],
    morph_scores: Optional[np.ndarray],
    morphology_weight: float,
    *,
    spot_proportions: pd.DataFrame,
    cell_ids: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Assign every cell to a cell type via per-spot Hungarian matching.

    Groups cells by spot, reconciles integer counts to the actual number of
    cells observed in each spot, then calls :func:`hungarian_assign_spot`.

    Args:
        integer_counts: (I, T) DataFrame of integer cell counts per spot.
            Index must be spot IDs; columns must be cell type names.
        cell_to_spot: (C,) int array mapping each cell to its spot index
            (0-based position into ``integer_counts.index``).
        cell_types: Ordered list of T cell type names (same as columns of
            ``integer_counts``).
        morph_scores: (C, T) morphology score matrix for all cells, or None.
        morphology_weight: Weight for morphology scores in [0, 1].
        spot_proportions: (I, T) DataFrame of continuous proportions; used
            as confidence fallback when morphology is unavailable.
        cell_ids: Optional (C,) array of cell identifiers.  Defaults to
            0-based integer indices.

    Returns:
        DataFrame with columns: ``spot_id``, ``cell_id``, one column per
        cell type (score), ``assigned_type``, ``confidence``.
    """
    spot_index = integer_counts.index.tolist()  # ordered list of spot IDs
    n_cells_total = len(cell_to_spot)

    if cell_ids is None:
        cell_ids = np.arange(n_cells_total)

    rows = []
    for spot_pos, spot_id in enumerate(spot_index):
        # Mask of cells belonging to this spot
        mask = cell_to_spot == spot_pos
        cell_indices = np.where(mask)[0]
        n_actual = len(cell_indices)

        if n_actual == 0:
            continue

        raw_counts = integer_counts.loc[spot_id, cell_types].values.astype(np.float64)

        # Reconcile: re-round to actual cell count if totals differ
        if int(raw_counts.sum()) != n_actual:
            n_expected = int(raw_counts.sum())
            if n_expected > 0 and abs(n_actual - n_expected) / n_expected > 0.20:
                logger.warning(
                    "Spot %s: nuclei_counts=%d, cells_in_spot=%d (%.0f%% disagreement)",
                    spot_id,
                    n_expected,
                    n_actual,
                    100 * abs(n_actual - n_expected) / n_expected,
                )
            raw_counts = round_counts_largest_remainder(
                raw_counts / max(raw_counts.sum(), 1.0) * n_actual,
                n_actual,
            ).astype(np.float64)

        # Per-spot morphology scores (subset of rows)
        spot_morph = morph_scores[cell_indices] if morph_scores is not None else None

        assignments = hungarian_assign_spot(
            n_cells=n_actual,
            integer_counts=raw_counts,
            cell_types=cell_types,
            morph_scores=spot_morph,
            morphology_weight=morphology_weight,
        )

        # Build score matrix for confidence reporting
        if spot_morph is not None and morphology_weight > 0.0:
            score_matrix = spot_morph  # (n_actual, T)
        else:
            # Broadcast spot proportions to all cells in spot
            prop_row = spot_proportions.loc[spot_id, cell_types].values  # (T,)
            score_matrix = np.tile(prop_row, (n_actual, 1))

        for local_i, (cell_idx, assigned) in enumerate(zip(cell_indices, assignments)):
            t_idx = cell_types.index(assigned)
            confidence = float(score_matrix[local_i, t_idx])
            row = {
                "spot_id": spot_id,
                "cell_id": cell_ids[cell_idx],
            }
            for t_j, ct in enumerate(cell_types):
                row[ct] = float(score_matrix[local_i, t_j])
            row["assigned_type"] = assigned
            row["confidence"] = confidence
            rows.append(row)

    if not rows:
        cols = ["spot_id", "cell_id"] + list(cell_types) + ["assigned_type", "confidence"]
        return pd.DataFrame(columns=cols)

    return pd.DataFrame(rows)


def bayesian_assign_cells(
    morph_scores: np.ndarray,
    cell_to_spot: np.ndarray,
    proportion_prior: np.ndarray,
    detection_mask: np.ndarray,
    *,
    cell_ids: Optional[np.ndarray] = None,
    spot_ids: Optional[list] = None,
    cell_types: Optional[list] = None,
    eps: float = 1e-10,
) -> pd.DataFrame:
    """Assign cells to types via Bayesian posterior with protein detection mask.

    For each cell c in spot s:
        P(type_t | c, s) ∝ morph[c,t] × prior[s,t] × detected[s,t]

    No count constraints — cells are assigned independently.
    Absent types (detected=False) get exactly 0 probability.
    Uses p_base (protein-only) as prior to avoid morphology double-counting.

    Args:
        morph_scores: (C, T) per-cell morphology scores from LLP.
        cell_to_spot: (C,) int spot index per cell.
        proportion_prior: (I, T) spot-level proportions (p_base, protein-only).
        detection_mask: (I, T) boolean mask. True = type present.
        cell_ids: (C,) cell identifiers. Defaults to integer indices.
        spot_ids: List of I spot identifiers.
        cell_types: List of T type names.
        eps: Minimum value floor for prior.

    Returns:
        DataFrame (C rows): spot_id, cell_id, per-type posterior scores,
        assigned_type, confidence.
    """
    C, T = morph_scores.shape
    detection = detection_mask.astype(np.float64)

    # Floor the prior to prevent near-zero suppression of detected rare types
    prior = np.clip(proportion_prior, eps, None)

    posterior = np.zeros((C, T))
    for c in range(C):
        s = cell_to_spot[c]
        posterior[c] = morph_scores[c] * prior[s] * detection[s]
        row_sum = posterior[c].sum()
        if row_sum > eps:
            posterior[c] /= row_sum
        elif detection[s].sum() > 0:
            # Morphology uninformative — fall back to prior over detected types
            fallback = prior[s] * detection[s]
            fallback_sum = fallback.sum()
            posterior[c] = fallback / fallback_sum if fallback_sum > eps else detection[s] / detection[s].sum()
        else:
            # No types detected — uniform over all types
            posterior[c] = 1.0 / T

    assigned_idx = posterior.argmax(axis=1)
    confidence = posterior.max(axis=1)

    if cell_ids is None:
        cell_ids = np.arange(C)
    if cell_types is None:
        cell_types = [f"type_{i}" for i in range(T)]
    if spot_ids is None:
        spot_ids = [f"spot_{i}" for i in range(int(cell_to_spot.max()) + 1)]

    result = pd.DataFrame(posterior, columns=cell_types)
    result.insert(0, "spot_id", [spot_ids[cell_to_spot[c]] for c in range(C)])
    result.insert(1, "cell_id", cell_ids)
    result["assigned_type"] = [cell_types[i] for i in assigned_idx]
    result["confidence"] = confidence

    return result


def assign_cells(
    cell_prop: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_to_spot: np.ndarray,
    *,
    cell_ids: Optional[np.ndarray] = None,
    morphology_weight: float = 0.5,
    existing_integer_counts: Optional[pd.DataFrame] = None,
    output_folder: Optional[str] = None,
    sample_name: str = "sample",
    device: str = "cpu",
    random_state: int = 42,
    assignment_method: str = "hungarian",
    detection_mask: Optional[np.ndarray] = None,
    proportion_prior: Optional[np.ndarray] = None,
    morph_scores_precomputed: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Assign individual cells to types using proportions + optional morphology.

    Three-stage pipeline:
      Stage 1: Discretize proportions -> integer counts per spot
      Stage 2: (bayesian only) morph_scores_precomputed supplies morphology scores
      Stage 3: Hungarian assignment of cells to types within each spot

    Args:
        cell_prop: (I, T) spot-level proportions.
        nuclei_counts: (I,) nuclei count per spot.
        cell_to_spot: (C,) int spot index per cell (positional into cell_prop).
        cell_ids: (C,) cell/nucleus identifiers, or None.
        morphology_weight: Nudge strength (0=none, 1=full). Default 0.5.
        existing_integer_counts: (I, T) pre-computed integer counts to reuse.
        output_folder: Directory for embedding cache.
        sample_name: Sample name for cache key.
        device: Torch device for morphology scoring. Default "cpu".
        random_state: Seed for deterministic assignment. Default 42.
        assignment_method: "hungarian" or "bayesian". Default "hungarian".
        detection_mask: (I, T) boolean mask required for bayesian assignment.
        proportion_prior: (I, T) alternative prior for bayesian assignment.
        morph_scores_precomputed: (C, T) precomputed morphology scores for bayesian.

    Returns:
        DataFrame (C rows): spot_id, cell_id, per-type scores, assigned_type, confidence.

    Raises:
        ValueError: If cell_to_spot indices exceed the number of spots in cell_prop.
    """
    cell_types = cell_prop.columns.tolist()
    n_cells = len(cell_to_spot)

    # Guard against empty cell_to_spot
    if len(cell_to_spot) == 0:
        cols = ["spot_id", "cell_id"] + cell_types + ["assigned_type", "confidence"]
        return pd.DataFrame(columns=cols)

    # Validate alignment
    max_spot_idx = int(cell_to_spot.max())
    if max_spot_idx >= len(cell_prop):
        raise ValueError(f"cell_to_spot max index {max_spot_idx} >= number of spots {len(cell_prop)}")

    logger.info("assign_cells: %d cells, %d spots, %d types", n_cells, len(cell_prop), len(cell_types))

    # --- Stage 1: Discretize ---
    integer_counts = discretize_proportions(cell_prop, nuclei_counts, existing_integer_counts)
    logger.info("Stage 1: Discretized proportions to integer counts")

    # --- Stage 2: Morphology scores (optional, bayesian only) ---
    morph_scores = None

    # --- Bayesian dispatch ---
    if assignment_method == "bayesian":
        effective_morph = morph_scores_precomputed if morph_scores_precomputed is not None else morph_scores
        if effective_morph is None:
            raise ValueError("Provide morph_scores_precomputed (a (C, n_types) array) for bayesian assignment.")
        if detection_mask is None:
            raise ValueError("detection_mask required for bayesian assignment.")
        effective_prior = proportion_prior if proportion_prior is not None else cell_prop.values

        return bayesian_assign_cells(
            morph_scores=effective_morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=effective_prior,
            detection_mask=detection_mask,
            cell_ids=cell_ids,
            spot_ids=cell_prop.index.tolist(),
            cell_types=cell_prop.columns.tolist(),
        )

    # --- Stage 3: Hungarian assignment ---
    result = assign_cells_to_types(
        integer_counts=integer_counts,
        cell_to_spot=cell_to_spot,
        cell_types=cell_types,
        morph_scores=morph_scores,
        morphology_weight=morphology_weight,
        spot_proportions=cell_prop,
        cell_ids=cell_ids,
    )
    logger.info("Stage 3: Assigned %d cells to types", len(result))

    return result
