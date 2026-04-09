"""Gate identifiability audit, subtype proportion splitting, and per-cell protein-gate splitting."""

import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


@dataclass
class IdentifiabilityReport:
    gate_variance: Dict[str, float]
    condition_number: float
    effective_rank: int
    collinear_pairs: List[Tuple[str, str]]
    passed: bool


DEFAULT_SUBTYPE_CONFIG = {
    "CD8+ T cells": {
        "gate_marker": "PD-1",
        "subtypes": ["CD8T_exhausted", "CD8T_regular"],
        "min_gate_variance": 0.01,
        "gate_mode": "continuous",
    },
    "Macrophages": {
        "gate_marker": "VISTA",
        "subtypes": ["Mac_VISTA+", "Mac_VISTA-"],
        "min_gate_variance": 0.01,
        "gate_mode": "continuous",
    },
    "Endothelial": {
        "gate_marker": "PTEN",
        "subtypes": ["Endo_PTEN+", "Endo_PTEN-"],
        "min_gate_variance": 0.01,
        "gate_mode": "continuous",
    },
    "Fibroblasts": {
        "gate_marker": "PCNA",
        "subtypes": ["Fib_proliferating", "Fib_quiescent"],
        "min_gate_variance": 0.01,
        "gate_mode": "continuous",
    },
    "Epithelial": {
        "gate_marker": "PCNA",
        "subtypes": ["Epi_proliferating", "Epi_quiescent"],
        "min_gate_variance": 0.01,
        "gate_mode": "continuous",
    },
}


def audit_gate_identifiability(
    parent_props: pd.DataFrame,
    gates: np.ndarray,            # (I, T, M_func)
    type_names: list,
    func_marker_names: list,
    config: dict,
    collinearity_threshold: float = 0.8,
    min_effective_rank: int = 8,
) -> IdentifiabilityReport:
    """Run Phase 0 identifiability checks on gate-split subtype proportions.

    Args:
        parent_props: (n_spots, n_parent_types) DataFrame of parent proportions.
        gates: (n_spots, n_types, n_func_markers) ndarray of gate values in [0,1].
        type_names: List of parent type names matching gates axis 1.
        func_marker_names: List of functional marker names matching gates axis 2.
        config: Subtype config dict mapping parent_type -> gate_marker, subtypes, etc.
        collinearity_threshold: Spearman |r| above which gates are flagged collinear.
        min_effective_rank: Minimum effective rank for go/no-go.

    Returns:
        IdentifiabilityReport with gate variance, condition number, etc.
    """
    # 1. Gate variance per configured type
    gate_variance = {}
    valid_splits = {}  # parent_type -> (type_idx, func_idx)
    for parent_type, cfg in config.items():
        marker = cfg["gate_marker"]
        if parent_type not in type_names or marker not in func_marker_names:
            continue
        t_idx = type_names.index(parent_type)
        m_idx = func_marker_names.index(marker)
        g = gates[:, t_idx, m_idx]
        var_g = float(np.var(g))
        gate_variance[parent_type] = var_g
        if var_g >= cfg.get("min_gate_variance", 0.01):
            valid_splits[parent_type] = (t_idx, m_idx)

    # 2. Build subtype proportion matrix
    subtype_cols = []
    parent_props_np = parent_props.values
    col_list = list(parent_props.columns)

    for col_name in col_list:
        col_idx = col_list.index(col_name)
        if col_name in valid_splits:
            t_idx, m_idx = valid_splits[col_name]
            g = gates[:, t_idx, m_idx]
            p = parent_props_np[:, col_idx]
            subtype_cols.append(p * g)
            subtype_cols.append(p * (1 - g))
        else:
            subtype_cols.append(parent_props_np[:, col_idx])

    subtype_matrix = np.column_stack(subtype_cols)

    # 3. Condition number and effective rank
    sv = np.linalg.svd(subtype_matrix, compute_uv=False)
    condition_number = float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf")
    effective_rank = int(np.sum(sv > 0.01 * sv[0]))

    # 4. Pairwise gate correlation for shared markers (same functional marker index)
    collinear_pairs = []
    split_types = list(valid_splits.keys())
    for i_a in range(len(split_types)):
        for i_b in range(i_a + 1, len(split_types)):
            type_a, type_b = split_types[i_a], split_types[i_b]
            t_a, m_a = valid_splits[type_a]
            t_b, m_b = valid_splits[type_b]
            # Only check types sharing the same functional marker
            if m_a != m_b:
                continue
            g_a = gates[:, t_a, m_a]
            g_b = gates[:, t_b, m_b]
            rho, _ = stats.spearmanr(g_a, g_b)
            if abs(rho) > collinearity_threshold:
                collinear_pairs.append((type_a, type_b))

    passed = (
        effective_rank >= min_effective_rank
        and condition_number < 100
        and len(collinear_pairs) == 0
    )
    return IdentifiabilityReport(
        gate_variance=gate_variance,
        condition_number=condition_number,
        effective_rank=effective_rank,
        collinear_pairs=collinear_pairs,
        passed=passed,
    )


def build_subtype_proportions(
    parent_props: pd.DataFrame,
    gates: np.ndarray,
    type_names: list,
    func_marker_names: list,
    config: dict,
) -> pd.DataFrame:
    """Split parent proportions into subtypes using functional gates.

    Args:
        parent_props: (n_spots, n_parent_types) parent cell type proportions.
        gates: (n_spots, n_types, n_func) gate values in [0, 1].
        type_names: Parent type names (axis 1 of gates).
        func_marker_names: Functional marker names (axis 2 of gates).
        config: Subtype config dict.

    Returns:
        DataFrame (n_spots, n_subtypes) with row sums preserved.
    """
    columns = []
    col_names = []

    for col_name in parent_props.columns:
        p = parent_props[col_name].values

        if col_name not in config:
            columns.append(p)
            col_names.append(col_name)
            continue

        cfg = config[col_name]
        marker = cfg["gate_marker"]
        if col_name not in type_names or marker not in func_marker_names:
            columns.append(p)
            col_names.append(col_name)
            continue

        t_idx = type_names.index(col_name)
        m_idx = func_marker_names.index(marker)
        g = gates[:, t_idx, m_idx].copy()

        # Check gate variance
        if np.var(g) < cfg.get("min_gate_variance", 0.01):
            columns.append(p)
            col_names.append(col_name)
            continue

        # Apply gate mode
        mode = cfg.get("gate_mode", "continuous")
        if mode == "binary_median":
            g = (g >= np.median(g)).astype(float)
        elif mode == "binary_0.5":
            g = (g >= 0.5).astype(float)

        columns.append(p * g)
        col_names.append(cfg["subtypes"][0])
        columns.append(p * (1 - g))
        col_names.append(cfg["subtypes"][1])

    return pd.DataFrame(
        np.column_stack(columns), index=parent_props.index, columns=col_names,
    )


def permute_gates_within_type(gates: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle gate values across spots independently per (type, marker) pair.

    Preserves per-type marginal distributions but destroys spatial structure.
    Used for permutation null in Phase 1 evaluation.
    """
    shuffled = gates.copy()
    _, T, M = shuffled.shape
    for t in range(T):
        for m in range(M):
            rng.shuffle(shuffled[:, t, m])
    return shuffled


# ---------------------------------------------------------------------------
# Protein-gate-based per-cell subtype splitting (Phase 2)
# ---------------------------------------------------------------------------

def split_by_protein_gates(
    cell_assignments: Dict[str, str],
    protein_gates_df: pd.DataFrame,
    proportions: pd.DataFrame,
    cell_spot_map: pd.DataFrame,
    validated_pairs: List[Tuple[str, str]],
    min_subtype_cells: int = 50,
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Split cell types into functional subtypes using per-cell protein gates.

    For each validated (type, marker) pair:
    - Cells of the parent type are re-assigned to ``{type}_{marker}_pos`` or
      ``{type}_{marker}_neg`` based on their binary gate call.
    - Spot-level proportions are recomputed using cell-count ratios within each
      spot.  When a spot has zero cells of type T but non-zero proportion
      (soft QP assignment), the global positive fraction is used as a fallback.

    Pairs where either resulting subtype would have fewer than ``min_subtype_cells``
    cells across the whole region are skipped (the parent type is kept as-is).

    Args:
        cell_assignments: Dict mapping cell_id → parent type name.
        protein_gates_df: DataFrame indexed by cell_id with columns
            ``func_{marker}_{type}_gate`` (0/1 int).  Produced by
            ``run_sace_protein()``.
        proportions: (n_spots, n_types) DataFrame of spot-level proportions.
            Rows = spot barcodes, columns = parent type names.
        cell_spot_map: DataFrame with columns ``cell_id`` and ``spot_id``
            mapping each cell to its Visium spot.
        validated_pairs: List of ``(type_name, marker_name)`` tuples to split.
            Only pairs with a matching gate column in ``protein_gates_df`` and
            sufficient cells in each subtype are applied.
        min_subtype_cells: Minimum number of cells required in both pos and neg
            subtypes to proceed with the split (default 50).

    Returns:
        Tuple of:
            updated_assignments: Dict[cell_id → updated type name].
                Cells from split types are re-labelled; others unchanged.
            updated_proportions: DataFrame (n_spots, n_updated_types) with the
                same row index as ``proportions``.  Columns for split types are
                replaced by two subtype columns; non-split types are preserved.
    """
    # Build a Series for fast spot lookup: cell_id → spot_id
    if "cell_id" in cell_spot_map.columns and "spot_id" in cell_spot_map.columns:
        spot_lookup = cell_spot_map.set_index("cell_id")["spot_id"]
    else:
        raise ValueError("cell_spot_map must have 'cell_id' and 'spot_id' columns")

    updated_assignments = dict(cell_assignments)
    updated_props = proportions.copy()

    for parent_type, marker in validated_pairs:
        gate_col = f"func_{marker}_{parent_type}_gate"
        if gate_col not in protein_gates_df.columns:
            logger.debug("split_by_protein_gates: gate column %s not found; skipping", gate_col)
            continue
        if parent_type not in proportions.columns:
            logger.debug("split_by_protein_gates: parent type %s not in proportions; skipping", parent_type)
            continue

        # Identify cells of this type
        type_cell_ids = [
            cid for cid, ct in cell_assignments.items() if ct == parent_type
        ]
        if not type_cell_ids:
            logger.debug("split_by_protein_gates: no cells of type %s; skipping", parent_type)
            continue

        type_gates = protein_gates_df.loc[
            protein_gates_df.index.intersection(type_cell_ids), gate_col
        ]

        n_pos = int((type_gates == 1).sum())
        n_neg = int((type_gates == 0).sum())
        n_total = n_pos + n_neg

        if n_pos < min_subtype_cells or n_neg < min_subtype_cells:
            logger.info(
                "split_by_protein_gates: skipping (%s, %s) — n_pos=%d, n_neg=%d "
                "(min_subtype_cells=%d)",
                parent_type, marker, n_pos, n_neg, min_subtype_cells,
            )
            continue

        global_pos_frac = n_pos / n_total if n_total > 0 else 0.5
        pos_label = f"{parent_type}_{marker}_pos"
        neg_label = f"{parent_type}_{marker}_neg"

        logger.info(
            "split_by_protein_gates: splitting %s by %s → %s (n=%d) / %s (n=%d)",
            parent_type, marker, pos_label, n_pos, neg_label, n_neg,
        )

        # Update cell assignments
        for cid in type_gates.index:
            updated_assignments[cid] = pos_label if type_gates[cid] == 1 else neg_label

        # Build per-spot cell counts for the split
        cell_ids_series = pd.Series(type_gates.values, index=type_gates.index)
        cell_spots = spot_lookup.reindex(type_gates.index)
        combined = pd.DataFrame({"gate": cell_ids_series, "spot_id": cell_spots})
        combined = combined.dropna(subset=["spot_id"])

        spot_counts = combined.groupby("spot_id")["gate"].agg(
            n_pos_in_spot=lambda x: (x == 1).sum(),
            n_T_in_spot="count",
        )
        spot_counts["pos_frac"] = spot_counts["n_pos_in_spot"] / spot_counts["n_T_in_spot"]

        # Recompute proportions
        parent_col = updated_props[parent_type]
        pos_col = pd.Series(0.0, index=updated_props.index)
        neg_col = pd.Series(0.0, index=updated_props.index)

        for spot_id in updated_props.index:
            p_t = parent_col[spot_id]
            if p_t == 0.0:
                continue
            if spot_id in spot_counts.index:
                pos_frac = spot_counts.loc[spot_id, "pos_frac"]
            else:
                # No cells of this type mapped to this spot — use global fraction
                pos_frac = global_pos_frac

            pos_col[spot_id] = p_t * pos_frac
            neg_col[spot_id] = p_t * (1.0 - pos_frac)

        # Replace parent column with two subtype columns
        insert_pos = updated_props.columns.get_loc(parent_type)
        updated_props = updated_props.drop(columns=[parent_type])
        updated_props.insert(insert_pos, neg_label, neg_col)
        updated_props.insert(insert_pos, pos_label, pos_col)

    return updated_assignments, updated_props
