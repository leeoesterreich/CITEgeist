"""Distribute deconvolved GEX to individual cells."""
import pandas as pd
from typing import Dict


def distribute_gex_to_cells(
    deconvolved_gex: pd.DataFrame,
    assignments: Dict[int, str],
    nucleus_spot_map: pd.DataFrame,
) -> pd.DataFrame:
    """
    Distribute spot-level deconvolved GEX to individual cells.

    Each cell of a given type in a spot receives an equal share of that
    spot-type's deconvolved expression.

    Args:
        deconvolved_gex: DataFrame indexed by 'spot_id:::cell_type' with genes as columns
        assignments: Dict mapping nucleus_id -> cell_type
        nucleus_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns

    Returns:
        DataFrame indexed by nucleus_id with genes as columns
    """
    # Build nucleus info
    nucleus_info = nucleus_spot_map.copy()
    nucleus_info['cell_type'] = nucleus_info['nucleus_id'].map(assignments)

    # Initialize output
    genes = deconvolved_gex.columns
    cell_gex = pd.DataFrame(
        index=nucleus_info['nucleus_id'],
        columns=genes,
        dtype=float
    )
    cell_gex[:] = 0.0

    # Group by spot and cell type
    for (spot_id, cell_type), group in nucleus_info.groupby(['spot_id', 'cell_type']):
        layer_key = f"{spot_id}:::{cell_type}"

        if layer_key not in deconvolved_gex.index:
            continue

        # Get expression for this spot-type
        layer_expr = deconvolved_gex.loc[layer_key]

        # Equal split among cells
        n_cells = len(group)
        per_cell_expr = layer_expr / n_cells

        # Assign to each cell
        for nid in group['nucleus_id']:
            cell_gex.loc[nid] = per_cell_expr.values

    return cell_gex
