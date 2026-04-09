"""Unit tests for subtype splitting via protein gates."""
import pytest
import pandas as pd
import numpy as np
from CITEgeist.model.annotation.subtype_splitting import split_by_protein_gates


def _make_inputs(n_cells=200, n_spots=10):
    """Build minimal synthetic inputs for split_by_protein_gates."""
    rng = np.random.default_rng(42)
    cell_ids = [f"cell_{i:04d}" for i in range(n_cells)]
    spot_ids = [f"spot_{i:02d}" for i in range(n_spots)]

    # All cells are Macrophages
    cell_assignments = {c: "Macrophage" for c in cell_ids}

    # Each cell gets a binary PD-L1 gate: 50% positive
    gates = rng.integers(0, 2, size=n_cells)
    protein_gates_df = pd.DataFrame(
        {"func_PD-L1_Macrophage_gate": gates},
        index=cell_ids,
    )

    # Uniform spot proportions
    proportions = pd.DataFrame(
        {"Macrophage": np.ones(n_spots) * 0.3},
        index=spot_ids,
    )

    # Map cells to spots evenly
    cell_spot_map = pd.DataFrame({
        "cell_id": cell_ids,
        "spot_id": [spot_ids[i % n_spots] for i in range(n_cells)],
    })

    validated_pairs = [("Macrophage", "PD-L1")]
    return cell_assignments, protein_gates_df, proportions, cell_spot_map, validated_pairs


def test_split_creates_subtypes():
    assignments, gates_df, props, csmap, pairs = _make_inputs(n_cells=200)
    new_assignments, new_props = split_by_protein_gates(
        assignments, gates_df, props, csmap, pairs, min_subtype_cells=10
    )
    # Should produce _pos and _neg variants
    unique_types = set(new_assignments.values())
    assert "Macrophage_PD-L1_pos" in unique_types or "Macrophage_PD-L1_neg" in unique_types, (
        f"No subtypes created; got: {unique_types}"
    )


def test_split_new_props_shape_matches_original_spots():
    assignments, gates_df, props, csmap, pairs = _make_inputs(n_cells=200)
    _, new_props = split_by_protein_gates(
        assignments, gates_df, props, csmap, pairs, min_subtype_cells=10
    )
    # New props should have same number of rows (spots)
    assert new_props.shape[0] == props.shape[0]


def test_split_skips_when_too_few_cells():
    """min_subtype_cells=10000 forces skip — parent type is kept unchanged."""
    assignments, gates_df, props, csmap, pairs = _make_inputs(n_cells=200)
    new_assignments, _ = split_by_protein_gates(
        assignments, gates_df, props, csmap, pairs, min_subtype_cells=10000
    )
    assert set(new_assignments.values()) == {"Macrophage"}, (
        "Expected parent type preserved when min_subtype_cells not met"
    )


def test_split_with_no_validated_pairs_is_noop():
    assignments, gates_df, props, csmap, _ = _make_inputs(n_cells=200)
    new_assignments, new_props = split_by_protein_gates(
        assignments, gates_df, props, csmap, validated_pairs=[]
    )
    assert new_assignments == assignments
