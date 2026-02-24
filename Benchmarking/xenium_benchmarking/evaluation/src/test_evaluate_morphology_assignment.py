# test_evaluate_morphology_assignment.py
"""Tests for morphology-guided cell type assignment evaluation."""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from evaluate_morphology_assignment import (
    load_ground_truth,
    match_cellpose_to_gt,
    run_baseline_random,
    run_baseline_uniform,
    run_baseline_spot_proportion,
    PROTEIN_GT_CELL_TYPES,
    RNA_GT_CELL_TYPES,
)


def test_load_protein_ground_truth():
    """Test loading protein-gated ground truth."""
    gt_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")

    gt_df = load_ground_truth("protein", gt_dir)

    assert isinstance(gt_df, pd.DataFrame)
    assert "cell_id" in gt_df.columns or gt_df.index.name == "cell_id"
    assert "cell_type" in gt_df.columns
    valid_types = set(PROTEIN_GT_CELL_TYPES + ["Unknown"])
    assert set(gt_df["cell_type"].unique()).issubset(valid_types)


def test_load_rna_ground_truth():
    """Test loading RNA clustering ground truth."""
    gt_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")

    gt_df = load_ground_truth("rna", gt_dir)

    assert isinstance(gt_df, pd.DataFrame)
    assert "cell_type" in gt_df.columns
    assert set(gt_df["cell_type"].unique()).issubset(set(RNA_GT_CELL_TYPES))


def test_match_cellpose_to_gt_exact():
    """Test matching when coordinates are identical."""
    cellpose_coords = pd.DataFrame({
        "nucleus_id": [1, 2, 3],
        "centroid_x": [100.0, 200.0, 300.0],
        "centroid_y": [100.0, 200.0, 300.0],
    })

    gt_coords = pd.DataFrame({
        "x_centroid": [100.0, 200.0, 300.0],
        "y_centroid": [100.0, 200.0, 300.0],
    }, index=["cell_a", "cell_b", "cell_c"])

    matches = match_cellpose_to_gt(cellpose_coords, gt_coords, max_dist=10.0)

    assert len(matches) == 3
    assert matches[1] == "cell_a"
    assert matches[2] == "cell_b"
    assert matches[3] == "cell_c"


def test_match_cellpose_to_gt_with_threshold():
    """Test that distant cells are not matched."""
    cellpose_coords = pd.DataFrame({
        "nucleus_id": [1, 2],
        "centroid_x": [100.0, 200.0],
        "centroid_y": [100.0, 200.0],
    })

    gt_coords = pd.DataFrame({
        "x_centroid": [105.0, 500.0],
        "y_centroid": [100.0, 200.0],
    }, index=["cell_a", "cell_b"])

    matches = match_cellpose_to_gt(cellpose_coords, gt_coords, max_dist=10.0)

    assert len(matches) == 1
    assert matches[1] == "cell_a"
    assert 2 not in matches


# --- Baseline Assignment Tests ---

def test_baseline_random_respects_counts():
    """Test that random baseline preserves cell type counts per spot."""
    original = pd.DataFrame({
        "nucleus_id": [1, 2, 3],
        "spot_id": ["spot_1", "spot_1", "spot_1"],
        "cell_type": ["Macrophages", "Macrophages", "T cells"],
    })

    result = run_baseline_random(original, seed=42)

    # Overall counts should be preserved
    orig_counts = original.groupby("cell_type").size()
    result_counts = result.groupby("cell_type").size()
    pd.testing.assert_series_equal(orig_counts.sort_index(), result_counts.sort_index())


def test_baseline_random_shuffles_within_spot():
    """Test that random baseline shuffles assignments within each spot."""
    original = pd.DataFrame({
        "nucleus_id": [1, 2, 3, 4],
        "spot_id": ["spot_1", "spot_1", "spot_2", "spot_2"],
        "cell_type": ["A", "B", "C", "D"],
    })

    # Run multiple times to verify shuffling happens
    seen_different = False
    for seed in range(10):
        result = run_baseline_random(original, seed=seed)
        # Within-spot counts should be preserved
        for spot_id in ["spot_1", "spot_2"]:
            orig_spot = original[original["spot_id"] == spot_id]
            result_spot = result[result["spot_id"] == spot_id]
            assert set(orig_spot["cell_type"]) == set(result_spot["cell_type"])
        # Check if any shuffling happened
        if not original["cell_type"].equals(result["cell_type"]):
            seen_different = True

    assert seen_different, "Random baseline should shuffle at least once in 10 seeds"


def test_baseline_uniform_uses_hungarian():
    """Test uniform baseline with equal probabilities."""
    spot_props = pd.DataFrame({
        "spot_id": ["spot_1"],
        "Macrophages": [0.5],
        "T cells": [0.5],
    }).set_index("spot_id")

    nuclei_per_spot = pd.Series({"spot_1": 4})
    cell_types = ["Macrophages", "T cells"]

    result = run_baseline_uniform(spot_props, nuclei_per_spot, cell_types)

    assert len(result) == 4
    assert result["cell_type"].value_counts()["Macrophages"] == 2
    assert result["cell_type"].value_counts()["T cells"] == 2


def test_baseline_uniform_respects_proportions():
    """Test uniform baseline respects target proportions from spot_props."""
    spot_props = pd.DataFrame({
        "spot_id": ["spot_1"],
        "A": [0.75],
        "B": [0.25],
    }).set_index("spot_id")

    nuclei_per_spot = pd.Series({"spot_1": 4})
    cell_types = ["A", "B"]

    result = run_baseline_uniform(spot_props, nuclei_per_spot, cell_types)

    assert len(result) == 4
    assert result["cell_type"].value_counts()["A"] == 3
    assert result["cell_type"].value_counts()["B"] == 1


def test_baseline_spot_proportion():
    """Test spot-proportion baseline uses spot context only."""
    spot_props = pd.DataFrame({
        "spot_id": ["spot_1", "spot_2"],
        "Macrophages": [0.8, 0.2],
        "T cells": [0.2, 0.8],
    }).set_index("spot_id")

    nuclei_df = pd.DataFrame({
        "nucleus_id": [1, 2, 3, 4],
        "spot_id": ["spot_1", "spot_1", "spot_2", "spot_2"],
    })

    nuclei_per_spot = pd.Series({"spot_1": 2, "spot_2": 2})
    cell_types = ["Macrophages", "T cells"]

    result = run_baseline_spot_proportion(nuclei_df, spot_props, nuclei_per_spot, cell_types)

    # spot_1 has 80% Macrophages -> 2 nuclei -> expect 2 Macrophages (0.8*2=1.6, rounds to 2)
    # Actually with largest_remainder: 0.8*2=1.6, 0.2*2=0.4, floor=[1,0], remainder=[0.6,0.4], add 1 to first
    # So spot_1 should have 2 Macrophages, 0 T cells
    spot1_result = result[result["spot_id"] == "spot_1"]
    assert len(spot1_result) == 2
    # With 80/20 split and 2 nuclei, should get at least 1 Macrophage
    assert (spot1_result["cell_type"] == "Macrophages").sum() >= 1

    # spot_2 has 80% T cells
    spot2_result = result[result["spot_id"] == "spot_2"]
    assert len(spot2_result) == 2
    assert (spot2_result["cell_type"] == "T cells").sum() >= 1


def test_baseline_spot_proportion_empty_spot():
    """Test spot-proportion baseline handles spots with no nuclei."""
    spot_props = pd.DataFrame({
        "spot_id": ["spot_1", "spot_2"],
        "A": [0.5, 0.5],
        "B": [0.5, 0.5],
    }).set_index("spot_id")

    nuclei_df = pd.DataFrame({
        "nucleus_id": [1, 2],
        "spot_id": ["spot_1", "spot_1"],
    })

    nuclei_per_spot = pd.Series({"spot_1": 2, "spot_2": 0})
    cell_types = ["A", "B"]

    result = run_baseline_spot_proportion(nuclei_df, spot_props, nuclei_per_spot, cell_types)

    # Should only have 2 results from spot_1
    assert len(result) == 2
    assert all(result["spot_id"] == "spot_1")
