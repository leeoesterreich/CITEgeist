# test_evaluate_morphology_assignment.py
"""Tests for morphology-guided cell type assignment evaluation."""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from evaluate_morphology_assignment import (
    load_ground_truth,
    match_cellpose_to_gt,
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
