# test_evaluate_morphology_assignment.py
"""Tests for morphology-guided cell type assignment evaluation."""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from evaluate_morphology_assignment import load_ground_truth, PROTEIN_GT_CELL_TYPES, RNA_GT_CELL_TYPES


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
