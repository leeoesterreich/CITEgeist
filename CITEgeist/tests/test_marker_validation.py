"""Tests for single-cell marker gene validation."""
import numpy as np
import pandas as pd
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from model.marker_validation import compute_marker_scores, summarize_validation


def test_compute_marker_scores_positive_case():
    gex_data = {
        "spot_A:::Macrophages": {"CD68": 5.0, "EPCAM": 0.1, "VIM": 0.2},
        "spot_A:::Cancer": {"CD68": 0.1, "EPCAM": 4.0, "VIM": 0.3},
    }
    gex_df = pd.DataFrame(gex_data).T

    assignments = pd.DataFrame({
        "nucleus_id": ["n1", "n2"],
        "barcode": ["spot_A", "spot_A"],
        "cell_type": ["Macrophages", "Cancer"],
    })

    rna_markers = {"Macrophages": ["CD68"], "Cancer": ["EPCAM"]}
    scores = compute_marker_scores(assignments, gex_df, rna_markers)

    assert len(scores) == 2
    assert scores.loc[scores["nucleus_id"] == "n1", "marker_score"].values[0] > 0
    assert scores.loc[scores["nucleus_id"] == "n2", "marker_score"].values[0] > 0


def test_compute_marker_scores_negative_case():
    """Test that spots with low assigned-type markers score negatively.

    Simulates a spot where the deconvolved Cancer row has very low EPCAM
    (assigned type marker) but high CD68 (other type marker), as would happen
    when the assigned type has minimal presence in that spot.
    """
    gex_data = {
        # Cancer deconvolved row: low EPCAM (Cancer marker), high CD68 (Macrophage marker)
        "spot_A:::Cancer": {"CD68": 5.0, "EPCAM": 0.1},
        # Macrophages deconvolved row: high CD68 (Macrophage marker), low EPCAM
        "spot_A:::Macrophages": {"CD68": 4.0, "EPCAM": 0.2},
    }
    gex_df = pd.DataFrame(gex_data).T

    # Nucleus assigned to Cancer, but Cancer row shows low EPCAM -> should score negative
    assignments = pd.DataFrame({
        "nucleus_id": ["n1"],
        "barcode": ["spot_A"],
        "cell_type": ["Cancer"],
    })

    rna_markers = {"Macrophages": ["CD68"], "Cancer": ["EPCAM"]}
    scores = compute_marker_scores(assignments, gex_df, rna_markers)

    assert len(scores) == 1
    # Cancer marker (EPCAM=0.1) < other marker mean (CD68=5.0 from Macrophages row via Cancer GEX row)
    assert scores.loc[scores["nucleus_id"] == "n1", "marker_score"].values[0] < 0
    assert not scores.loc[scores["nucleus_id"] == "n1", "markers_above_others"].values[0]


def test_summarize_validation():
    scores_df = pd.DataFrame({
        "nucleus_id": ["n1", "n2", "n3", "n4"],
        "barcode": ["s1", "s1", "s2", "s2"],
        "cell_type": ["Cancer", "Cancer", "Macrophages", "Macrophages"],
        "marker_score": [2.0, 1.5, 3.0, -0.5],
        "markers_above_others": [True, True, True, False],
    })

    summary = summarize_validation(scores_df)
    assert "per_type" in summary
    assert "overall" in summary
    assert summary["per_type"]["Cancer"]["fraction_correct"] == 1.0
    assert summary["per_type"]["Macrophages"]["fraction_correct"] == 0.5
