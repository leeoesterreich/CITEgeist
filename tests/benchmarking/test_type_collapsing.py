# tests/benchmarking/test_type_collapsing.py
"""Tests for deterministic 9→6 cell type collapsing."""
import numpy as np
import pandas as pd
import pytest


def test_collapse_maps_nk_to_tcell():
    """NK_Cell should merge into T_Cell."""
    from Benchmarking.morphology_audit.src.type_collapsing import collapse_cell_types

    labels = pd.Series(["NK_Cell", "T_Cell", "Cancer_Epithelial"])
    result = collapse_cell_types(labels)
    assert result.iloc[0] == "T_Cell"
    assert result.iloc[1] == "T_Cell"
    assert result.iloc[2] == "Cancer_Epithelial"


def test_collapse_maps_bcell_to_plasma():
    """B_Cell should merge into Plasma_Cell."""
    from Benchmarking.morphology_audit.src.type_collapsing import collapse_cell_types

    labels = pd.Series(["B_Cell", "Plasma_Cell"])
    result = collapse_cell_types(labels)
    assert result.iloc[0] == "Plasma_Cell"
    assert result.iloc[1] == "Plasma_Cell"


def test_collapse_drops_mast_cell():
    """Mast_Cell should be dropped (returns NaN or filtered)."""
    from Benchmarking.morphology_audit.src.type_collapsing import collapse_cell_types

    labels = pd.Series(["Mast_Cell", "Cancer_Epithelial", "Fibroblast"])
    result = collapse_cell_types(labels)
    # Mast_Cell maps to None/NaN
    assert pd.isna(result.iloc[0])
    assert result.iloc[1] == "Cancer_Epithelial"
    assert result.iloc[2] == "Fibroblast"


def test_collapse_preserves_major_types():
    """Major types should pass through unchanged."""
    from Benchmarking.morphology_audit.src.type_collapsing import collapse_cell_types

    major = ["Cancer_Epithelial", "Fibroblast", "Plasma_Cell",
             "Macrophage", "T_Cell", "Endothelial"]
    labels = pd.Series(major)
    result = collapse_cell_types(labels)
    for i, t in enumerate(major):
        assert result.iloc[i] == t


def test_collapse_unknown_becomes_nan():
    """Unknown cell types should become NaN."""
    from Benchmarking.morphology_audit.src.type_collapsing import collapse_cell_types

    labels = pd.Series(["Unknown", "Cancer_Epithelial"])
    result = collapse_cell_types(labels)
    assert pd.isna(result.iloc[0])


def test_get_final_type_list():
    """Should return the 6 canonical types in fixed order."""
    from Benchmarking.morphology_audit.src.type_collapsing import FINAL_TYPES

    assert len(FINAL_TYPES) == 6
    assert "Cancer_Epithelial" in FINAL_TYPES
    assert "T_Cell" in FINAL_TYPES
    assert "NK_Cell" not in FINAL_TYPES
    assert "Mast_Cell" not in FINAL_TYPES
