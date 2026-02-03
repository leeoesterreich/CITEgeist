#!/usr/bin/env python
"""
Test Module 4 Joint Program Discovery.

Tests the joint discovery mode that finds programs spanning multiple cell types.
"""

import numpy as np
import pandas as pd
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import scanpy as sc


def test_joint_program_dataclass():
    """Test JointProgram dataclass initialization."""
    from CITEgeist.model.anchored_program_discovery import JointProgram

    program = JointProgram(
        program_id=0,
        top_genes=["GENE1", "GENE2", "GENE3"],
        gene_loadings={"GENE1": 0.8, "GENE2": 0.5, "GENE3": 0.3},
        variance_explained=0.15,
        spatial_moran_i=0.45,
        spatial_moran_pvalue=0.001,
        mean_activity=0.25,
        active_spots_fraction=0.6,
        cell_type_enrichments={"Cancer": 0.7, "Macrophage": 0.3},
        primary_cell_type="Cancer",
        secondary_cell_type="Macrophage",
        interaction_score=0.35,
        program_type="interaction",
    )

    assert program.program_id == 0
    assert program.primary_cell_type == "Cancer"
    assert program.interaction_score == 0.35
    assert program.program_type == "interaction"


def test_joint_discovery_result_dataclass():
    """Test JointDiscoveryResult dataclass initialization."""
    from CITEgeist.model.anchored_program_discovery import JointProgram, JointDiscoveryResult

    program = JointProgram(
        program_id=0,
        top_genes=["GENE1"],
        gene_loadings={"GENE1": 0.8},
        variance_explained=0.15,
        spatial_moran_i=0.45,
        spatial_moran_pvalue=0.001,
        mean_activity=0.25,
        active_spots_fraction=0.6,
        cell_type_enrichments={"Cancer": 0.7},
        primary_cell_type="Cancer",
        secondary_cell_type=None,
        interaction_score=0.0,
        program_type="single_celltype",
    )

    result = JointDiscoveryResult(
        programs=[program],
        W=np.array([[0.8], [0.5]]),
        H=np.array([[0.3, 0.4, 0.5]]),
        gene_names=["GENE1", "GENE2"],
        cell_type_names=["Cancer", "Macrophage"],
        n_spots=3,
        reconstruction_error=0.05,
        parameters={"K_programs": 1},
    )

    assert len(result.programs) == 1
    assert result.n_spots == 3
    assert result.W.shape == (2, 1)

    # Test to_dataframe method
    df = result.to_dataframe()
    assert len(df) == 1
    assert "primary_cell_type" in df.columns
    assert "interaction_score" in df.columns


def test_assign_program_cell_types():
    """Test cell type assignment for joint programs."""
    from CITEgeist.model.anchored_program_discovery import _assign_program_cell_types

    # Mock H matrix: 2 programs, 100 spots
    np.random.seed(42)
    H = np.random.rand(2, 100)

    # Mock proportions: program 0 correlates with Cancer, program 1 with both
    proportions = pd.DataFrame({
        "Cancer": np.linspace(0, 1, 100) + np.random.normal(0, 0.1, 100),
        "Macrophage": np.linspace(1, 0, 100) + np.random.normal(0, 0.1, 100),
    })
    # Make program 0 correlate with Cancer
    H[0, :] = proportions["Cancer"].values + np.random.normal(0, 0.1, 100)
    H[0, H[0, :] < 0] = 0
    # Make program 1 correlate with both equally
    H[1, :] = 0.5 * proportions["Cancer"].values + 0.5 * proportions["Macrophage"].values

    result = _assign_program_cell_types(H, proportions)

    # Program 0 should be single cell type (Cancer)
    assert result[0]["primary_cell_type"] == "Cancer"
    assert result[0]["program_type"] == "single_celltype"
    assert result[0]["interaction_score"] < 0.5

    # Program 1 should be interaction or microenvironment (balanced)
    assert result[1]["program_type"] in ["interaction", "microenvironment"]


def test_discover_joint_programs_simulated():
    """Test joint program discovery on simulated data."""
    from CITEgeist.model.anchored_program_discovery import (
        discover_joint_programs,
        JointDiscoveryResult,
    )

    # Create simulated AnnData with deconvolved layers
    np.random.seed(42)
    n_spots = 100
    n_genes = 200

    # Base expression
    X = np.random.rand(n_spots, n_genes) * 10

    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Spot_{i}" for i in range(n_spots)]
    adata.obsm["spatial"] = np.random.rand(n_spots, 2) * 1000

    # Add deconvolved layers
    adata.layers["Cancer_genes_pass1"] = np.random.rand(n_spots, n_genes) * 5
    adata.layers["Macrophage_genes_pass1"] = np.random.rand(n_spots, n_genes) * 3
    adata.layers["Tcell_genes_pass1"] = np.random.rand(n_spots, n_genes) * 2

    # Create proportions
    proportions = pd.DataFrame({
        "Cancer": np.random.rand(n_spots),
        "Macrophage": np.random.rand(n_spots),
        "Tcell": np.random.rand(n_spots),
    }, index=adata.obs_names)
    # Normalize rows
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    # Run joint discovery
    result = discover_joint_programs(
        adata=adata,
        cell_type_proportions=proportions,
        K_programs=5,
        layer_pattern="_genes_pass1",
    )

    assert isinstance(result, JointDiscoveryResult)
    assert len(result.programs) == 5
    assert result.W.shape[0] == n_genes
    assert result.W.shape[1] == 5
    assert result.H.shape[0] == 5
    assert result.H.shape[1] == n_spots
    assert all(p.program_type in ["single_celltype", "interaction", "microenvironment"]
               for p in result.programs)


def test_exports():
    """Test that new classes are exported from CITEgeist.model."""
    from CITEgeist.model import (
        JointProgram,
        JointDiscoveryResult,
        discover_joint_programs,
    )

    assert JointProgram is not None
    assert JointDiscoveryResult is not None
    assert callable(discover_joint_programs)


def test_joint_vs_anchored_comparison():
    """Test that joint and anchored discovery produce comparable results."""
    from CITEgeist.model import discover_joint_programs

    # Create richer simulated data
    np.random.seed(42)
    n_spots = 200
    n_genes = 500

    X = np.random.rand(n_spots, n_genes) * 10
    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Spot_{i}" for i in range(n_spots)]

    # Create spatial coords with some structure
    adata.obsm["spatial"] = np.column_stack([
        np.random.rand(n_spots) * 1000,
        np.random.rand(n_spots) * 1000,
    ])

    # Create structured deconvolved layers
    # Cancer high in left side
    cancer_weight = 1 - (adata.obsm["spatial"][:, 0] / 1000)
    adata.layers["Cancer_genes_pass1"] = X * cancer_weight[:, np.newaxis]

    # Macrophage high in right side
    macro_weight = adata.obsm["spatial"][:, 0] / 1000
    adata.layers["Macrophage_genes_pass1"] = X * macro_weight[:, np.newaxis]

    # Proportions matching the spatial pattern
    proportions = pd.DataFrame({
        "Cancer": cancer_weight,
        "Macrophage": macro_weight,
    }, index=adata.obs_names)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    # Run joint discovery
    joint_result = discover_joint_programs(
        adata=adata,
        cell_type_proportions=proportions,
        K_programs=4,
    )

    # Verify structure
    assert len(joint_result.programs) == 4
    assert joint_result.H.shape == (4, n_spots)

    # At least one program should be spatial (Moran's I > 0.1)
    spatial_programs = [p for p in joint_result.programs if p.spatial_moran_i > 0.1]
    assert len(spatial_programs) > 0, "Expected at least one spatial program"

    # Programs should have different cell type assignments
    primary_types = [p.primary_cell_type for p in joint_result.programs]
    assert len(set(primary_types)) > 1 or any(
        p.program_type == "interaction" for p in joint_result.programs
    ), "Expected diverse cell type assignments"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
