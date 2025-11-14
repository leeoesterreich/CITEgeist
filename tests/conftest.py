"""
Pytest configuration and fixtures for CITEgeist tests.

This module provides reusable fixtures and configuration for testing
the CITEgeist spatial transcriptomics deconvolution framework.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

import pytest
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

# Add the parent directory to the system path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# ==================== Directory Fixtures ====================

@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory that is cleaned up after tests."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Cleanup
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


@pytest.fixture
def mock_gurobi_license(temp_output_dir):
    """Create a mock Gurobi license file for testing."""
    license_path = os.path.join(temp_output_dir, "gurobi.lic")
    with open(license_path, 'w') as f:
        f.write("# Mock Gurobi license for testing\n")
    return license_path


# ==================== Data Generation Fixtures ====================

@pytest.fixture
def n_spots():
    """Number of spots for test data."""
    return 100


@pytest.fixture
def n_genes():
    """Number of genes for test data."""
    return 200


@pytest.fixture
def n_proteins():
    """Number of proteins for test data."""
    return 18


@pytest.fixture
def n_celltypes():
    """Number of cell types for test data."""
    return 9


@pytest.fixture
def cell_type_names():
    """Standard cell type names used in CITEgeist."""
    return [
        "B-cells", "CAFs", "Cancer Epithelial", "Endothelial",
        "Myeloid", "Normal Epithelial", "PVL", "Plasmablasts", "T-cells"
    ]


@pytest.fixture
def mock_spatial_coords(n_spots):
    """Generate mock spatial coordinates in a hexagonal grid pattern."""
    np.random.seed(42)
    # Create a rough grid with some noise
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    x = np.repeat(np.arange(grid_size), grid_size)[:n_spots]
    y = np.tile(np.arange(grid_size), grid_size)[:n_spots]

    # Add small random noise
    coords = np.column_stack([
        x * 10 + np.random.randn(n_spots) * 0.5,
        y * 10 + np.random.randn(n_spots) * 0.5
    ])
    return coords


@pytest.fixture
def mock_gene_expression(n_spots, n_genes):
    """Generate mock gene expression count matrix (sparse)."""
    np.random.seed(42)
    # Create sparse count data with ~20% sparsity
    density = 0.8
    data = sp.random(n_spots, n_genes, density=density, format='csr', random_state=42)
    # Scale to count-like integers
    data.data = np.round(data.data * 100).astype(int)
    return data


@pytest.fixture
def mock_protein_expression(n_spots, n_proteins):
    """Generate mock protein/antibody expression matrix."""
    np.random.seed(42)
    # Proteins are typically less sparse than genes
    data = np.random.lognormal(mean=2, sigma=1, size=(n_spots, n_proteins))
    return data


@pytest.fixture
def mock_gene_names(n_genes):
    """Generate mock gene names."""
    return [f"GENE_{i:04d}" for i in range(n_genes)]


@pytest.fixture
def mock_protein_names(cell_type_names):
    """Generate mock protein names matching cell types."""
    proteins = []
    for cell_type in cell_type_names:
        proteins.append(f"{cell_type}_Protein_1")
        proteins.append(f"{cell_type}_Protein_2")
    return proteins


@pytest.fixture
def mock_spot_names(n_spots):
    """Generate mock spot names."""
    return [f"spot_{i:04d}" for i in range(n_spots)]


# ==================== AnnData Fixtures ====================

@pytest.fixture
def mock_gex_adata(mock_gene_expression, mock_gene_names, mock_spot_names, mock_spatial_coords):
    """Create a mock gene expression AnnData object."""
    adata = AnnData(
        X=mock_gene_expression,
        obs=pd.DataFrame(index=mock_spot_names),
        var=pd.DataFrame(index=mock_gene_names)
    )
    adata.var['feature_types'] = 'Gene Expression'
    adata.obsm['spatial'] = mock_spatial_coords
    return adata


@pytest.fixture
def mock_protein_adata(mock_protein_expression, mock_protein_names, mock_spot_names, mock_spatial_coords):
    """Create a mock protein/antibody capture AnnData object."""
    adata = AnnData(
        X=mock_protein_expression,
        obs=pd.DataFrame(index=mock_spot_names),
        var=pd.DataFrame(index=mock_protein_names)
    )
    adata.var['feature_types'] = 'Antibody Capture'
    adata.obsm['spatial'] = mock_spatial_coords
    return adata


@pytest.fixture
def mock_combined_adata(mock_gex_adata, mock_protein_adata):
    """Create a combined AnnData with both gene expression and protein data."""
    import scanpy as sc
    # Concatenate along variable axis
    combined = sc.concat([mock_gex_adata, mock_protein_adata], axis=1)
    return combined


# ==================== Cell Profile Fixtures ====================

@pytest.fixture
def mock_cell_profile_dict(cell_type_names, mock_protein_names):
    """Create a mock cell profile dictionary matching the expected format."""
    profiles = {}
    protein_idx = 0
    for cell_type in cell_type_names:
        profiles[cell_type] = {
            "Major": [mock_protein_names[protein_idx], mock_protein_names[protein_idx + 1]]
        }
        protein_idx += 2
    return profiles


@pytest.fixture
def mock_cell_proportions(n_spots, cell_type_names):
    """Generate mock ground truth cell proportions."""
    np.random.seed(42)
    # Generate random proportions that sum to 1
    props = np.random.dirichlet(np.ones(len(cell_type_names)), size=n_spots)
    df = pd.DataFrame(props, columns=cell_type_names)
    df.index = [f"spot_{i:04d}" for i in range(n_spots)]
    return df


# ==================== Ground Truth Fixtures ====================

@pytest.fixture
def mock_ground_truth_layers(temp_output_dir, n_genes, n_spots, cell_type_names):
    """Create mock ground truth gene expression layers for benchmarking."""
    gt_dir = os.path.join(temp_output_dir, "ground_truth")
    os.makedirs(gt_dir, exist_ok=True)

    np.random.seed(42)
    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    spot_names = [f"spot_{i:04d}" for i in range(n_spots)]

    for cell_type in cell_type_names:
        # Generate random expression data
        data = np.random.lognormal(mean=1, sigma=1, size=(n_genes, n_spots))
        df = pd.DataFrame(data, index=gene_names, columns=spot_names)

        # Save as ground truth
        filepath = os.path.join(gt_dir, f"{cell_type}_GT.csv")
        df.to_csv(filepath)

    return gt_dir


# ==================== Utility Fixtures ====================

@pytest.fixture
def sample_name():
    """Standard sample name for tests."""
    return "test_sample"


@pytest.fixture
def mock_antibody_mapping(mock_protein_names):
    """Create mock antibody to profile mapping."""
    mapping = {}
    for i, protein in enumerate(mock_protein_names):
        mapping[protein] = i
    return mapping


# ==================== Skip Markers ====================

def pytest_configure(config):
    """Configure custom markers."""
    config.addinivalue_line(
        "markers", "unit: Unit tests for individual functions"
    )
    config.addinivalue_line(
        "markers", "integration: Integration tests for full workflows"
    )
    config.addinivalue_line(
        "markers", "slow: Tests that take a long time to run"
    )
    config.addinivalue_line(
        "markers", "requires_gurobi: Tests that require Gurobi license"
    )
    config.addinivalue_line(
        "markers", "requires_data: Tests that require actual data files"
    )


def pytest_collection_modifyitems(config, items):
    """Automatically mark tests based on naming conventions."""
    for item in items:
        # Mark tests in test_integration_*.py as integration tests
        if "integration" in item.nodeid:
            item.add_marker(pytest.mark.integration)

        # Mark tests with 'gurobi' in name as requiring Gurobi
        if "gurobi" in item.nodeid.lower():
            item.add_marker(pytest.mark.requires_gurobi)

        # Mark optimization tests as slow
        if "optimize" in item.nodeid.lower() or "optimization" in item.nodeid.lower():
            item.add_marker(pytest.mark.slow)
