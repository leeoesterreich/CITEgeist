"""
QC test fixtures.

Provides mock single-cell AnnData, proportions, and ground truth
for testing all QC modules.
"""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData


# ==================== Constants ====================

QC_CELL_TYPES = [
    "Macrophages", "CD8_T_Cells", "CD4_T_Cells", "B_Cells",
    "Epithelial", "Endothelial", "Fibroblasts", "Monocytes",
    "Dendritic_Cells",
]

N_CELLS = 500
N_SPOTS = 100
N_GENES = 300


# ==================== Core Fixtures ====================

@pytest.fixture
def cell_types():
    """Standard QC cell type names."""
    return QC_CELL_TYPES


@pytest.fixture
def mock_sc_adata():
    """Mock single-cell AnnData with raw counts, cell types, and spot IDs.

    Creates 500 cells across 100 spots with 300 genes including
    MT- genes for mitochondrial QC and known canonical markers.
    """
    rng = np.random.default_rng(42)

    # Gene names: include MT genes and some canonical markers
    gene_names = [f"GENE_{i:04d}" for i in range(N_GENES - 10)]
    gene_names += ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-ATP6"]
    gene_names += ["CD68", "CD8A", "EPCAM", "PECAM1", "COL1A1"]
    assert len(gene_names) == N_GENES

    # Assign cells to types and spots
    cell_type_assignments = rng.choice(QC_CELL_TYPES, size=N_CELLS)
    spot_ids = [f"spot_{i % N_SPOTS:04d}" for i in range(N_CELLS)]
    spatial_coords = rng.random((N_CELLS, 2)) * 100

    # Raw count matrix — sparse, ~20% nonzero
    X = sp.random(N_CELLS, N_GENES, density=0.2, format="csr", random_state=42)
    X.data = np.round(X.data * 50).astype(np.float32)

    # Boost MT genes slightly for some cells to create variation
    X_dense = X.toarray()
    mt_idx = [gene_names.index(g) for g in gene_names if g.startswith("MT-")]
    X_dense[:50, mt_idx] *= 5  # First 50 cells have high MT%
    X = sp.csr_matrix(X_dense)

    adata = AnnData(
        X=X,
        obs=pd.DataFrame(
            {"cell_type": cell_type_assignments, "spot_id": spot_ids},
            index=[f"cell_{i:04d}" for i in range(N_CELLS)],
        ),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm["spatial"] = spatial_coords
    return adata


@pytest.fixture
def mock_proportions():
    """Mock spot-level proportions (100 spots × 9 types), sum to 1."""
    rng = np.random.default_rng(42)
    props = rng.dirichlet(np.ones(len(QC_CELL_TYPES)), size=N_SPOTS)
    return pd.DataFrame(
        props,
        columns=QC_CELL_TYPES,
        index=[f"spot_{i:04d}" for i in range(N_SPOTS)],
    )


@pytest.fixture
def mock_gt_proportions():
    """Mock ground truth proportions (slightly different from predicted)."""
    rng = np.random.default_rng(123)
    props = rng.dirichlet(np.ones(len(QC_CELL_TYPES)), size=N_SPOTS)
    return pd.DataFrame(
        props,
        columns=QC_CELL_TYPES,
        index=[f"spot_{i:04d}" for i in range(N_SPOTS)],
    )


@pytest.fixture
def mock_gt_gex_layers():
    """Mock ground truth GEX layers: dict of type → DataFrame (spots × genes)."""
    rng = np.random.default_rng(456)
    gene_names = [f"GENE_{i:04d}" for i in range(50)]  # Smaller gene set for GT
    spot_names = [f"spot_{i:04d}" for i in range(N_SPOTS)]
    layers = {}
    for ct in QC_CELL_TYPES:
        data = rng.lognormal(mean=1, sigma=1, size=(N_SPOTS, len(gene_names)))
        layers[ct] = pd.DataFrame(data, index=spot_names, columns=gene_names)
    return layers
