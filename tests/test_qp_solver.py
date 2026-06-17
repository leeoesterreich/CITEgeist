"""Smoke test for cuOPT QP proportion solver.

Marked requires_cuopt: skipped unless GPU node with cuOPT available.
Tests that the solver runs end-to-end on a tiny synthetic problem and
returns valid proportions (non-negative, sum ≤ 1).
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# Guard: skip entire module if cuOPT is not importable or CUDA is unavailable.
# pytest.importorskip alone is not enough because cuOPT's rmm dependency calls
# cudaGetDevice at import time and raises CUDARuntimeError on CPU-only nodes.
try:
    import cuopt  # noqa: F401

    _CUOPT_AVAILABLE = True
except Exception:
    _CUOPT_AVAILABLE = False

if not _CUOPT_AVAILABLE:
    pytest.skip("cuOPT not installed or CUDA unavailable — requires GPU node", allow_module_level=True)

pytestmark = pytest.mark.requires_cuopt


def _make_synthetic_adata(n_spots=10, n_markers=7, n_types=4, seed=0):
    """Build a minimal AnnData with antibody expression for QP testing."""
    rng = np.random.default_rng(seed)
    marker_names = [f"CD{i+1}-1" for i in range(n_markers)]
    spot_names = [f"spot_{i:02d}" for i in range(n_spots)]
    X = rng.lognormal(0, 1, size=(n_spots, n_markers)).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index(spot_names)
    adata.var_names = pd.Index(marker_names)
    adata.obsm["spatial"] = rng.uniform(0, 100, size=(n_spots, 2))
    return adata, marker_names


def _make_profiles(marker_names, n_types=4):
    """Build a minimal cell_profile_dict in the expected format."""
    type_names = [f"Type{i}" for i in range(n_types)]
    rng = np.random.default_rng(1)
    profiles = {}
    for t in type_names:
        profiles[t] = {"Major": rng.choice(marker_names, size=3, replace=False).tolist()}
    return profiles


def test_map_antibodies_v2_returns_correct_shape():
    from CITEgeist.model.deconvolution.qp_solver import map_antibodies_to_profiles_v2

    adata, marker_names = _make_synthetic_adata()
    profiles = _make_profiles(marker_names)
    marker_data, marker_names_out, assignment_matrix, cell_type_names = map_antibodies_to_profiles_v2(
        adata, profiles
    )
    n_types = len(profiles)
    assert marker_data.shape[0] == adata.n_obs
    assert assignment_matrix.shape == (len(marker_names_out), n_types)
    assert len(cell_type_names) == n_types


def test_cuopt_proportions_non_negative_and_bounded():
    """Run CitegeistModel.run_cell_proportion_model on tiny synthetic data."""
    from CITEgeist.model.citegeist_model import CitegeistModel

    adata, marker_names = _make_synthetic_adata(n_spots=5)
    adata.layers["antibody"] = adata.X.copy()
    adata.layers["counts"] = (adata.X * 10).astype(int)
    profiles = _make_profiles(marker_names, n_types=3)
    model = CitegeistModel(adata, simulation=True)
    result = model.run_cell_proportion_model(profiles=profiles)
    prop = result.proportions.values
    assert np.all(prop >= -1e-6), "Negative proportions"
    assert np.all(prop.sum(axis=1) <= 1.05), "Proportions sum > 1.05"
