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


def _make_sim_pair(n_side=6, n_genes=40, seed=0):
    """Build separate GEX and antibody AnnDatas on a spot grid, as a simulation run supplies them."""
    rng = np.random.default_rng(seed)
    n_spots = n_side * n_side
    spot_names = pd.Index([f"spot_{i:03d}" for i in range(n_spots)])
    xx, yy = np.meshgrid(np.arange(n_side), np.arange(n_side))
    coords = np.column_stack([xx.ravel(), yy.ravel()]).astype(float) * 100.0

    marker_names = [f"Marker{i}" for i in range(6)]  # simulation data carries no "-1" suffix
    ab = rng.poisson(rng.uniform(5, 50, size=(n_spots, len(marker_names)))).astype(np.float32)
    adata_cite = ad.AnnData(X=ab, obs=pd.DataFrame(index=spot_names.copy()))
    adata_cite.var_names = pd.Index(marker_names)
    adata_cite.obsm["spatial"] = coords

    counts = rng.poisson(rng.uniform(0.5, 20, size=(n_spots, n_genes))).astype(np.float32)
    adata_gex = ad.AnnData(X=counts, obs=pd.DataFrame(index=spot_names.copy()))
    adata_gex.var_names = pd.Index([f"Gene{i}" for i in range(n_genes)])
    adata_gex.obsm["spatial"] = coords.copy()
    return adata_gex, adata_cite, marker_names


def test_cuopt_proportions_non_negative_and_bounded(tmp_path):
    """Run the simulation-mode QP path (CitegeistModel + cuOPT) on tiny synthetic data."""
    from CITEgeist.model.citegeist_model import CitegeistModel

    adata_gex, adata_cite, marker_names = _make_sim_pair()
    profiles = _make_profiles(marker_names, n_types=3)
    model = CitegeistModel(
        sample_name="qp_smoke",
        output_folder=tmp_path,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    model.cell_profile_dict = profiles
    global_props, _finetuned = model.run_cell_proportion_model(method="qp")
    prop = global_props.drop(columns=["recon_error"], errors="ignore")
    assert prop.shape == (adata_cite.n_obs, len(profiles))
    assert list(prop.columns) == list(profiles)
    values = prop.to_numpy(dtype=float)
    assert np.all(values >= -1e-6), "Negative proportions"
    assert np.all(values.sum(axis=1) <= 1.05), "Proportions sum > 1.05"
