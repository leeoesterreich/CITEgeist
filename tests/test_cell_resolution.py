"""Tests for cell-resolution mode in CITEgeist."""
import numpy as np
import pytest
import scanpy as sc


def _make_test_adata(n_obs=100, n_genes=50, n_proteins=10):
    """Create minimal test AnnData with spatial coords."""
    rng = np.random.default_rng(42)
    adata_gex = sc.AnnData(X=rng.poisson(5, (n_obs, n_genes)).astype(np.float64))
    adata_gex.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata_gex.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata_gex.obsm["spatial"] = rng.uniform(0, 1000, (n_obs, 2))

    adata_protein = sc.AnnData(X=rng.exponential(2, (n_obs, n_proteins)).astype(np.float64))
    adata_protein.var_names = [f"Protein_{i}" for i in range(n_proteins)]
    adata_protein.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata_protein.obsm["spatial"] = adata_gex.obsm["spatial"].copy()

    return adata_gex, adata_protein


class TestResolutionPresets:
    """Test resolution parameter and preset loading."""

    def test_default_resolution_is_spot(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        model = CitegeistModel(
            sample_name="test",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
        )
        assert model.resolution == "spot"
        assert model.resolution_params["lambda_sparse"] == 0.0
        assert model.resolution_params["laplacian_k"] == 8

    def test_cell_resolution_loads_preset(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        model = CitegeistModel(
            sample_name="test",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="cell",
        )
        assert model.resolution == "cell"
        assert model.resolution_params["lambda_sparse"] == 0.1
        assert model.resolution_params["laplacian_k"] == 50
        assert model.resolution_params["pass2_library_slack"] == 1.5

    def test_invalid_resolution_raises(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        with pytest.raises(ValueError, match="resolution must be"):
            CitegeistModel(
                sample_name="test",
                output_folder=str(tmp_path),
                simulation=True,
                gene_expression_adata=adata_gex,
                antibody_capture_adata=adata_protein,
                resolution="invalid",
            )

    def test_user_can_override_preset(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        model = CitegeistModel(
            sample_name="test",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="cell",
            resolution_overrides={"lambda_sparse": 0.2, "laplacian_k": 30},
        )
        assert model.resolution_params["lambda_sparse"] == 0.2
        assert model.resolution_params["laplacian_k"] == 30
        assert model.resolution_params["pass2_library_slack"] == 1.5
