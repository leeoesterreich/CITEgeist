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


class TestPass1Sparsity:
    """Test that lambda_sparse produces near-one-hot assignments on synthetic single-cell data."""

    def test_sparsity_produces_near_one_hot(self):
        """With lambda_sparse > 0, cells should be assigned predominantly to one type."""
        from CITEgeist.model.gurobi_impl import optimize_cell_proportions_per_marker

        rng = np.random.default_rng(42)
        n_cells = 20
        n_markers = 6
        n_types = 3

        # Simulate 3 cell types with 2 markers each
        # Type 0: markers 0,1 high; Type 1: markers 2,3 high; Type 2: markers 4,5 high
        X = rng.exponential(0.5, (n_cells, n_markers))
        for i in range(n_cells):
            true_type = i % n_types
            X[i, true_type * 2] += 5.0
            X[i, true_type * 2 + 1] += 5.0

        assignment = np.zeros((n_markers, n_types))
        assignment[0, 0] = assignment[1, 0] = 1
        assignment[2, 1] = assignment[3, 1] = 1
        assignment[4, 2] = assignment[5, 2] = 1

        marker_names = [f"M{i}" for i in range(n_markers)]
        type_names = [f"Type_{i}" for i in range(n_types)]

        # With sparsity
        Y_sparse, _, _ = optimize_cell_proportions_per_marker(
            marker_level_data=X,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=type_names,
            lambda_sparse=0.1,
            lambda_laplacian=0.0,
            max_iterations=5,
        )

        # Check near-one-hot: max proportion per cell should be > 0.7
        max_props = Y_sparse.max(axis=1)
        assert np.mean(max_props > 0.7) >= 0.8, (
            f"Expected >=80% cells with max proportion >0.7, got {np.mean(max_props > 0.7):.2f}"
        )

    def test_zero_sparsity_matches_original(self):
        """lambda_sparse=0 should give same results as before."""
        from CITEgeist.model.gurobi_impl import optimize_cell_proportions_per_marker

        rng = np.random.default_rng(42)
        n_spots = 15
        n_markers = 4
        n_types = 2
        X = rng.exponential(2, (n_spots, n_markers))
        assignment = np.zeros((n_markers, n_types))
        assignment[0, 0] = assignment[1, 0] = 1
        assignment[2, 1] = assignment[3, 1] = 1

        Y_no_sparse, _, _ = optimize_cell_proportions_per_marker(
            marker_level_data=X,
            marker_names=[f"M{i}" for i in range(n_markers)],
            assignment_matrix=assignment,
            cell_type_names=["A", "B"],
            lambda_sparse=0.0,
            lambda_laplacian=0.0,
            max_iterations=5,
        )

        # Should be identical to calling without the parameter (default=0)
        Y_default, _, _ = optimize_cell_proportions_per_marker(
            marker_level_data=X,
            marker_names=[f"M{i}" for i in range(n_markers)],
            assignment_matrix=assignment,
            cell_type_names=["A", "B"],
            lambda_laplacian=0.0,
            max_iterations=5,
        )

        np.testing.assert_allclose(Y_no_sparse, Y_default, atol=1e-4)


class TestResolutionDispatch:
    """Test that CitegeistModel dispatches to correct Pass 1/2 based on resolution."""

    def test_cell_mode_passes_sparsity_to_pass1(self, tmp_path):
        """In cell mode, run_cell_proportion_model should pass lambda_sparse from preset."""
        adata_gex, adata_protein = _make_test_adata(n_obs=30, n_genes=20, n_proteins=6)
        from CITEgeist.model import CitegeistModel

        model = CitegeistModel(
            sample_name="test_cell",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="cell",
        )

        # Verify the model stores cell-level params
        assert model.resolution == "cell"
        assert model.resolution_params["lambda_sparse"] > 0

    def test_spot_mode_no_sparsity(self, tmp_path):
        """In spot mode, lambda_sparse should be 0."""
        adata_gex, adata_protein = _make_test_adata(n_obs=30, n_genes=20, n_proteins=6)
        from CITEgeist.model import CitegeistModel

        model = CitegeistModel(
            sample_name="test_spot",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="spot",
        )

        assert model.resolution == "spot"
        assert model.resolution_params["lambda_sparse"] == 0.0


class TestPass2CellLevel:
    """Test cell-level gene expression optimization."""

    def test_basic_true_count_estimation(self):
        """Verify true count estimation recovers dropout for expected genes."""
        from CITEgeist.model.gurobi_impl import estimate_true_expression_cell

        rng = np.random.default_rng(42)
        n_cells = 30
        n_genes = 20
        n_types = 2

        # Simulate cells: type 0 expresses genes 0-9, type 1 expresses genes 10-19
        X_obs = np.zeros((n_cells, n_genes), dtype=np.float64)
        cell_types = np.zeros(n_cells, dtype=int)
        Y_assignments = np.zeros((n_cells, n_types))

        for i in range(n_cells):
            ct = i % n_types
            cell_types[i] = ct
            Y_assignments[i, ct] = 0.95
            Y_assignments[i, 1 - ct] = 0.05

            # Express type-specific genes with some dropout
            type_genes = range(ct * 10, ct * 10 + 10)
            for g in type_genes:
                if rng.random() > 0.3:  # 30% dropout
                    X_obs[i, g] = rng.poisson(10)

        coords = rng.uniform(0, 100, (n_cells, 2))

        # Build enrichment weights: gene g is enriched in its type
        enrichment = np.ones((n_types, n_genes)) * 0.1
        enrichment[0, :10] = 1.0
        enrichment[1, 10:] = 1.0

        X_true = estimate_true_expression_cell(
            X_obs=X_obs,
            Y_assignments=Y_assignments,
            coords=coords,
            enrichment_weights=enrichment,
            library_slack=1.5,
            lambda_spatial=0.01,
            spatial_k=10,
        )

        assert X_true.shape == X_obs.shape

        # X_true should be non-negative
        assert np.all(X_true >= 0)

        # Dropout genes should have some recovery (not all zeros where expected)
        for ct in range(n_types):
            ct_cells = np.where(cell_types == ct)[0]
            type_genes = range(ct * 10, ct * 10 + 10)
            obs_zeros = np.sum(X_obs[np.ix_(ct_cells, list(type_genes))] == 0)
            true_zeros = np.sum(X_true[np.ix_(ct_cells, list(type_genes))] == 0)
            # Should recover at least some dropouts
            assert true_zeros < obs_zeros, (
                f"Type {ct}: expected fewer zeros in X_true ({true_zeros}) than X_obs ({obs_zeros})"
            )

    def test_library_size_bounded(self):
        """X_true total counts should not exceed library_slack * observed."""
        from CITEgeist.model.gurobi_impl import estimate_true_expression_cell

        rng = np.random.default_rng(42)
        n_cells = 10
        n_genes = 10
        X_obs = rng.poisson(5, (n_cells, n_genes)).astype(np.float64)
        Y_assignments = np.zeros((n_cells, 1))
        Y_assignments[:, 0] = 1.0
        coords = rng.uniform(0, 100, (n_cells, 2))
        enrichment = np.ones((1, n_genes))

        slack = 1.5
        X_true = estimate_true_expression_cell(
            X_obs=X_obs,
            Y_assignments=Y_assignments,
            coords=coords,
            enrichment_weights=enrichment,
            library_slack=slack,
            lambda_spatial=0.0,
            spatial_k=5,
        )

        obs_lib = X_obs.sum(axis=1)
        true_lib = X_true.sum(axis=1)
        for i in range(n_cells):
            assert true_lib[i] <= slack * obs_lib[i] + 1.0, (
                f"Cell {i}: true lib {true_lib[i]:.1f} > slack * obs lib {slack * obs_lib[i]:.1f}"
            )
