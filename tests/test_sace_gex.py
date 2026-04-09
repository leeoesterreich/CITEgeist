"""Tests for SACE per-cell GEX deconvolution."""
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse import issparse


class TestBuildKernelMatrix:
    """Test spatial kernel matrix construction."""

    @pytest.fixture
    def grid_coords(self):
        """4 spots on a 2x2 grid, spacing=100um."""
        return np.array([[0, 0], [100, 0], [0, 100], [100, 100]], dtype=float)

    def test_returns_sparse_csr(self, grid_coords):
        from CITEgeist.model.gex.sace_gex import build_kernel_matrix
        K = build_kernel_matrix(grid_coords, bandwidth=200.0)
        assert issparse(K)
        assert K.format == "csr"
        assert K.shape == (4, 4)

    def test_rows_sum_to_one(self, grid_coords):
        from CITEgeist.model.gex.sace_gex import build_kernel_matrix
        K = build_kernel_matrix(grid_coords, bandwidth=200.0)
        row_sums = np.array(K.sum(axis=1)).ravel()
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_all_weights_nonnegative(self, grid_coords):
        from CITEgeist.model.gex.sace_gex import build_kernel_matrix
        K = build_kernel_matrix(grid_coords, bandwidth=200.0)
        assert K.min() >= 0.0

    def test_diagonal_is_largest_per_row(self, grid_coords):
        """Self-weight should be >= any neighbor weight."""
        from CITEgeist.model.gex.sace_gex import build_kernel_matrix
        K = build_kernel_matrix(grid_coords, bandwidth=200.0)
        dense = K.toarray()
        for i in range(dense.shape[0]):
            assert dense[i, i] >= dense[i].max() - 1e-10

    def test_auto_bandwidth(self, grid_coords):
        from CITEgeist.model.gex.sace_gex import build_kernel_matrix
        K = build_kernel_matrix(grid_coords, bandwidth=None)
        assert K.shape == (4, 4)
        assert K.nnz > 4

    def test_small_bandwidth_is_sparse(self):
        """Very small bandwidth should mostly keep self-weights."""
        from CITEgeist.model.gex.sace_gex import build_kernel_matrix
        coords = np.array([[0, 0], [1000, 0], [0, 1000]], dtype=float)
        K = build_kernel_matrix(coords, bandwidth=1.0)
        dense = K.toarray()
        np.testing.assert_allclose(dense[0, 1], 0.0, atol=1e-6)


class TestSaceEmInternals:
    """Test E-step allocation and M-step updates."""

    @pytest.fixture
    def simple_problem(self):
        """2 spots, 2 types, 3 genes. Known analytical solution."""
        N, T, G = 2, 2, 3
        Y = np.array([[100, 50, 30],
                       [60, 80, 40]], dtype=float)
        proportions = np.array([[0.8, 0.2],
                                 [0.3, 0.7]])
        mu = np.array([[[0.5, 0.3, 0.2],
                         [0.2, 0.4, 0.4]],
                        [[0.5, 0.3, 0.2],
                         [0.2, 0.4, 0.4]]])
        lib_sizes = Y.sum(axis=1)
        B = proportions * lib_sizes[:, None]
        return {"Y": Y, "proportions": proportions, "mu": mu, "B": B,
                "N": N, "T": T, "G": G}

    def test_e_step_conserves_counts(self, simple_problem):
        from CITEgeist.model.gex.sace_gex import _e_step
        p = simple_problem
        E_x = _e_step(p["Y"], p["B"], p["mu"])
        assert E_x.shape == (2, 2, 3)
        reconstructed = E_x.sum(axis=1)
        np.testing.assert_allclose(reconstructed, p["Y"], atol=1e-10)

    def test_e_step_absent_type_gets_zero(self, simple_problem):
        from CITEgeist.model.gex.sace_gex import _e_step
        p = simple_problem
        B = p["B"].copy()
        B[0, 1] = 0.0
        E_x = _e_step(p["Y"], B, p["mu"])
        np.testing.assert_allclose(E_x[0, 1, :], 0.0, atol=1e-10)
        np.testing.assert_allclose(E_x[0, 0, :], p["Y"][0, :], atol=1e-10)

    def test_m_step_B_equals_allocated_total(self, simple_problem):
        from CITEgeist.model.gex.sace_gex import _e_step, _m_step_B
        p = simple_problem
        E_x = _e_step(p["Y"], p["B"], p["mu"])
        B_new = _m_step_B(E_x)
        expected = E_x.sum(axis=2)
        np.testing.assert_allclose(B_new, expected, atol=1e-10)

    def test_m_step_mu_is_normalized(self, simple_problem):
        from CITEgeist.model.gex.sace_gex import _e_step, _m_step_mu
        p = simple_problem
        E_x = _e_step(p["Y"], p["B"], p["mu"])
        K = sp.eye(2, format="csr")
        mu_new = _m_step_mu(E_x, K, n_0=10, eps=1e-6)
        profile_sums = mu_new.sum(axis=2)
        np.testing.assert_allclose(profile_sums, 1.0, atol=1e-6)

    def test_m_step_mu_global_shrinkage(self):
        """With n_0 very large, local profiles should equal global."""
        from CITEgeist.model.gex.sace_gex import _e_step, _m_step_mu
        N, T, G = 3, 2, 4
        Y = np.random.RandomState(42).poisson(50, size=(N, G)).astype(float)
        mu = np.ones((N, T, G)) / G
        B = np.ones((N, T)) * 50.0
        E_x = _e_step(Y, B, mu)
        K = sp.eye(N, format="csr")
        mu_new = _m_step_mu(E_x, K, n_0=1e8, eps=1e-6)
        for t in range(T):
            for s in range(1, N):
                np.testing.assert_allclose(
                    mu_new[s, t, :], mu_new[0, t, :], atol=1e-4
                )


class TestRunSace:
    """Integration tests for the full SACE EM loop."""

    @pytest.fixture
    def xenium_like_problem(self):
        """Simulated problem: 10 spots, 3 types, 20 genes, ~5 cells/spot."""
        rng = np.random.RandomState(42)
        N, T, G = 10, 3, 20
        n_cells_per_spot = rng.randint(3, 8, size=N)

        Y = rng.poisson(30, size=(N, G)).astype(float)
        raw_props = rng.dirichlet([2, 3, 5], size=N)
        proportions = pd.DataFrame(
            raw_props,
            index=[f"spot_{i}" for i in range(N)],
            columns=["TypeA", "TypeB", "TypeC"],
        )

        cell_ids = []
        cell_spot_indices = []
        cell_types = []
        cell_coords = []
        spot_coords = rng.uniform(0, 1000, size=(N, 2))

        for s in range(N):
            nc = n_cells_per_spot[s]
            counts_per_type = np.round(raw_props[s] * nc).astype(int)
            diff = nc - counts_per_type.sum()
            counts_per_type[np.argmax(counts_per_type)] += diff

            for t_idx, t_name in enumerate(["TypeA", "TypeB", "TypeC"]):
                for _ in range(max(0, counts_per_type[t_idx])):
                    cid = f"cell_{len(cell_ids)}"
                    cell_ids.append(cid)
                    cell_spot_indices.append(s)
                    cell_types.append(t_name)
                    jitter = rng.normal(0, 10, size=2)
                    cell_coords.append(spot_coords[s] + jitter)

        cell_assignments = dict(zip(cell_ids, cell_types))
        cell_spot_map = pd.DataFrame({
            "cell_id": cell_ids,
            "spot_barcode": [f"spot_{i}" for i in cell_spot_indices],
            "spot_idx": cell_spot_indices,
            "x": [c[0] for c in cell_coords],
            "y": [c[1] for c in cell_coords],
        })
        gene_names = [f"gene_{g}" for g in range(G)]

        return {
            "Y": Y, "proportions": proportions,
            "cell_assignments": cell_assignments,
            "cell_spot_map": cell_spot_map,
            "spot_coords": spot_coords,
            "gene_names": gene_names,
            "N": N, "T": T, "G": G,
            "total_cells": len(cell_ids),
        }

    def test_returns_correct_types(self, xenium_like_problem):
        from CITEgeist.model.gex.sace_gex import run_sace
        p = xenium_like_problem
        spot_type_gex, cell_adata, diagnostics = run_sace(
            spot_counts=p["Y"], proportions=p["proportions"],
            cell_assignments=p["cell_assignments"],
            cell_spot_map=p["cell_spot_map"],
            spot_coords=p["spot_coords"], gene_names=p["gene_names"], max_iter=3)
        assert isinstance(spot_type_gex, dict)
        assert isinstance(cell_adata, sc.AnnData)
        assert isinstance(diagnostics, dict)

    def test_spot_type_gex_shapes(self, xenium_like_problem):
        from CITEgeist.model.gex.sace_gex import run_sace
        p = xenium_like_problem
        spot_type_gex, _, _ = run_sace(
            spot_counts=p["Y"], proportions=p["proportions"],
            cell_assignments=p["cell_assignments"],
            cell_spot_map=p["cell_spot_map"],
            spot_coords=p["spot_coords"], gene_names=p["gene_names"], max_iter=3)
        for s_idx, profile in spot_type_gex.items():
            assert profile.shape == (p["T"], p["G"])

    def test_count_conservation(self, xenium_like_problem):
        """Sum of per-cell GEX per spot must equal spot total."""
        from CITEgeist.model.gex.sace_gex import run_sace
        p = xenium_like_problem
        _, cell_adata, _ = run_sace(
            spot_counts=p["Y"], proportions=p["proportions"],
            cell_assignments=p["cell_assignments"],
            cell_spot_map=p["cell_spot_map"],
            spot_coords=p["spot_coords"], gene_names=p["gene_names"], max_iter=3)
        spot_barcodes = cell_adata.obs["spot_barcode"].values
        for s in range(p["N"]):
            spot_name = f"spot_{s}"
            mask = spot_barcodes == spot_name
            if mask.sum() == 0:
                continue
            cell_sum = cell_adata.X[mask].sum(axis=0)
            np.testing.assert_allclose(
                np.asarray(cell_sum).ravel(), p["Y"][s], atol=1e-6,
                err_msg=f"Count conservation violated for {spot_name}")

    def test_cell_adata_has_diagnostics(self, xenium_like_problem):
        from CITEgeist.model.gex.sace_gex import run_sace
        p = xenium_like_problem
        _, cell_adata, _ = run_sace(
            spot_counts=p["Y"], proportions=p["proportions"],
            cell_assignments=p["cell_assignments"],
            cell_spot_map=p["cell_spot_map"],
            spot_coords=p["spot_coords"], gene_names=p["gene_names"], max_iter=3)
        for col in ["allocation_entropy", "n_eff_type", "shrinkage_alpha", "B_st"]:
            assert col in cell_adata.obs.columns, f"Missing diagnostic: {col}"

    def test_convergence_diagnostics(self, xenium_like_problem):
        from CITEgeist.model.gex.sace_gex import run_sace
        p = xenium_like_problem
        _, _, diagnostics = run_sace(
            spot_counts=p["Y"], proportions=p["proportions"],
            cell_assignments=p["cell_assignments"],
            cell_spot_map=p["cell_spot_map"],
            spot_coords=p["spot_coords"], gene_names=p["gene_names"], max_iter=5)
        assert "log_likelihood_trace" in diagnostics
        assert "allocation_change" in diagnostics
        assert len(diagnostics["log_likelihood_trace"]) <= 5

    def test_single_cell_spot(self, xenium_like_problem):
        """Spots with 1 cell should give that cell all counts."""
        from CITEgeist.model.gex.sace_gex import run_sace
        p = xenium_like_problem
        mask = p["cell_spot_map"]["spot_barcode"] == "spot_0"
        first_cell = p["cell_spot_map"][mask].iloc[0]["cell_id"]
        keep = ~mask | (p["cell_spot_map"]["cell_id"] == first_cell)
        csm = p["cell_spot_map"][keep].copy()
        assignments = {k: v for k, v in p["cell_assignments"].items() if k in csm["cell_id"].values}

        _, cell_adata, _ = run_sace(
            spot_counts=p["Y"], proportions=p["proportions"],
            cell_assignments=assignments, cell_spot_map=csm,
            spot_coords=p["spot_coords"], gene_names=p["gene_names"], max_iter=3)
        solo_mask = cell_adata.obs["spot_barcode"] == "spot_0"
        assert solo_mask.sum() == 1
        solo_expr = np.asarray(cell_adata.X[solo_mask]).ravel()
        np.testing.assert_allclose(solo_expr, p["Y"][0], atol=1e-6)


class TestMarkerGuidedInit:
    """Tests for the marker-guided initialization function."""

    def _make_structured_problem(self, rng, N=50, T=3, G=100, M=6):
        """Create structured data where each type has distinct marker-gene correlations."""
        # True gene profiles: each type has a distinct set of marker genes
        true_profiles = np.zeros((T, G))
        true_profiles[0, :30] = 1.0   # Type A -> genes 0-29
        true_profiles[1, 30:60] = 1.0  # Type B -> genes 30-59
        true_profiles[2, 60:90] = 1.0  # Type C -> genes 60-89
        true_profiles += 0.01
        true_profiles /= true_profiles.sum(axis=1, keepdims=True)

        # True antibody profiles: each type has 2 distinct markers
        props = rng.dirichlet([5, 5, 5], size=N)
        lib_sizes = rng.poisson(1000, size=N).astype(float)

        Y = np.zeros((N, G))
        for s in range(N):
            rate = np.zeros(G)
            for t in range(T):
                rate += props[s, t] * lib_sizes[s] * true_profiles[t]
            Y[s] = rng.poisson(rate)

        # Antibody signals correlated with type proportions
        # TypeA -> markers 0,1; TypeB -> markers 2,3; TypeC -> markers 4,5
        antibody_data = np.zeros((N, M))
        antibody_data[:, 0:2] = props[:, 0:1] * 100 + rng.normal(0, 5, (N, 2))
        antibody_data[:, 2:4] = props[:, 1:2] * 100 + rng.normal(0, 5, (N, 2))
        antibody_data[:, 4:6] = props[:, 2:3] * 100 + rng.normal(0, 5, (N, 2))
        antibody_data = np.maximum(antibody_data, 0)

        antibody_names = ["markerA1", "markerA2", "markerB1", "markerB2", "markerC1", "markerC2"]
        cell_profile_dict = {
            "TypeA": {"Major": ["markerA1", "markerA2"]},
            "TypeB": {"Major": ["markerB1", "markerB2"]},
            "TypeC": {"Major": ["markerC1", "markerC2"]},
        }
        type_names = ["TypeA", "TypeB", "TypeC"]

        return {
            "Y": Y, "props": props, "lib_sizes": lib_sizes,
            "antibody_data": antibody_data, "antibody_names": antibody_names,
            "cell_profile_dict": cell_profile_dict, "type_names": type_names,
            "true_profiles": true_profiles, "N": N, "T": T, "G": G,
        }

    def test_returns_correct_shapes(self):
        from CITEgeist.model.gex.sace_gex import _marker_guided_init
        rng = np.random.RandomState(42)
        p = self._make_structured_problem(rng)

        mu_global, diag = _marker_guided_init(
            p["Y"], p["antibody_data"], p["antibody_names"],
            p["cell_profile_dict"], p["type_names"],
        )

        assert mu_global is not None
        assert mu_global.shape == (p["T"], p["G"])
        assert diag["status"] == "ok"

    def test_profiles_are_normalized(self):
        from CITEgeist.model.gex.sace_gex import _marker_guided_init
        rng = np.random.RandomState(42)
        p = self._make_structured_problem(rng)

        mu_global, _ = _marker_guided_init(
            p["Y"], p["antibody_data"], p["antibody_names"],
            p["cell_profile_dict"], p["type_names"],
        )

        row_sums = mu_global.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-6)

    def test_profiles_differ_across_types(self):
        """Marker-guided init should produce distinguishable type profiles."""
        from CITEgeist.model.gex.sace_gex import _marker_guided_init
        rng = np.random.RandomState(42)
        p = self._make_structured_problem(rng)

        mu_global, diag = _marker_guided_init(
            p["Y"], p["antibody_data"], p["antibody_names"],
            p["cell_profile_dict"], p["type_names"],
        )

        # Profiles should have low pairwise correlation (types are distinct)
        assert diag["mean_row_corr"] < 0.95

        # Each type should have highest weight on its marker genes
        assert mu_global[0, :30].mean() > mu_global[0, 30:60].mean()
        assert mu_global[1, 30:60].mean() > mu_global[1, :30].mean()
        assert mu_global[2, 60:90].mean() > mu_global[2, :30].mean()

    def test_missing_markers_falls_back_to_uniform(self):
        """Type with no matching markers should get uniform profile."""
        from CITEgeist.model.gex.sace_gex import _marker_guided_init
        rng = np.random.RandomState(42)
        p = self._make_structured_problem(rng)

        # Remove markers for TypeA
        cell_profile_dict_partial = {
            "TypeA": {"Major": ["nonexistent_marker"]},
            "TypeB": {"Major": ["markerB1", "markerB2"]},
            "TypeC": {"Major": ["markerC1", "markerC2"]},
        }

        mu_global, diag = _marker_guided_init(
            p["Y"], p["antibody_data"], p["antibody_names"],
            cell_profile_dict_partial, p["type_names"],
        )

        # Should not raise; TypeA gets uniform profile but result is valid
        assert mu_global is not None
        assert mu_global.shape == (p["T"], p["G"])

    def test_zero_genes_get_eps_floor(self):
        """Genes with zero counts across all spots should get eps profile."""
        from CITEgeist.model.gex.sace_gex import _marker_guided_init
        rng = np.random.RandomState(42)
        p = self._make_structured_problem(rng)

        Y_sparse = p["Y"].copy()
        Y_sparse[:, 90:] = 0  # last 10 genes zero

        mu_global, diag = _marker_guided_init(
            Y_sparse, p["antibody_data"], p["antibody_names"],
            p["cell_profile_dict"], p["type_names"],
        )

        assert diag["G_nonzero"] == 90
        assert mu_global[:, 90:].max() < 1e-5  # zero genes get eps floor

    def test_run_sace_uses_marker_guided_init(self):
        """run_sace with antibody_data should use marker-guided init path."""
        from CITEgeist.model.gex.sace_gex import run_sace
        rng = np.random.RandomState(42)
        N, T, G, M = 10, 3, 20, 6
        Y = rng.poisson(30, size=(N, G)).astype(float)
        props = pd.DataFrame(
            rng.dirichlet([2, 3, 5], size=N),
            index=[f"spot_{i}" for i in range(N)],
            columns=["TypeA", "TypeB", "TypeC"],
        )
        spot_coords = rng.uniform(0, 1000, size=(N, 2))
        cell_assignments = {f"cell_{i}": "TypeA" for i in range(N)}
        cell_spot_map = pd.DataFrame({
            "cell_id": [f"cell_{i}" for i in range(N)],
            "spot_barcode": [f"spot_{i}" for i in range(N)],
            "spot_idx": list(range(N)),
            "x": spot_coords[:, 0],
            "y": spot_coords[:, 1],
        })
        antibody_data = rng.poisson(10, size=(N, M)).astype(float)
        antibody_names = ["mA1", "mA2", "mB1", "mB2", "mC1", "mC2"]
        cell_profile_dict = {
            "TypeA": {"Major": ["mA1", "mA2"]},
            "TypeB": {"Major": ["mB1", "mB2"]},
            "TypeC": {"Major": ["mC1", "mC2"]},
        }

        # Should not raise -- marker-guided path is exercised
        spot_type_gex, cell_adata, diagnostics = run_sace(
            Y, props, cell_assignments, cell_spot_map,
            spot_coords, [f"gene_{g}" for g in range(G)],
            spotwise_profiles_init=None, max_iter=1,
            antibody_data=antibody_data,
            antibody_names=antibody_names,
            cell_profile_dict=cell_profile_dict,
        )
        assert len(spot_type_gex) == N

    def test_run_sace_falls_back_without_antibody_data(self):
        """run_sace without antibody_data uses confounded proportional init."""
        from CITEgeist.model.gex.sace_gex import run_sace
        rng = np.random.RandomState(42)
        N, T, G = 10, 3, 20
        Y = rng.poisson(30, size=(N, G)).astype(float)
        props = pd.DataFrame(
            rng.dirichlet([2, 3, 5], size=N),
            index=[f"spot_{i}" for i in range(N)],
            columns=["TypeA", "TypeB", "TypeC"],
        )
        spot_coords = rng.uniform(0, 1000, size=(N, 2))
        cell_assignments = {f"cell_{i}": "TypeA" for i in range(N)}
        cell_spot_map = pd.DataFrame({
            "cell_id": [f"cell_{i}" for i in range(N)],
            "spot_barcode": [f"spot_{i}" for i in range(N)],
            "spot_idx": list(range(N)),
            "x": spot_coords[:, 0],
            "y": spot_coords[:, 1],
        })

        # Should not raise -- confounded init fallback is used
        spot_type_gex, cell_adata, diagnostics = run_sace(
            Y, props, cell_assignments, cell_spot_map,
            spot_coords, [f"gene_{g}" for g in range(G)],
            spotwise_profiles_init=None, max_iter=1,
        )
        assert len(spot_type_gex) == N
