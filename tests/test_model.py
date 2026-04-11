"""
Unit tests for CITEgeist model preprocessing and core functionality.

Tests initialization, data splitting, preprocessing, and utility functions.
"""

import os
import sys
import tempfile

import pytest
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from CITEgeist.model.citegeist_model import CitegeistModel


# ==================== Model Initialization Tests ====================

@pytest.mark.unit
class TestModelInitialization:
    """Tests for CitegeistModel initialization."""

    def test_init_simulation_mode(self, temp_output_dir, sample_name, mock_gex_adata, mock_protein_adata):
        """Test initialization in simulation mode."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        assert model.gene_expression_adata is not None
        assert model.antibody_capture_adata is not None
        assert model.adata is None
        assert model.sample_name == sample_name
        assert model.output_folder == temp_output_dir
        assert os.path.exists(temp_output_dir)

    def test_init_non_simulation_mode(self, temp_output_dir, sample_name, mock_combined_adata):
        """Test initialization in non-simulation mode."""
        model = CitegeistModel(
            sample_name=sample_name,
            adata=mock_combined_adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        assert model.adata is not None
        assert model.gene_expression_adata is None
        assert model.antibody_capture_adata is None

    def test_init_simulation_missing_data(self, temp_output_dir, sample_name, mock_gex_adata):
        """Test initialization fails when simulation data is incomplete."""
        with pytest.raises(ValueError, match="both .* must be provided"):
            CitegeistModel(
                sample_name=sample_name,
                output_folder=temp_output_dir,
                simulation=True,
                gene_expression_adata=mock_gex_adata,
                antibody_capture_adata=None
            )

    def test_init_non_simulation_missing_adata(self, temp_output_dir, sample_name):
        """Test initialization fails without adata in non-simulation mode."""
        with pytest.raises(ValueError, match=r"`?adata`?.*must be provided"):
            CitegeistModel(
                sample_name=sample_name,
                output_folder=temp_output_dir,
                simulation=False
            )

    def test_init_missing_output_folder(self, sample_name, mock_gex_adata, mock_protein_adata):
        """Test initialization fails without output folder."""
        with pytest.raises(ValueError, match="output_folder must be provided"):
            CitegeistModel(
                sample_name=sample_name,
                simulation=True,
                gene_expression_adata=mock_gex_adata,
                antibody_capture_adata=mock_protein_adata,
                output_folder=None
            )

    def test_init_creates_output_directory(self, sample_name, mock_gex_adata, mock_protein_adata):
        """Test that initialization creates output directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = os.path.join(temp_dir, "nonexistent", "path")

            model = CitegeistModel(
                sample_name=sample_name,
                output_folder=output_path,
                simulation=True,
                gene_expression_adata=mock_gex_adata,
                antibody_capture_adata=mock_protein_adata
            )

            assert os.path.exists(output_path)

    def test_repr(self, temp_output_dir, sample_name, mock_gex_adata, mock_protein_adata):
        """Test __repr__ method."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        repr_str = repr(model)
        assert "CitegeistModel" in repr_str
        assert "Loaded" in repr_str

    def test_str(self, temp_output_dir, sample_name, mock_gex_adata, mock_protein_adata):
        """Test __str__ method."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        str_repr = str(model)
        assert "CitegeistModel Summary" in str_repr
        assert temp_output_dir in str_repr


# ==================== Data Splitting Tests ====================

@pytest.mark.unit
class TestDataSplitting:
    """Tests for data splitting functionality."""

    def test_split_adata_basic(self, temp_output_dir, sample_name, mock_combined_adata):
        """Test basic AnnData splitting."""
        model = CitegeistModel(
            sample_name=sample_name,
            adata=mock_combined_adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        model.split_adata()

        assert model.gene_expression_adata is not None
        assert model.antibody_capture_adata is not None
        assert model.gene_expression_adata.n_vars > 0
        assert model.antibody_capture_adata.n_vars > 0

    def test_split_adata_correct_features(self, temp_output_dir, sample_name, n_genes, n_proteins):
        """Test that splitting correctly separates gene and protein features."""
        # Create combined data with known features
        n_spots = 50
        gene_data = sp.random(n_spots, n_genes, density=0.8, format='csr')
        protein_data = np.random.rand(n_spots, n_proteins)

        combined_data = np.hstack([gene_data.toarray(), protein_data])
        gene_names = [f"GENE_{i}" for i in range(n_genes)]
        protein_names = [f"PROT_{i}" for i in range(n_proteins)]

        adata = AnnData(X=combined_data)
        adata.var_names = gene_names + protein_names
        adata.var['feature_types'] = ['Gene Expression'] * n_genes + ['Antibody Capture'] * n_proteins

        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        model.split_adata()

        assert model.gene_expression_adata.n_vars == n_genes
        assert model.antibody_capture_adata.n_vars == n_proteins

    def test_split_adata_missing_feature_types(self, temp_output_dir, sample_name):
        """Test splitting fails without feature_types column."""
        adata = AnnData(X=np.random.rand(10, 20))

        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        with pytest.raises(ValueError, match="feature_types.*missing"):
            model.split_adata()

    def test_split_adata_no_gene_expression(self, temp_output_dir, sample_name):
        """Test splitting fails when no Gene Expression features exist."""
        adata = AnnData(X=np.random.rand(10, 20))
        adata.var['feature_types'] = ['Antibody Capture'] * 20

        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        with pytest.raises(ValueError, match="No 'Gene Expression' features"):
            model.split_adata()

    def test_split_adata_no_antibody_capture(self, temp_output_dir, sample_name):
        """Test splitting fails when no Antibody Capture features exist."""
        adata = AnnData(X=np.random.rand(10, 20))
        adata.var['feature_types'] = ['Gene Expression'] * 20

        model = CitegeistModel(
            sample_name=sample_name,
            adata=adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        with pytest.raises(ValueError, match="No 'Antibody Capture' features"):
            model.split_adata()

    def test_split_adata_already_split(self, temp_output_dir, sample_name, mock_combined_adata):
        """Test that splitting fails if data is already split."""
        model = CitegeistModel(
            sample_name=sample_name,
            adata=mock_combined_adata,
            output_folder=temp_output_dir,
            simulation=False
        )

        model.split_adata()

        # Try to split again
        with pytest.raises(ValueError, match="already be split"):
            model.split_adata()


# ==================== Cell Profile Management Tests ====================

@pytest.mark.unit
class TestCellProfileManagement:
    """Tests for cell profile dictionary management."""

    def test_load_cell_profile_dict_valid(self, temp_output_dir, sample_name,
                                          mock_gex_adata, mock_protein_adata,
                                          mock_cell_profile_dict):
        """Test loading valid cell profile dictionary."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        model.load_cell_profile_dict(mock_cell_profile_dict)

        assert model.cell_profile_dict is not None
        assert len(model.cell_profile_dict) > 0

    def test_load_cell_profile_dict_invalid(self, temp_output_dir, sample_name,
                                           mock_gex_adata, mock_protein_adata):
        """Test loading invalid cell profile dictionary."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        with pytest.raises((ValueError, AssertionError)):
            model.load_cell_profile_dict("not a dict")


# ==================== Utility Functions Tests ====================

@pytest.mark.unit
class TestUtilityFunctions:
    """Tests for static utility functions."""

    def test_winsorize_basic(self):
        """Test basic winsorization."""
        matrix = np.array([[1, 2, 3, 4, 100],
                          [1, 2, 3, 4, 5]])

        winsorized = CitegeistModel.winsorize(matrix, lower_percentile=10, upper_percentile=90)

        # Extreme values should be clipped
        assert np.max(winsorized) < 100
        assert winsorized.shape == matrix.shape

    def test_winsorize_no_outliers(self):
        """Test winsorization with no outliers."""
        matrix = np.random.rand(10, 10)

        winsorized = CitegeistModel.winsorize(matrix, lower_percentile=5, upper_percentile=95)

        # Should be similar to original
        assert np.allclose(winsorized, matrix, atol=0.1)

    def test_row_normalize_basic(self):
        """Test basic row normalization."""
        matrix = np.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]], dtype=float)

        normalized = CitegeistModel.row_normalize(matrix, target_sum=100)

        # Each row should sum to target_sum
        row_sums = normalized.sum(axis=1)
        assert np.allclose(row_sums, 100)

    def test_row_normalize_different_targets(self):
        """Test row normalization with different target sums."""
        matrix = np.random.rand(20, 30)

        for target in [1, 10, 1000, 10000]:
            normalized = CitegeistModel.row_normalize(matrix, target_sum=target)
            row_sums = normalized.sum(axis=1)
            assert np.allclose(row_sums, target)

    def test_row_normalize_sparse(self):
        """Test row normalization with sparse matrix."""
        matrix = sp.random(20, 30, density=0.5, format='csr').toarray()

        normalized = CitegeistModel.row_normalize(matrix, target_sum=1000)

        row_sums = normalized.sum(axis=1)
        assert np.allclose(row_sums, 1000)

    def test_global_clr_basic(self):
        """Test global CLR transformation."""
        matrix = np.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]], dtype=float)

        clr_matrix = CitegeistModel.global_clr(matrix)

        # CLR should maintain shape
        assert clr_matrix.shape == matrix.shape
        # Values should be positive
        assert np.all(clr_matrix > 0)

    def test_global_clr_with_zeros(self):
        """Test CLR transformation with zeros."""
        matrix = np.array([[0, 1, 2],
                          [3, 0, 5],
                          [6, 7, 0]], dtype=float)

        # Should handle zeros with epsilon
        clr_matrix = CitegeistModel.global_clr(matrix, epsilon=1e-6)

        assert clr_matrix.shape == matrix.shape
        assert np.all(clr_matrix > 0)

    def test_global_clr_epsilon_effect(self):
        """Test that epsilon prevents division by zero."""
        matrix = np.zeros((5, 5))

        # Should not raise error due to epsilon
        clr_matrix = CitegeistModel.global_clr(matrix, epsilon=1.0)

        assert clr_matrix.shape == matrix.shape
        assert not np.any(np.isnan(clr_matrix))
        assert not np.any(np.isinf(clr_matrix))


# ==================== Gurobi Registration Tests ====================

@pytest.mark.unit
class TestGurobiRegistration:
    """Tests for Gurobi license registration."""

    def test_register_gurobi_valid_file(self, temp_output_dir, sample_name,
                                        mock_gex_adata, mock_protein_adata,
                                        mock_gurobi_license):
        """Test registering valid Gurobi license file."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        # Should not raise error
        model.register_gurobi(mock_gurobi_license)

        assert "GRB_LICENSE_FILE" in os.environ
        assert os.environ["GRB_LICENSE_FILE"] == mock_gurobi_license

    def test_register_gurobi_missing_file(self, temp_output_dir, sample_name,
                                         mock_gex_adata, mock_protein_adata):
        """Test registering nonexistent license file."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        fake_path = "/nonexistent/path/to/license.lic"

        with pytest.raises(FileNotFoundError, match="License file not found"):
            model.register_gurobi(fake_path)


# ==================== Edge Cases ====================

@pytest.mark.unit
class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_model_with_minimal_data(self, temp_output_dir, sample_name):
        """Test model with minimal data (1 spot, 1 gene, 1 protein)."""
        gex_adata = AnnData(X=sp.csr_matrix([[1]]))
        gex_adata.var['feature_types'] = ['Gene Expression']
        gex_adata.obsm['spatial'] = np.array([[0, 0]])

        protein_adata = AnnData(X=np.array([[1.0]]))
        protein_adata.var['feature_types'] = ['Antibody Capture']
        protein_adata.obsm['spatial'] = np.array([[0, 0]])

        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=gex_adata,
            antibody_capture_adata=protein_adata
        )

        assert model is not None

    def test_winsorize_single_value(self):
        """Test winsorization with matrix of single value."""
        matrix = np.ones((10, 10))

        winsorized = CitegeistModel.winsorize(matrix)

        assert np.allclose(winsorized, matrix)

    def test_row_normalize_zero_sum_row(self):
        """Test row normalization with zero-sum row."""
        matrix = np.array([[1, 2, 3],
                          [0, 0, 0],
                          [4, 5, 6]], dtype=float)

        # Should handle zero-sum rows (though might produce nan/inf)
        normalized = CitegeistModel.row_normalize(matrix, target_sum=100)

        # Non-zero rows should normalize correctly
        assert np.isclose(normalized[0].sum(), 100)
        assert np.isclose(normalized[2].sum(), 100)

    def test_model_state_tracking(self, temp_output_dir, sample_name,
                                  mock_gex_adata, mock_protein_adata):
        """Test that model correctly tracks preprocessing state."""
        model = CitegeistModel(
            sample_name=sample_name,
            output_folder=temp_output_dir,
            simulation=True,
            gene_expression_adata=mock_gex_adata,
            antibody_capture_adata=mock_protein_adata
        )

        assert model.preprocessed_gex is False
        assert model.preprocessed_antibody is False
        assert model.results == {}
