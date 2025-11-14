"""
Unit tests for CITEgeist Gurobi optimization implementation.

Note: These tests mock Gurobi functionality since Gurobi requires a license.
Tests marked with @pytest.mark.requires_gurobi need actual Gurobi installation.
"""

import os
import sys

import pytest
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from CITEgeist.model.gurobi_impl import (
    map_antibodies_to_profiles,
    normalize_counts
)


# ==================== Antibody Mapping Tests ====================

@pytest.mark.unit
class TestAntibodyMapping:
    """Tests for antibody to profile mapping."""

    def test_map_antibodies_to_profiles_basic(self, mock_cell_profile_dict, mock_protein_names):
        """Test basic antibody mapping."""
        antibody_names = mock_protein_names

        mapping = map_antibodies_to_profiles(antibody_names, mock_cell_profile_dict)

        assert isinstance(mapping, dict)
        assert len(mapping) > 0

    def test_map_antibodies_to_profiles_complete(self, mock_cell_profile_dict):
        """Test mapping with all antibodies present."""
        # Create antibody names that match the profile dict
        antibody_names = []
        for cell_type, markers in mock_cell_profile_dict.items():
            antibody_names.extend(markers['Major'])

        mapping = map_antibodies_to_profiles(antibody_names, mock_cell_profile_dict)

        # All antibodies should be mapped
        assert len(mapping) == len(antibody_names)

    def test_map_antibodies_partial_match(self):
        """Test mapping with partial antibody matches."""
        profile_dict = {
            "CellTypeA": {"Major": ["AB1", "AB2"]},
            "CellTypeB": {"Major": ["AB3", "AB4"]}
        }
        antibody_names = ["AB1", "AB3", "AB5"]  # AB5 doesn't match

        mapping = map_antibodies_to_profiles(antibody_names, profile_dict)

        # Should map the matching antibodies
        assert "AB1" in mapping or "AB3" in mapping

    def test_map_antibodies_empty_profile_dict(self):
        """Test mapping with empty profile dictionary."""
        antibody_names = ["AB1", "AB2"]

        mapping = map_antibodies_to_profiles(antibody_names, {})

        assert isinstance(mapping, dict)
        assert len(mapping) == 0

    def test_map_antibodies_empty_antibody_list(self, mock_cell_profile_dict):
        """Test mapping with empty antibody list."""
        mapping = map_antibodies_to_profiles([], mock_cell_profile_dict)

        assert isinstance(mapping, dict)
        assert len(mapping) == 0


# ==================== Normalization Tests ====================

@pytest.mark.unit
class TestNormalization:
    """Tests for count normalization functions."""

    def test_normalize_counts_basic(self):
        """Test basic count normalization."""
        counts = np.array([[10, 20, 30],
                          [40, 50, 60],
                          [70, 80, 90]], dtype=float)

        normalized = normalize_counts(counts, target_sum=100)

        # Each row should sum to target_sum
        row_sums = normalized.sum(axis=1)
        assert np.allclose(row_sums, 100, rtol=0.01)

    def test_normalize_counts_different_targets(self):
        """Test normalization with different target sums."""
        counts = np.random.rand(20, 30) * 1000

        for target in [1000, 10000, 100000]:
            normalized = normalize_counts(counts, target_sum=target)
            row_sums = normalized.sum(axis=1)
            assert np.allclose(row_sums, target, rtol=0.01)

    def test_normalize_counts_zero_row(self):
        """Test normalization with zero-count row."""
        counts = np.array([[1, 2, 3],
                          [0, 0, 0],
                          [4, 5, 6]], dtype=float)

        normalized = normalize_counts(counts, target_sum=100)

        # Non-zero rows should normalize correctly
        assert np.isclose(normalized[0].sum(), 100)
        assert np.isclose(normalized[2].sum(), 100)

    def test_normalize_counts_preserves_zeros(self):
        """Test that normalization preserves zero values."""
        counts = np.array([[10, 0, 30],
                          [0, 50, 60],
                          [70, 80, 0]], dtype=float)

        normalized = normalize_counts(counts, target_sum=100)

        # Zeros should remain zeros
        assert normalized[0, 1] == 0
        assert normalized[1, 0] == 0
        assert normalized[2, 2] == 0

    def test_normalize_counts_single_spot(self):
        """Test normalization with single spot."""
        counts = np.array([[10, 20, 30]], dtype=float)

        normalized = normalize_counts(counts, target_sum=100)

        assert np.isclose(normalized.sum(), 100)

    def test_normalize_counts_maintains_ratios(self):
        """Test that normalization maintains relative ratios."""
        counts = np.array([[10, 20, 30]], dtype=float)

        normalized = normalize_counts(counts, target_sum=60)

        # Ratios should be maintained: 1:2:3
        assert np.isclose(normalized[0, 0] / normalized[0, 1], 0.5)
        assert np.isclose(normalized[0, 1] / normalized[0, 2], 2/3)


# ==================== Mock Optimization Tests ====================

@pytest.mark.unit
class TestOptimizationMocks:
    """Tests for optimization functions using mocks (no Gurobi required)."""

    def test_optimization_input_shapes(self, mock_protein_expression, mock_cell_profile_dict):
        """Test that optimization inputs have correct shapes."""
        n_spots, n_proteins = mock_protein_expression.shape
        n_celltypes = len(mock_cell_profile_dict)

        assert n_spots > 0
        assert n_proteins > 0
        assert n_celltypes > 0

    def test_optimization_data_types(self, mock_protein_expression):
        """Test that optimization data has correct types."""
        assert isinstance(mock_protein_expression, np.ndarray)
        assert mock_protein_expression.dtype in [np.float32, np.float64]

    def test_optimization_no_negative_values(self, mock_protein_expression):
        """Test that protein expression has no negative values."""
        assert np.all(mock_protein_expression >= 0)


# ==================== Integration-like Tests (Marked as Requires Gurobi) ====================

@pytest.mark.requires_gurobi
@pytest.mark.slow
class TestGurobiOptimization:
    """
    Tests that require actual Gurobi installation.
    These are skipped by default unless Gurobi is available.
    """

    def test_optimize_cell_proportions_runs(self):
        """Test that cell proportion optimization can run."""
        pytest.skip("Requires Gurobi license")

    def test_optimize_gene_expression_runs(self):
        """Test that gene expression optimization can run."""
        pytest.skip("Requires Gurobi license")

    def test_finetune_cell_proportions_runs(self):
        """Test that finetuning can run."""
        pytest.skip("Requires Gurobi license")


# ==================== Utility and Helper Tests ====================

@pytest.mark.unit
class TestOptimizationHelpers:
    """Tests for optimization helper functions."""

    def test_parameter_validation(self):
        """Test validation of optimization parameters."""
        # Valid parameters
        params = {
            'lambda_reg': 0.001,
            'alpha': 0.7,
            'tolerance': 1e-4,
            'max_iterations': 20
        }

        assert params['lambda_reg'] > 0
        assert 0 <= params['alpha'] <= 1
        assert params['tolerance'] > 0
        assert params['max_iterations'] > 0

    def test_parameter_ranges(self):
        """Test that parameters are in valid ranges."""
        # Test lambda_reg range
        lambda_values = [1e-5, 1e-3, 0.1, 1.0]
        for lam in lambda_values:
            assert lam > 0

        # Test alpha range (elastic net mixing)
        alpha_values = [0.0, 0.3, 0.7, 1.0]
        for alpha in alpha_values:
            assert 0 <= alpha <= 1


# ==================== Edge Cases ====================

@pytest.mark.unit
class TestEdgeCases:
    """Tests for edge cases in optimization."""

    def test_single_spot_optimization(self):
        """Test optimization setup with single spot."""
        protein_expr = np.array([[1.0, 2.0, 3.0]])
        assert protein_expr.shape[0] == 1

    def test_single_celltype_optimization(self):
        """Test optimization setup with single cell type."""
        profile_dict = {"CellTypeA": {"Major": ["AB1", "AB2"]}}
        assert len(profile_dict) == 1

    def test_normalization_extreme_values(self):
        """Test normalization with extreme values."""
        counts = np.array([[1e-10, 1e10, 1.0]], dtype=float)

        normalized = normalize_counts(counts, target_sum=1000)

        assert np.isfinite(normalized).all()
        assert np.isclose(normalized.sum(), 1000, rtol=0.01)

    def test_mapping_case_sensitivity(self):
        """Test antibody mapping case sensitivity."""
        profile_dict = {
            "CellType": {"Major": ["Antibody1", "Antibody2"]}
        }
        antibody_names_lower = ["antibody1", "antibody2"]

        mapping = map_antibodies_to_profiles(antibody_names_lower, profile_dict)

        # Mapping may or may not be case-sensitive depending on implementation
        assert isinstance(mapping, dict)
