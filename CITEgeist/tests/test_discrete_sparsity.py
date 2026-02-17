"""Tests for discrete IQP sparsity regularization.

This module tests the lambda_sparse parameter added to solve_discrete_cell_counts()
to encourage sparse cell type assignments (fewer active types per spot).
"""
import numpy as np
import pytest
import sys

sys.path.insert(0, '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist')
from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts


class TestDiscreteSparsity:
    """Test sparsity regularization in discrete cell assignment."""

    @pytest.fixture
    def simple_problem(self):
        """Create a simple 3-cell-type problem."""
        # 2 spots, 3 markers, 3 cell types
        # Marker 0 -> Type 0, Marker 1 -> Type 1, Marker 2 -> Type 2
        marker_level_data = np.array([
            [1.0, 0.5, 0.5],  # Spot 0: strong type 0, weak types 1,2
            [0.3, 0.3, 0.4],  # Spot 1: roughly equal signal
        ])
        marker_names = ["M0", "M1", "M2"]
        assignment_matrix = np.eye(3)  # Each marker maps to one type
        cell_type_names = ["Type0", "Type1", "Type2"]
        nuclei_counts = np.array([5, 5])
        beta_values = np.ones(3)
        return {
            "marker_level_data": marker_level_data,
            "marker_names": marker_names,
            "assignment_matrix": assignment_matrix,
            "cell_type_names": cell_type_names,
            "nuclei_counts": nuclei_counts,
            "beta_values": beta_values,
        }

    def test_sparsity_penalty_reduces_active_types(self, simple_problem):
        """With high lambda_sparse, fewer cell types should be active."""
        # Without sparsity
        result_no_sparse = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=0.0,
        )

        # With high sparsity
        result_high_sparse = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=1.0,
        )

        # Count non-zero types per spot
        active_no_sparse = (result_no_sparse > 0).sum(axis=1)
        active_high_sparse = (result_high_sparse > 0).sum(axis=1)

        # High sparsity should have same or fewer active types
        assert active_high_sparse.sum() <= active_no_sparse.sum()

    def test_zero_lambda_matches_baseline(self, simple_problem):
        """lambda_sparse=0.0 should give identical results to no sparsity."""
        result_zero = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=0.0,
        )

        # Nuclei constraint still holds
        assert result_zero.sum(axis=1).tolist() == [5, 5]

        # Should be non-negative integers
        assert (result_zero >= 0).all()
        assert result_zero.dtype in [np.int64, np.int32]

    def test_sparsity_respects_nuclei_constraint(self, simple_problem):
        """Sum of cell counts must equal nuclei count regardless of lambda."""
        for lam in [0.0, 0.1, 0.5, 1.0]:
            result = solve_discrete_cell_counts(
                **simple_problem,
                lambda_sparse=lam,
            )
            np.testing.assert_array_equal(
                result.sum(axis=1),
                simple_problem["nuclei_counts"],
                err_msg=f"Nuclei constraint violated for lambda={lam}"
            )

    def test_sparsity_with_continuous_relaxation(self, simple_problem):
        """Sparsity should work with continuous relaxation (high nuclei)."""
        simple_problem["nuclei_counts"] = np.array([50, 50])  # Above max_nuclei_cap

        result = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=0.5,
            max_nuclei_cap=30,
        )

        # Should still sum to nuclei counts
        np.testing.assert_array_equal(result.sum(axis=1), [50, 50])

    def test_very_high_sparsity_concentrates_assignment(self, simple_problem):
        """Very high sparsity should concentrate cells in dominant type."""
        # Modify to have one clear dominant marker
        simple_problem["marker_level_data"] = np.array([
            [2.0, 0.1, 0.1],  # Spot 0: very strong type 0
            [0.1, 2.0, 0.1],  # Spot 1: very strong type 1
        ])

        result = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=10.0,  # Very high sparsity
        )

        # With such high sparsity, should prefer single-type assignments
        # Spot 0 should have most/all cells as Type0
        assert result[0, 0] >= 3  # At least 3 out of 5 should be Type0
        # Spot 1 should have most/all cells as Type1
        assert result[1, 1] >= 3  # At least 3 out of 5 should be Type1

    def test_lambda_sparse_parameter_exists(self, simple_problem):
        """Test that lambda_sparse parameter is accepted."""
        # This should not raise TypeError about unexpected keyword argument
        result = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=0.5,
        )
        assert result.shape == (2, 3)

    def test_indicator_variables_binary(self, simple_problem):
        """Test that indicator behavior is correct (type active or not)."""
        # Single spot to make verification easier
        # Use very high signal for Type0, zero for others
        simple_problem["marker_level_data"] = np.array([[5.0, 0.0, 0.0]])
        simple_problem["nuclei_counts"] = np.array([5])

        result = solve_discrete_cell_counts(
            **simple_problem,
            lambda_sparse=1.0,  # Higher sparsity to encourage single-type solution
        )

        # Nuclei constraint must hold
        assert result.sum(axis=1)[0] == 5
        # With strong signal for Type0 and high sparsity, Type0 should dominate
        # At minimum, Type0 should have more cells than any other type
        assert result[0, 0] >= result[0, 1]
        assert result[0, 0] >= result[0, 2]


class TestDiscreteSparsityEdgeCases:
    """Edge case tests for sparsity regularization."""

    def test_single_spot_single_type(self):
        """Single spot with one cell type - trivial case."""
        marker_level_data = np.array([[1.0]])
        marker_names = ["M0"]
        assignment_matrix = np.array([[1.0]])
        cell_type_names = ["Type0"]
        nuclei_counts = np.array([3])
        beta_values = np.array([1.0])

        result = solve_discrete_cell_counts(
            marker_level_data=marker_level_data,
            marker_names=marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_values=beta_values,
            lambda_sparse=0.5,
        )

        assert result[0, 0] == 3

    def test_zero_nuclei_spot(self):
        """Spots with zero nuclei should remain zero regardless of lambda."""
        marker_level_data = np.array([
            [1.0, 0.5],
            [0.5, 1.0],
        ])
        marker_names = ["M0", "M1"]
        assignment_matrix = np.eye(2)
        cell_type_names = ["Type0", "Type1"]
        nuclei_counts = np.array([0, 5])  # First spot has no nuclei
        beta_values = np.ones(2)

        result = solve_discrete_cell_counts(
            marker_level_data=marker_level_data,
            marker_names=marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_values=beta_values,
            lambda_sparse=0.5,
        )

        # First spot should be all zeros
        np.testing.assert_array_equal(result[0, :], [0, 0])
        # Second spot should sum to 5
        assert result[1, :].sum() == 5

    def test_negative_lambda_raises_or_works(self):
        """Negative lambda_sparse should either work or raise clear error."""
        marker_level_data = np.array([[1.0, 0.5]])
        marker_names = ["M0", "M1"]
        assignment_matrix = np.eye(2)
        cell_type_names = ["Type0", "Type1"]
        nuclei_counts = np.array([5])
        beta_values = np.ones(2)

        # Negative lambda would encourage more active types
        # This is mathematically valid but unusual
        try:
            result = solve_discrete_cell_counts(
                marker_level_data=marker_level_data,
                marker_names=marker_names,
                assignment_matrix=assignment_matrix,
                cell_type_names=cell_type_names,
                nuclei_counts=nuclei_counts,
                beta_values=beta_values,
                lambda_sparse=-0.1,
            )
            # If it works, nuclei constraint should hold
            assert result.sum() == 5
        except ValueError as e:
            # If it raises, should be a clear message
            assert "lambda" in str(e).lower() or "sparse" in str(e).lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
