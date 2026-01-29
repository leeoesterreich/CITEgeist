"""Tests for marker exclusivity weighting."""
import inspect

import numpy as np
import pytest

from CITEgeist.model.gurobi_impl import (
    compute_marker_exclusivity,
    deconvolute_local_cell_proportions_per_marker,
    finetune_cell_proportions_per_marker,
)


class TestComputeMarkerExclusivity:
    """Test exclusivity score computation."""

    def _make_data(self, N=200, seed=42):
        """Create synthetic marker data with known specificity patterns."""
        rng = np.random.RandomState(seed)

        # 3 cell types: A, B, C
        Y = np.zeros((N, 3))
        Y[:70, 0] = 0.8   # Region 1: mostly A
        Y[:70, 1] = 0.1
        Y[:70, 2] = 0.1
        Y[70:140, 0] = 0.1
        Y[70:140, 1] = 0.8  # Region 2: mostly B
        Y[70:140, 2] = 0.1
        Y[140:, 0] = 0.1
        Y[140:, 1] = 0.1
        Y[140:, 2] = 0.8   # Region 3: mostly C

        # Marker 0: specific to A
        # Marker 1: specific to B
        # Marker 2: ubiquitous (high everywhere, assigned to C)
        # Marker 3: specific to C
        S = np.zeros((N, 4))
        S[:, 0] = Y[:, 0] * 0.9 + rng.normal(0, 0.05, N)
        S[:, 1] = Y[:, 1] * 0.9 + rng.normal(0, 0.05, N)
        S[:, 2] = 0.5 + rng.normal(0, 0.05, N)
        S[:, 3] = Y[:, 2] * 0.9 + rng.normal(0, 0.05, N)
        S = np.clip(S, 0, 1)

        assignment = np.zeros((4, 3))
        assignment[0, 0] = 1
        assignment[1, 1] = 1
        assignment[2, 2] = 1
        assignment[3, 2] = 1

        marker_owners = [[0], [1], [2], [2]]

        return S, Y, assignment, marker_owners

    def test_specific_markers_score_higher(self):
        """Specific markers should have higher exclusivity than ubiquitous ones."""
        S, Y, assignment, marker_owners = self._make_data()
        exclusivity = compute_marker_exclusivity(S, Y, marker_owners, assignment)

        assert exclusivity.shape == (4,)
        assert exclusivity[0] > exclusivity[2], f"Specific marker 0 ({exclusivity[0]:.3f}) should > ubiquitous ({exclusivity[2]:.3f})"
        assert exclusivity[1] > exclusivity[2], f"Specific marker 1 ({exclusivity[1]:.3f}) should > ubiquitous ({exclusivity[2]:.3f})"
        assert exclusivity[3] > exclusivity[2], f"Specific marker 3 ({exclusivity[3]:.3f}) should > ubiquitous ({exclusivity[2]:.3f})"

    def test_exclusivity_range(self):
        """Exclusivity scores should be in [0.3, 1.0] (floored at 0.3)."""
        S, Y, assignment, marker_owners = self._make_data()
        exclusivity = compute_marker_exclusivity(S, Y, marker_owners, assignment)

        assert np.all(exclusivity >= 0.3), f"Min exclusivity {exclusivity.min():.3f} < 0.3 floor"
        assert np.all(exclusivity <= 1.0), f"Max exclusivity {exclusivity.max():.3f} > 1.0"

    def test_unowned_markers_get_default(self):
        """Markers with no owner should get exclusivity = 1.0 (neutral)."""
        S, Y, assignment, marker_owners = self._make_data()
        S_extended = np.column_stack([S, np.random.rand(S.shape[0])])
        assignment_extended = np.vstack([assignment, np.zeros((1, 3))])
        marker_owners_extended = marker_owners + [[]]

        exclusivity = compute_marker_exclusivity(S_extended, Y, marker_owners_extended, assignment_extended)
        assert exclusivity[4] == 1.0, f"Unowned marker should get 1.0, got {exclusivity[4]:.3f}"

    def test_shared_marker_uses_combined_owners(self):
        """Shared markers should correlate with combined owner proportions."""
        rng = np.random.RandomState(42)
        N = 200
        Y = np.zeros((N, 3))
        Y[:100, 0] = 0.4
        Y[:100, 1] = 0.4
        Y[:100, 2] = 0.2
        Y[100:, 0] = 0.1
        Y[100:, 1] = 0.1
        Y[100:, 2] = 0.8

        S = np.zeros((N, 2))
        S[:, 0] = (Y[:, 0] + Y[:, 1]) * 0.5 + rng.normal(0, 0.05, N)
        S[:, 1] = Y[:, 2] * 0.9 + rng.normal(0, 0.05, N)
        S = np.clip(S, 0, 1)

        assignment = np.zeros((2, 3))
        assignment[0, 0] = 1
        assignment[0, 1] = 1
        assignment[1, 2] = 1

        marker_owners = [[0, 1], [2]]

        exclusivity = compute_marker_exclusivity(S, Y, marker_owners, assignment)
        assert exclusivity[0] > 0.5, f"Shared marker should be discriminative, got {exclusivity[0]:.3f}"


class TestExclusivityInFinetuning:
    """Test that exclusivity weights are accepted by finetuning functions."""

    def test_deconvolute_accepts_marker_exclusivity(self):
        """deconvolute_local_cell_proportions_per_marker should accept marker_exclusivity."""
        sig = inspect.signature(deconvolute_local_cell_proportions_per_marker)
        assert "marker_exclusivity" in sig.parameters, \
            "deconvolute_local_cell_proportions_per_marker missing marker_exclusivity parameter"
        param = sig.parameters["marker_exclusivity"]
        assert param.default is None, f"Default should be None, got {param.default}"

    def test_finetune_accepts_marker_exclusivity(self):
        """finetune_cell_proportions_per_marker should accept marker_exclusivity."""
        sig = inspect.signature(finetune_cell_proportions_per_marker)
        assert "marker_exclusivity" in sig.parameters, \
            "finetune_cell_proportions_per_marker missing marker_exclusivity parameter"
        param = sig.parameters["marker_exclusivity"]
        assert param.default is None, f"Default should be None, got {param.default}"
