"""Tests for adaptive per-marker baseline (alpha/intercept)."""
import inspect

import numpy as np
import pytest

from CITEgeist.model.gurobi_impl import (
    deconvolute_local_cell_proportions_per_marker,
    finetune_cell_proportions_per_marker,
    optimize_cell_proportions_per_marker,
)


class TestMarkerBaseline:
    """Test that per-marker alpha (baseline) is learned correctly."""

    def _make_data(self, N=200, seed=42):
        """Create synthetic data with one ubiquitous marker and one specific marker.

        Markers:
        - Marker 0: specific to cell type 0 (signal ~ Y[:,0], baseline ~0)
        - Marker 1: ubiquitous with variation (baseline=0.6, signal += 0.3*Y[:,1])
        - Marker 2: specific to cell type 2 (signal ~ Y[:,2], baseline ~0)
        """
        rng = np.random.RandomState(seed)
        T = 3

        # Ground truth proportions
        Y_true = rng.dirichlet([1, 1, 1], size=N)

        # Marker signals
        M = 3
        S = np.zeros((N, M))
        S[:, 0] = Y_true[:, 0] * 0.8 + rng.normal(0, 0.02, N)  # specific
        S[:, 1] = 0.6 + Y_true[:, 1] * 0.3 + rng.normal(0, 0.02, N)  # ubiquitous
        S[:, 2] = Y_true[:, 2] * 0.8 + rng.normal(0, 0.02, N)  # specific
        S = np.clip(S, 0, 1)

        # Normalize per column (like map_antibodies_to_profiles_v2)
        col_max = np.max(S, axis=0)
        S = S / col_max

        assignment = np.eye(M, T)
        cell_type_names = ["TypeA", "TypeB", "TypeC"]
        marker_names = ["specific_A", "ubiquitous_B", "specific_C"]

        return S, assignment, cell_type_names, marker_names, Y_true

    def test_returns_alpha_values(self):
        """optimize_cell_proportions_per_marker should return alpha array."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        result = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=5, warn_only=True,
        )
        assert len(result) == 4, f"Expected 4 return values, got {len(result)}"
        alpha_values = result[3]
        assert alpha_values.shape == (3,), f"Alpha shape should be (3,), got {alpha_values.shape}"

    def test_ubiquitous_marker_learns_nonzero_alpha(self):
        """Ubiquitous marker should learn alpha > 0."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        _, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=10, warn_only=True,
        )
        assert alpha_values[1] > 0.1, (
            f"Ubiquitous marker alpha={alpha_values[1]:.3f}, expected > 0.1"
        )

    def test_specific_markers_learn_near_zero_alpha(self):
        """Specific markers should learn alpha near 0."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        _, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=10, warn_only=True,
        )
        assert alpha_values[0] < 0.15, (
            f"Specific marker 0 alpha={alpha_values[0]:.3f}, expected < 0.15"
        )
        assert alpha_values[2] < 0.15, (
            f"Specific marker 2 alpha={alpha_values[2]:.3f}, expected < 0.15"
        )

    def test_alpha_clipped_to_range(self):
        """Alpha should be in [0, alpha_max]."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        _, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=10, warn_only=True,
            alpha_max=0.8,
        )
        assert np.all(alpha_values >= 0), f"Alpha min={alpha_values.min():.3f} < 0"
        assert np.all(alpha_values <= 0.8), f"Alpha max={alpha_values.max():.3f} > 0.8"

    def test_proportions_improve_with_baseline(self):
        """With baseline, proportions for ubiquitous marker's cell type should correlate with GT."""
        S, assignment, ct_names, m_names, Y_true = self._make_data(N=300)
        Y_values, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=15, warn_only=True,
        )
        r = np.corrcoef(Y_true[:, 1], Y_values[:, 1])[0, 1]
        assert r > 0.3, f"TypeB correlation with GT = {r:.3f}, expected > 0.3"


class TestBaselineInFinetuning:
    """Test that finetuning functions accept marker_alpha."""

    def test_deconvolute_accepts_marker_alpha(self):
        """deconvolute_local_cell_proportions_per_marker should accept marker_alpha."""
        sig = inspect.signature(deconvolute_local_cell_proportions_per_marker)
        assert "marker_alpha" in sig.parameters, \
            "deconvolute_local_cell_proportions_per_marker missing marker_alpha parameter"
        assert sig.parameters["marker_alpha"].default is None

    def test_finetune_accepts_marker_alpha(self):
        """finetune_cell_proportions_per_marker should accept marker_alpha."""
        sig = inspect.signature(finetune_cell_proportions_per_marker)
        assert "marker_alpha" in sig.parameters, \
            "finetune_cell_proportions_per_marker missing marker_alpha parameter"
        assert sig.parameters["marker_alpha"].default is None


class TestOrchestration:
    """Test that orchestration passes alpha through."""

    def test_model_passes_alpha_to_finetuning(self):
        """run_cell_proportion_model should handle and pass alpha_values."""
        source = inspect.getsource(
            __import__("CITEgeist.model.citegeist_model", fromlist=["CitegeistModel"]).CitegeistModel.run_cell_proportion_model
        )
        assert "alpha_values" in source, \
            "run_cell_proportion_model should handle alpha_values"
        assert "marker_alpha" in source, \
            "run_cell_proportion_model should pass marker_alpha to finetuning"
