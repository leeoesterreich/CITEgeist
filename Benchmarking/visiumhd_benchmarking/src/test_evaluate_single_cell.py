"""Tests for single-cell evaluation.

TDD tests for Task 8: Evaluation pipeline for single-cell type assignment
using Hungarian assignment against ground truth.
"""
import sys
from pathlib import Path
import numpy as np
import pytest

# Add paths
_src_dir = Path(__file__).parent
_repo_root = _src_dir.parent.parent.parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from evaluate_single_cell import (
    evaluate_spot_assignment,
    compute_accuracy_metrics,
    hungarian_assign,
    evaluate_batch,
)


class TestHungarianAssign:
    """Test Hungarian assignment with count constraints."""

    def test_perfect_assignment(self):
        """Test when attention strongly indicates correct types."""
        # 5 nuclei, 3 types
        # Attention strongly suggests: [0, 0, 1, 2, 2]
        cost_matrix = -np.array([
            [0.9, 0.05, 0.05],  # Nucleus 0 -> type 0
            [0.85, 0.1, 0.05],  # Nucleus 1 -> type 0
            [0.05, 0.9, 0.05],  # Nucleus 2 -> type 1
            [0.1, 0.05, 0.85],  # Nucleus 3 -> type 2
            [0.05, 0.1, 0.85],  # Nucleus 4 -> type 2
        ])
        target_counts = np.array([2, 1, 2])  # 2 of type 0, 1 of type 1, 2 of type 2

        assignments = hungarian_assign(cost_matrix, target_counts)

        assert len(assignments) == 5
        assert (assignments == 0).sum() == 2  # Exactly 2 assigned to type 0
        assert (assignments == 1).sum() == 1  # Exactly 1 assigned to type 1
        assert (assignments == 2).sum() == 2  # Exactly 2 assigned to type 2

    def test_ambiguous_assignment(self):
        """Test when attention is less clear but counts constrain result."""
        # All nuclei have similar attention - Hungarian should still respect counts
        cost_matrix = -np.array([
            [0.4, 0.3, 0.3],
            [0.35, 0.35, 0.3],
            [0.3, 0.4, 0.3],
            [0.3, 0.3, 0.4],
            [0.33, 0.33, 0.34],
        ])
        target_counts = np.array([2, 2, 1])

        assignments = hungarian_assign(cost_matrix, target_counts)

        # Should still respect count constraints
        assert (assignments == 0).sum() == 2
        assert (assignments == 1).sum() == 2
        assert (assignments == 2).sum() == 1

    def test_zero_count_type(self):
        """Test when one type has zero expected count."""
        cost_matrix = -np.array([
            [0.8, 0.1, 0.1],
            [0.7, 0.2, 0.1],
            [0.6, 0.3, 0.1],
        ])
        target_counts = np.array([2, 1, 0])  # No type 2

        assignments = hungarian_assign(cost_matrix, target_counts)

        assert len(assignments) == 3
        assert (assignments == 0).sum() == 2
        assert (assignments == 1).sum() == 1
        assert (assignments == 2).sum() == 0

    def test_single_nucleus(self):
        """Test edge case with single nucleus."""
        cost_matrix = -np.array([[0.8, 0.1, 0.1]])
        target_counts = np.array([1, 0, 0])

        assignments = hungarian_assign(cost_matrix, target_counts)

        assert len(assignments) == 1
        assert assignments[0] == 0


class TestEvaluateSpotAssignment:
    """Test single spot evaluation."""

    def test_evaluate_spot_assignment(self):
        """Test single spot evaluation with Hungarian assignment."""
        # Mock attention weights (5 nuclei, 3 types)
        attention = np.array([
            [0.8, 0.1, 0.1],  # Likely type 0
            [0.7, 0.2, 0.1],  # Likely type 0
            [0.1, 0.8, 0.1],  # Likely type 1
            [0.2, 0.1, 0.7],  # Likely type 2
            [0.1, 0.1, 0.8],  # Likely type 2
        ])

        gt_proportions = np.array([0.4, 0.2, 0.4])  # 2, 1, 2
        gt_types = np.array([0, 0, 1, 2, 2])

        results = evaluate_spot_assignment(
            attention, gt_proportions, gt_types,
            use_hungarian=True
        )

        assert 'accuracy' in results
        assert 'assignments' in results
        assert 0 <= results['accuracy'] <= 1
        # With clear attention, should get high accuracy
        assert results['accuracy'] >= 0.8

    def test_evaluate_with_argmax(self):
        """Test evaluation using simple argmax (no constraints)."""
        attention = np.array([
            [0.8, 0.1, 0.1],
            [0.7, 0.2, 0.1],
            [0.1, 0.8, 0.1],
            [0.2, 0.1, 0.7],
            [0.1, 0.1, 0.8],
        ])

        gt_proportions = np.array([0.4, 0.2, 0.4])
        gt_types = np.array([0, 0, 1, 2, 2])

        results = evaluate_spot_assignment(
            attention, gt_proportions, gt_types,
            use_hungarian=False
        )

        assert 'accuracy' in results
        assert 'assignments' in results
        # Argmax should also work well with clear attention
        assert results['accuracy'] >= 0.8

    def test_evaluate_returns_assignments(self):
        """Test that assignments are returned correctly."""
        attention = np.array([
            [0.9, 0.05, 0.05],
            [0.05, 0.9, 0.05],
            [0.05, 0.05, 0.9],
        ])
        gt_proportions = np.array([0.33, 0.33, 0.34])
        gt_types = np.array([0, 1, 2])

        results = evaluate_spot_assignment(
            attention, gt_proportions, gt_types,
            use_hungarian=True
        )

        assert len(results['assignments']) == 3
        assert results['assignments'].dtype in [np.int32, np.int64, int]


class TestComputeAccuracyMetrics:
    """Test accuracy metric computation."""

    def test_compute_accuracy_metrics(self):
        """Test accuracy metric computation."""
        predictions = np.array([0, 0, 1, 2, 2])
        ground_truth = np.array([0, 0, 1, 2, 1])  # One mistake (last is 1, predicted 2)

        metrics = compute_accuracy_metrics(predictions, ground_truth, n_types=3)

        assert 'overall_accuracy' in metrics
        assert 'per_type_f1' in metrics
        assert metrics['overall_accuracy'] == 0.8  # 4/5 correct

    def test_perfect_accuracy(self):
        """Test with perfect predictions."""
        predictions = np.array([0, 1, 2, 0, 1])
        ground_truth = np.array([0, 1, 2, 0, 1])

        metrics = compute_accuracy_metrics(predictions, ground_truth, n_types=3)

        assert metrics['overall_accuracy'] == 1.0
        assert metrics['macro_f1'] == 1.0

    def test_type_names(self):
        """Test with custom type names."""
        predictions = np.array([0, 0, 1])
        ground_truth = np.array([0, 1, 1])
        type_names = ['TypeA', 'TypeB', 'TypeC']

        metrics = compute_accuracy_metrics(
            predictions, ground_truth, n_types=3, type_names=type_names
        )

        assert 'TypeA' in metrics['per_type_f1']
        assert 'TypeB' in metrics['per_type_f1']
        assert 'TypeC' in metrics['per_type_f1']

    def test_confusion_matrix(self):
        """Test confusion matrix generation."""
        predictions = np.array([0, 0, 1, 1])
        ground_truth = np.array([0, 1, 0, 1])

        metrics = compute_accuracy_metrics(predictions, ground_truth, n_types=2)

        cm = metrics['confusion_matrix']
        assert cm.shape == (2, 2)
        # Check diagonal (correct predictions)
        assert cm[0, 0] == 1  # True type 0, predicted 0
        assert cm[1, 1] == 1  # True type 1, predicted 1
        # Check off-diagonal (errors)
        assert cm[0, 1] == 1  # True type 0, predicted 1
        assert cm[1, 0] == 1  # True type 1, predicted 0

    def test_per_type_accuracy(self):
        """Test per-type accuracy calculation."""
        predictions = np.array([0, 0, 1, 1, 2])
        ground_truth = np.array([0, 0, 1, 2, 2])

        metrics = compute_accuracy_metrics(predictions, ground_truth, n_types=3)

        # Type 0: 2/2 correct = 100%
        assert metrics['per_type_accuracy']['Type_0'] == 1.0
        # Type 1: 1/1 correct = 100%
        assert metrics['per_type_accuracy']['Type_1'] == 1.0
        # Type 2: 1/2 correct = 50%
        assert metrics['per_type_accuracy']['Type_2'] == 0.5


class TestEvaluateBatch:
    """Test batch evaluation across multiple spots."""

    def test_evaluate_batch_basic(self):
        """Test batch evaluation with multiple spots."""
        # Create mock data for 3 spots
        spot_data = []

        # Spot 1: 3 nuclei, perfect attention
        spot_data.append({
            'attention': np.array([
                [0.9, 0.05, 0.05],
                [0.05, 0.9, 0.05],
                [0.05, 0.05, 0.9],
            ]),
            'gt_proportions': np.array([0.33, 0.33, 0.34]),
            'gt_types': np.array([0, 1, 2]),
        })

        # Spot 2: 4 nuclei
        spot_data.append({
            'attention': np.array([
                [0.8, 0.1, 0.1],
                [0.8, 0.1, 0.1],
                [0.1, 0.8, 0.1],
                [0.1, 0.1, 0.8],
            ]),
            'gt_proportions': np.array([0.5, 0.25, 0.25]),
            'gt_types': np.array([0, 0, 1, 2]),
        })

        # Spot 3: 5 nuclei with some ambiguity
        spot_data.append({
            'attention': np.array([
                [0.6, 0.3, 0.1],
                [0.5, 0.4, 0.1],
                [0.1, 0.7, 0.2],
                [0.2, 0.2, 0.6],
                [0.1, 0.1, 0.8],
            ]),
            'gt_proportions': np.array([0.4, 0.2, 0.4]),
            'gt_types': np.array([0, 0, 1, 2, 2]),
        })

        results = evaluate_batch(spot_data, n_types=3)

        assert 'overall_accuracy' in results
        assert 'per_spot_accuracy' in results
        assert 'per_type_metrics' in results
        assert len(results['per_spot_accuracy']) == 3
        assert 0 <= results['overall_accuracy'] <= 1

    def test_evaluate_batch_empty(self):
        """Test batch evaluation with empty list."""
        results = evaluate_batch([], n_types=3)

        assert results['overall_accuracy'] == 0.0
        assert results['n_spots'] == 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
