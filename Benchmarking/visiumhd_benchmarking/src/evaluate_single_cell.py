"""Single-cell assignment evaluation.

Evaluates MIL attention weights against ground truth cell types
using Hungarian assignment for constrained evaluation.

This module provides tools to:
1. Convert MIL attention weights to discrete cell assignments
2. Use Hungarian algorithm for count-constrained assignment
3. Compute accuracy, F1, and confusion metrics
4. Evaluate across batches of spots
"""
import sys
from pathlib import Path
import numpy as np
from typing import Dict, Optional, List, Any
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix

# Add repo root
REPO_ROOT = Path(__file__).parent.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def hungarian_assign(
    cost_matrix: np.ndarray,
    target_counts: np.ndarray,
) -> np.ndarray:
    """Perform Hungarian assignment with count constraints.

    Given a cost matrix and target counts per type, finds optimal
    assignment of nuclei to types that minimizes total cost while
    respecting the count constraints.

    Args:
        cost_matrix: (N, K) cost for assigning nucleus i to type k
                     (use negative attention as cost to maximize attention)
        target_counts: (K,) number of nuclei to assign to each type

    Returns:
        (N,) assignments: type index for each nucleus
    """
    n_nuclei = cost_matrix.shape[0]
    n_types = len(target_counts)

    # Expand cost matrix to handle counts
    # Create one column per slot (so if type k has count c, create c columns)
    expanded_cols = []
    col_to_type = []
    for k, count in enumerate(target_counts):
        for _ in range(int(count)):
            expanded_cols.append(cost_matrix[:, k])
            col_to_type.append(k)

    if len(expanded_cols) == 0:
        # No assignments to make - return zeros
        return np.zeros(n_nuclei, dtype=int)

    expanded_cost = np.column_stack(expanded_cols)
    n_slots = expanded_cost.shape[1]

    # Handle size mismatches with padding
    if n_nuclei > n_slots:
        # More nuclei than slots - add dummy columns with high cost
        padding = np.full((n_nuclei, n_nuclei - n_slots), 1e6)
        expanded_cost = np.hstack([expanded_cost, padding])
        col_to_type.extend([-1] * (n_nuclei - n_slots))
    elif n_slots > n_nuclei:
        # More slots than nuclei - add dummy rows with high cost
        padding = np.full((n_slots - n_nuclei, expanded_cost.shape[1]), 1e6)
        expanded_cost = np.vstack([expanded_cost, padding])

    # Run Hungarian algorithm
    row_ind, col_ind = linear_sum_assignment(expanded_cost)

    # Map back to type assignments
    assignments = np.zeros(n_nuclei, dtype=int)
    for r, c in zip(row_ind, col_ind):
        if r < n_nuclei and c < len(col_to_type):
            typ = col_to_type[c]
            if typ >= 0:
                assignments[r] = typ

    return assignments


def evaluate_spot_assignment(
    attention: np.ndarray,
    gt_proportions: np.ndarray,
    gt_types: np.ndarray,
    use_hungarian: bool = True,
) -> Dict[str, Any]:
    """Evaluate single-cell assignment for one spot.

    Converts MIL attention weights to discrete cell type assignments
    and compares against ground truth.

    Args:
        attention: (N, K) attention weights from MIL model
                   Higher values indicate higher likelihood for that type
        gt_proportions: (K,) ground truth proportions (used for count constraints)
        gt_types: (N,) ground truth cell type indices for each nucleus
        use_hungarian: If True, use Hungarian assignment with count constraints
                       If False, use simple argmax

    Returns:
        Dictionary containing:
            - accuracy: Overall assignment accuracy
            - assignments: (N,) predicted type for each nucleus
            - gt_types: (N,) ground truth types (for reference)
    """
    n_nuclei, n_types = attention.shape

    if use_hungarian:
        # Convert proportions to integer counts
        counts = np.round(gt_proportions * n_nuclei).astype(int)

        # Adjust for rounding errors to ensure sum equals n_nuclei
        diff = n_nuclei - counts.sum()
        if diff > 0:
            # Add to largest proportion type
            counts[np.argmax(gt_proportions)] += diff
        elif diff < 0:
            # Remove from largest count type
            counts[np.argmax(counts)] += diff

        # Use negative attention as cost (minimize cost = maximize attention)
        cost = -attention
        assignments = hungarian_assign(cost, counts)
    else:
        # Simple argmax assignment
        assignments = np.argmax(attention, axis=1)

    # Compute accuracy
    accuracy = accuracy_score(gt_types, assignments)

    return {
        'accuracy': accuracy,
        'assignments': assignments,
        'gt_types': gt_types,
    }


def compute_accuracy_metrics(
    predictions: np.ndarray,
    ground_truth: np.ndarray,
    n_types: int,
    type_names: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Compute comprehensive accuracy metrics.

    Args:
        predictions: (N,) predicted cell type indices
        ground_truth: (N,) ground truth cell type indices
        n_types: Number of cell types
        type_names: Optional list of cell type names

    Returns:
        Dictionary containing:
            - overall_accuracy: Fraction of correct predictions
            - macro_f1: Macro-averaged F1 score
            - weighted_f1: Weighted F1 score
            - per_type_f1: Dict mapping type names to F1 scores
            - per_type_accuracy: Dict mapping type names to per-type accuracy
            - confusion_matrix: (n_types, n_types) confusion matrix
    """
    if type_names is None:
        type_names = [f"Type_{i}" for i in range(n_types)]

    # Overall accuracy
    overall_acc = accuracy_score(ground_truth, predictions)

    # F1 scores
    f1_per_type = f1_score(
        ground_truth, predictions,
        labels=list(range(n_types)),
        average=None,
        zero_division=0,
    )
    macro_f1 = f1_score(
        ground_truth, predictions,
        average='macro',
        zero_division=0,
    )
    weighted_f1 = f1_score(
        ground_truth, predictions,
        average='weighted',
        zero_division=0,
    )

    # Confusion matrix
    cm = confusion_matrix(
        ground_truth, predictions,
        labels=list(range(n_types)),
    )

    # Per-type accuracy (recall for each type)
    per_type_acc = {}
    for i, name in enumerate(type_names):
        mask = ground_truth == i
        if mask.sum() > 0:
            per_type_acc[name] = accuracy_score(
                ground_truth[mask], predictions[mask]
            )
        else:
            per_type_acc[name] = np.nan

    return {
        'overall_accuracy': overall_acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'per_type_f1': dict(zip(type_names, f1_per_type)),
        'per_type_accuracy': per_type_acc,
        'confusion_matrix': cm,
    }


def evaluate_batch(
    spot_data: List[Dict[str, np.ndarray]],
    n_types: int,
    type_names: Optional[List[str]] = None,
    use_hungarian: bool = True,
) -> Dict[str, Any]:
    """Evaluate single-cell assignment across multiple spots.

    Args:
        spot_data: List of dicts, each containing:
            - attention: (N_i, K) attention weights
            - gt_proportions: (K,) ground truth proportions
            - gt_types: (N_i,) ground truth cell types
        n_types: Number of cell types
        type_names: Optional list of cell type names
        use_hungarian: Whether to use Hungarian assignment

    Returns:
        Dictionary containing:
            - overall_accuracy: Accuracy across all nuclei
            - per_spot_accuracy: List of per-spot accuracies
            - per_type_metrics: Detailed per-type metrics
            - n_spots: Number of spots evaluated
            - n_nuclei: Total number of nuclei
    """
    if len(spot_data) == 0:
        return {
            'overall_accuracy': 0.0,
            'per_spot_accuracy': [],
            'per_type_metrics': {},
            'n_spots': 0,
            'n_nuclei': 0,
        }

    if type_names is None:
        type_names = [f"Type_{i}" for i in range(n_types)]

    all_predictions = []
    all_ground_truth = []
    per_spot_accuracy = []

    for spot in spot_data:
        attention = spot['attention']
        gt_proportions = spot['gt_proportions']
        gt_types = spot['gt_types']

        result = evaluate_spot_assignment(
            attention, gt_proportions, gt_types,
            use_hungarian=use_hungarian,
        )

        all_predictions.append(result['assignments'])
        all_ground_truth.append(gt_types)
        per_spot_accuracy.append(result['accuracy'])

    # Concatenate all predictions and ground truth
    all_predictions = np.concatenate(all_predictions)
    all_ground_truth = np.concatenate(all_ground_truth)

    # Compute overall metrics
    overall_metrics = compute_accuracy_metrics(
        all_predictions, all_ground_truth,
        n_types=n_types,
        type_names=type_names,
    )

    return {
        'overall_accuracy': overall_metrics['overall_accuracy'],
        'macro_f1': overall_metrics['macro_f1'],
        'weighted_f1': overall_metrics['weighted_f1'],
        'per_spot_accuracy': per_spot_accuracy,
        'per_type_metrics': {
            'f1': overall_metrics['per_type_f1'],
            'accuracy': overall_metrics['per_type_accuracy'],
        },
        'confusion_matrix': overall_metrics['confusion_matrix'],
        'n_spots': len(spot_data),
        'n_nuclei': len(all_predictions),
    }


def compare_with_random_baseline(
    spot_data: List[Dict[str, np.ndarray]],
    n_types: int,
    n_trials: int = 100,
    seed: int = 42,
) -> Dict[str, Any]:
    """Compare Hungarian assignment against random baseline.

    Args:
        spot_data: List of spot data dicts
        n_types: Number of cell types
        n_trials: Number of random trials for baseline
        seed: Random seed

    Returns:
        Dictionary with Hungarian vs random comparison
    """
    np.random.seed(seed)

    # Evaluate Hungarian
    hungarian_results = evaluate_batch(
        spot_data, n_types, use_hungarian=True
    )

    # Evaluate random baseline (multiple trials)
    random_accuracies = []
    for trial in range(n_trials):
        trial_predictions = []
        trial_ground_truth = []

        for spot in spot_data:
            n_nuclei = len(spot['gt_types'])
            gt_proportions = spot['gt_proportions']

            # Random assignment respecting counts
            counts = np.round(gt_proportions * n_nuclei).astype(int)
            diff = n_nuclei - counts.sum()
            if diff > 0:
                counts[np.argmax(gt_proportions)] += diff
            elif diff < 0:
                counts[np.argmax(counts)] += diff

            # Create random assignment
            assignments = []
            for k, c in enumerate(counts):
                assignments.extend([k] * int(c))
            np.random.shuffle(assignments)
            assignments = np.array(assignments[:n_nuclei])

            trial_predictions.append(assignments)
            trial_ground_truth.append(spot['gt_types'])

        # Compute accuracy
        all_pred = np.concatenate(trial_predictions)
        all_gt = np.concatenate(trial_ground_truth)
        random_accuracies.append(accuracy_score(all_gt, all_pred))

    return {
        'hungarian_accuracy': hungarian_results['overall_accuracy'],
        'random_mean': np.mean(random_accuracies),
        'random_std': np.std(random_accuracies),
        'improvement': hungarian_results['overall_accuracy'] - np.mean(random_accuracies),
        'improvement_ratio': (
            hungarian_results['overall_accuracy'] / np.mean(random_accuracies)
            if np.mean(random_accuracies) > 0 else float('inf')
        ),
    }


if __name__ == '__main__':
    # Run simple test
    import json

    print("Running simple evaluation test...")

    # Create mock data
    spot_data = [
        {
            'attention': np.array([
                [0.9, 0.05, 0.05],
                [0.05, 0.9, 0.05],
                [0.05, 0.05, 0.9],
            ]),
            'gt_proportions': np.array([0.33, 0.33, 0.34]),
            'gt_types': np.array([0, 1, 2]),
        },
        {
            'attention': np.array([
                [0.8, 0.1, 0.1],
                [0.7, 0.2, 0.1],
                [0.1, 0.8, 0.1],
                [0.1, 0.1, 0.8],
            ]),
            'gt_proportions': np.array([0.5, 0.25, 0.25]),
            'gt_types': np.array([0, 0, 1, 2]),
        },
    ]

    results = evaluate_batch(spot_data, n_types=3)
    print(f"Overall accuracy: {results['overall_accuracy']:.4f}")
    print(f"Macro F1: {results['macro_f1']:.4f}")
    print(f"Per-spot accuracy: {results['per_spot_accuracy']}")

    # Compare with random
    comparison = compare_with_random_baseline(spot_data, n_types=3)
    print(f"\nHungarian vs Random:")
    print(f"  Hungarian: {comparison['hungarian_accuracy']:.4f}")
    print(f"  Random: {comparison['random_mean']:.4f} +/- {comparison['random_std']:.4f}")
    print(f"  Improvement: {comparison['improvement']:.4f}")
