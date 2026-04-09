"""Module 3b: Per-Nucleus Cell Type Assignment.

Assigns individual nuclei to cell types using spot-level proportions.
Nuclei are detected via StarDist segmentation (see segmentation.py).

Two assignment methods are available:
1. Random (default): Random assignment within spots, respecting cell type counts.
2. Morphology-guided: Uses nucleus/cell morphology features with soft-label classifier
   and Hungarian matching. Provides minimal improvement (~1% accuracy gain).
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional

from ..morphology.morphology_features import extract_nucleus_features, extract_cell_features, largest_remainder_discretize
from ..morphology.soft_label_classifier import SoftLabelClassifier
from .hungarian_assignment import assign_nuclei_to_types


@dataclass
class NucleusAssignmentResult:
    """Result of nucleus assignment."""
    assignments: Dict[int, str]  # nucleus_id -> cell_type name
    morphology_features: Optional[pd.DataFrame]  # features for all nuclei (None if random)
    classifier: Optional[SoftLabelClassifier]  # trained classifier (None if random)
    assignment_probs: Optional[pd.DataFrame]  # probability matrix (None if random)
    method: str = "random"  # "random" or "morphology"


def random_assign_nuclei_to_types(
    nucleus_ids: np.ndarray,
    counts: np.ndarray,
    rng: Optional[np.random.Generator] = None,
) -> Dict[int, int]:
    """
    Randomly assign nuclei to cell types respecting count constraints.

    Args:
        nucleus_ids: (n_nuclei,) unique identifier for each nucleus
        counts: (n_types,) integer counts per cell type
        rng: Optional random generator for reproducibility

    Returns:
        Dict mapping nucleus_id -> cell_type_index
    """
    if rng is None:
        rng = np.random.default_rng()

    n_nuclei = len(nucleus_ids)
    n_slots = int(counts.sum())

    # Build slot list from counts
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * int(count))

    # Shuffle slots randomly
    rng.shuffle(slots)

    # Assign nuclei to shuffled slots
    assignments = {}
    for i, nid in enumerate(nucleus_ids):
        if i < len(slots):
            assignments[int(nid)] = slots[i]
        else:
            # More nuclei than slots - assign uniformly at random
            assignments[int(nid)] = rng.integers(0, len(counts))

    return assignments


def run_nucleus_assignment(
    mask: np.ndarray,
    nuclei_spot_map: pd.DataFrame,
    proportions: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_types: List[str],
    cell_mask: Optional[np.ndarray] = None,
    use_morphology: bool = False,
    random_seed: Optional[int] = None,
) -> NucleusAssignmentResult:
    """
    Run full nucleus assignment pipeline.

    Args:
        mask: Cellpose label mask (H, W) with nucleus labels
        nuclei_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns
        proportions: DataFrame with 'spot_id' and one column per cell type
        nuclei_counts: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column names in proportions)
        cell_mask: Optional cell mask from watershed (same labels as nucleus mask).
                   If provided and use_morphology=True, extracts cell-level features.
        use_morphology: If True, use morphology-guided Hungarian assignment.
                        If False (default), use random assignment within spots.
                        Note: Morphology provides only ~1% accuracy improvement
                        over random in benchmarks.
        random_seed: Random seed for reproducibility (only used when use_morphology=False)

    Returns:
        NucleusAssignmentResult with assignments and metadata
    """
    # Set up proportions lookup
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Random assignment path (default)
    if not use_morphology:
        rng = np.random.default_rng(random_seed)
        all_assignments = {}

        for spot_id in nuclei_spot_map['spot_id'].unique():
            spot_nuclei = nuclei_spot_map[nuclei_spot_map['spot_id'] == spot_id]
            nucleus_ids = spot_nuclei['nucleus_id'].values

            # Get proportions and nuclei count for this spot
            spot_proportions = spot_props.loc[spot_id].values
            n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

            # Convert proportions to counts
            counts = largest_remainder_discretize(spot_proportions, n_nuclei)

            # Random assignment
            type_assignments = random_assign_nuclei_to_types(nucleus_ids, counts, rng)

            # Convert type indices to names
            for nid, type_idx in type_assignments.items():
                all_assignments[nid] = cell_types[type_idx]

        return NucleusAssignmentResult(
            assignments=all_assignments,
            morphology_features=None,
            classifier=None,
            assignment_probs=None,
            method="random",
        )

    # Morphology-guided assignment path
    # Step 1: Extract morphology features
    if cell_mask is not None:
        # Use cell-level features (12 features)
        morph_df = extract_cell_features(mask, cell_mask)
        morph_df = morph_df.rename(columns={'cell_id': 'nucleus_id'})
        feature_cols = [
            'nucleus_area', 'nucleus_circularity', 'nucleus_eccentricity',
            'nucleus_solidity', 'nucleus_aspect_ratio',
            'cell_area', 'cell_circularity', 'cell_eccentricity',
            'cell_solidity', 'cell_aspect_ratio',
            'nc_ratio', 'cytoplasm_area',
        ]
    else:
        # Use nuclear features only (7 features)
        morph_df = extract_nucleus_features(mask)
        feature_cols = ['area', 'circularity', 'eccentricity', 'solidity',
                        'major_axis_length', 'minor_axis_length', 'aspect_ratio']

    # Merge with spot assignments
    morph_df = morph_df.merge(nuclei_spot_map, on='nucleus_id', how='inner')

    # Step 2: Build training data (soft labels from spot proportions)
    # Get soft labels for each nucleus
    y_soft = morph_df['spot_id'].map(lambda s: spot_props.loc[s].values).values
    y_soft = np.vstack(y_soft)  # (n_nuclei, n_types)

    # Feature matrix
    feature_cols = [c for c in feature_cols if c in morph_df.columns]
    X = morph_df[feature_cols].values

    # Step 3: Train classifier
    clf = SoftLabelClassifier(n_cell_types=len(cell_types))
    clf.fit(X, y_soft)

    # Step 4: Predict probabilities
    probs = clf.predict_proba(X)
    probs_df = pd.DataFrame(probs, columns=cell_types)
    probs_df['nucleus_id'] = morph_df['nucleus_id'].values
    probs_df['spot_id'] = morph_df['spot_id'].values

    # Step 5: Per-spot Hungarian assignment
    all_assignments = {}

    for spot_id in morph_df['spot_id'].unique():
        spot_mask = morph_df['spot_id'] == spot_id
        spot_nuclei = morph_df[spot_mask]
        spot_probs = probs[spot_mask]

        # Get proportions and nuclei count for this spot
        spot_proportions = spot_props.loc[spot_id].values
        n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

        # Convert proportions to counts
        counts = largest_remainder_discretize(spot_proportions, n_nuclei)

        # Run Hungarian assignment
        nucleus_ids = spot_nuclei['nucleus_id'].values
        type_assignments = assign_nuclei_to_types(spot_probs, counts, nucleus_ids)

        # Convert type indices to names
        for nid, type_idx in type_assignments.items():
            all_assignments[nid] = cell_types[type_idx]

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=morph_df,
        classifier=clf,
        assignment_probs=probs_df,
        method="morphology",
    )


