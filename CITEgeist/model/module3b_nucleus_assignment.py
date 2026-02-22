"""Module 3b: Per-Nucleus Cell Type Assignment.

Assigns individual nuclei to cell types using:
1. Morphology features from Cellpose masks
2. Soft-label classifier trained on spot proportions
3. Hungarian algorithm for optimal assignment
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional

from .morphology_features import extract_nucleus_features, largest_remainder_discretize
from .soft_label_classifier import SoftLabelClassifier
from .hungarian_assignment import assign_nuclei_to_types


@dataclass
class NucleusAssignmentResult:
    """Result of nucleus assignment."""
    assignments: Dict[int, str]  # nucleus_id -> cell_type name
    morphology_features: pd.DataFrame  # features for all nuclei
    classifier: SoftLabelClassifier  # trained classifier
    assignment_probs: pd.DataFrame  # probability matrix used for assignment


def run_nucleus_assignment(
    mask: np.ndarray,
    nuclei_spot_map: pd.DataFrame,
    proportions: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_types: List[str],
) -> NucleusAssignmentResult:
    """
    Run full nucleus assignment pipeline.

    Args:
        mask: Cellpose label mask (H, W) with nucleus labels
        nuclei_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns
        proportions: DataFrame with 'spot_id' and one column per cell type
        nuclei_counts: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column names in proportions)

    Returns:
        NucleusAssignmentResult with assignments and metadata
    """
    # Step 1: Extract morphology features
    morph_df = extract_nucleus_features(mask)

    # Merge with spot assignments
    morph_df = morph_df.merge(nuclei_spot_map, on='nucleus_id', how='inner')

    # Step 2: Build training data (soft labels from spot proportions)
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Get soft labels for each nucleus
    y_soft = morph_df['spot_id'].map(lambda s: spot_props.loc[s].values).values
    y_soft = np.vstack(y_soft)  # (n_nuclei, n_types)

    # Feature matrix
    feature_cols = ['area', 'circularity', 'eccentricity', 'solidity',
                    'major_axis_length', 'minor_axis_length', 'aspect_ratio']
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
    )
