"""Nuclear morphology feature extraction from Cellpose masks."""
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table


def extract_nucleus_features(mask: np.ndarray) -> pd.DataFrame:
    """
    Extract morphology features from a labeled nucleus mask.

    Args:
        mask: 2D numpy array where each unique non-zero integer labels a nucleus.
              Background should be 0.

    Returns:
        DataFrame with columns:
            - nucleus_id: label from mask
            - centroid_x, centroid_y: nucleus center coordinates
            - area: pixel count
            - perimeter: boundary length
            - circularity: 4*pi*area / perimeter^2 (1.0 = perfect circle)
            - eccentricity: 0 = circle, 1 = line
            - solidity: area / convex_hull_area
            - major_axis_length, minor_axis_length: fitted ellipse axes
            - aspect_ratio: major / minor axis
    """
    if mask.max() == 0:
        # Return empty DataFrame with correct columns
        return pd.DataFrame(columns=[
            'nucleus_id', 'centroid_x', 'centroid_y', 'area', 'perimeter',
            'circularity', 'eccentricity', 'solidity',
            'major_axis_length', 'minor_axis_length', 'aspect_ratio'
        ])

    # Extract properties using skimage
    props = regionprops_table(
        mask,
        properties=[
            'label',
            'centroid',
            'area',
            'perimeter',
            'eccentricity',
            'solidity',
            'major_axis_length',
            'minor_axis_length',
        ]
    )

    df = pd.DataFrame(props)

    # Rename columns for clarity
    df = df.rename(columns={
        'label': 'nucleus_id',
        'centroid-0': 'centroid_y',  # skimage uses row, col order
        'centroid-1': 'centroid_x',
    })

    # Compute derived features
    # Circularity: 4*pi*area / perimeter^2
    # Protect against zero perimeter (single pixel)
    perimeter = df['perimeter'].replace(0, np.nan)
    df['circularity'] = (4 * np.pi * df['area']) / (perimeter ** 2)
    df['circularity'] = df['circularity'].fillna(1.0).clip(0, 1)

    # Aspect ratio: major / minor (protect against zero)
    minor = df['minor_axis_length'].replace(0, np.nan)
    df['aspect_ratio'] = df['major_axis_length'] / minor
    df['aspect_ratio'] = df['aspect_ratio'].fillna(1.0)

    # Reorder columns
    columns = [
        'nucleus_id', 'centroid_x', 'centroid_y', 'area', 'perimeter',
        'circularity', 'eccentricity', 'solidity',
        'major_axis_length', 'minor_axis_length', 'aspect_ratio'
    ]

    return df[columns]


def largest_remainder_discretize(proportions: np.ndarray, n_total: int) -> np.ndarray:
    """
    Convert proportions to integer counts using largest remainder method.

    Ensures counts sum exactly to n_total while respecting proportions.

    Args:
        proportions: Array of proportions (should sum to ~1.0)
        n_total: Total count to distribute

    Returns:
        Integer array of counts summing to n_total
    """
    if n_total == 0:
        return np.zeros(len(proportions), dtype=int)

    # Scale proportions to target total
    scaled = proportions * n_total

    # Take floor as initial allocation
    counts = np.floor(scaled).astype(int)

    # Compute remainders
    remainders = scaled - counts

    # Number of additional units to allocate
    n_remaining = n_total - counts.sum()

    # Allocate to types with largest remainders
    if n_remaining > 0:
        # Get indices sorted by remainder (descending)
        indices = np.argsort(-remainders)
        # Give one extra count to top n_remaining types
        for i in range(n_remaining):
            counts[indices[i]] += 1

    return counts
