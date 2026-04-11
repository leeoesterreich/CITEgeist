"""Nuclear morphology feature extraction from segmentation masks and image patches."""

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
        return pd.DataFrame(
            columns=[
                "nucleus_id",
                "centroid_x",
                "centroid_y",
                "area",
                "perimeter",
                "circularity",
                "eccentricity",
                "solidity",
                "major_axis_length",
                "minor_axis_length",
                "aspect_ratio",
            ]
        )

    # Extract properties using skimage
    props = regionprops_table(
        mask,
        properties=[
            "label",
            "centroid",
            "area",
            "perimeter",
            "eccentricity",
            "solidity",
            "major_axis_length",
            "minor_axis_length",
        ],
    )

    df = pd.DataFrame(props)

    # Rename columns for clarity
    df = df.rename(
        columns={
            "label": "nucleus_id",
            "centroid-0": "centroid_y",  # skimage uses row, col order
            "centroid-1": "centroid_x",
        }
    )

    # Compute derived features
    # Circularity: 4*pi*area / perimeter^2
    # Protect against zero perimeter (single pixel)
    perimeter = df["perimeter"].replace(0, np.nan)
    df["circularity"] = (4 * np.pi * df["area"]) / (perimeter**2)
    df["circularity"] = df["circularity"].fillna(1.0).clip(0, 1)

    # Aspect ratio: major / minor (protect against zero)
    minor = df["minor_axis_length"].replace(0, np.nan)
    df["aspect_ratio"] = df["major_axis_length"] / minor
    df["aspect_ratio"] = df["aspect_ratio"].fillna(1.0)

    # Reorder columns
    columns = [
        "nucleus_id",
        "centroid_x",
        "centroid_y",
        "area",
        "perimeter",
        "circularity",
        "eccentricity",
        "solidity",
        "major_axis_length",
        "minor_axis_length",
        "aspect_ratio",
    ]

    return df[columns]


def extract_cell_features(
    nucleus_mask: np.ndarray,
    cell_mask: np.ndarray,
) -> pd.DataFrame:
    """
    Extract morphology features from paired nucleus and cell masks.

    Assumes nucleus_mask and cell_mask use the same integer labels,
    where each label corresponds to a nucleus-cell pair.

    Args:
        nucleus_mask: 2D labeled array of nuclei (from segmentation)
        cell_mask: 2D labeled array of cells (from watershed)

    Returns:
        DataFrame with columns:
            - cell_id: label from masks
            - centroid_x, centroid_y: nucleus center coordinates
            Nuclear features (prefix 'nucleus_'):
            - nucleus_area, nucleus_circularity, nucleus_eccentricity,
              nucleus_solidity, nucleus_aspect_ratio
            Cell features (prefix 'cell_'):
            - cell_area, cell_circularity, cell_eccentricity,
              cell_solidity, cell_aspect_ratio
            Ratio features:
            - nc_ratio: nucleus_area / cell_area
            - cytoplasm_area: cell_area - nucleus_area
    """
    if nucleus_mask.shape != cell_mask.shape:
        raise ValueError(f"Shape mismatch: nucleus_mask {nucleus_mask.shape} vs " f"cell_mask {cell_mask.shape}")

    # Get nucleus features
    nuc_df = extract_nucleus_features(nucleus_mask)
    if len(nuc_df) == 0:
        return pd.DataFrame()

    nuc_df = nuc_df.rename(
        columns={
            "nucleus_id": "cell_id",
            "area": "nucleus_area",
            "circularity": "nucleus_circularity",
            "eccentricity": "nucleus_eccentricity",
            "solidity": "nucleus_solidity",
            "aspect_ratio": "nucleus_aspect_ratio",
        }
    )
    # Drop columns we'll rename differently
    nuc_df = nuc_df.drop(columns=["perimeter", "major_axis_length", "minor_axis_length"], errors="ignore")

    # Get cell features using regionprops
    if cell_mask.max() == 0:
        return pd.DataFrame()

    props = regionprops_table(
        cell_mask,
        properties=[
            "label",
            "area",
            "perimeter",
            "eccentricity",
            "solidity",
            "major_axis_length",
            "minor_axis_length",
        ],
    )
    cell_df = pd.DataFrame(props)
    cell_df = cell_df.rename(columns={"label": "cell_id"})

    # Compute cell circularity
    perimeter = cell_df["perimeter"].replace(0, np.nan)
    cell_df["cell_circularity"] = (4 * np.pi * cell_df["area"]) / (perimeter**2)
    cell_df["cell_circularity"] = cell_df["cell_circularity"].fillna(1.0).clip(0, 1)

    # Compute cell aspect ratio
    minor = cell_df["minor_axis_length"].replace(0, np.nan)
    cell_df["cell_aspect_ratio"] = cell_df["major_axis_length"] / minor
    cell_df["cell_aspect_ratio"] = cell_df["cell_aspect_ratio"].fillna(1.0)

    # Rename and select columns
    cell_df = cell_df.rename(
        columns={
            "area": "cell_area",
            "eccentricity": "cell_eccentricity",
            "solidity": "cell_solidity",
        }
    )
    cell_df = cell_df[
        ["cell_id", "cell_area", "cell_circularity", "cell_eccentricity", "cell_solidity", "cell_aspect_ratio"]
    ]

    # Merge nucleus and cell features
    merged = nuc_df.merge(cell_df, on="cell_id", how="inner")

    # Compute ratio features
    merged["nc_ratio"] = merged["nucleus_area"] / merged["cell_area"].replace(0, np.nan)
    merged["nc_ratio"] = merged["nc_ratio"].fillna(1.0).clip(0, 1)
    merged["cytoplasm_area"] = merged["cell_area"] - merged["nucleus_area"]
    merged["cytoplasm_area"] = merged["cytoplasm_area"].clip(lower=0)

    return merged


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
        n_types = len(proportions)
        # Give one extra count to top n_remaining types (wrap around if needed)
        for i in range(n_remaining):
            counts[indices[i % n_types]] += 1

    return counts
