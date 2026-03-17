"""Nuclear morphology feature extraction from segmentation masks and image patches."""
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops, label
from skimage.feature import graycomatrix, graycoprops
from scipy.ndimage import sobel


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
        raise ValueError(
            f"Shape mismatch: nucleus_mask {nucleus_mask.shape} vs "
            f"cell_mask {cell_mask.shape}"
        )

    # Get nucleus features
    nuc_df = extract_nucleus_features(nucleus_mask)
    if len(nuc_df) == 0:
        return pd.DataFrame()

    nuc_df = nuc_df.rename(columns={
        'nucleus_id': 'cell_id',
        'area': 'nucleus_area',
        'circularity': 'nucleus_circularity',
        'eccentricity': 'nucleus_eccentricity',
        'solidity': 'nucleus_solidity',
        'aspect_ratio': 'nucleus_aspect_ratio',
    })
    # Drop columns we'll rename differently
    nuc_df = nuc_df.drop(columns=[
        'perimeter', 'major_axis_length', 'minor_axis_length'
    ], errors='ignore')

    # Get cell features using regionprops
    if cell_mask.max() == 0:
        return pd.DataFrame()

    props = regionprops_table(
        cell_mask,
        properties=[
            'label',
            'area',
            'perimeter',
            'eccentricity',
            'solidity',
            'major_axis_length',
            'minor_axis_length',
        ]
    )
    cell_df = pd.DataFrame(props)
    cell_df = cell_df.rename(columns={'label': 'cell_id'})

    # Compute cell circularity
    perimeter = cell_df['perimeter'].replace(0, np.nan)
    cell_df['cell_circularity'] = (4 * np.pi * cell_df['area']) / (perimeter ** 2)
    cell_df['cell_circularity'] = cell_df['cell_circularity'].fillna(1.0).clip(0, 1)

    # Compute cell aspect ratio
    minor = cell_df['minor_axis_length'].replace(0, np.nan)
    cell_df['cell_aspect_ratio'] = cell_df['major_axis_length'] / minor
    cell_df['cell_aspect_ratio'] = cell_df['cell_aspect_ratio'].fillna(1.0)

    # Rename and select columns
    cell_df = cell_df.rename(columns={
        'area': 'cell_area',
        'eccentricity': 'cell_eccentricity',
        'solidity': 'cell_solidity',
    })
    cell_df = cell_df[['cell_id', 'cell_area', 'cell_circularity',
                       'cell_eccentricity', 'cell_solidity', 'cell_aspect_ratio']]

    # Merge nucleus and cell features
    merged = nuc_df.merge(cell_df, on='cell_id', how='inner')

    # Compute ratio features
    merged['nc_ratio'] = merged['nucleus_area'] / merged['cell_area'].replace(0, np.nan)
    merged['nc_ratio'] = merged['nc_ratio'].fillna(1.0).clip(0, 1)
    merged['cytoplasm_area'] = merged['cell_area'] - merged['nucleus_area']
    merged['cytoplasm_area'] = merged['cytoplasm_area'].clip(lower=0)

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


def extract_extended_features(patch: np.ndarray) -> np.ndarray:
    """
    Extract 12 morphology features from a 2-channel patch.

    Args:
        patch: (2, H, W) array with DAPI (channel 0) and boundary (channel 1)

    Returns:
        (12,) feature vector:
        - Basic (6): dapi_mean, dapi_std, dapi_area, boundary_mean, boundary_std, channel_corr
        - Shape (3): circularity, eccentricity, solidity
        - Texture (3): dapi_entropy, dapi_contrast, boundary_gradient_mag
    """
    dapi = patch[0].astype(np.float32)
    boundary = patch[1].astype(np.float32)

    # Basic features (6)
    dapi_mean = float(dapi.mean())
    dapi_std = float(dapi.std())
    dapi_area = float((dapi > dapi.mean()).sum()) if dapi.mean() > 0 else 0.0
    boundary_mean = float(boundary.mean())
    boundary_std = float(boundary.std())

    # Channel correlation
    if dapi.std() > 1e-6 and boundary.std() > 1e-6:
        corr = np.corrcoef(dapi.flatten(), boundary.flatten())[0, 1]
        channel_corr = float(corr) if not np.isnan(corr) else 0.0
    else:
        channel_corr = 0.0

    # Shape features (3) - from thresholded DAPI
    threshold = dapi.mean() + 0.5 * dapi.std() if dapi.std() > 0 else dapi.mean()
    binary = dapi > threshold
    labeled = label(binary)

    if labeled.max() > 0:
        props = regionprops(labeled)
        largest = max(props, key=lambda p: p.area)
        circularity = 4 * np.pi * largest.area / (largest.perimeter ** 2) if largest.perimeter > 0 else 0.0
        eccentricity = float(largest.eccentricity)
        solidity = float(largest.solidity)
    else:
        circularity, eccentricity, solidity = 0.0, 0.0, 0.0

    # Texture features (3)
    # GLCM entropy and contrast
    dapi_uint8 = (dapi / (dapi.max() + 1e-6) * 255).astype(np.uint8)
    if dapi_uint8.max() > dapi_uint8.min():
        glcm = graycomatrix(dapi_uint8, distances=[1], angles=[0], levels=256, symmetric=True, normed=True)
        dapi_contrast = float(graycoprops(glcm, 'contrast')[0, 0])
        # Entropy from normalized GLCM
        glcm_norm = glcm[:, :, 0, 0]
        glcm_norm = glcm_norm / (glcm_norm.sum() + 1e-10)
        dapi_entropy = float(-np.sum(glcm_norm * np.log2(glcm_norm + 1e-10)))
    else:
        dapi_contrast, dapi_entropy = 0.0, 0.0

    # Boundary gradient magnitude
    grad_x = sobel(boundary, axis=0)
    grad_y = sobel(boundary, axis=1)
    boundary_gradient_mag = float(np.sqrt(grad_x**2 + grad_y**2).mean())

    features = np.array([
        dapi_mean, dapi_std, dapi_area, boundary_mean, boundary_std, channel_corr,
        circularity, eccentricity, solidity,
        dapi_entropy, dapi_contrast, boundary_gradient_mag
    ], dtype=np.float32)

    # Replace any remaining NaN with 0
    features = np.nan_to_num(features, nan=0.0)

    return features


def extract_patch(image: np.ndarray, x: float, y: float, size: int = 64) -> np.ndarray:
    """
    Extract a patch centered at (x, y) from a multi-channel image.

    Args:
        image: (H, W, C) or (C, H, W) image
        x, y: center coordinates (in pixels)
        size: patch size (square)

    Returns:
        (C, size, size) patch array
    """
    # Ensure channel-first format
    if image.ndim == 3 and image.shape[2] <= 4:  # (H, W, C)
        image = np.transpose(image, (2, 0, 1))

    C, H, W = image.shape
    half = size // 2

    x, y = int(round(x)), int(round(y))

    # Clamp to image bounds
    x0, x1 = max(0, x - half), min(W, x + half)
    y0, y1 = max(0, y - half), min(H, y + half)

    patch = np.zeros((C, size, size), dtype=image.dtype)

    # Compute offsets for centering
    px0 = half - (x - x0)
    py0 = half - (y - y0)
    px1 = px0 + (x1 - x0)
    py1 = py0 + (y1 - y0)

    patch[:, py0:py1, px0:px1] = image[:, y0:y1, x0:x1]

    return patch
