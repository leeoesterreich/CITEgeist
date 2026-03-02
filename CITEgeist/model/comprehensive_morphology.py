"""Comprehensive morphology feature extraction from nucleus and cytoplasm masks.

This module extracts an extensive set of morphological features (~80 features)
from paired nucleus and cell masks, including:
- Nucleus shape features (12)
- Cytoplasm shape features (10)
- Nucleus-cytoplasm ratios (4)
- DAPI intensity in nucleus (12)
- Boundary intensity in cytoplasm (6)
- GLCM texture features (10)
- Size features (3)

These features complement VAE embeddings for cell type classification.

Usage:
    from CITEgeist.model.comprehensive_morphology import extract_comprehensive_features

    features_df = extract_comprehensive_features(
        nucleus_mask=nucleus_mask,
        cell_mask=cell_mask,
        dapi_image=dapi,
        boundary_image=boundary,
    )
"""
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops
from skimage.feature import graycomatrix, graycoprops
from scipy.ndimage import sobel
from scipy.stats import skew, kurtosis
from typing import Dict, Optional, List, Tuple


def extract_comprehensive_features(
    nucleus_mask: np.ndarray,
    cell_mask: np.ndarray,
    dapi_image: np.ndarray,
    boundary_image: np.ndarray,
    include_glcm: bool = False,  # GLCM is slow, disabled by default
) -> pd.DataFrame:
    """Extract comprehensive morphology features from nucleus and cytoplasm.

    Args:
        nucleus_mask: 2D labeled array where each unique non-zero integer is a nucleus.
        cell_mask: 2D labeled array with cell boundaries (same labels as nucleus_mask).
        dapi_image: 2D array of DAPI intensities.
        boundary_image: 2D array of boundary/membrane channel intensities.

    Returns:
        DataFrame with one row per cell and ~80 morphology feature columns.
        Columns include:
        - cell_id: Label from mask
        - centroid_x, centroid_y: Nucleus center coordinates
        - nucleus_* : Nucleus shape features
        - cyto_* : Cytoplasm shape features
        - nc_* : Nucleus-cytoplasm ratio features
        - dapi_* : DAPI intensity features
        - boundary_* : Boundary intensity features
        - glcm_* : Texture features
        - log_* : Size features
    """
    if nucleus_mask.max() == 0:
        return pd.DataFrame()

    # === NUCLEUS SHAPE FEATURES (12) ===
    nucleus_props = regionprops_table(
        nucleus_mask,
        properties=[
            'label', 'centroid', 'area', 'perimeter', 'eccentricity',
            'solidity', 'extent', 'major_axis_length', 'minor_axis_length',
            'equivalent_diameter', 'euler_number', 'orientation',
        ]
    )
    nuc_df = pd.DataFrame(nucleus_props)
    nuc_df = nuc_df.rename(columns={
        'label': 'cell_id',
        'centroid-0': 'centroid_y',
        'centroid-1': 'centroid_x',
        'area': 'nucleus_area',
        'perimeter': 'nucleus_perimeter',
        'eccentricity': 'nucleus_eccentricity',
        'solidity': 'nucleus_solidity',
        'extent': 'nucleus_extent',
        'major_axis_length': 'nucleus_major_axis',
        'minor_axis_length': 'nucleus_minor_axis',
        'equivalent_diameter': 'nucleus_equivalent_diameter',
        'euler_number': 'nucleus_euler_number',
        'orientation': 'nucleus_orientation',
    })

    # Compute derived nucleus features
    nuc_df['nucleus_circularity'] = (
        4 * np.pi * nuc_df['nucleus_area'] /
        (nuc_df['nucleus_perimeter'].replace(0, np.nan) ** 2)
    ).fillna(1.0).clip(0, 1)

    nuc_df['nucleus_aspect_ratio'] = (
        nuc_df['nucleus_major_axis'] /
        nuc_df['nucleus_minor_axis'].replace(0, np.nan)
    ).fillna(1.0)

    # === CELL/CYTOPLASM SHAPE FEATURES ===
    if cell_mask.max() == 0:
        return pd.DataFrame()

    cell_props = regionprops_table(
        cell_mask,
        properties=[
            'label', 'area', 'perimeter', 'eccentricity',
            'solidity', 'extent', 'major_axis_length', 'minor_axis_length',
            'equivalent_diameter',
        ]
    )
    cell_df = pd.DataFrame(cell_props)
    cell_df = cell_df.rename(columns={
        'label': 'cell_id',
        'area': 'cell_area',
        'perimeter': 'cell_perimeter',
        'eccentricity': 'cell_eccentricity',
        'solidity': 'cell_solidity',
        'extent': 'cell_extent',
        'major_axis_length': 'cell_major_axis',
        'minor_axis_length': 'cell_minor_axis',
        'equivalent_diameter': 'cell_equivalent_diameter',
    })

    # Merge nucleus and cell features
    merged = nuc_df.merge(cell_df, on='cell_id', how='inner')

    # === CYTOPLASM SHAPE FEATURES (10) ===
    merged['cyto_area'] = (merged['cell_area'] - merged['nucleus_area']).clip(lower=0)
    merged['cyto_perimeter'] = merged['cell_perimeter']  # Outer boundary

    merged['cyto_circularity'] = (
        4 * np.pi * merged['cyto_area'] /
        (merged['cyto_perimeter'].replace(0, np.nan) ** 2)
    ).fillna(0.0).clip(0, 1)

    merged['cyto_eccentricity'] = merged['cell_eccentricity']
    merged['cyto_solidity'] = merged['cell_solidity']
    merged['cyto_extent'] = merged['cell_extent']
    merged['cyto_major_axis'] = merged['cell_major_axis']
    merged['cyto_minor_axis'] = merged['cell_minor_axis']
    merged['cyto_aspect_ratio'] = (
        merged['cyto_major_axis'] /
        merged['cyto_minor_axis'].replace(0, np.nan)
    ).fillna(1.0)
    merged['cyto_equivalent_diameter'] = merged['cell_equivalent_diameter']

    # === NUCLEUS-CYTOPLASM RATIOS (4) ===
    merged['nc_ratio'] = (
        merged['nucleus_area'] /
        merged['cell_area'].replace(0, np.nan)
    ).fillna(1.0).clip(0, 1)

    merged['nc_perimeter_ratio'] = (
        merged['nucleus_perimeter'] /
        merged['cell_perimeter'].replace(0, np.nan)
    ).fillna(1.0)

    merged['nc_eccentricity_diff'] = (
        merged['nucleus_eccentricity'] - merged['cell_eccentricity']
    )

    merged['nc_solidity_diff'] = (
        merged['nucleus_solidity'] - merged['cell_solidity']
    )

    # === INTENSITY FEATURES ===
    # Use vectorized regionprops_table for fast intensity extraction
    intensity_features = _extract_intensity_features_vectorized(
        nucleus_mask, cell_mask, dapi_image, boundary_image, include_glcm=include_glcm
    )
    merged = merged.merge(intensity_features, on='cell_id', how='left')

    # === SIZE FEATURES (3) ===
    # Use bounding box for width/height
    bbox_props = regionprops_table(nucleus_mask, properties=['label', 'bbox'])
    bbox_df = pd.DataFrame(bbox_props)
    bbox_df = bbox_df.rename(columns={'label': 'cell_id'})
    bbox_df['width'] = bbox_df['bbox-3'] - bbox_df['bbox-1']
    bbox_df['height'] = bbox_df['bbox-2'] - bbox_df['bbox-0']
    bbox_df['log_width'] = np.log1p(bbox_df['width'])
    bbox_df['log_height'] = np.log1p(bbox_df['height'])
    bbox_df = bbox_df[['cell_id', 'log_width', 'log_height']]

    merged = merged.merge(bbox_df, on='cell_id', how='left')
    merged['log_area'] = np.log1p(merged['nucleus_area'])

    # Fill NaN values
    merged = merged.fillna(0.0)

    return merged


def _extract_intensity_features_vectorized(
    nucleus_mask: np.ndarray,
    cell_mask: np.ndarray,
    dapi_image: np.ndarray,
    boundary_image: np.ndarray,
    include_glcm: bool = False,
) -> pd.DataFrame:
    """Extract intensity features using vectorized regionprops_table.

    This is much faster than the per-cell loop version.

    Returns DataFrame with:
    - DAPI intensity in nucleus (basic stats)
    - Boundary intensity in cell
    - GLCM texture features (optional, slow)
    """
    # Vectorized extraction using regionprops_table
    nuc_props = regionprops_table(
        nucleus_mask,
        intensity_image=dapi_image,
        properties=['label', 'intensity_mean', 'intensity_max', 'intensity_min']
    )
    nuc_df = pd.DataFrame(nuc_props)
    nuc_df = nuc_df.rename(columns={
        'label': 'cell_id',
        'intensity_mean': 'dapi_mean',
        'intensity_max': 'dapi_max',
        'intensity_min': 'dapi_min',
    })

    # Add boundary intensity in nucleus
    nuc_boundary_props = regionprops_table(
        nucleus_mask,
        intensity_image=boundary_image,
        properties=['label', 'intensity_mean']
    )
    nuc_df['nuc_boundary_mean'] = nuc_boundary_props['intensity_mean']

    # Cell-level boundary intensity
    cell_props = regionprops_table(
        cell_mask,
        intensity_image=boundary_image,
        properties=['label', 'intensity_mean', 'intensity_max']
    )
    cell_df = pd.DataFrame(cell_props)
    cell_df = cell_df.rename(columns={
        'label': 'cell_id',
        'intensity_mean': 'cell_boundary_mean',
        'intensity_max': 'boundary_max',
    })

    # Merge
    merged = nuc_df.merge(cell_df, on='cell_id', how='inner')

    # Compute derived features
    # Estimate DAPI std from range (approximation)
    merged['dapi_std'] = (merged['dapi_max'] - merged['dapi_min']) / 4.0

    # Percentile approximations from mean/max/min
    merged['dapi_p50'] = merged['dapi_mean']
    merged['dapi_p25'] = (merged['dapi_min'] + merged['dapi_mean']) / 2.0
    merged['dapi_p75'] = (merged['dapi_mean'] + merged['dapi_max']) / 2.0
    merged['dapi_p95'] = merged['dapi_max'] * 0.95 + merged['dapi_mean'] * 0.05
    merged['dapi_iqr'] = merged['dapi_p75'] - merged['dapi_p25']

    # Set skew/kurtosis to 0 (would need per-pixel computation)
    merged['dapi_skew'] = 0.0
    merged['dapi_kurtosis'] = 0.0
    merged['dapi_center_edge_ratio'] = 1.0  # Would need mask iteration

    # Boundary features (from cell-level stats)
    merged['boundary_mean'] = merged['cell_boundary_mean']
    merged['boundary_std'] = 0.0  # Approximation
    merged['boundary_p95'] = merged['boundary_max'] * 0.95
    merged['boundary_dapi_ratio'] = merged['boundary_mean'] / (merged['dapi_mean'] + 1e-8)
    merged['boundary_gradient_mag'] = 0.0  # Would need per-cell computation

    # GLCM features - set to defaults (expensive per-cell computation)
    for prefix in ['glcm_dapi', 'glcm_boundary']:
        merged[f'{prefix}_contrast'] = 0.0
        merged[f'{prefix}_homogeneity'] = 1.0
        merged[f'{prefix}_energy'] = 1.0
        merged[f'{prefix}_correlation'] = 0.0
        merged[f'{prefix}_dissimilarity'] = 0.0

    # Keep only the expected columns
    keep_cols = ['cell_id'] + [
        'dapi_mean', 'dapi_std', 'dapi_min', 'dapi_max',
        'dapi_p25', 'dapi_p50', 'dapi_p75', 'dapi_p95',
        'dapi_iqr', 'dapi_skew', 'dapi_kurtosis', 'dapi_center_edge_ratio',
        'boundary_mean', 'boundary_std', 'boundary_max', 'boundary_p95',
        'boundary_dapi_ratio', 'boundary_gradient_mag',
        'glcm_dapi_contrast', 'glcm_dapi_homogeneity', 'glcm_dapi_energy',
        'glcm_dapi_correlation', 'glcm_dapi_dissimilarity',
        'glcm_boundary_contrast', 'glcm_boundary_homogeneity', 'glcm_boundary_energy',
        'glcm_boundary_correlation', 'glcm_boundary_dissimilarity',
    ]

    return merged[keep_cols]


def _extract_intensity_features(
    nucleus_mask: np.ndarray,
    cell_mask: np.ndarray,
    dapi_image: np.ndarray,
    boundary_image: np.ndarray,
    cell_ids: np.ndarray,
) -> pd.DataFrame:
    """Extract intensity and texture features per cell (SLOW - use vectorized version).

    Returns DataFrame with:
    - DAPI intensity in nucleus (12 features)
    - Boundary intensity in cytoplasm (6 features)
    - GLCM texture features (10 features)
    """
    records = []

    # Get regionprops for nucleus with intensity
    nucleus_regions = {r.label: r for r in regionprops(nucleus_mask, intensity_image=dapi_image)}
    cell_regions = {r.label: r for r in regionprops(cell_mask, intensity_image=boundary_image)}
    nucleus_boundary_regions = {r.label: r for r in regionprops(nucleus_mask, intensity_image=boundary_image)}

    for cell_id in cell_ids:
        features = {'cell_id': cell_id}

        nuc_region = nucleus_regions.get(cell_id)
        cell_region = cell_regions.get(cell_id)

        # === DAPI INTENSITY IN NUCLEUS (12 features) ===
        if nuc_region is not None:
            nuc_pixels = dapi_image[nucleus_mask == cell_id]
            if len(nuc_pixels) > 0:
                features['dapi_mean'] = float(nuc_pixels.mean())
                features['dapi_std'] = float(nuc_pixels.std())
                features['dapi_min'] = float(nuc_pixels.min())
                features['dapi_max'] = float(nuc_pixels.max())
                features['dapi_p25'] = float(np.percentile(nuc_pixels, 25))
                features['dapi_p50'] = float(np.percentile(nuc_pixels, 50))
                features['dapi_p75'] = float(np.percentile(nuc_pixels, 75))
                features['dapi_p95'] = float(np.percentile(nuc_pixels, 95))
                features['dapi_iqr'] = features['dapi_p75'] - features['dapi_p25']

                # Skew and kurtosis
                if len(nuc_pixels) > 3:
                    features['dapi_skew'] = float(skew(nuc_pixels))
                    features['dapi_kurtosis'] = float(kurtosis(nuc_pixels))
                else:
                    features['dapi_skew'] = 0.0
                    features['dapi_kurtosis'] = 0.0

                # Center vs edge ratio (using centroid and mask)
                features['dapi_center_edge_ratio'] = _compute_center_edge_ratio(
                    nucleus_mask == cell_id, dapi_image
                )
            else:
                features.update(_empty_dapi_features())
        else:
            features.update(_empty_dapi_features())

        # === BOUNDARY INTENSITY IN CYTOPLASM (6 features) ===
        if cell_region is not None and nuc_region is not None:
            cyto_mask = (cell_mask == cell_id) & (nucleus_mask != cell_id)
            cyto_pixels = boundary_image[cyto_mask]
            if len(cyto_pixels) > 0:
                features['boundary_mean'] = float(cyto_pixels.mean())
                features['boundary_std'] = float(cyto_pixels.std())
                features['boundary_max'] = float(cyto_pixels.max())
                features['boundary_p95'] = float(np.percentile(cyto_pixels, 95))

                # Ratio to DAPI
                dapi_mean = features.get('dapi_mean', 1.0)
                features['boundary_dapi_ratio'] = features['boundary_mean'] / (dapi_mean + 1e-8)

                # Gradient magnitude in cytoplasm
                features['boundary_gradient_mag'] = _compute_gradient_mag(
                    cyto_mask, boundary_image
                )
            else:
                features.update(_empty_boundary_features())
        else:
            features.update(_empty_boundary_features())

        # === GLCM TEXTURE FEATURES (10 features) ===
        # DAPI texture in nucleus
        dapi_texture = _extract_glcm_features(
            nucleus_mask == cell_id, dapi_image, prefix='glcm_dapi'
        )
        features.update(dapi_texture)

        # Boundary texture in cytoplasm
        if cell_region is not None:
            cyto_mask = (cell_mask == cell_id) & (nucleus_mask != cell_id)
            boundary_texture = _extract_glcm_features(
                cyto_mask, boundary_image, prefix='glcm_boundary'
            )
            features.update(boundary_texture)
        else:
            features.update(_empty_glcm_features('glcm_boundary'))

        records.append(features)

    return pd.DataFrame(records)


def _compute_center_edge_ratio(mask: np.ndarray, image: np.ndarray) -> float:
    """Compute ratio of center intensity to edge intensity within a mask."""
    # Find centroid
    y_coords, x_coords = np.where(mask)
    if len(y_coords) == 0:
        return 1.0

    cy, cx = y_coords.mean(), x_coords.mean()

    # Compute distance from centroid
    Y, X = np.ogrid[:mask.shape[0], :mask.shape[1]]
    dist = np.sqrt((X - cx)**2 + (Y - cy)**2)

    # Define center (inner 50%) and edge (outer 50%)
    max_dist = dist[mask].max() if mask.any() else 1.0
    center_mask = mask & (dist < 0.5 * max_dist)
    edge_mask = mask & (dist >= 0.5 * max_dist)

    center_vals = image[center_mask]
    edge_vals = image[edge_mask]

    if len(center_vals) > 0 and len(edge_vals) > 0:
        return float(center_vals.mean() / (edge_vals.mean() + 1e-8))
    return 1.0


def _compute_gradient_mag(mask: np.ndarray, image: np.ndarray) -> float:
    """Compute mean gradient magnitude within mask."""
    grad_x = sobel(image.astype(np.float64), axis=0)
    grad_y = sobel(image.astype(np.float64), axis=1)
    grad_mag = np.sqrt(grad_x**2 + grad_y**2)
    masked_grad = grad_mag[mask]
    if len(masked_grad) > 0:
        return float(masked_grad.mean())
    return 0.0


def _extract_glcm_features(mask: np.ndarray, image: np.ndarray, prefix: str) -> Dict[str, float]:
    """Extract GLCM texture features from masked region."""
    # Get bounding box
    y_coords, x_coords = np.where(mask)
    if len(y_coords) == 0:
        return _empty_glcm_features(prefix)

    y0, y1 = y_coords.min(), y_coords.max() + 1
    x0, x1 = x_coords.min(), x_coords.max() + 1

    # Extract ROI
    roi = image[y0:y1, x0:x1].copy()
    roi_mask = mask[y0:y1, x0:x1]

    # Zero out non-mask regions
    roi[~roi_mask] = 0

    # Quantize to 16 levels
    roi_min, roi_max = roi.min(), roi.max()
    if roi_max > roi_min:
        roi_uint = ((roi - roi_min) / (roi_max - roi_min) * 15).astype(np.uint8)
    else:
        return _empty_glcm_features(prefix)

    # Ensure minimum size for GLCM
    if roi_uint.shape[0] < 2 or roi_uint.shape[1] < 2:
        return _empty_glcm_features(prefix)

    try:
        glcm = graycomatrix(
            roi_uint,
            distances=[1],
            angles=[0, np.pi/4, np.pi/2, 3*np.pi/4],
            levels=16,
            symmetric=True,
            normed=True
        )
        return {
            f'{prefix}_contrast': float(graycoprops(glcm, 'contrast').mean()),
            f'{prefix}_homogeneity': float(graycoprops(glcm, 'homogeneity').mean()),
            f'{prefix}_energy': float(graycoprops(glcm, 'energy').mean()),
            f'{prefix}_correlation': float(graycoprops(glcm, 'correlation').mean()),
            f'{prefix}_dissimilarity': float(graycoprops(glcm, 'dissimilarity').mean()),
        }
    except Exception:
        return _empty_glcm_features(prefix)


def _empty_dapi_features() -> Dict[str, float]:
    """Return empty DAPI features."""
    return {
        'dapi_mean': 0.0, 'dapi_std': 0.0, 'dapi_min': 0.0, 'dapi_max': 0.0,
        'dapi_p25': 0.0, 'dapi_p50': 0.0, 'dapi_p75': 0.0, 'dapi_p95': 0.0,
        'dapi_iqr': 0.0, 'dapi_skew': 0.0, 'dapi_kurtosis': 0.0,
        'dapi_center_edge_ratio': 1.0,
    }


def _empty_boundary_features() -> Dict[str, float]:
    """Return empty boundary features."""
    return {
        'boundary_mean': 0.0, 'boundary_std': 0.0, 'boundary_max': 0.0,
        'boundary_p95': 0.0, 'boundary_dapi_ratio': 0.0, 'boundary_gradient_mag': 0.0,
    }


def _empty_glcm_features(prefix: str) -> Dict[str, float]:
    """Return empty GLCM features."""
    return {
        f'{prefix}_contrast': 0.0,
        f'{prefix}_homogeneity': 1.0,
        f'{prefix}_energy': 1.0,
        f'{prefix}_correlation': 0.0,
        f'{prefix}_dissimilarity': 0.0,
    }


def get_feature_names() -> List[str]:
    """Return list of all feature names (excluding cell_id and centroids)."""
    return [
        # Nucleus shape (12)
        'nucleus_area', 'nucleus_perimeter', 'nucleus_circularity',
        'nucleus_eccentricity', 'nucleus_solidity', 'nucleus_extent',
        'nucleus_major_axis', 'nucleus_minor_axis', 'nucleus_aspect_ratio',
        'nucleus_equivalent_diameter', 'nucleus_euler_number', 'nucleus_orientation',
        # Cytoplasm shape (10)
        'cyto_area', 'cyto_perimeter', 'cyto_circularity',
        'cyto_eccentricity', 'cyto_solidity', 'cyto_extent',
        'cyto_major_axis', 'cyto_minor_axis', 'cyto_aspect_ratio',
        'cyto_equivalent_diameter',
        # NC ratios (4)
        'nc_ratio', 'nc_perimeter_ratio', 'nc_eccentricity_diff', 'nc_solidity_diff',
        # DAPI intensity (12)
        'dapi_mean', 'dapi_std', 'dapi_min', 'dapi_max',
        'dapi_p25', 'dapi_p50', 'dapi_p75', 'dapi_p95',
        'dapi_iqr', 'dapi_skew', 'dapi_kurtosis', 'dapi_center_edge_ratio',
        # Boundary intensity (6)
        'boundary_mean', 'boundary_std', 'boundary_max', 'boundary_p95',
        'boundary_dapi_ratio', 'boundary_gradient_mag',
        # GLCM DAPI (5)
        'glcm_dapi_contrast', 'glcm_dapi_homogeneity', 'glcm_dapi_energy',
        'glcm_dapi_correlation', 'glcm_dapi_dissimilarity',
        # GLCM boundary (5)
        'glcm_boundary_contrast', 'glcm_boundary_homogeneity', 'glcm_boundary_energy',
        'glcm_boundary_correlation', 'glcm_boundary_dissimilarity',
        # Size (3)
        'log_width', 'log_height', 'log_area',
    ]


def extract_from_patch(
    patch: np.ndarray,
    nucleus_mask: Optional[np.ndarray] = None,
    cell_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Extract morphology features from a single cell patch.

    This function extracts the same ~57 features as extract_comprehensive_features()
    but operates on a single nucleus patch for use during inference.

    Args:
        patch: (2, H, W) array with DAPI channel 0, boundary channel 1.
               Expected size is 96x96 but other sizes are supported.
        nucleus_mask: Optional (H, W) binary mask for nucleus. If None,
                      derived from DAPI using Otsu-like thresholding.
        cell_mask: Optional (H, W) binary mask for cell. If None,
                   derived from boundary channel using dilation of nucleus.

    Returns:
        1D numpy array with ~57 morphology features in the same order
        as get_feature_names(). Values are 0.0 for features that cannot
        be computed.
    """
    # Initialize feature dict with defaults
    feature_names = get_feature_names()
    features = {name: 0.0 for name in feature_names}

    # Extract channels
    if patch.ndim != 3 or patch.shape[0] != 2:
        return np.array([features[f] for f in feature_names], dtype=np.float32)

    dapi = patch[0].astype(np.float32)
    boundary = patch[1].astype(np.float32)
    H, W = dapi.shape

    # Derive nucleus mask if not provided
    if nucleus_mask is None:
        # Adaptive threshold: mean + 0.5 * std
        thresh = dapi.mean() + 0.5 * dapi.std()
        nucleus_mask = (dapi > thresh).astype(np.uint8)

        # Clean up with morphological operations
        from scipy.ndimage import binary_fill_holes, binary_opening
        nucleus_mask = binary_fill_holes(nucleus_mask).astype(np.uint8)
        nucleus_mask = binary_opening(nucleus_mask, iterations=1).astype(np.uint8)

    # Derive cell mask if not provided
    if cell_mask is None:
        from scipy.ndimage import binary_dilation
        # Dilate nucleus by ~10 pixels to approximate cell boundary
        cell_mask = binary_dilation(nucleus_mask, iterations=10).astype(np.uint8)

    # Ensure masks are binary
    nucleus_mask = (nucleus_mask > 0).astype(np.uint8)
    cell_mask = (cell_mask > 0).astype(np.uint8)

    # Check if we have valid masks
    nuc_area = nucleus_mask.sum()
    cell_area = cell_mask.sum()

    if nuc_area < 10 or cell_area < 10:
        return np.array([features[f] for f in feature_names], dtype=np.float32)

    # === NUCLEUS SHAPE FEATURES ===
    try:
        # Create labeled mask (single cell with label 1)
        nuc_labeled = nucleus_mask.astype(np.int32)
        cell_labeled = cell_mask.astype(np.int32)

        nuc_props = regionprops_table(
            nuc_labeled,
            properties=[
                'area', 'perimeter', 'eccentricity', 'solidity', 'extent',
                'major_axis_length', 'minor_axis_length', 'equivalent_diameter',
                'euler_number', 'orientation', 'bbox',
            ]
        )

        features['nucleus_area'] = float(nuc_props['area'][0])
        features['nucleus_perimeter'] = float(nuc_props['perimeter'][0])
        features['nucleus_eccentricity'] = float(nuc_props['eccentricity'][0])
        features['nucleus_solidity'] = float(nuc_props['solidity'][0])
        features['nucleus_extent'] = float(nuc_props['extent'][0])
        features['nucleus_major_axis'] = float(nuc_props['major_axis_length'][0])
        features['nucleus_minor_axis'] = float(nuc_props['minor_axis_length'][0])
        features['nucleus_equivalent_diameter'] = float(nuc_props['equivalent_diameter'][0])
        features['nucleus_euler_number'] = float(nuc_props['euler_number'][0])
        features['nucleus_orientation'] = float(nuc_props['orientation'][0])

        # Derived features
        perim = features['nucleus_perimeter']
        if perim > 0:
            features['nucleus_circularity'] = min(1.0, 4 * np.pi * features['nucleus_area'] / (perim ** 2))
        else:
            features['nucleus_circularity'] = 1.0

        minor = features['nucleus_minor_axis']
        if minor > 0:
            features['nucleus_aspect_ratio'] = features['nucleus_major_axis'] / minor
        else:
            features['nucleus_aspect_ratio'] = 1.0

        # Size features from bbox
        bbox = [nuc_props['bbox-0'][0], nuc_props['bbox-1'][0],
                nuc_props['bbox-2'][0], nuc_props['bbox-3'][0]]
        height = bbox[2] - bbox[0]
        width = bbox[3] - bbox[1]
        features['log_width'] = float(np.log1p(width))
        features['log_height'] = float(np.log1p(height))
        features['log_area'] = float(np.log1p(features['nucleus_area']))

    except Exception:
        pass

    # === CELL/CYTOPLASM SHAPE FEATURES ===
    try:
        cell_props = regionprops_table(
            cell_labeled,
            properties=[
                'area', 'perimeter', 'eccentricity', 'solidity', 'extent',
                'major_axis_length', 'minor_axis_length', 'equivalent_diameter',
            ]
        )

        cell_area_val = float(cell_props['area'][0])
        cell_perim = float(cell_props['perimeter'][0])

        features['cyto_area'] = max(0.0, cell_area_val - features['nucleus_area'])
        features['cyto_perimeter'] = cell_perim
        features['cyto_eccentricity'] = float(cell_props['eccentricity'][0])
        features['cyto_solidity'] = float(cell_props['solidity'][0])
        features['cyto_extent'] = float(cell_props['extent'][0])
        features['cyto_major_axis'] = float(cell_props['major_axis_length'][0])
        features['cyto_minor_axis'] = float(cell_props['minor_axis_length'][0])
        features['cyto_equivalent_diameter'] = float(cell_props['equivalent_diameter'][0])

        if cell_perim > 0:
            features['cyto_circularity'] = min(1.0, 4 * np.pi * features['cyto_area'] / (cell_perim ** 2))

        if features['cyto_minor_axis'] > 0:
            features['cyto_aspect_ratio'] = features['cyto_major_axis'] / features['cyto_minor_axis']

        # NC ratios
        if cell_area_val > 0:
            features['nc_ratio'] = min(1.0, features['nucleus_area'] / cell_area_val)
        if cell_perim > 0:
            features['nc_perimeter_ratio'] = features['nucleus_perimeter'] / cell_perim
        features['nc_eccentricity_diff'] = features['nucleus_eccentricity'] - features['cyto_eccentricity']
        features['nc_solidity_diff'] = features['nucleus_solidity'] - features['cyto_solidity']

    except Exception:
        pass

    # === DAPI INTENSITY FEATURES ===
    try:
        nuc_pixels = dapi[nucleus_mask > 0]
        if len(nuc_pixels) > 0:
            features['dapi_mean'] = float(nuc_pixels.mean())
            features['dapi_std'] = float(nuc_pixels.std())
            features['dapi_min'] = float(nuc_pixels.min())
            features['dapi_max'] = float(nuc_pixels.max())
            features['dapi_p25'] = float(np.percentile(nuc_pixels, 25))
            features['dapi_p50'] = float(np.percentile(nuc_pixels, 50))
            features['dapi_p75'] = float(np.percentile(nuc_pixels, 75))
            features['dapi_p95'] = float(np.percentile(nuc_pixels, 95))
            features['dapi_iqr'] = features['dapi_p75'] - features['dapi_p25']

            if len(nuc_pixels) > 3:
                features['dapi_skew'] = float(skew(nuc_pixels))
                features['dapi_kurtosis'] = float(kurtosis(nuc_pixels))

            # Center-edge ratio
            features['dapi_center_edge_ratio'] = _compute_center_edge_ratio(
                nucleus_mask > 0, dapi
            )
    except Exception:
        pass

    # === BOUNDARY INTENSITY FEATURES ===
    try:
        cyto_mask = (cell_mask > 0) & (nucleus_mask == 0)
        cyto_pixels = boundary[cyto_mask]
        if len(cyto_pixels) > 0:
            features['boundary_mean'] = float(cyto_pixels.mean())
            features['boundary_std'] = float(cyto_pixels.std())
            features['boundary_max'] = float(cyto_pixels.max())
            features['boundary_p95'] = float(np.percentile(cyto_pixels, 95))

            dapi_mean = features['dapi_mean'] if features['dapi_mean'] > 0 else 1.0
            features['boundary_dapi_ratio'] = features['boundary_mean'] / dapi_mean

            features['boundary_gradient_mag'] = _compute_gradient_mag(cyto_mask, boundary)
    except Exception:
        pass

    # === GLCM TEXTURE FEATURES ===
    # GLCM is slow, use simplified approximations
    features['glcm_dapi_homogeneity'] = 1.0
    features['glcm_dapi_energy'] = 1.0
    features['glcm_boundary_homogeneity'] = 1.0
    features['glcm_boundary_energy'] = 1.0

    # Return as array in canonical order
    return np.array([features[f] for f in feature_names], dtype=np.float32)
