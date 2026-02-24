"""Watershed cell segmentation using nuclei as seeds."""
import numpy as np
from skimage.segmentation import watershed
from skimage.filters import sobel


def watershed_from_nuclei(
    nucleus_mask: np.ndarray,
    boundary_img: np.ndarray,
    use_gradient: bool = True,
) -> np.ndarray:
    """
    Expand nuclear masks to cell boundaries using watershed.

    Args:
        nucleus_mask: 2D labeled array where each nucleus has unique integer ID.
                      Background = 0.
        boundary_img: 2D float array of boundary channel (membrane stain).
                      Higher values = membrane signal.
        use_gradient: If True, use Sobel gradient of boundary_img as elevation.
                      If False, use inverted boundary_img directly.

    Returns:
        Cell mask with same labels as nucleus_mask, expanded to cell boundaries.
    """
    if nucleus_mask.shape != boundary_img.shape:
        raise ValueError(
            f"Shape mismatch: nucleus_mask {nucleus_mask.shape} vs "
            f"boundary_img {boundary_img.shape}"
        )

    # Compute elevation map for watershed
    if use_gradient:
        # Sobel gradient - edges become ridges (high values)
        elevation = sobel(boundary_img.astype(np.float64))
    else:
        # Invert so membrane (high) becomes ridge
        elevation = -boundary_img.astype(np.float64)

    # Generate markers from nucleus centroids
    # Each nucleus label becomes a seed
    markers = nucleus_mask.copy()

    # Run watershed
    cell_mask = watershed(elevation, markers=markers)

    return cell_mask.astype(np.int32)
