"""Extract nucleus patches from H&E images.

Extracts 224x224 RGB patches around detected nuclei for ViT input.
"""
import numpy as np
from typing import Tuple, Optional
from skimage.transform import resize


def extract_nucleus_patch(
    image: np.ndarray,
    mask: np.ndarray,
    nucleus_id: int,
    output_size: int = 224,
    expansion: float = 0.75,
) -> np.ndarray:
    """Extract a patch around a nucleus.

    Args:
        image: (H, W, 3) RGB image
        mask: (H, W) segmentation mask with nucleus IDs
        nucleus_id: ID of nucleus to extract
        output_size: Output patch size (square)
        expansion: Fraction to expand bounding box

    Returns:
        (output_size, output_size, 3) RGB patch
    """
    # Get nucleus bounding box
    nucleus_mask = mask == nucleus_id
    if not nucleus_mask.any():
        raise ValueError(f"Nucleus {nucleus_id} not found in mask")

    rows = np.where(nucleus_mask.any(axis=1))[0]
    cols = np.where(nucleus_mask.any(axis=0))[0]

    r_min, r_max = rows.min(), rows.max()
    c_min, c_max = cols.min(), cols.max()

    # Expand bounding box
    height = r_max - r_min
    width = c_max - c_min

    r_expand = int(height * expansion / 2)
    c_expand = int(width * expansion / 2)

    r_min = max(0, r_min - r_expand)
    r_max = min(image.shape[0], r_max + r_expand)
    c_min = max(0, c_min - c_expand)
    c_max = min(image.shape[1], c_max + c_expand)

    # Extract and resize
    crop = image[r_min:r_max, c_min:c_max]

    # Resize to output_size
    resized = resize(crop, (output_size, output_size),
                     preserve_range=True, anti_aliasing=True)

    return resized.astype(np.uint8)


def normalize_he_patch(patch: np.ndarray) -> np.ndarray:
    """Normalize H&E patch for ViT input.

    Args:
        patch: (H, W, 3) RGB patch, uint8

    Returns:
        (3, H, W) normalized float32 tensor
    """
    # Convert to float [0, 1]
    normalized = patch.astype(np.float32) / 255.0

    # HWC to CHW
    normalized = np.transpose(normalized, (2, 0, 1))

    return normalized


def extract_patches_for_spot(
    image: np.ndarray,
    mask: np.ndarray,
    nucleus_ids: np.ndarray,
    output_size: int = 224,
    expansion: float = 0.75,
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract patches for all nuclei in a spot.

    Args:
        image: (H, W, 3) RGB image
        mask: (H, W) segmentation mask
        nucleus_ids: IDs of nuclei in this spot
        output_size: Output patch size
        expansion: Bounding box expansion

    Returns:
        patches: (N, 3, H, W) normalized patches
        valid_ids: IDs of successfully extracted nuclei
    """
    patches = []
    valid_ids = []

    for nid in nucleus_ids:
        try:
            patch = extract_nucleus_patch(
                image, mask, int(nid),
                output_size=output_size,
                expansion=expansion
            )
            normalized = normalize_he_patch(patch)
            patches.append(normalized)
            valid_ids.append(nid)
        except (ValueError, IndexError):
            continue

    if not patches:
        return np.empty((0, 3, output_size, output_size), dtype=np.float32), np.array([])

    return np.stack(patches), np.array(valid_ids)
