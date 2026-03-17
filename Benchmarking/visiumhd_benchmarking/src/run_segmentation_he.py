"""StarDist segmentation for H&E images.

Segments nuclei from H&E WSIs using StarDist, processing in tiles
to handle large images.
"""
import numpy as np
from typing import List, Tuple, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Lazy import StarDist segmenter (heavy dependency)
_segmenter = None


def get_segmenter(modality: str = 'he'):
    """Get or create StarDist segmenter (singleton).

    Args:
        modality: Segmentation modality ('he' for H&E, 'dapi' for fluorescence)

    Returns:
        StarDistSegmenter instance
    """
    global _segmenter
    if _segmenter is None:
        import importlib.util
        # Import from CITEgeist package
        _src_dir = Path(__file__).parent
        REPO_ROOT = _src_dir.parent.parent.parent
        spec = importlib.util.spec_from_file_location(
            'segmentation',
            REPO_ROOT / 'CITEgeist/model/segmentation.py'
        )
        seg_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(seg_module)
        _segmenter = seg_module.StarDistSegmenter(modality=modality)
    return _segmenter


def preprocess_he_for_stardist(image: np.ndarray) -> np.ndarray:
    """Preprocess H&E image for StarDist nuclei detection.

    StarDist H&E model expects RGB uint8 input. Normalization is handled
    internally by csbdeep. This function ensures correct dtype/shape.

    Args:
        image: (H, W, 3) RGB H&E image or (H, W) grayscale, uint8

    Returns:
        Preprocessed image for StarDist (RGB uint8)
    """
    # Convert grayscale to RGB if needed
    if image.ndim == 2:
        image = np.stack([image] * 3, axis=-1)

    # Ensure uint8
    if image.dtype != np.uint8:
        if image.max() <= 1.0:
            image = (image * 255).astype(np.uint8)
        else:
            image = image.astype(np.uint8)

    return image


def segment_tile(
    tile: np.ndarray,
    modality: str = 'he',
    prob_thresh: Optional[float] = None,
    nms_thresh: Optional[float] = None,
    scale: Optional[float] = None,
) -> np.ndarray:
    """Segment nuclei in a single tile using StarDist.

    Args:
        tile: (H, W, 3) RGB tile or (H, W) grayscale
        modality: Segmentation modality ('he' or 'dapi')
        prob_thresh: Probability threshold for detection (StarDist default if None)
        nms_thresh: Non-maximum suppression threshold (StarDist default if None)
        scale: Rescale factor for image (e.g., 0.5 to downscale 2x)

    Returns:
        (H, W) integer mask with unique ID per nucleus
    """
    segmenter = get_segmenter(modality)

    # Preprocess
    processed = preprocess_he_for_stardist(tile)

    # Build kwargs for StarDist predict_instances
    kwargs = {}
    if prob_thresh is not None:
        kwargs['prob_thresh'] = prob_thresh
    if nms_thresh is not None:
        kwargs['nms_thresh'] = nms_thresh
    if scale is not None:
        kwargs['scale'] = scale

    # Run StarDist segmentation
    masks, _centroids_df = segmenter.segment(processed, **kwargs)

    return masks.astype(np.int32)


def stitch_masks(
    tiles: List[Tuple[np.ndarray, Tuple[int, int]]],
    output_shape: Tuple[int, int],
    overlap: int = 0,
) -> np.ndarray:
    """Stitch tile masks into full image mask.

    Handles overlapping regions by keeping existing labels (first-come-first-served)
    and relabeling new masks to avoid ID conflicts.

    Args:
        tiles: List of (mask, (row_offset, col_offset)) tuples
        output_shape: (H, W) of full output
        overlap: Overlap between tiles in pixels

    Returns:
        Stitched mask with relabeled unique IDs
    """
    output = np.zeros(output_shape, dtype=np.int32)
    max_label = 0

    for mask, (row_off, col_off) in tiles:
        h, w = mask.shape

        # Relabel mask to avoid conflicts with existing labels
        relabeled = mask.copy()
        relabeled[relabeled > 0] += max_label

        # Determine the region to place (handle overlap by taking center)
        # For overlap, we skip the first overlap//2 pixels except at edges
        row_start = row_off + overlap // 2 if row_off > 0 else row_off
        col_start = col_off + overlap // 2 if col_off > 0 else col_off

        # Calculate end positions, clamping to output shape
        row_end = min(row_off + h, output_shape[0])
        col_end = min(col_off + w, output_shape[1])

        # Calculate corresponding region in the mask
        mask_row_start = row_start - row_off
        mask_col_start = col_start - col_off
        mask_row_end = mask_row_start + (row_end - row_start)
        mask_col_end = mask_col_start + (col_end - col_start)

        # Extract regions
        output_region = output[row_start:row_end, col_start:col_end]
        mask_region = relabeled[mask_row_start:mask_row_end, mask_col_start:mask_col_end]

        # Only overwrite where output is zero (first-come-first-served)
        zero_mask = output_region == 0
        output_region[zero_mask] = mask_region[zero_mask]

        # Update max label for next tile
        if output.max() > max_label:
            max_label = output.max()

    return output


def segment_wsi(
    wsi_path: Path,
    output_path: Path,
    tile_size: int = 2048,
    overlap: int = 128,
    modality: str = 'he',
    prob_thresh: Optional[float] = None,
    nms_thresh: Optional[float] = None,
) -> np.ndarray:
    """Segment entire WSI in tiles using StarDist.

    Processes large whole slide images by tiling, segmenting each tile,
    and stitching results together.

    Args:
        wsi_path: Path to WSI TIFF
        output_path: Path to save mask .npy
        tile_size: Tile size in pixels
        overlap: Overlap between tiles
        modality: Segmentation modality ('he' or 'dapi')
        prob_thresh: Probability threshold for detection
        nms_thresh: Non-maximum suppression threshold

    Returns:
        Full segmentation mask
    """
    from PIL import Image
    from tqdm import tqdm

    # Allow large images
    Image.MAX_IMAGE_PIXELS = None

    wsi_path = Path(wsi_path)
    output_path = Path(output_path)

    logger.info(f"Loading WSI from {wsi_path}")
    wsi = np.array(Image.open(wsi_path))
    h, w = wsi.shape[:2]

    logger.info(f"WSI shape: {wsi.shape}")

    # Calculate tile positions
    tiles = []
    step = tile_size - overlap

    # Generate all positions
    positions = []
    for r in range(0, h, step):
        for c in range(0, w, step):
            positions.append((r, c))

    logger.info(f"Processing {len(positions)} tiles")

    for row, col in tqdm(positions, desc="Segmenting tiles"):
        row_end = min(row + tile_size, h)
        col_end = min(col + tile_size, w)

        tile = wsi[row:row_end, col:col_end]

        # Segment this tile
        mask = segment_tile(
            tile,
            modality=modality,
            prob_thresh=prob_thresh,
            nms_thresh=nms_thresh,
        )
        tiles.append((mask, (row, col)))

    # Stitch all masks together
    logger.info("Stitching masks")
    full_mask = stitch_masks(tiles, (h, w), overlap=overlap)

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save mask
    np.save(output_path, full_mask)
    logger.info(f"Saved mask to {output_path}, {full_mask.max()} nuclei detected")

    return full_mask


def extract_centroids(mask: np.ndarray) -> np.ndarray:
    """Extract nucleus centroids from segmentation mask.

    Args:
        mask: (H, W) integer mask with unique ID per nucleus

    Returns:
        (N, 2) array of (y, x) centroid coordinates
    """
    from scipy import ndimage

    labels = np.unique(mask)
    labels = labels[labels > 0]  # Exclude background

    if len(labels) == 0:
        return np.array([]).reshape(0, 2)

    centroids = ndimage.center_of_mass(
        mask > 0,
        mask,
        labels
    )

    return np.array(centroids)


def main():
    """Command-line interface for WSI segmentation."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Segment nuclei from H&E WSI using StarDist'
    )
    parser.add_argument(
        'wsi_path',
        type=Path,
        help='Path to WSI TIFF file'
    )
    parser.add_argument(
        'output_path',
        type=Path,
        help='Path to save output mask (.npy)'
    )
    parser.add_argument(
        '--tile-size',
        type=int,
        default=2048,
        help='Tile size in pixels (default: 2048)'
    )
    parser.add_argument(
        '--overlap',
        type=int,
        default=128,
        help='Overlap between tiles (default: 128)'
    )
    parser.add_argument(
        '--modality',
        type=str,
        default='he',
        choices=['he', 'dapi'],
        help='Segmentation modality (default: he)'
    )
    parser.add_argument(
        '--prob-thresh',
        type=float,
        default=None,
        help='StarDist probability threshold (default: model default)'
    )
    parser.add_argument(
        '--nms-thresh',
        type=float,
        default=None,
        help='StarDist NMS threshold (default: model default)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Run segmentation
    segment_wsi(
        wsi_path=args.wsi_path,
        output_path=args.output_path,
        tile_size=args.tile_size,
        overlap=args.overlap,
        modality=args.modality,
        prob_thresh=args.prob_thresh,
        nms_thresh=args.nms_thresh,
    )


if __name__ == '__main__':
    main()
