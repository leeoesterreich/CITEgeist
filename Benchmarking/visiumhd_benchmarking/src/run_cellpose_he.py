"""Cellpose segmentation for H&E images.

Segments nuclei from H&E WSIs using Cellpose, processing in tiles
to handle large images.
"""
import numpy as np
from typing import List, Tuple, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Lazy import Cellpose (heavy dependency)
_cellpose_model = None


def get_cellpose_model(model_type: str = 'nuclei', gpu: bool = True):
    """Get or create Cellpose model (singleton).

    Args:
        model_type: Cellpose model type ('nuclei', 'cyto', etc.)
        gpu: Use GPU if available

    Returns:
        Cellpose model instance
    """
    global _cellpose_model
    if _cellpose_model is None:
        from cellpose import models
        _cellpose_model = models.Cellpose(model_type=model_type, gpu=gpu)
    return _cellpose_model


def preprocess_he_for_cellpose(image: np.ndarray) -> np.ndarray:
    """Preprocess H&E image for Cellpose nuclei detection.

    Cellpose nuclei model works best on images where nuclei are bright.
    For H&E, hematoxylin stains nuclei blue/purple (dark), so we invert.

    Args:
        image: (H, W, 3) RGB H&E image or (H, W) grayscale, uint8

    Returns:
        Preprocessed image for Cellpose (grayscale, nuclei bright), uint8
    """
    # Convert to grayscale - nuclei appear darker in H&E
    if image.ndim == 3:
        # Standard luminance conversion
        gray = 0.299 * image[..., 0] + 0.587 * image[..., 1] + 0.114 * image[..., 2]
        gray = gray.astype(np.float32)
    else:
        gray = image.astype(np.float32)

    # Invert so nuclei are bright (Cellpose expects bright objects on dark background)
    gray = 255.0 - gray

    # Normalize to [0, 255] range
    min_val = gray.min()
    max_val = gray.max()
    if max_val - min_val > 1e-8:
        gray = (gray - min_val) / (max_val - min_val) * 255.0
    else:
        gray = np.zeros_like(gray)

    return gray.astype(np.uint8)


def segment_tile(
    tile: np.ndarray,
    diameter: float = 30,
    model_type: str = 'nuclei',
    gpu: bool = True,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
) -> np.ndarray:
    """Segment nuclei in a single tile.

    Args:
        tile: (H, W, 3) RGB tile or (H, W) grayscale
        diameter: Expected nucleus diameter in pixels
        model_type: Cellpose model type
        gpu: Use GPU if available
        flow_threshold: Flow error threshold for Cellpose
        cellprob_threshold: Cell probability threshold

    Returns:
        (H, W) integer mask with unique ID per nucleus
    """
    model = get_cellpose_model(model_type, gpu)

    # Preprocess
    processed = preprocess_he_for_cellpose(tile)

    # Run Cellpose - channels=[0,0] for grayscale
    masks, flows, styles = model.eval(
        processed,
        diameter=diameter,
        channels=[0, 0],
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )

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
    diameter: float = 30,
    gpu: bool = True,
    model_type: str = 'nuclei',
) -> np.ndarray:
    """Segment entire WSI in tiles.

    Processes large whole slide images by tiling, segmenting each tile,
    and stitching results together.

    Args:
        wsi_path: Path to WSI TIFF
        output_path: Path to save mask .npy
        tile_size: Tile size in pixels
        overlap: Overlap between tiles
        diameter: Expected nucleus diameter
        gpu: Use GPU
        model_type: Cellpose model type

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
            diameter=diameter,
            model_type=model_type,
            gpu=gpu,
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
        description='Segment nuclei from H&E WSI using Cellpose'
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
        '--diameter',
        type=float,
        default=30,
        help='Expected nucleus diameter (default: 30)'
    )
    parser.add_argument(
        '--no-gpu',
        action='store_true',
        help='Disable GPU usage'
    )
    parser.add_argument(
        '--model-type',
        type=str,
        default='nuclei',
        help='Cellpose model type (default: nuclei)'
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
        diameter=args.diameter,
        gpu=not args.no_gpu,
        model_type=args.model_type,
    )


if __name__ == '__main__':
    main()
