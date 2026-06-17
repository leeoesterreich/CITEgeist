"""Segmentation QC PDF report for visual sniff-testing."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle

from .segmentation import SegmentationQC

logger = logging.getLogger(__name__)

_TILES_PER_PAGE = 24  # 6 columns × 4 rows
_TILE_SAMPLES = 24  # 8 low + 8 median + 8 high


def generate_segmentation_qc_pdf(
    output_folder: str,
    sample_name: str,
    image: np.ndarray,
    masks: np.ndarray,
    centroids_xy: np.ndarray,
    spot_centers_xy: np.ndarray,
    spot_radius_px: float,
    nuclei_count_raw: pd.Series,
    qc: Optional[SegmentationQC] = None,
) -> str:
    """Generate a multi-page segmentation QC PDF.

    Page 1: Full tissue overview with mask overlay + count histogram.
    Pages 2+: Spot tile gallery sampling low/median/high count spots.

    Args:
        output_folder: Directory to write the PDF.
        sample_name: Sample identifier for the title.
        image: (H, W, 3) RGB uint8 histology image.
        masks: (H, W) int32 StarDist label array (0 = background).
            All-zero masks (patchwise dummy) trigger centroid-dot fallback.
        centroids_xy: (n_nuclei, 2) nucleus centroids in image frame.
        spot_centers_xy: (n_spots, 2) spot centers in image frame.
        spot_radius_px: Spot radius in image-frame pixels.
        nuclei_count_raw: Per-spot nuclei counts indexed by spot name.
        qc: Optional SegmentationQC metrics for the title bar.

    Returns:
        Path to the written PDF file.
    """
    out_path = str(Path(output_folder) / f"{sample_name}_segmentation_qc.pdf")
    has_masks = masks.max() > 0

    with PdfPages(out_path) as pdf:
        # --- Page 1: Overview ---
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))

        # Left: tissue + overlay
        ax = axes[0]
        ax.imshow(image)
        if has_masks:
            boundary = _extract_boundaries(masks)
            overlay = np.zeros((*masks.shape, 4), dtype=np.float32)
            overlay[boundary, :] = [1.0, 0.2, 0.2, 0.6]
            ax.imshow(overlay)
        elif len(centroids_xy) > 0:
            ax.scatter(
                centroids_xy[:, 0],
                centroids_xy[:, 1],
                s=2,
                c="red",
                alpha=0.5,
                linewidths=0,
            )
            ax.text(
                0.02,
                0.02,
                "Centroid-only (no full masks)",
                transform=ax.transAxes,
                fontsize=8,
                color="yellow",
                backgroundcolor="black",
            )

        for i in range(len(spot_centers_xy)):
            circ = Circle(
                spot_centers_xy[i],
                spot_radius_px,
                fill=False,
                edgecolor="cyan",
                linewidth=0.5,
                alpha=0.4,
            )
            ax.add_patch(circ)
        ax.set_axis_off()

        title_parts = [f"{sample_name}"]
        if qc is not None:
            title_parts.append(
                f"n={qc.n_after_area_filter} nuclei | "
                f"{qc.density_per_mm2:.0f}/mm² ({qc.density_flag}) | "
                f"px={qc.pixel_size_um:.3f} µm/px"
            )
        ax.set_title("\n".join(title_parts), fontsize=10)

        # Right: histogram
        ax2 = axes[1]
        counts = nuclei_count_raw.values
        ax2.hist(counts, bins=max(10, int(counts.max()) + 1), edgecolor="black")
        ax2.set_xlabel("Nuclei per spot")
        ax2.set_ylabel("Number of spots")
        ax2.set_title(f"Distribution (n={len(counts)} spots)")
        mean_c = counts.mean()
        ax2.axvline(mean_c, color="red", linestyle="--", label=f"mean={mean_c:.1f}")
        ax2.legend(fontsize=8)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # --- Pages 2+: Spot tile gallery ---
        selected_indices = _select_representative_spots(
            nuclei_count_raw,
            n_total=_TILE_SAMPLES,
        )

        for page_start in range(0, len(selected_indices), _TILES_PER_PAGE):
            page_indices = selected_indices[page_start : page_start + _TILES_PER_PAGE]
            n_cols = 6
            n_rows = 4
            fig, axes_grid = plt.subplots(
                n_rows,
                n_cols,
                figsize=(18, 12),
            )
            axes_flat = axes_grid.flatten()

            for ax_i, spot_i in enumerate(page_indices):
                ax = axes_flat[ax_i]
                cx, cy = spot_centers_xy[spot_i]
                r = spot_radius_px
                y0 = max(0, int(cy - r * 1.2))
                y1 = min(image.shape[0], int(cy + r * 1.2))
                x0 = max(0, int(cx - r * 1.2))
                x1 = min(image.shape[1], int(cx + r * 1.2))
                tile = image[y0:y1, x0:x1]

                if tile.size == 0:
                    ax.set_visible(False)
                    continue

                ax.imshow(tile)

                if has_masks:
                    mask_tile = masks[y0:y1, x0:x1]
                    boundary_tile = _extract_boundaries(mask_tile)
                    overlay_tile = np.zeros((*mask_tile.shape, 4), dtype=np.float32)
                    overlay_tile[boundary_tile, :] = [1.0, 0.2, 0.2, 0.8]
                    ax.imshow(overlay_tile)
                else:
                    in_tile = (
                        (centroids_xy[:, 0] >= x0)
                        & (centroids_xy[:, 0] < x1)
                        & (centroids_xy[:, 1] >= y0)
                        & (centroids_xy[:, 1] < y1)
                    )
                    if in_tile.any():
                        ax.scatter(
                            centroids_xy[in_tile, 0] - x0,
                            centroids_xy[in_tile, 1] - y0,
                            s=8,
                            c="red",
                            alpha=0.7,
                            linewidths=0,
                        )

                count = int(nuclei_count_raw.iloc[spot_i])
                mean_val = counts.mean()
                std_val = counts.std()

                border_color = "black"
                if count == 0:
                    border_color = "red"
                elif std_val > 0 and count > mean_val + 3 * std_val:
                    border_color = "orange"

                for spine in ax.spines.values():
                    spine.set_edgecolor(border_color)
                    spine.set_linewidth(2 if border_color != "black" else 0.5)

                ax.set_title(f"n={count}", fontsize=8, color=border_color)
                ax.set_xticks([])
                ax.set_yticks([])

            for ax_i in range(len(page_indices), len(axes_flat)):
                axes_flat[ax_i].set_visible(False)

            plt.suptitle(f"{sample_name} — Spot Tiles (sorted by count)", fontsize=11)
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig)
            plt.close(fig)

    logger.info("Segmentation QC PDF: %s", out_path)
    return out_path


def _extract_boundaries(label_image: np.ndarray) -> np.ndarray:
    """Return boolean mask of label boundaries (1-pixel erosion)."""
    from scipy.ndimage import maximum_filter

    dilated = maximum_filter(label_image, size=3)
    return (dilated != label_image) & (label_image > 0) | (dilated > 0) & (label_image == 0)


def _select_representative_spots(
    nuclei_count_raw: pd.Series,
    n_total: int = 24,
) -> list:
    """Select spot indices covering low/median/high count range."""
    n_per_group = n_total // 3
    sorted_idx = np.argsort(nuclei_count_raw.values)
    n_spots = len(sorted_idx)

    if n_spots <= n_total:
        return list(sorted_idx)

    low = sorted_idx[:n_per_group].tolist()
    mid_start = max(0, n_spots // 2 - n_per_group // 2)
    mid = sorted_idx[mid_start : mid_start + n_per_group].tolist()
    high = sorted_idx[-n_per_group:].tolist()

    selected = []
    seen = set()
    for idx in low + mid + high:
        if idx not in seen:
            selected.append(idx)
            seen.add(idx)

    return selected
