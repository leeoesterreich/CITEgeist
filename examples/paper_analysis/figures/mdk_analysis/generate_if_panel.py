#!/usr/bin/env python -u
"""Generate publication-quality IF panel for Figure 4F.

Reads 60x z-stack TIFFs for MCF7 WT and D538G (5-day -E2),
creates max-intensity projections, applies per-channel contrast,
and produces a side-by-side composite with scale bar and legend.

Output: panels/panel-E.png (+ panel-E.svg wrapper for assembly)

Channel mapping: ch0=Hoechst342 (Blue/nuclei), ch1=AF555 (Green/E-cadherin),
ch2=AF647 (Magenta/MDK).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tifffile
from matplotlib import pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))
from _shared.generate_utils import PanelContext, load_panel_sizes, panel_argparser

DATA_ROOT = Path(
    "/path/to/imaging_data"
    "/A.V. Lee Lab's files - Alex Chang/1. Projects/CITEgeist_wetlab"
    "/HW120 T47D MCF7 ESR1mut E2 Depravation MDK IF"
    "/HW120 Midkine Experimental MCF7 T47D ESR1m IF"
)
# FOV choice rationale (2026-05-18 visual review against all 6 WT and 5 D538G
# 60x z-stacks; contact sheets in handoff note
# vault://CITEgeist/handoffs/2026-05-18-fov-review-fig4f.md):
# - WT z-stack 1: cleanest E-cadherin membrane outlines, low/sparse MDK
#   consistent with WT baseline. Stacks 2/3 denser but blurrier E-cad; 4 sparse;
#   5 defocused; 6 tissue edge.
# - D538G z-stack 1: uniform dense monolayer, widespread MDK signal, good
#   channel balance. Stack 2 brighter MDK but overexposed; 3 off-centre;
#   4/5 lower contrast.
WT_TIFF = DATA_ROOT / "MCF7 WT 5 Day" / "60x z stack 1.tif"
D538G_TIFF = DATA_ROOT / "MCF7 D538G 5 Day" / "60x z stack 1.tif"

PANELS_DIR = SCRIPT_DIR / "panels"

CH_HOECHST = 0
CH_AF555 = 1  # E-cadherin
CH_AF647 = 2  # MDK

LO_PCT, HI_PCT = 1.0, 99.5

# 60x Nikon: ~0.1 µm/px typical
PIXEL_UM = 0.11
SCALE_BAR_UM = 10


def max_project(tiff_path: Path) -> np.ndarray:
    """Read z-stack TIFF and return max-intensity projection (Y, X, C) as float32."""
    data = tifffile.imread(str(tiff_path))  # (Z, Y, X, C) uint16
    mip = data.max(axis=0).astype(np.float32)  # (Y, X, C)
    return mip


def contrast_stretch(img: np.ndarray) -> np.ndarray:
    """Percentile-based contrast stretch to [0, 1]."""
    lo = np.percentile(img, LO_PCT)
    hi = np.percentile(img, HI_PCT)
    if hi <= lo:
        return np.zeros_like(img, dtype=np.float32)
    return np.clip((img - lo) / (hi - lo), 0, 1)


def make_composite(mip: np.ndarray) -> np.ndarray:
    """Create RGB composite from (Y, X, C) MIP.

    Channel mapping matches Hunter's Fiji display:
      Blue  = Hoechst342 (nuclei)
      Green = AF555 (E-cadherin)
      Red   = AF647 (MDK)
    Pure RGB (no cross-channel mixing) so each colour reads cleanly at print size.
    """
    hoechst = contrast_stretch(mip[:, :, CH_HOECHST])
    ecad = contrast_stretch(mip[:, :, CH_AF555])
    mdk = contrast_stretch(mip[:, :, CH_AF647])

    r = mdk
    g = ecad
    b = hoechst

    return np.stack([r, g, b], axis=-1)


def add_scale_bar(ax, pixel_um: float, bar_um: float, color="white"):
    """Add scale bar to bottom-right of axis. Label uses Nature's >=10pt minimum.

    Position the bar with enough right-edge clearance that the "µm" label is
    not truncated when the image is rendered into a constrained panel area.
    """
    bar_px = bar_um / pixel_um
    img_h, img_w = ax.get_images()[0].get_array().shape[:2]
    # Right-edge clearance accounts for the "10 µm" label being wider than the bar.
    right_margin = max(img_w * 0.10, bar_px * 0.6)
    x_start = img_w - bar_px - right_margin
    y_pos = img_h - img_h * 0.07
    ax.plot([x_start, x_start + bar_px], [y_pos, y_pos], color=color, linewidth=3)
    ax.text(
        x_start + bar_px / 2,
        y_pos - img_h * 0.025,
        f"{bar_um} µm",
        ha="center",
        va="bottom",
        color=color,
        fontsize=10,
        fontfamily="Arial",
        fontweight="bold",
    )


def draw_panel_F(ax):
    """Draw the IF comparison panel into the given axes (MCF7 WT vs D538G, 5-day -E2).

    Note: the function name is historical; this renders Figure 4 panel **E**.
    """
    print("  Reading WT z-stack...")
    mip_wt = max_project(WT_TIFF)
    print(f"    Shape: {mip_wt.shape}")

    print("  Reading D538G z-stack...")
    mip_d538g = max_project(D538G_TIFF)
    print(f"    Shape: {mip_d538g.shape}")

    rgb_wt = make_composite(mip_wt)
    rgb_d538g = make_composite(mip_d538g)

    # Split the single axis into two side-by-side images
    fig = ax.get_figure()
    ax.set_visible(False)
    pos = ax.get_position()

    # Reserve bottom ~12% of panel for the channel legend. Bumped from 7%
    # to add ~10% panel-height gap between the legend strip and the bottom
    # edge of the micrographs (legend itself sits at y=0.02). Manifest
    # aspect E was tuned for 7% — at 12% the images become very slightly
    # height-limited inside their axes (~3mm horizontal whitespace each
    # side); acceptable for the breathing room win.
    legend_band = 0.12
    # Column headers placed INSIDE the upper image region (with white
    # bbox background). Drawing above the image causes the title text
    # to overflow the panel slot top after compose-stage scaling, even
    # with a reserved top band — matplotlib's bbox_inches="tight" still
    # tightens the panel SVG bounding box around the text.
    img_h_frac = 1.0 - legend_band
    img_y0 = pos.y0 + pos.height * legend_band
    img_h = pos.height * img_h_frac

    # Tighter split: each image 0.495 wide with a 1% gap (was 0.48 / 4% gap).
    ax1 = fig.add_axes([pos.x0, img_y0, pos.width * 0.495, img_h])
    ax2 = fig.add_axes([pos.x0 + pos.width * 0.505, img_y0, pos.width * 0.495, img_h])

    ax1.imshow(rgb_wt)
    ax1.axis("off")
    add_scale_bar(ax1, PIXEL_UM, SCALE_BAR_UM)

    ax2.imshow(rgb_d538g)
    ax2.axis("off")

    title_bbox = dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.92)
    ax1.text(
        0.5,
        0.90,
        "MCF7 WT",
        transform=ax1.transAxes,
        fontsize=10,
        fontfamily="Arial",
        ha="center",
        va="top",
        color="black",
        bbox=title_bbox,
    )
    ax2.text(
        0.5,
        0.90,
        "MCF7 D538G",
        transform=ax2.transAxes,
        fontsize=10,
        fontfamily="Arial",
        ha="center",
        va="top",
        color="black",
        bbox=title_bbox,
    )

    # Channel legend along the reserved bottom band, in Fiji-display colours.
    # Short labels (full fluorophore detail lives in the figure caption).
    legend_labels = [
        ("Nuclei", "#3D9BFF"),
        ("E-cadherin", "#39B54A"),
        ("MDK", "#E83A3A"),
    ]
    legend_y = pos.y0 + pos.height * 0.02
    centres = [0.18, 0.50, 0.82]
    for (name, color), cx_frac in zip(legend_labels, centres):
        cx = pos.x0 + pos.width * cx_frac
        fig.text(
            cx,
            legend_y,
            f"■ {name}",
            fontsize=10,
            fontfamily="Arial",
            color=color,
            transform=fig.transFigure,
            va="bottom",
            ha="center",
        )


def main():
    parser = panel_argparser("Generate IF panel from raw z-stacks.")
    parser.add_argument("--label", default="E")
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()
    panel_sizes = load_panel_sizes(args)

    sizes_dict = None
    if panel_sizes:
        sizes_dict = {
            label: {"w_mm": v[0], "h_mm": v[1]} if isinstance(v, (list, tuple)) else v
            for label, v in panel_sizes.items()
        }

    from pathlib import Path

    out_dir = Path(args.output_dir) if args.output_dir else PANELS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Generating IF panel {args.label}...")
    with PanelContext(args.label, sizes_dict, out_dir) as (fig, ax):
        draw_panel_F(ax)
    print("Done.")


if __name__ == "__main__":
    main()
