"""Shared spatial plotting utilities for manuscript figures.

Usage (from any figure generate.py, after sys.path.insert for figures/):
    from _shared.spatial_utils import compute_spot_radius, draw_pie_spots
    from _shared.spatial_utils import load_visium_he, draw_he_background
    from _shared.spatial_utils import inject_spatial_metadata
    from _shared.spatial_utils import compute_tissue_crop_from_coords, compute_tissue_crop
    from _shared.spatial_utils import crop_he_to_tissue
"""
from __future__ import annotations

import json
import os

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from scipy.spatial import cKDTree

DEFAULT_ROOT = "/path/to/CITEgeist_public_data/processed_files"


def compute_spot_radius(coords_xy: np.ndarray) -> float:
    """Return 0.45 x median nearest-neighbour distance in data units.

    Parameters
    ----------
    coords_xy : (N, 2) array of spot centre coordinates
    """
    tree = cKDTree(coords_xy)
    dists, _ = tree.query(coords_xy, k=2)  # k=2: self + nearest neighbour
    return 0.45 * float(np.median(dists[:, 1]))


def draw_pie_spots(
    ax: plt.Axes,
    coords_xy: np.ndarray,
    proportions: np.ndarray,
    colors: list[str],
    spot_radius: float,
) -> None:
    """Draw each spatial spot as pie wedges using a single PatchCollection.

    Uses Wedge + PatchCollection for O(1) figure overhead regardless of N.
    Do NOT use ax.pie() per spot -- it creates a new figure internally each call.

    Parameters
    ----------
    ax          : target matplotlib Axes
    coords_xy   : (N, 2) float array of spot centres in data units
    proportions : (N, K) float array; each row must sum to <=1
                  (caller should zero out fracs < 0.02 and renormalise first)
    colors      : list of K hex/named colours, one per cell type (same order
                  as proportions columns)
    spot_radius : outer radius of each wedge in data units; use
                  compute_spot_radius() to derive this from the coordinate array
    """
    all_wedges: list[mpatches.Wedge] = []
    all_facecolors: list[str] = []

    for i in range(len(coords_xy)):
        x, y = float(coords_xy[i, 0]), float(coords_xy[i, 1])
        theta1 = 0.0
        for k, frac in enumerate(proportions[i]):
            if frac < 1e-9:
                continue
            theta2 = theta1 + float(frac) * 360.0
            all_wedges.append(mpatches.Wedge((x, y), spot_radius, theta1, theta2))
            all_facecolors.append(colors[k])
            theta1 = theta2

    if all_wedges:
        pc = PatchCollection(all_wedges, match_original=False)
        pc.set_facecolor(all_facecolors)
        pc.set_linewidth(0)
        ax.add_collection(pc)

    pad = spot_radius * 2
    ax.set_xlim(float(coords_xy[:, 0].min()) - pad, float(coords_xy[:, 0].max()) + pad)
    ax.set_ylim(float(coords_xy[:, 1].min()) - pad, float(coords_xy[:, 1].max()) + pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def load_visium_he(
    specimen: str,
    processed_root: str = DEFAULT_ROOT,
) -> "tuple[np.ndarray, float, pd.DataFrame]":
    """Load a Visium SpaceRanger output and return H&E image + scaled spot coordinates.

    Parameters
    ----------
    specimen : str
        Specimen directory name under ``processed_root`` (e.g. ``"HCC22-088-P1-S1"``).
    processed_root : str
        Root directory containing per-specimen SpaceRanger ``outs/`` trees.

    Returns
    -------
    hires_img : np.ndarray, shape (H, W, 3) or (H, W, 4)
        H&E hires image loaded as float32 in [0, 1].
    scalef : float
        ``tissue_hires_scalef`` from ``scalefactors_json.json``; typically ~0.15.
    positions : pd.DataFrame
        Indexed by barcode with columns ``x_hires`` and ``y_hires`` giving each
        spot's centre in hires-image pixel space
        (``pxl_col_in_fullres * scalef`` and ``pxl_row_in_fullres * scalef``).

    Notes
    -----
    Handles both SpaceRanger v1 (no header) and v2 (with header) variants of
    ``tissue_positions.csv``. Spots must be plotted in hires-pixel coordinates
    to align with :func:`draw_he_background`.
    """
    spatial_dir = os.path.join(processed_root, specimen, "outs", "spatial")

    # --- H&E image ---
    hires_img = mpimg.imread(os.path.join(spatial_dir, "tissue_hires_image.png"))

    # --- Scale factor ---
    with open(os.path.join(spatial_dir, "scalefactors_json.json")) as fh:
        scalef = float(json.load(fh)["tissue_hires_scalef"])

    # --- Spot positions ---
    csv_path = os.path.join(spatial_dir, "tissue_positions.csv")
    # Detect header: SpaceRanger v2 starts with "barcode,..."
    with open(csv_path) as fh:
        first_line = fh.readline()
    has_header = first_line.startswith("barcode")

    _COLS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if has_header:
        positions = pd.read_csv(csv_path, index_col=0)
    else:
        positions = pd.read_csv(csv_path, header=None, names=_COLS, index_col=0)

    positions["x_hires"] = positions["pxl_col_in_fullres"] * scalef
    positions["y_hires"] = positions["pxl_row_in_fullres"] * scalef

    return hires_img, scalef, positions


def inject_spatial_metadata(
    adata: "anndata.AnnData",
    specimen: str,
    processed_root: str = DEFAULT_ROOT,
) -> str:
    """Load H&E image and scale factors into ``adata.uns['spatial']``.

    If ``adata.uns['spatial']`` already contains an entry (e.g. from
    ``sq.read.visium`` which uses the SpaceRanger slide ID as the key),
    that entry is reused instead of adding a duplicate.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object to inject metadata into (modified in-place).
    specimen : str
        Specimen directory name under ``processed_root``.
    processed_root : str
        Root directory containing per-specimen SpaceRanger ``outs/`` trees.

    Returns
    -------
    library_id : str
        The library ID key under which the spatial metadata reside — either
        the pre-existing key or ``specimen`` if a new entry was created.
    """
    if "spatial" in adata.uns and adata.uns["spatial"]:
        # SpaceRanger already populated uns['spatial'] via sq.read.visium —
        # reuse the existing entry rather than adding a duplicate.
        existing_id = list(adata.uns["spatial"].keys())[0]
        return existing_id

    spatial_dir = os.path.join(processed_root, specimen, "outs", "spatial")

    # Load hires image
    hires_img = mpimg.imread(os.path.join(spatial_dir, "tissue_hires_image.png"))

    # Load scale factors
    with open(os.path.join(spatial_dir, "scalefactors_json.json")) as fh:
        scalefactors = json.load(fh)

    if "spatial" not in adata.uns:
        adata.uns["spatial"] = {}

    adata.uns["spatial"][specimen] = {
        "images": {"hires": hires_img},
        "scalefactors": scalefactors,
    }

    return specimen


def compute_tissue_crop_from_coords(
    coords_fullres: np.ndarray,
    pad_frac: float = 0.05,
) -> "tuple[int, int, int, int]":
    """Compute a tissue bounding box from full-resolution spot coordinates.

    Filters out NaN/Inf rows, then returns a padded bounding box matching
    squidpy's ``CropCoords(x0, y0, x1, y1)`` convention.

    Parameters
    ----------
    coords_fullres : (N, 2) array
        Full-resolution pixel coordinates; column 0 = x (col), column 1 = y (row).
    pad_frac : float
        Fractional padding added to each side of the bounding box (default 0.05).

    Returns
    -------
    (left, top, right, bottom) : tuple of int
        Bounding box as ``(x0, y0, x1, y1)`` in full-resolution pixel space.

    Raises
    ------
    ValueError
        If no finite coordinates remain after filtering.
    """
    coords = np.asarray(coords_fullres, dtype=float)
    finite_mask = np.all(np.isfinite(coords), axis=1)
    coords = coords[finite_mask]

    if len(coords) == 0:
        raise ValueError("No finite coordinates found in coords_fullres")

    x_min, y_min = coords[:, 0].min(), coords[:, 1].min()
    x_max, y_max = coords[:, 0].max(), coords[:, 1].max()

    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac

    left = int(x_min - x_pad)
    top = int(y_min - y_pad)
    right = int(x_max + x_pad)
    bottom = int(y_max + y_pad)

    return left, top, right, bottom


def compute_tissue_crop(
    adata: "anndata.AnnData",
    pad_frac: float = 0.05,
) -> "tuple[int, int, int, int]":
    """Convenience wrapper: extract ``obsm['spatial']`` and call :func:`compute_tissue_crop_from_coords`.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have ``obsm['spatial']`` with full-resolution coordinates
        (column 0 = x/col, column 1 = y/row).
    pad_frac : float
        Fractional padding (default 0.05).

    Returns
    -------
    (left, top, right, bottom) : tuple of int
    """
    return compute_tissue_crop_from_coords(adata.obsm["spatial"], pad_frac=pad_frac)


def crop_he_to_tissue(
    hires_img: np.ndarray,
    scalef: float,
    crop_fullres: "tuple[int, int, int, int]",
) -> "tuple[np.ndarray, int, int]":
    """Crop an H&E hires image to the tissue bounding box.

    Parameters
    ----------
    hires_img : np.ndarray, shape (H, W, C)
        H&E hires image as returned by :func:`load_visium_he`.
    scalef : float
        ``tissue_hires_scalef`` from SpaceRanger (converts fullres -> hires pixels).
    crop_fullres : (left, top, right, bottom)
        Full-resolution bounding box as ``(x0, y0, x1, y1)`` such as returned by
        :func:`compute_tissue_crop_from_coords`.

    Returns
    -------
    cropped_img : np.ndarray
        Cropped hires image slice.
    x_offset : int
        Hires-pixel x offset of the crop origin (subtract from ``x_hires`` coords
        to get coordinates relative to ``cropped_img``).
    y_offset : int
        Hires-pixel y offset of the crop origin.
    """
    left, top, right, bottom = crop_fullres
    h, w = hires_img.shape[:2]

    x0 = max(0, int(left * scalef))
    y0 = max(0, int(top * scalef))
    x1 = min(w, int(right * scalef))
    y1 = min(h, int(bottom * scalef))

    cropped_img = hires_img[y0:y1, x0:x1]
    return cropped_img, x0, y0


def draw_he_background(ax: "plt.Axes", hires_img: np.ndarray) -> None:
    """Display H&E hires image as the bottom layer of *ax* in hires-pixel coordinates.

    After calling this function, spots and other overlays drawn on *ax* using
    coordinates produced by :func:`load_visium_he` (``x_hires``, ``y_hires``)
    will be correctly registered on the image.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes; should be freshly created or have no prior data limits set.
    hires_img : np.ndarray
        H&E image array as returned by :func:`load_visium_he` or
        :func:`crop_he_to_tissue`.

    Notes
    -----
    * Image is rendered at ``zorder=0`` so all subsequent artists appear on top.
    * Axis limits use pixel-center convention: ``xlim=(-0.5, w-0.5)``,
      ``ylim=(h-0.5, -0.5)`` (inverted, top-left origin).
    * Tick marks and spines are hidden for a clean publication panel.
    """
    h, w = hires_img.shape[:2]
    ax.imshow(hires_img, zorder=0)
    ax.set_xlim(-0.5, w - 0.5)
    ax.set_ylim(h - 0.5, -0.5)  # image convention: y increases downward
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def tighten_axes_to_spots(
    ax,
    coords_plot,
    spot_radius,
):
    """Tighten axis limits to the spot bounding box, clipping excess H&E.

    Call AFTER draw_he_background + scatter/PatchCollection to restrict
    the visible area to where spots actually exist.  H&E tissue that extends
    beyond the Visium capture grid is clipped.
    """
    pad = spot_radius * 1.5
    x_min, x_max = float(coords_plot[:, 0].min()), float(coords_plot[:, 0].max())
    y_min, y_max = float(coords_plot[:, 1].min()), float(coords_plot[:, 1].max())
    ax.set_xlim(x_min - pad, x_max + pad)
    ax.set_ylim(y_max + pad, y_min - pad)  # inverted y (image convention)
