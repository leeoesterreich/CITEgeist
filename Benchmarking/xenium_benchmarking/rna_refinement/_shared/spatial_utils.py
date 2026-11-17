"""Pie-chart spot plotting used by rmse_split_check.part_c.

compute_spot_radius and draw_pie_spots are copied verbatim from the manuscript figure helpers (_shared/spatial_utils.py),
so the public Benchmarking tree does not depend on the manuscript tree.
"""

from __future__ import annotations

import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection
from scipy.spatial import cKDTree


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


