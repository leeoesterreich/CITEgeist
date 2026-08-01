#!/usr/bin/env python -u
"""
Generate breast cancer benchmark figures for CITEgeist manuscript.

Panels:
  A — Overall proportion Pearson r across 5 methods (bar + individual regions)
  B — Per-cell-type proportion heatmap (methods x cell types)
  C — GEX correlation comparison: Tangram vs Cell2Location per cell type
  D — GEX RMSE comparison: Tangram vs Cell2Location per cell type
"""

import json
import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR.parent))
from _shared.generate_utils import PanelContext, load_panel_sizes, panel_argparser, place_legend
from figure_style import METHOD_COLORS, apply_style

apply_style()
plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
        "figure.constrained_layout.use": False,
        "svg.fonttype": "none",
    }
)

BASE_DIR = Path("/path/to/CITEgeist_analysis/benchmarks/xenium_benchmarking")
RESULTS_DIR = BASE_DIR / "evaluation" / "results_breast"
OUTPUT_DIR = SCRIPT_DIR / "panels"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
COMPOSITE_DIR = SCRIPT_DIR / "composite"
COMPOSITE_DIR.mkdir(parents=True, exist_ok=True)

# Wong (2011) colorblind-safe palette — Nature recommended
# Ordered so adjacent bars are maximally distinguishable
METHOD_COLORS_BREAST = {
    "Cell2Location": "#0072B2",  # blue
    "Seurat": "#009E73",  # bluish green
    "CARD": "#CC79A7",  # reddish purple
    "RCTD": "#E69F00",  # orange
    "Tangram": "#D55E00",  # vermillion
}
HEATMAP_CMAP = "YlGnBu"

BREAST_CELL_TYPES = [
    "Cancer Epithelial",
    "Normal Epithelial",
    "Macrophages",
    "Fibroblasts",
    "T cells",
    "Endothelial",
    "Perivascular",
    "Plasma cells",
]

CT_SHORT = {
    "Cancer Epithelial": "Cancer Epi.",
    "Normal Epithelial": "Normal Epi.",
    "Macrophages": "Macrophages",
    "Fibroblasts": "Fibroblasts",
    "T cells": "T cells",
    "Endothelial": "Endothelial",
    "Perivascular": "Perivascular",
    "Plasma cells": "Plasma",
}

METHOD_ORDER = ["Cell2Location", "Seurat", "CARD", "RCTD", "Tangram"]
METHOD_SHORT = {
    "Cell2Location": "C2L",
    "Seurat": "Seurat",
    "CARD": "CARD",
    "RCTD": "RCTD",
    "Tangram": "Tangram",
}
N_REGIONS = 5


def _text_color_for_bg(val, vmin=0, vmax=1):
    norm = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
    # YlGnBu: light at low values, dark at high values
    return "white" if norm > 0.6 else "black"


def load_proportion_results():
    data = {}
    for method in METHOD_ORDER:
        path = RESULTS_DIR / f"{method}_results.json"
        if path.exists():
            with open(path) as f:
                data[method] = json.load(f)
    return data


def load_gex_results():
    path = RESULTS_DIR / "gex_comparison.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def panel_a_overall_proportions(prop_data, ax):
    methods = [m for m in METHOD_ORDER if m in prop_data]
    mean_rs = [prop_data[m]["summary"]["overall_mean_pearson_r"] for m in methods]
    region_rs = [[r["overall_pearson_r"] for r in prop_data[m]["regions"]] for m in methods]
    sem_rs = [float(np.nanstd(r, ddof=1) / np.sqrt(np.sum(~np.isnan(r)))) if len(r) > 1 else 0.0 for r in region_rs]

    x = np.arange(len(methods))
    colors = [METHOD_COLORS_BREAST[m] for m in methods]

    ax.bar(
        x,
        mean_rs,
        width=0.6,
        color=colors,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.9,
        zorder=2,
        yerr=sem_rs,
        error_kw=dict(ecolor="#333333", elinewidth=0.9, capsize=3, capthick=0.9, zorder=4),
    )

    for i, regions in enumerate(region_rs):
        jitter = np.random.default_rng(42).uniform(-0.18, 0.18, len(regions))
        ax.scatter(x[i] + jitter, regions, color="black", s=24, alpha=0.6, zorder=3, edgecolors="white", linewidths=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=30, ha="right")
    ax.set_ylabel("Pearson r", labelpad=8)
    ax.set_ylim(0, 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.02,
        0.97,
        f"n = {N_REGIONS} regions",
        transform=ax.transAxes,
        va="top",
        ha="left",
        color="#666666",
    )

    for i, mr in enumerate(mean_rs):
        ax.text(i, mr + sem_rs[i] + 0.04, f"{mr:.2f}", ha="center", va="bottom", fontweight="bold")


def panel_b_pertype_heatmap(prop_data, ax=None, cbar_ax=None):
    """Draw the per-type heatmap into *ax*.  If ax is None, used standalone (composite only)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 3), dpi=300)
        cbar_ax = None

    methods = [m for m in METHOD_ORDER if m in prop_data]
    ct_labels = [CT_SHORT[ct] for ct in BREAST_CELL_TYPES]

    matrix = np.full((len(methods), len(BREAST_CELL_TYPES)), np.nan)
    for i, m in enumerate(methods):
        s = prop_data[m]["summary"]
        for j, ct in enumerate(BREAST_CELL_TYPES):
            matrix[i, j] = s.get(f"{ct}_mean_r", np.nan)

    im = ax.imshow(matrix, cmap=HEATMAP_CMAP, aspect="auto", vmin=0, vmax=1, rasterized=True)
    ax.set_xticks(np.arange(len(ct_labels)))
    ax.set_xticklabels(ct_labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels([{"Cell2Location": "C2L"}.get(m, m) for m in methods])

    for i in range(len(methods)):
        for j in range(len(BREAST_CELL_TYPES)):
            val = matrix[i, j]
            if np.isnan(val):
                ax.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        fill=True,
                        facecolor="#f0f0f0",
                        edgecolor="white",
                        lw=0.5,
                        hatch="//",
                        linewidth=0.3,
                    )
                )
                ax.text(j, i, "n/a", ha="center", va="center", color="#999999", style="italic")
            else:
                color = _text_color_for_bg(val)
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontweight="medium")

    if standalone:
        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    else:
        cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label("Pearson r", labelpad=4)
    cbar.ax.tick_params(labelsize=9)

    if standalone:
        plt.close(fig)


def draw_panel_B_into_fig(fig, prop_data):
    """Draw Panel B heatmap with colorbar into *fig* (PanelContext figonly=True)."""
    methods = [m for m in METHOD_ORDER if m in prop_data]
    ct_labels = [CT_SHORT[ct] for ct in BREAST_CELL_TYPES]

    matrix = np.full((len(methods), len(BREAST_CELL_TYPES)), np.nan)
    for i, m in enumerate(methods):
        s = prop_data[m]["summary"]
        for j, ct in enumerate(BREAST_CELL_TYPES):
            matrix[i, j] = s.get(f"{ct}_mean_r", np.nan)

    # Use GridSpec inside the figure to place heatmap + narrow colorbar
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 0.04], wspace=0.08)
    ax_hm = fig.add_subplot(gs[0, 0])
    ax_cb = fig.add_subplot(gs[0, 1])

    im = ax_hm.imshow(matrix, cmap=HEATMAP_CMAP, aspect="auto", vmin=0, vmax=1, rasterized=True)
    ax_hm.set_xticks(np.arange(len(ct_labels)))
    ax_hm.set_xticklabels(ct_labels, rotation=45, ha="right")
    ax_hm.set_yticks(np.arange(len(methods)))
    _METHOD_ABBREV = {"Cell2Location": "C2L"}
    ax_hm.set_yticklabels([_METHOD_ABBREV.get(m, m) for m in methods], rotation=0, ha="right")
    ax_hm.tick_params(axis="y", pad=2)

    for i in range(len(methods)):
        for j in range(len(BREAST_CELL_TYPES)):
            val = matrix[i, j]
            if np.isnan(val):
                ax_hm.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        fill=True,
                        facecolor="#f0f0f0",
                        edgecolor="white",
                        lw=0.5,
                        hatch="//",
                        linewidth=0.3,
                    )
                )
                ax_hm.text(j, i, "n/a", ha="center", va="center", color="#999999", style="italic")
            else:
                color = _text_color_for_bg(val)
                ax_hm.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontweight="medium")

    cbar = fig.colorbar(im, cax=ax_cb)
    cbar.set_label("Pearson r", labelpad=4)
    cbar.ax.tick_params(labelsize=9)

    fig.subplots_adjust(bottom=0.28)


def _gex_data_for_metric(gex_data, metric):
    """Per-cell-type means and per-region arrays for *metric* (e.g. 'rmse', 'pearson_r').

    Returns (tangram_means, c2l_means, tangram_per_region, c2l_per_region) where
    *_per_region is a list (one entry per cell type) of length-n_regions arrays.
    """
    tangram_vals, c2l_vals = [], []
    tangram_per_region, c2l_per_region = [], []
    for ct in BREAST_CELL_TYPES:
        t_regions = np.array(
            [r["per_celltype"].get(ct, {}).get(metric, np.nan) for r in gex_data["Tangram"]["regions"]],
            dtype=float,
        )
        c_regions = np.array(
            [r["per_celltype"].get(ct, {}).get(metric, np.nan) for r in gex_data["Cell2Location"]["regions"]],
            dtype=float,
        )
        tangram_vals.append(np.nanmean(t_regions))
        c2l_vals.append(np.nanmean(c_regions))
        tangram_per_region.append(t_regions)
        c2l_per_region.append(c_regions)
    return tangram_vals, c2l_vals, tangram_per_region, c2l_per_region


def panel_c_gex_correlation(gex_data, ax):
    if "Tangram" not in gex_data or "Cell2Location" not in gex_data:
        print("  Panel C skipped")
        return

    tangram_rs, c2l_rs, tangram_per_region, c2l_per_region = _gex_data_for_metric(gex_data, "pearson_r")
    tangram_sem = [
        float(np.nanstd(a, ddof=1) / np.sqrt(np.sum(~np.isnan(a)))) if np.sum(~np.isnan(a)) > 1 else 0.0
        for a in tangram_per_region
    ]
    c2l_sem = [
        float(np.nanstd(a, ddof=1) / np.sqrt(np.sum(~np.isnan(a)))) if np.sum(~np.isnan(a)) > 1 else 0.0
        for a in c2l_per_region
    ]
    ct_labels = [CT_SHORT[ct] for ct in BREAST_CELL_TYPES]
    x = np.arange(len(ct_labels))
    width = 0.35
    ekw = dict(ecolor="#333333", elinewidth=0.9, capsize=2.5, capthick=0.9, zorder=4)

    ax.bar(
        x - width / 2,
        tangram_rs,
        width,
        label="Tangram",
        color=METHOD_COLORS_BREAST["Tangram"],
        alpha=0.9,
        edgecolor="white",
        yerr=tangram_sem,
        error_kw=ekw,
    )
    ax.bar(
        x + width / 2,
        c2l_rs,
        width,
        label="C2L",
        color=METHOD_COLORS_BREAST["Cell2Location"],
        alpha=0.9,
        edgecolor="white",
        yerr=c2l_sem,
        error_kw=ekw,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(ct_labels, rotation=45, ha="right")
    ax.set_ylabel("Pearson r (log1p GEX)", labelpad=8)
    ax.set_ylim(0, 0.65)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # pad=0.10 pushes legend offset to ~-0.55 to clear 45°-rotated x-tick labels
    place_legend(ax, position="bottom", ncol=2, borderpad=0.10)


def panel_d_gex_rmse(gex_data, ax):
    if "Tangram" not in gex_data or "Cell2Location" not in gex_data:
        print("  Panel D skipped")
        return

    tangram_rmse, c2l_rmse, tangram_per_region, c2l_per_region = _gex_data_for_metric(gex_data, "rmse")
    tangram_sem = [
        float(np.nanstd(a, ddof=1) / np.sqrt(np.sum(~np.isnan(a)))) if np.sum(~np.isnan(a)) > 1 else 0.0
        for a in tangram_per_region
    ]
    c2l_sem = [
        float(np.nanstd(a, ddof=1) / np.sqrt(np.sum(~np.isnan(a)))) if np.sum(~np.isnan(a)) > 1 else 0.0
        for a in c2l_per_region
    ]
    ct_labels = [CT_SHORT[ct] for ct in BREAST_CELL_TYPES]
    x = np.arange(len(ct_labels))
    width = 0.35
    ekw = dict(ecolor="#333333", elinewidth=0.9, capsize=2.5, capthick=0.9, zorder=4)

    ax.bar(
        x - width / 2,
        tangram_rmse,
        width,
        label="Tangram",
        color=METHOD_COLORS_BREAST["Tangram"],
        alpha=0.9,
        edgecolor="white",
        yerr=tangram_sem,
        error_kw=ekw,
    )
    ax.bar(
        x + width / 2,
        c2l_rmse,
        width,
        label="C2L",
        color=METHOD_COLORS_BREAST["Cell2Location"],
        alpha=0.9,
        edgecolor="white",
        yerr=c2l_sem,
        error_kw=ekw,
    )

    # ── Reference lines: per-gene σ on log1p scale of the breast Xenium GT.
    # Values precomputed (5 regions × 8 cell types × ≈280 genes, n=11,160
    # gene-by-celltype-by-region observations) — see commit message.
    # Above the dashed σ_p90 line, RMSE exceeds the spread of 90% of genes
    # in the panel; above the solid σ_median line, the error already exceeds
    # the spread of the typical gene.
    SIGMA_MEAN = 0.314
    SIGMA_P90 = 0.651
    xmin, xmax = x.min() - 0.5, x.max() + 0.5
    ax.hlines(
        SIGMA_MEAN,
        xmin,
        xmax,
        colors="#444444",
        linestyles="-",
        linewidth=0.9,
        zorder=5,
    )
    ax.hlines(
        SIGMA_P90,
        xmin,
        xmax,
        colors="#444444",
        linestyles="--",
        linewidth=0.9,
        zorder=5,
    )
    # Anchor labels at LEFT edge with ha="left" so the validator's
    # text-overflow check (which assumes text-anchor=start regardless of
    # the actual SVG anchor) doesn't trigger composite shrinking.
    # fontsize=9 → SVG "9px" → validator measures ~6.75pt (matplotlib serializes
    # the user-units fontsize as px; validator converts px×0.75 to pt). Anything
    # ≤7 here would render as ≤5.25pt and trip the 6.0pt print-minimum floor.
    # (no white bbox needed — transparent fill used instead)
    # 50% white plate restores the legibility halo behind σ-line labels while
    # letting bars show through (Issue 1/7 swung opaque → fully transparent;
    # this is the compromise: text reads cleanly, bars are not occluded).
    _label_bbox = dict(facecolor="white", alpha=0.5, edgecolor="none", pad=1.0)
    ax.text(
        xmin + 0.05,
        SIGMA_MEAN + 0.012,
        f"mean per-gene σ = {SIGMA_MEAN:.2f}",
        ha="left",
        va="bottom",
        color="#444444",
        fontsize=9,
        bbox=_label_bbox,
        zorder=6,
    )
    ax.text(
        xmin + 0.05,
        SIGMA_P90 + 0.012,
        f"90th-pct per-gene σ = {SIGMA_P90:.2f}",
        ha="left",
        va="bottom",
        color="#444444",
        fontsize=9,
        bbox=_label_bbox,
        zorder=6,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(ct_labels, rotation=45, ha="right")
    ax.set_ylabel("RMSE (log1p GEX)", labelpad=8)
    ax.set_xlim(xmin, xmax)
    # Headroom above the dashed line so its label doesn't collide with the top spine
    ax.set_ylim(0, max(1.0, SIGMA_P90 + 0.18))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # pad=0.10 pushes legend offset to ~-0.55 to clear 45°-rotated x-tick labels
    place_legend(ax, position="bottom", ncol=2, borderpad=0.10)


def composite_figure(prop_data, gex_data):
    fig = plt.figure(figsize=(11, 8.5))
    gs = gridspec.GridSpec(
        2,
        3,
        figure=fig,
        width_ratios=[1, 1.3, 0.05],
        hspace=0.50,
        wspace=0.40,
        left=0.07,
        right=0.93,
        top=0.95,
        bottom=0.10,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    cbar_ax = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a_overall_proportions(prop_data, ax=ax_a)
    panel_b_pertype_heatmap(prop_data, ax=ax_b, cbar_ax=cbar_ax)
    panel_c_gex_correlation(gex_data, ax=ax_c)
    panel_d_gex_rmse(gex_data, ax=ax_d)

    fig.savefig(COMPOSITE_DIR / "breast_benchmark_composite.pdf")
    fig.savefig(COMPOSITE_DIR / "breast_benchmark_composite.png", dpi=300)
    plt.close(fig)
    print("  Composite figure saved")


def save_panels(panel_sizes: dict | None = None) -> None:
    prop_data = load_proportion_results()
    gex_data = load_gex_results()

    # Convert list format [w_mm, h_mm] → dict format {w_mm: ..., h_mm: ...} for PanelContext
    sizes_dict: dict | None = None
    if panel_sizes:
        sizes_dict = {
            label: {"w_mm": v[0], "h_mm": v[1]} if isinstance(v, (list, tuple)) else v
            for label, v in panel_sizes.items()
        }

    with PanelContext("A", sizes_dict, OUTPUT_DIR, default_figsize=(4.5, 3.5)) as (fig, ax):
        panel_a_overall_proportions(prop_data, ax)

    with PanelContext(
        "B", sizes_dict, OUTPUT_DIR, default_figsize=(6.5, 3), figonly=True, constrained_layout=False
    ) as fig:
        draw_panel_B_into_fig(fig, prop_data)

    with PanelContext("C", sizes_dict, OUTPUT_DIR, default_figsize=(6, 3.5)) as (fig, ax):
        panel_c_gex_correlation(gex_data, ax)

    with PanelContext("D", sizes_dict, OUTPUT_DIR, default_figsize=(6, 3.5)) as (fig, ax):
        panel_d_gex_rmse(gex_data, ax)


def main() -> None:
    parser = panel_argparser(description=__doc__)
    args = parser.parse_args()
    panel_sizes = load_panel_sizes(args)

    print("Loading results...")
    prop_data = load_proportion_results()
    gex_data = load_gex_results()
    print(f"  Proportion methods: {list(prop_data.keys())}")
    print(f"  GEX methods: {list(gex_data.keys())}")

    print("\nGenerating individual panels...")
    save_panels(panel_sizes=panel_sizes)

    print("\nGenerating composite figure...")
    composite_figure(prop_data, gex_data)

    print(f"\nAll figures saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
