#!/usr/bin/env python3
"""
Simulated Benchmarking Figure — CITEgeist Patterns v10.

Panels:
  A — Proportion benchmark (high_seg + mixed), all methods, Pearson r ± std
  B — GEX RMSE bar chart, methods with GEX output
  C — Tissue robustness scatter: high_seg r vs mixed r per method
  D — Spatial pie maps: high_seg condition (GT + 4 methods incl. Seurat)
  E — Spatial pie maps: mixed condition (GT + 4 methods incl. Seurat)

Usage:
    python generate.py               # panels + composite
    python generate.py --panels-only # panels only (fast)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(SCRIPT_DIR.parent))
import matplotlib.patheffects as pe
from _shared.generate_utils import PanelContext, load_panel_sizes, panel_argparser, place_legend
from _shared.style import METHOD_COLORS, METHOD_DISPLAY, PALETTE, apply_style
from adjustText import adjust_text

apply_style()
plt.rcParams.update({"svg.fonttype": "none"})

import matplotlib.gridspec as gridspec
from _shared.spatial_utils import compute_spot_radius, draw_pie_spots

SIM_ROOT = PROJECT / "benchmarks" / "simulation_benchmarking"

# Canonical cell type order for simulation data (matches GT CSV column order)
CELL_TYPES_ORDERED = [
    "B-cells",
    "CAFs",
    "Cancer Epithelial",
    "Endothelial",
    "Myeloid",
    "Normal Epithelial",
    "PVL",
    "Plasmablasts",
    "T-cells",
]

# Colors for simulation cell types (mapped to closest biological equivalent)
SIM_CELL_COLORS = {
    "B-cells": "#6baed6",
    "CAFs": "#fe9929",
    "Cancer Epithelial": "#2171b5",
    "Endothelial": "#636363",
    "Myeloid": "#41ab5d",
    "Normal Epithelial": "#7fb8d4",
    "PVL": "#fd8d3c",
    "Plasmablasts": "#5ba3cf",
    "T-cells": "#c51b8a",
}

BENCH = PROJECT / "benchmarks" / "simulation_benchmarking" / "Figures"
PANELS_DIR = SCRIPT_DIR / "panels"
COMPOSITE_DIR = SCRIPT_DIR / "composite"
PANELS_DIR.mkdir(parents=True, exist_ok=True)
COMPOSITE_DIR.mkdir(parents=True, exist_ok=True)

# Methods to show (CITEgeist_NB excluded — only cuOPT QP)
PROP_METHODS = ["CITEgeist", "Cell2Location", "RCTD", "CARD", "Seurat", "Tangram"]
GEX_METHODS = ["CITEgeist", "Cell2Location", "Tangram"]

CITEGEIST_COLOR = METHOD_COLORS["CITEgeist"]


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def load_prop(csv_path: Path) -> pd.DataFrame:
    """Per-replicate mean Pearson r (one row per Method × Replicate)."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing benchmark data: {csv_path}")
    df = pd.read_csv(csv_path)
    df = df[df["Method"].isin(PROP_METHODS)]
    return df.groupby(["Method", "Replicate"], sort=False)["corr"].first().reset_index()


def load_gex(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing benchmark data: {csv_path}")
    df = pd.read_csv(csv_path)
    df = df[df["Method"].isin(GEX_METHODS)]
    return df[["Method", "Replicate", "Average RMSE"]].copy()


def _threshold_renorm(df: pd.DataFrame) -> pd.DataFrame:
    """Zero proportions < 0.02 and renormalise each row to sum to 1."""
    ct_cols = [c for c in CELL_TYPES_ORDERED if c in df.columns]
    d = df[ct_cols].copy()
    d[d < 0.02] = 0.0
    row_sums = d.sum(axis=1)
    return d.div(row_sums.where(row_sums > 0, 1.0), axis=0)


def _load_sim_spatial(condition: str) -> dict[str, pd.DataFrame]:
    """Load spot proportions + coordinates for four methods, one simulation condition.

    Returns a dict mapping method display name → DataFrame with columns:
        spot_x, spot_y, B-cells, CAFs, Cancer Epithelial, ...

    condition : "high_seg" or "mixed"
    Uses replicate Wu_rep_0 / Wu_ST_0 (fixed for reproducibility).

    Raises FileNotFoundError if any expected data file is missing.
    """
    # ── Ground truth (has coordinates) ────────────────────────────────────────
    gt_csv = SIM_ROOT / "replicates" / condition / "ST_sim" / "Wu_ST_0_prop.csv"
    if not gt_csv.exists():
        raise FileNotFoundError(f"GT prop CSV not found: {gt_csv}")
    gt = pd.read_csv(gt_csv, index_col="spot")
    coords = gt[["spot_x", "spot_y"]].copy()
    gt_props = _threshold_renorm(gt)

    # ── CITEgeist QP ──────────────────────────────────────────────────────────
    cg_csv = (
        SIM_ROOT
        / "CITEgeist"
        / "results_qp"
        / condition
        / condition
        / "Wu_rep_0"
        / "citegeist"
        / "Wu_rep_0_cell_prop_finetuned_results.csv"
    )
    if not cg_csv.exists():
        raise FileNotFoundError(f"CITEgeist QP CSV not found: {cg_csv}")
    cg = pd.read_csv(cg_csv, index_col=0)

    # ── Cell2Location ──────────────────────────────────────────────────────────
    c2l_csv = SIM_ROOT / "Cell2Location" / condition / "cell2location_map_0" / "cell2loc_deconv_predictions.csv"
    if not c2l_csv.exists():
        raise FileNotFoundError(f"C2L CSV not found: {c2l_csv}")
    c2l = pd.read_csv(c2l_csv, index_col="spot")

    # ── RCTD — CSV sits directly at condition level ────────────────────────────
    rctd_csv = SIM_ROOT / "RCTD" / condition / "Wu_ST_0_RCTD_deconv_predictions.csv"
    if not rctd_csv.exists():
        raise FileNotFoundError(f"RCTD CSV not found: {rctd_csv}")
    rctd = pd.read_csv(rctd_csv, index_col=0)

    # ── Seurat ─────────────────────────────────────────────────────────────────
    seurat_csv = SIM_ROOT / "Seurat" / condition / "output" / "Wu_rep_0_Seurat_deconv_predictions.csv"
    if not seurat_csv.exists():
        raise FileNotFoundError(f"Seurat CSV not found: {seurat_csv}")
    seurat = pd.read_csv(seurat_csv, index_col=0)

    # ── Tangram ───────────────────────────────────────────────────────────────
    tangram_csv = SIM_ROOT / "Tangram" / condition / "tangram_map_0" / "cell_type_proportions.csv"
    if not tangram_csv.exists():
        raise FileNotFoundError(f"Tangram CSV not found: {tangram_csv}")
    tangram = pd.read_csv(tangram_csv, index_col=0)

    # ── CARD ──────────────────────────────────────────────────────────────────
    card_csv = SIM_ROOT / "CARD" / condition / "Wu_ST_0_CARD_deconv_predictions.csv"
    if not card_csv.exists():
        raise FileNotFoundError(f"CARD CSV not found: {card_csv}")
    card = pd.read_csv(card_csv, index_col=0)

    # ── Align to common spot index ─────────────────────────────────────────────
    common = (
        gt_props.index.intersection(cg.index)
        .intersection(c2l.index)
        .intersection(rctd.index)
        .intersection(seurat.index)
        .intersection(tangram.index)
        .intersection(card.index)
    )

    def _fmt(df: pd.DataFrame) -> pd.DataFrame:
        props = _threshold_renorm(df.loc[common])
        props = props.reindex(columns=CELL_TYPES_ORDERED, fill_value=0.0)
        return pd.concat([coords.loc[common], props], axis=1)

    return {
        "Ground Truth": _fmt(gt.loc[common]),
        "CITEgeist": _fmt(cg.loc[common]),
        "Cell2Location": _fmt(c2l.loc[common]),
        "RCTD": _fmt(rctd.loc[common]),
        "CARD": _fmt(card.loc[common]),
        "Seurat": _fmt(seurat.loc[common]),
        "Tangram": _fmt(tangram.loc[common]),
    }


def _draw_spatial_row(
    fig: plt.Figure,
    subplot_spec,
    condition: str,
) -> None:
    """Draw 7 pie-map subplots (GT | CITEgeist | C2L | RCTD | CARD | Seurat | Tangram) into subplot_spec."""
    data = _load_sim_spatial(condition)
    method_keys = ["Ground Truth", "CITEgeist", "Cell2Location", "RCTD", "CARD", "Seurat", "Tangram"]
    method_titles = ["Ground Truth", "CITEgeist", "Cell2Location", "RCTD", "CARD", "Seurat", "Tangram"]

    gs_sub = gridspec.GridSpecFromSubplotSpec(
        1,
        7,
        subplot_spec=subplot_spec,
        wspace=0.04,
    )

    gt_df = data["Ground Truth"]
    coords = gt_df[["spot_x", "spot_y"]].values
    spot_radius = compute_spot_radius(coords)
    ct_colors = [SIM_CELL_COLORS[ct] for ct in CELL_TYPES_ORDERED]

    for col_i, (key, title) in enumerate(zip(method_keys, method_titles)):
        ax = fig.add_subplot(gs_sub[0, col_i])
        df = data[key]
        props = df[CELL_TYPES_ORDERED].values
        draw_pie_spots(ax, coords, props, ct_colors, spot_radius)
        ax.set_title(title, fontweight="semibold", pad=4)


# ---------------------------------------------------------------------------
# Panel draw functions (called for both standalone and composite)
# ---------------------------------------------------------------------------


def draw_panel_A(ax: plt.Axes) -> None:
    """Proportion benchmark: grouped bar chart, high_seg (solid) vs mixed (hatched)."""
    hs = load_prop(BENCH / "prop_all_metrics_highseg_combined.csv")
    mx = load_prop(BENCH / "prop_all_metrics_mixed_combined.csv")

    hs_s = hs.groupby("Method")["corr"].agg(mean="mean", std="std")
    mx_s = mx.groupby("Method")["corr"].agg(mean="mean", std="std")

    n = len(PROP_METHODS)
    bar_w = 0.35
    x = np.arange(n)
    ekw = dict(ecolor="#444", elinewidth=0.8, capsize=3, capthick=0.8)

    for i, method in enumerate(PROP_METHODS):
        col = METHOD_COLORS.get(method, PALETTE["neutral"])
        hs_mean = hs_s.loc[method, "mean"] if method in hs_s.index else 0
        hs_std = hs_s.loc[method, "std"] if method in hs_s.index else 0
        mx_mean = mx_s.loc[method, "mean"] if method in mx_s.index else 0
        mx_std = mx_s.loc[method, "std"] if method in mx_s.index else 0
        ax.bar(
            i - bar_w / 2, hs_mean, bar_w, color=col, alpha=0.92, edgecolor="#555", yerr=hs_std, error_kw=ekw, zorder=3
        )
        ax.bar(
            i + bar_w / 2,
            mx_mean,
            bar_w,
            color=col,
            alpha=0.75,
            edgecolor="#555",
            yerr=mx_std,
            error_kw=ekw,
            zorder=3,
            hatch="///",
        )

    ax.set_xticks(x)
    # Panel A is narrow (~40 mm); abbreviate the longest method name (C2L) so the
    # x-tick labels render at >=7 pt (print minimum) without colliding (R2).
    _abbr = {"Cell2Location": "C2L"}
    ax.set_xticklabels(
        [_abbr.get(METHOD_DISPLAY.get(m, m), METHOD_DISPLAY.get(m, m)) for m in PROP_METHODS],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=7.5,
    )
    ax.set_ylabel("Pearson r", labelpad=12)
    ax.set_ylim(0.45, 1.05)
    ax.set_title("Proportion benchmark", fontsize=9)
    from matplotlib.patches import Patch

    # Constrained-layout-managed bottom legend: reserves space below the
    # rotated x-tick labels and shrinks the axes within the fixed panel canvas,
    # so the legend never clips or floats far below (no manual bbox_to_anchor).
    place_legend(
        ax,
        handles=[
            Patch(facecolor="#888", alpha=0.92, label="High seg"),
            Patch(facecolor="#888", alpha=0.45, hatch="///", label="Mixed"),
        ],
        position="bottom",
        fontsize=8,
        ncol=2,
        handlelength=1.4,
        handletextpad=0.4,
        borderpad=0.3,
    )
    ax.yaxis.grid(True, linestyle="--", alpha=0.35, zorder=0)
    ax.set_axisbelow(True)


def draw_panel_B(ax: plt.Axes) -> None:
    """GEX RMSE bar chart — methods with GEX output, high_seg + mixed."""
    hs = load_gex(BENCH / "all_GEX_metrics_high_seg.csv")
    mx = load_gex(BENCH / "all_GEX_metrics_mixed.csv")

    hs_s = hs.groupby("Method")["Average RMSE"].agg(mean="mean", std="std")
    mx_s = mx.groupby("Method")["Average RMSE"].agg(mean="mean", std="std")

    n = len(GEX_METHODS)
    bar_w = 0.35
    x = np.arange(n)
    ekw = dict(ecolor="#444", elinewidth=0.8, capsize=3, capthick=0.8)

    for i, method in enumerate(GEX_METHODS):
        col = METHOD_COLORS.get(method, PALETTE["neutral"])
        hs_mean = hs_s.loc[method, "mean"] if method in hs_s.index else 0
        hs_std = hs_s.loc[method, "std"] if method in hs_s.index else 0
        mx_mean = mx_s.loc[method, "mean"] if method in mx_s.index else 0
        mx_std = mx_s.loc[method, "std"] if method in mx_s.index else 0
        ax.bar(
            i - bar_w / 2, hs_mean, bar_w, color=col, alpha=0.92, edgecolor="#555", yerr=hs_std, error_kw=ekw, zorder=3
        )
        ax.bar(
            i + bar_w / 2,
            mx_mean,
            bar_w,
            color=col,
            alpha=0.75,
            edgecolor="#555",
            yerr=mx_std,
            error_kw=ekw,
            zorder=3,
            hatch="///",
        )

    # ── Reference lines: per-gene σ on the log1p scale of the simulated GT.
    # Computed by compute_sim_gene_sigma.py.
    # σ is computed over GT-EXPRESSED genes only (sum_arr > 0, n=8,013 of 29,733),
    # matching Fig1D/3E's expressed-panel convention.  Unexpressed genes (all-zero GT)
    # are excluded because they contribute no biological variability.
    # SOLID line = mean per-gene σ over expressed genes.
    # DASHED line = 90th-pct per-gene σ over expressed genes.
    SIGMA_MEAN = 0.214  # SIGMA_MEAN_EXPRESSED (mean-of-expressed; supersedes job 9763090 median)
    SIGMA_P90 = 0.414  # SIGMA_P90_EXPRESSED from job 9763090
    xmin, xmax = -0.5, n - 0.5
    ax.hlines(SIGMA_MEAN, xmin, xmax, colors="#444444", linestyles="-", linewidth=0.9, zorder=5)
    ax.hlines(SIGMA_P90, xmin, xmax, colors="#444444", linestyles="--", linewidth=0.9, zorder=5)
    # Anchor labels at LEFT edge with ha="left" so the validator's
    # text-overflow check (which assumes text-anchor=start regardless of
    # the actual SVG anchor) doesn't trigger composite shrinking.
    # fontsize=9 -> SVG "9px" -> validator measures ~6.75pt; below 7 trips the
    # 6.0pt print-minimum floor.
    # 50% white plate restores the legibility halo behind σ-line labels while
    # letting bars show through (Issue 1/7 swung opaque → fully transparent;
    # this is the compromise: text reads cleanly, bars are not occluded).
    _label_bbox = dict(facecolor="white", alpha=0.5, edgecolor="none", pad=1.0)
    ax.text(
        xmin + 0.05,
        SIGMA_MEAN + 0.004,
        f"mean σ = {SIGMA_MEAN:.2f}",
        ha="left",
        va="bottom",
        color="#444444",
        fontsize=9,
        bbox=_label_bbox,
        zorder=6,
    )
    ax.text(
        xmin + 0.05,
        SIGMA_P90 + 0.004,
        f"90th %σ = {SIGMA_P90:.2f}",
        ha="left",
        va="bottom",
        color="#444444",
        fontsize=9,
        bbox=_label_bbox,
        zorder=6,
    )

    ax.set_xticks(x)
    ax.set_xticklabels([METHOD_DISPLAY.get(m, m) for m in GEX_METHODS], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("RMSE (log1p GEX)", labelpad=8)
    ax.set_title("GEX accuracy")
    ax.yaxis.grid(True, linestyle="--", alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    # Headroom above the dashed line (P90=0.414) so its label clears the top spine.
    # Bars top ~0.17; dashed line at 0.414 is well above them.
    ax.set_ylim(0, max(0.50, SIGMA_P90 + 0.07))


def draw_panel_C(ax: plt.Axes) -> None:
    """Tissue robustness scatter: high_seg r vs mixed r per method."""
    hs = load_prop(BENCH / "prop_all_metrics_highseg_combined.csv")
    mx = load_prop(BENCH / "prop_all_metrics_mixed_combined.csv")

    hs_mean = hs.groupby("Method")["corr"].mean()
    mx_mean = mx.groupby("Method")["corr"].mean()

    abbreviations = {"Cell2Location": "C2L"}
    texts: list[plt.Text] = []
    point_x: list[float] = []
    point_y: list[float] = []

    for method in PROP_METHODS:
        if method not in hs_mean.index or method not in mx_mean.index:
            continue

        x = hs_mean[method]
        y = mx_mean[method]
        point_x.append(x)
        point_y.append(y)
        col = METHOD_COLORS.get(method, PALETTE["neutral"])
        # CARD light orange (#fdd0a2) is unreadable as text; darken for this panel.
        if method == "CARD":
            col = "#d49a6a"

        ax.scatter(
            x,
            y,
            s=50,
            marker="o",
            facecolors=col,
            edgecolors="#111111",
            linewidths=1.0,
            zorder=3,
        )

        display = abbreviations.get(
            METHOD_DISPLAY.get(method, method),
            METHOD_DISPLAY.get(method, method),
        )

        t = ax.text(
            x,
            y,
            display,
            fontsize=8,
            color=col,
            ha="left",
            va="center",
            path_effects=[pe.withStroke(linewidth=1.0, foreground="#111111")],
        )
        texts.append(t)

    # Set axes padding BEFORE drawing the ref line and adjustText.
    xvals = [hs_mean[m] for m in PROP_METHODS if m in hs_mean.index]
    yvals = [mx_mean[m] for m in PROP_METHODS if m in mx_mean.index]
    if xvals and yvals:
        lo = min(min(xvals), min(yvals)) - 0.05
        hi = max(max(xvals), max(yvals)) + 0.50
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)

    # y = x reference line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "--", color="#aaaaaa", linewidth=0.8, zorder=1)
    # labelpad=8 clears tick labels (framework was flagging axis_title_tick=7)
    ax.set_xlabel("High-seg Pearson r", labelpad=8)
    ax.set_ylabel("Mixed tissue Pearson r", labelpad=16)
    ax.set_title("Tissue robustness")

    texts, _ = adjust_text(
        texts,
        x=point_x,
        y=point_y,
        target_x=point_x,
        target_y=point_y,
        ax=ax,
        force_explode=(4.0, 4.0),
        force_static=(1.5, 1.5),
        expand=(1.2, 1.2),
        avoid_self=True,
        ensure_inside_axes=True,
        min_arrow_len=0,
        iter_lim=2000,
    )


# ---------------------------------------------------------------------------
# Save standalone panels
# ---------------------------------------------------------------------------


def save_panels(panel_sizes: dict | None = None) -> None:
    # Convert list format [w_mm, h_mm] → dict format {w_mm: ..., h_mm: ...} for PanelContext
    sizes_dict: dict | None = None
    if panel_sizes:
        sizes_dict = {
            label: {"w_mm": v[0], "h_mm": v[1]} if isinstance(v, (list, tuple)) else v
            for label, v in panel_sizes.items()
        }

    for label, fn in [("A", draw_panel_A), ("B", draw_panel_B), ("C", draw_panel_C)]:
        with PanelContext(label, sizes_dict, PANELS_DIR) as (fig, ax):
            fn(ax)

    # D and E: spatial pie map rows (standalone — need figonly + constrained_layout=False)
    for label, condition in [("D", "high_seg"), ("E", "mixed")]:
        with PanelContext(
            label,
            sizes_dict,
            PANELS_DIR,
            default_figsize=(16, 4),
            figonly=True,
            constrained_layout=False,
        ) as fig:
            fig.patch.set_facecolor("white")
            gs = gridspec.GridSpec(1, 1, figure=fig, left=0.02, right=0.98, top=0.88, bottom=0.04)
            _draw_spatial_row(fig, gs[0, 0], condition)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = panel_argparser("Generate simulated benchmarking figure panels.")
    args = parser.parse_args()
    panel_sizes = load_panel_sizes(args)

    print("Simulated benchmarking figure")
    save_panels(panel_sizes=panel_sizes)
    print("Done.")


if __name__ == "__main__":
    main()
