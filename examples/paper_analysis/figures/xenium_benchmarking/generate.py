#!/usr/bin/env python3
"""
Xenium Benchmarking Figure — CITEgeist Patterns v8.

Panels:
  A — Method comparison bar chart (6-type + 7-type), Pearson r ± s.e.m.
  B — Per-type Pearson r heatmap (6-type), all benchmarked methods
  C — Per-spot concordance: CITEgeist predicted vs ground truth proportions
  D — Spatial proportion pie maps: one map per method, each spot subdivided by
        cell-type proportion; GT / CITEgeist / Cell2Location / RCTD / Tangram / Seurat
  E — GEX deconvolution: per-type flat Pearson r, CITEgeist SACE vs Cell2Location
  F — Gene-set matched accuracy: CITEgeist vs C2L on own gene set and on
        Cell2Location's 180 training genes

Usage:
    python generate.py                  # all panels + composite
    python generate.py --panels-only    # panels only, skip composite

Data pipeline — run upstream scripts before regenerating:

  Panels A / B (proportion benchmark JSON):
    benchmarks/xenium_benchmarking/evaluation/src/compare_all_methods.py
    → results/method_comparison_singler_6type/*.json
    → results/method_comparison_singler_7type/*.json

  Panel E (GEX comparison CSV):
    benchmarks/xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py
    → results/gex_rerun_v8/xenium_gex_comparison_7type_summary.csv

  Panel D (proportion outputs per method):
    CITEgeist  → CITEgeist/output_qp_singler_6type/
    C2L        → Cell2Location/output_singler_6type/
    RCTD       → RCTD/output_singler_6type/
    Tangram    → Tangram/output_singler_6type/
    Seurat     → Seurat/output_singler_6type/

  Panel F (GEX layer CSVs for restricted-gene comparison):
    CITEgeist SACE layers → CITEgeist/results_sace_singler_7type/
    Cell2Location layers  → Cell2Location/output_singler_7type/
    GT GEX layers         → ground_truth_singler/ground_truth_7type/ground_truth_gex/
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

SCRIPT_DIR = Path(__file__).resolve().parent
# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(SCRIPT_DIR.parent))
from _shared.generate_utils import PanelContext, load_panel_sizes, panel_argparser, place_annotation, place_legend
from _shared.spatial_utils import compute_spot_radius, draw_pie_spots
from _shared.style import CELL_TYPE_COLORS, METHOD_COLORS, METHOD_DISPLAY, PALETTE, apply_style

apply_style()
plt.rcParams.update({"svg.fonttype": "none"})

# ── Data paths ────────────────────────────────────────────────────────────────
EVAL_DIR = PROJECT / "benchmarks" / "xenium_benchmarking" / "evaluation" / "results"
DIR_6TYPE = EVAL_DIR / "method_comparison_singler_6type"
DIR_7TYPE = EVAL_DIR / "method_comparison_singler_7type"

BENCH_XENIUM = PROJECT / "benchmarks" / "xenium_benchmarking"
GT_6TYPE_DIR = BENCH_XENIUM / "ground_truth_singler" / "ground_truth_6type" / "ground_truth"
GT_7TYPE_DIR = BENCH_XENIUM / "ground_truth_singler" / "ground_truth_7type" / "ground_truth"
CG_QP_DIR = BENCH_XENIUM / "CITEgeist" / "output_qp_singler_7type"
C2L_DIR = BENCH_XENIUM / "Cell2Location" / "output_singler_6type"
RCTD_DIR = BENCH_XENIUM / "RCTD" / "output_singler_6type"
TANGRAM_DIR = BENCH_XENIUM / "Tangram" / "output_singler_6type"
SEURAT_DIR = BENCH_XENIUM / "Seurat" / "output_singler_6type"
MORPH_DIR = BENCH_XENIUM / "scResolve" / "images" / "morphology"
GEX_RERUN_CSV = EVAL_DIR / "gex_rerun_v8" / "xenium_gex_comparison_7type_summary.csv"

# Paths for Panel F gene-coverage computation
GEX_CG_LAYERS_DIR = BENCH_XENIUM / "CITEgeist" / "results_sace_singler_7type"
GEX_C2L_LAYERS_DIR = BENCH_XENIUM / "Cell2Location" / "output_singler_7type"
GEX_GT_LAYERS_DIR = BENCH_XENIUM / "ground_truth_singler" / "ground_truth_7type" / "ground_truth_gex"

PANELS_DIR = SCRIPT_DIR / "panels"
COMPOSITE_DIR = SCRIPT_DIR / "composite"
PANELS_DIR.mkdir(parents=True, exist_ok=True)
COMPOSITE_DIR.mkdir(parents=True, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
METHOD_FILES_6TYPE: Dict[str, str] = {
    "CITEgeist_6type_results": "CITEgeist",
    "Cell2Location_6type_results": "Cell2Location",
    "RCTD_6type_results": "RCTD",
    "Tangram_6type_results": "Tangram",
    "Seurat_6type_results": "Seurat",
    "CARD_results": "CARD",
}

METHOD_FILES_7TYPE: Dict[str, str] = {
    "CITEgeist_7type_results": "CITEgeist",
    "Cell2Location_7type_results": "Cell2Location",
    "RCTD_7type_results": "RCTD",
    "Tangram_7type_results": "Tangram",
    "Seurat_7type_results": "Seurat",
    "CARD_results": "CARD",
}

METHOD_ORDER = ["CITEgeist", "Cell2Location", "RCTD", "Tangram", "Seurat", "CARD"]

# Cell types shown in spatial maps (6-type benchmark labels, in row order)
SPATIAL_CELL_TYPES = ["Epithelial", "Macrophages", "T cells"]
SPATIAL_REGION = 0

# Cell-type-specific colormaps: white → cell-type color
_CT_CMAPS: Dict[str, mcolors.LinearSegmentedColormap] = {}


def _ct_cmap(cell_type: str) -> mcolors.LinearSegmentedColormap:
    """Return a white→cell-type-color LinearSegmentedColormap."""
    if cell_type not in _CT_CMAPS:
        base_color = CELL_TYPE_COLORS.get(cell_type, PALETTE["neutral"])
        _CT_CMAPS[cell_type] = mcolors.LinearSegmentedColormap.from_list(f"cmap_{cell_type}", ["white", base_color])
    return _CT_CMAPS[cell_type]


# GEX panel: CSV method names → display
GEX_METHOD_MAP = {
    "CITEgeist_SACE": ("CITEgeist", METHOD_COLORS["CITEgeist"]),
    "Cell2Location": ("Cell2Location", METHOD_COLORS["Cell2Location"]),
    "Tangram": ("Tangram", METHOD_COLORS["Tangram"]),
}
GEX_CT_DISPLAY = {
    "B_cells": "B cells",
    "CD4pos_T_cells": "CD4+ T",
    "CD8pos_T_cells": "CD8+ T",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
    "Macrophages": "Macrophages",
}


# ── Data loading ──────────────────────────────────────────────────────────────


def load_method_results(
    result_dir: Path,
    method_files: Dict[str, str],
) -> Dict[str, Dict]:
    """Load per-region Pearson r from benchmark JSON files."""
    results: Dict[str, Dict] = {}
    for stem, label in method_files.items():
        path = result_dir / f"{stem}.json"
        if not path.exists():
            continue
        with open(path) as fh:
            d = json.load(fh)
        summary = d.get("summary", {})
        regions = d.get("regions", [])

        overall_mean = float(summary.get("overall_mean_pearson_r", np.nan))
        overall_std = float(summary.get("overall_std_pearson_r", np.nan))

        per_type: Dict[str, float] = {}
        cell_types: List[str] = summary.get("cell_types", [])
        for ct in cell_types:
            key = f"{ct}_mean_r"
            if key in summary:
                per_type[ct] = float(summary[key])

        region_r = [float(r.get("overall_pearson_r", np.nan)) for r in regions]

        results[label] = {
            "overall_mean": overall_mean,
            "overall_std": overall_std,
            "per_type": per_type,
            "region_r": region_r,
        }
    return results


def _normalise_col(name: str) -> str:
    """Normalise cell type column name: lowercase, replace dots/underscores with space."""
    return name.lower().replace(".", " ").replace("_", " ").strip()


def _align_prop_cols(df: pd.DataFrame, gt_cols: List[str]) -> pd.DataFrame:
    """
    Align a proportion DataFrame's columns to match gt_cols by normalised name.
    Returns DataFrame with columns renamed/reordered to gt_cols where possible.
    """
    norm_map = {_normalise_col(c): c for c in df.columns}
    rename = {}
    for gt_col in gt_cols:
        norm_gt = _normalise_col(gt_col)
        if norm_gt in norm_map:
            rename[norm_map[norm_gt]] = gt_col
    df = df.rename(columns=rename)
    # Keep only columns present in gt_cols
    keep = [c for c in gt_cols if c in df.columns]
    return df[keep] if keep else df


def load_spatial_data(
    region: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, pd.DataFrame],]:
    """
    Load spot proportions + spatial coordinates for one Xenium region.

    Returns:
        coords      — DataFrame with spot_x, spot_y (microns), indexed by spot_id
        gt_prop     — 6-type GT proportions (spots × cell types)
        method_props — dict: method_label → proportion DataFrame (spots × cell types)
    """
    sample = f"Xenium_region_{region}"

    # Coordinates come from 7-type GT (6-type region 0 has blank coords)
    gt7 = pd.read_csv(GT_7TYPE_DIR / f"{sample}_prop.csv", index_col=0)
    coords = gt7[["spot_x", "spot_y"]].copy()

    # 6-type GT proportions
    gt6 = pd.read_csv(GT_6TYPE_DIR / f"{sample}_prop.csv", index_col=0)
    meta_cols = {"n_cells", "spot_x", "spot_y", "x", "y"}
    gt_prop = gt6.drop(columns=[c for c in gt6.columns if c in meta_cols])
    gt_cols = list(gt_prop.columns)

    def _load_prop(base_dir: Path, fname_stem: str) -> Optional[pd.DataFrame]:
        for fname in (f"{fname_stem}_deconv_predictions.csv", f"{fname_stem}_cell_prop_finetuned_results.csv"):
            p = base_dir / sample / fname
            if p.exists():
                df = pd.read_csv(p, index_col=0)
                return _align_prop_cols(df, gt_cols)
        return None

    def _load_prop_merge_tcells(base_dir: Path, fname_stem: str) -> Optional[pd.DataFrame]:
        """Load props and merge CD4+/CD8+ T cells into 'T cells' before alignment."""
        for fname in (f"{fname_stem}_deconv_predictions.csv", f"{fname_stem}_cell_prop_finetuned_results.csv"):
            p = base_dir / sample / fname
            if p.exists():
                df = pd.read_csv(p, index_col=0)
                cd4 = [c for c in df.columns if "CD4" in c and "T" in c]
                cd8 = [c for c in df.columns if "CD8" in c and "T" in c]
                if cd4 and cd8:
                    df["T cells"] = df[cd4[0]] + df[cd8[0]]
                    df = df.drop(columns=cd4 + cd8)
                return _align_prop_cols(df, gt_cols)
        return None

    def _load_prop_remap(
        base_dir: Path, fname_stem: str, col_map: Dict[str, str], sum_cols: Dict[str, List[str]] = None
    ) -> Optional[pd.DataFrame]:
        """Load props, remap + sum columns, then align to GT."""
        for fname in (f"{fname_stem}_deconv_predictions.csv", f"{fname_stem}_cell_prop_finetuned_results.csv"):
            p = base_dir / sample / fname
            if p.exists():
                df = pd.read_csv(p, index_col=0)
                if sum_cols:
                    for target, sources in sum_cols.items():
                        present = [s for s in sources if s in df.columns]
                        if present:
                            df[target] = df[present].sum(axis=1)
                            df = df.drop(columns=[s for s in present if s != target])
                df = df.rename(columns=col_map)
                return _align_prop_cols(df, gt_cols)
        return None

    method_props: Dict[str, Optional[pd.DataFrame]] = {
        "CITEgeist": _load_prop_merge_tcells(CG_QP_DIR, sample),
        "Cell2Location": _load_prop(C2L_DIR, sample),
        "RCTD": _load_prop(RCTD_DIR, sample),
        "Tangram": _load_prop_merge_tcells(TANGRAM_DIR, sample),
        "Seurat": _load_prop(SEURAT_DIR, sample),
    }

    # Align all to common spots (GT × coords × all loaded methods)
    index_sets = [coords.index, gt_prop.index]
    for df in method_props.values():
        if df is not None:
            index_sets.append(df.index)
    common = index_sets[0]
    for idx in index_sets[1:]:
        common = common.intersection(idx)

    aligned: Dict[str, pd.DataFrame] = {}
    for label, df in method_props.items():
        if df is not None:
            aligned[label] = df.loc[common]

    return (
        coords.loc[common],
        gt_prop.loc[common],
        aligned,
    )


def load_gex_rerun_data() -> Optional[pd.DataFrame]:
    """Load GEX rerun summary CSV — per-cell-type rows only (excludes ALL)."""
    if not GEX_RERUN_CSV.exists():
        return None
    df = pd.read_csv(GEX_RERUN_CSV)
    return df[df["cell_type"] != "ALL"].copy()


def load_gex_all_summary() -> Dict[str, Dict]:
    """
    Load the ALL-row summary stats from GEX_RERUN_CSV.

    Returns dict keyed by CSV method name (e.g. 'CITEgeist_SACE', 'Cell2Location'):
        flat_r, flat_r_std, rmse, rmse_std, n_genes
    """
    if not GEX_RERUN_CSV.exists():
        raise FileNotFoundError(
            f"GEX summary CSV not found: {GEX_RERUN_CSV}\n"
            "Run: benchmarks/xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py"
        )
    df = pd.read_csv(GEX_RERUN_CSV)
    result: Dict[str, Dict] = {}
    for _, row in df[df["cell_type"] == "ALL"].iterrows():
        result[row["method"]] = {
            "flat_r": float(row["flattened_r_mean"]),
            "flat_r_std": float(row["flattened_r_std"]),
            "rmse": float(row["rmse_mean"]),
            "rmse_std": float(row["rmse_std"]),
            "n_genes": int(row["n_shared_genes_mean"]),
        }
    return result


# ── Helper ────────────────────────────────────────────────────────────────────


def _method_color(label: str) -> str:
    base = label.split("_")[0]
    return METHOD_COLORS.get(base, METHOD_COLORS.get(label, PALETTE["neutral"]))


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.10,
        label,
        transform=ax.transAxes,
        fontweight="bold",
        fontfamily="sans-serif",
        va="top",
        clip_on=False,
    )


# ── Panel A: method comparison bar chart ──────────────────────────────────────


def draw_panel_A(
    ax: plt.Axes,
    results_6: Dict,
    results_7: Dict,
) -> None:
    """Combined 6-type (solid) + 7-type (hatched) bar chart on a single axis."""
    methods = [m for m in METHOD_ORDER if m in results_6 or m in results_7]
    bar_w = 0.35
    x = np.arange(len(methods))
    ekw = dict(ecolor="#444", elinewidth=0.8, capsize=3, capthick=0.8)

    for i, method in enumerate(methods):
        col = _method_color(method)

        if method in results_6:
            mean6 = results_6[method]["overall_mean"]
            sem6 = results_6[method]["overall_std"] / np.sqrt(5)
            bar6 = ax.bar(x[i] - bar_w / 2, mean6, bar_w, color=col, alpha=0.92, yerr=sem6, error_kw=ekw, zorder=3)
            if method == "CITEgeist":
                bar6[0].set_edgecolor("black")
                bar6[0].set_linewidth(1.2)

        if method in results_7:
            mean7 = results_7[method]["overall_mean"]
            sem7 = results_7[method]["overall_std"] / np.sqrt(5)
            bar7 = ax.bar(
                x[i] + bar_w / 2, mean7, bar_w, color=col, alpha=0.45, hatch="///", yerr=sem7, error_kw=ekw, zorder=3
            )
            if method == "CITEgeist":
                bar7[0].set_edgecolor("black")
                bar7[0].set_linewidth(1.2)

    ax.set_xticks(x)
    ax.set_xticklabels([METHOD_DISPLAY.get(m, m) for m in methods], rotation=35, ha="right")
    ax.set_ylabel("Pearson r", labelpad=12)
    ax.set_ylim(0, 1.02)
    ax.set_title("Xenium Benchmark\n(SingleR GT)", pad=6)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch

    # Constrained-layout-managed bottom legend (auto-clears rotated x-ticks).
    place_legend(
        ax,
        handles=[
            Patch(facecolor="#888", alpha=0.92, label="6-type"),
            Patch(facecolor="#888", alpha=0.45, hatch="///", label="7-type"),
        ],
        position="bottom",
        ncol=2,
    )


# ── Panel B: per-type heatmap ─────────────────────────────────────────────────


def draw_panel_B(ax: plt.Axes, results_6: Dict) -> None:
    """Heatmap of per-type Pearson r for 6-type benchmark, all methods."""
    methods = [m for m in METHOD_ORDER if m in results_6]
    all_types: List[str] = []
    for m in methods:
        for ct in results_6[m]["per_type"]:
            if ct not in all_types:
                all_types.append(ct)

    if not all_types:
        ax.text(
            0.5,
            0.5,
            "No per-type data",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="gray",
        )
        ax.axis("off")
        return

    matrix = np.full((len(all_types), len(methods)), np.nan)
    for j, m in enumerate(methods):
        for i, ct in enumerate(all_types):
            matrix[i, j] = results_6[m]["per_type"].get(ct, np.nan)

    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-0.3, vmax=1.0, rasterized=True)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(
        [METHOD_DISPLAY.get(m, m) for m in methods],
        rotation=35,
        ha="right",
        fontsize=9,
    )
    ax.set_yticks(range(len(all_types)))
    ax.set_yticklabels(all_types, fontsize=9)
    ax.tick_params(axis="y", pad=2)
    ax.set_title("Per-type Pearson r\n(6-type)", pad=6)

    for i in range(len(all_types)):
        for j in range(len(methods)):
            v = matrix[i, j]
            if not np.isnan(v):
                # RdBu_r: extremes are dark (red/blue), middle is light
                midpoint = (-0.3 + 1.0) / 2  # 0.35
                text_color = "white" if abs(v - midpoint) > 0.35 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", color=text_color)

    cb = plt.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Pearson r", labelpad=4)  # was overlapping cbar ticks at labelpad=0


# ── Panel C: concordance scatter ──────────────────────────────────────────────


def draw_panel_C(ax: plt.Axes) -> None:
    """
    Per-spot predicted vs GT proportion scatter for CITEgeist (region 0).

    One scatter per cell type, all overlaid on a single axis, colored by
    cell type.  Per-type Pearson r values annotated in the legend.
    """
    try:
        coords, gt_prop, method_props = load_spatial_data(SPATIAL_REGION)
        cg_prop = method_props.get("CITEgeist", gt_prop)
    except Exception as exc:
        ax.text(
            0.5,
            0.5,
            f"Data load error:\n{exc}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="red",
            wrap=True,
        )
        ax.set_title("Concordance")
        return

    common_types = [ct for ct in gt_prop.columns if ct in cg_prop.columns]

    overall_gt, overall_pred = [], []

    for ct in common_types:
        gt_vals = gt_prop[ct].values
        pred_vals = cg_prop[ct].values
        mask = np.isfinite(gt_vals) & np.isfinite(pred_vals)
        if mask.sum() < 10:
            continue

        gt_v, pred_v = gt_vals[mask], pred_vals[mask]
        r, _ = pearsonr(gt_v, pred_v)
        color = CELL_TYPE_COLORS.get(ct, PALETTE["neutral"])
        ax.scatter(gt_v, pred_v, s=4, alpha=0.5, color=color, label=f"{ct}  r={r:.2f}", rasterized=True, zorder=3)

        # Per-type fit line at full opacity
        slope, intercept = np.polyfit(gt_v, pred_v, 1)
        xfit = np.linspace(gt_v.min(), gt_v.max(), 50)
        ax.plot(xfit, slope * xfit + intercept, color=color, alpha=1.0, linewidth=1.2, zorder=4)

        overall_gt.extend(gt_v.tolist())
        overall_pred.extend(pred_v.tolist())

    # Overall r + overall trend line in black
    if overall_gt:
        r_all, _ = pearsonr(overall_gt, overall_pred)
        og = np.asarray(overall_gt)
        op = np.asarray(overall_pred)
        s_all, i_all = np.polyfit(og, op, 1)
        xall = np.linspace(og.min(), og.max(), 50)
        ax.plot(xall, s_all * xall + i_all, color="black", linewidth=1.6, zorder=6)
        place_annotation(
            ax,
            0.50,
            0.05,
            f"Overall r = {r_all:.3f}",
            ha="center",
            va="bottom",
            fontweight="bold",
            color="black",
            fontsize=9,
        )

    # Identity line
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, zorder=2)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    ax.set_xlabel("GT proportion")
    ax.set_ylabel("CITEgeist predicted", labelpad=12)
    ax.set_title("Concordance", pad=6)

    # Compact in-axes legend: color → cell type + per-type r.
    # Placed in lower-right dead space (below the scatter cloud which is
    # densest in the lower-left near 0,0).  ncol=2 keeps height small.
    handles, labels_leg = ax.get_legend_handles_labels()
    if handles:
        legend = ax.legend(
            handles,
            labels_leg,
            loc="lower right",
            fontsize=9,
            ncol=2,
            framealpha=0.85,
            edgecolor="none",
            markerscale=1.5,
            handlelength=1.0,
            handletextpad=0.4,
            columnspacing=0.8,
            borderpad=0.5,
        )
        # Keep legend inside the axes bbox so it cannot clip past the panel edge.
        legend.set_clip_on(True)

    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)


# -- Panel D: spatial comparison maps -----------------------------------------------


def draw_panel_D_into_axes(
    fig: plt.Figure,
    outer_gs: gridspec.SubplotSpec,
    region: int = SPATIAL_REGION,
) -> None:
    """
    One spatial map per method: each spot rendered as a proportional pie chart.

    Cols: Ground Truth | CITEgeist | Cell2Location | RCTD | Tangram | Seurat | legend
    White background, black border per map.
    Cell types colored by CELL_TYPE_COLORS; shared legend on the right.
    """
    try:
        coords, gt_prop, method_props = load_spatial_data(region)
    except Exception as exc:
        ax = fig.add_subplot(outer_gs)
        ax.text(
            0.5,
            0.5,
            f"Spatial data unavailable:\n{exc}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="red",
        )
        ax.axis("off")
        return

    ORDERED_METHODS = ["CITEgeist", "Cell2Location", "RCTD", "Tangram", "Seurat"]
    method_labels = ["Ground Truth"] + [m for m in ORDERED_METHODS if m in method_props]
    method_dfs = [gt_prop] + [method_props[m] for m in ORDERED_METHODS if m in method_props]
    method_title_colors = ["#444444"] + [_method_color(m) for m in ORDERED_METHODS if m in method_props]

    x_vals = coords["spot_x"].values
    y_vals = coords["spot_y"].values
    x_min, x_max = float(x_vals.min()), float(x_vals.max())
    y_min, y_max = float(y_vals.min()), float(y_vals.max())
    pad = 55.0  # one spot radius of padding

    # Collect all cell type columns present in GT
    ct_cols = list(gt_prop.columns)
    ct_colors = [CELL_TYPE_COLORS.get(ct, PALETTE["neutral"]) for ct in ct_cols]

    n_maps = len(method_labels)  # expected 6 (GT + 5 methods)
    n_cols = 3
    n_rows = -(-n_maps // n_cols)  # ceil -> 2 rows for 6 maps
    # 3 map columns + 1 legend column (legend spans all rows)
    col_ratios = [1] * n_cols + [0.7]
    # hspace must clear the bottom row's method titles. These maps are tall and
    # narrow, so set_aspect("equal") shrinks each axes box HORIZONTALLY — the
    # axes fill the full cell height and the row-2 titles would otherwise sit
    # flush against the row-1 map borders. Tuned visually against the composite.
    gs_d = gridspec.GridSpecFromSubplotSpec(
        n_rows,
        n_cols + 1,
        subplot_spec=outer_gs,
        hspace=0.28,
        wspace=0.03,
        width_ratios=col_ratios,
    )

    spot_coords_gt = np.stack([x_vals, -y_vals], axis=1)
    spot_radius = compute_spot_radius(spot_coords_gt)

    for idx, (method_name, df, title_color) in enumerate(zip(method_labels, method_dfs, method_title_colors)):
        r, c = divmod(idx, n_cols)
        ax = fig.add_subplot(gs_d[r, c])
        ax.set_facecolor("white")
        spot_coords = np.stack([x_vals, -y_vals], axis=1)
        props_array = np.column_stack([df[ct].values if ct in df.columns else np.zeros(len(x_vals)) for ct in ct_cols])
        draw_pie_spots(ax, spot_coords, props_array, ct_colors, spot_radius)
        ax.set_xlim(x_min - pad, x_max + pad)
        ax.set_ylim(-y_max - pad, -y_min + pad)
        ax.set_aspect("equal")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(0.8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(method_name, fontweight="bold", color=title_color, pad=4, fontsize=8)

    # Legend axis spans all rows in the last column
    leg_ax = fig.add_subplot(gs_d[:, n_cols])
    leg_ax.axis("off")
    handles = [mpatches.Patch(facecolor=color, label=ct, linewidth=0) for ct, color in zip(ct_cols, ct_colors)]
    leg_ax.legend(
        handles=handles,
        loc="center left",
        frameon=False,
        handlelength=1.6,
        handleheight=1.6,
        borderpad=0,
        labelspacing=0.6,
        title="Cell type",
        fontsize=7,
        title_fontsize=7,
    )


# ── Panel E: GEX deconvolution accuracy ──────────────────────────────────────


def draw_panel_E(ax: plt.Axes, gex_df: Optional[pd.DataFrame]) -> None:
    """
    Grouped bar chart: per-type flat Pearson r for CITEgeist SACE vs Cell2Location.

    Data from gex_rerun_v8 (7-type SingleR GT, common gene set).
    """
    if gex_df is None:
        ax.text(
            0.5,
            0.5,
            "GEX data not found\n(run compare_xenium_gex_methods.py)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="gray",
            style="italic",
        )
        ax.set_title("GEX correlation")
        return

    # Cell types in display order (excluding ALL)
    ct_order = [k for k in GEX_CT_DISPLAY if k in gex_df["cell_type"].values]
    x = np.arange(len(ct_order))
    methods = [m for m in GEX_METHOD_MAP if m in gex_df["method"].values]
    n_methods = len(methods)
    width = 0.35
    offsets = np.linspace(
        -(n_methods - 1) * width / 2,
        (n_methods - 1) * width / 2,
        n_methods,
    )

    # Compact legend labels prevent the bottom legend from overflowing
    # the narrow panel-E width at 7pt minimum print font.
    _LEGEND_SHORT = {"Cell2Location": "C2L", "CITEgeist": "CITEgeist", "Tangram": "Tangram"}

    n_regions = 5
    for i, method in enumerate(methods):
        label, color = GEX_METHOD_MAP[method]
        legend_label = _LEGEND_SHORT.get(label, label)
        sub = gex_df[gex_df["method"] == method].set_index("cell_type")
        vals = [float(sub.loc[ct, "flattened_r_mean"]) if ct in sub.index else np.nan for ct in ct_order]
        # SEM across the 5 Xenium regions (CSV reports std per cell-type)
        errs = [
            float(sub.loc[ct, "flattened_r_std"]) / np.sqrt(n_regions) if ct in sub.index else np.nan for ct in ct_order
        ]

        bars = ax.bar(x + offsets[i], vals, width=width, color=color, label=legend_label, alpha=0.88, zorder=3)
        ax.errorbar(x + offsets[i], vals, yerr=errs, fmt="none", color="black", capsize=2.5, linewidth=0.8, zorder=4)

        # Highlight CITEgeist bars
        if method == "CITEgeist_SACE":
            for bar in bars:
                bar.set_edgecolor("black")
                bar.set_linewidth(0.8)

    # Annotate overall means in the upper-LEFT corner of the data area
    # (sparse — bars are skewed low on x-tick "ALL"-equivalent region absent).
    # Stack vertically so the two lines don't collide with each other.
    all_row = gex_df[gex_df["cell_type"] == "ALL"]
    for i, method in enumerate(methods):
        label, color = GEX_METHOD_MAP[method]
        row = all_row[all_row["method"] == method]
        if not row.empty:
            overall = float(row["flattened_r_mean"].values[0])
            ax.text(
                0.02,
                0.80 - i * 0.06,
                f"{label}: {overall:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontweight="bold",
                color=color,
                fontsize=7,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [GEX_CT_DISPLAY.get(ct, ct) for ct in ct_order],
        rotation=35,
        ha="right",
    )
    ax.set_ylabel("Flat Pearson r", labelpad=16)  # labelpad=12 caused R4 axis_title_y_tick collision (0.59 < 0.60)
    ax.set_title("GEX correlation", pad=4)
    # Constrained-layout-managed bottom legend: reserves space below the
    # 35deg-rotated cell-type x-tick labels (long names like "Macrophages")
    # within the fixed canvas, so the 3-method legend never clips off-panel.
    place_legend(
        ax,
        position="bottom",
        ncol=3,
        fontsize=7,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.8,
        borderpad=0.2,
    )
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", zorder=2)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, 1.05)


# ── Panel for new 3E: per-cell-type GEX RMSE with per-gene σ reference lines ──


def draw_panel_gex_rmse(ax: plt.Axes, gex_df: Optional[pd.DataFrame]) -> None:
    """
    Per-cell-type GEX RMSE (log1p): CITEgeist SACE vs Cell2Location vs Tangram.

    Bars: mean ± SEM across the 5 Xenium RCC regions (n=5 per cell type),
    SEM = rmse_std / sqrt(n_regions). Data from
    xenium_gex_comparison_7type_summary.csv.

    Reference lines mark per-gene σ on the log1p scale of the Xenium 7-type GT
    (n=14,175 gene-by-celltype-by-region observations):
      solid  — mean per-gene σ ≈ 0.35 ("typical gene" spread)
      dashed — 90th-pct per-gene σ ≈ 0.76 ("most-variable-gene" spread)

    The interpretation: bars above the mean line have deconvolution error
    larger than the typical gene's biological variation across spots.
    """
    if gex_df is None:
        ax.text(
            0.5,
            0.5,
            "GEX data not found\n(run compare_xenium_gex_methods.py)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="gray",
            style="italic",
        )
        ax.set_title("GEX RMSE")
        return

    ct_order = [k for k in GEX_CT_DISPLAY if k in gex_df["cell_type"].values]
    x = np.arange(len(ct_order))
    methods = [m for m in GEX_METHOD_MAP if m in gex_df["method"].values]
    n_methods = len(methods)
    width = 0.27
    offsets = np.linspace(
        -(n_methods - 1) * width / 2,
        (n_methods - 1) * width / 2,
        n_methods,
    )

    _LEGEND_SHORT = {"Cell2Location": "C2L", "CITEgeist": "CITEgeist", "Tangram": "Tangram"}

    n_regions = 5
    for i, method in enumerate(methods):
        label, color = GEX_METHOD_MAP[method]
        legend_label = _LEGEND_SHORT.get(label, label)
        sub = gex_df[gex_df["method"] == method].set_index("cell_type")
        vals = [float(sub.loc[ct, "rmse_mean"]) if ct in sub.index else np.nan for ct in ct_order]
        errs = [float(sub.loc[ct, "rmse_std"]) / np.sqrt(n_regions) if ct in sub.index else np.nan for ct in ct_order]

        bars = ax.bar(
            x + offsets[i],
            vals,
            width=width,
            color=color,
            label=legend_label,
            alpha=0.88,
            zorder=3,
        )
        ax.errorbar(
            x + offsets[i],
            vals,
            yerr=errs,
            fmt="none",
            color="#333333",
            capsize=2.5,
            linewidth=0.8,
            zorder=4,
        )

        if method == "CITEgeist_SACE":
            for bar in bars:
                bar.set_edgecolor("black")
                bar.set_linewidth(0.8)

    # Overall (ALL) RMSE annotation in upper-left
    all_row = gex_df[gex_df["cell_type"] == "ALL"]
    for i, method in enumerate(methods):
        label, color = GEX_METHOD_MAP[method]
        row = all_row[all_row["method"] == method]
        if not row.empty:
            overall = float(row["rmse_mean"].values[0])
            ax.text(
                0.02,
                0.95 - i * 0.06,
                f"{label}: {overall:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontweight="bold",
                color=color,
                fontsize=7,
            )

    # ── σ reference lines (Xenium RCC 7-type GT, log1p scale, precomputed)
    SIGMA_MEAN = 0.351
    SIGMA_P90 = 0.763
    xmin, xmax = x.min() - 0.45, x.max() + 0.45
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
    # text-overflow check (assumes text-anchor=start regardless of the
    # actual SVG anchor) doesn't trigger composite shrinking.
    # fontsize=9 → SVG "9px" → validator measures ~6.75pt > 6.0pt floor.
    # 50% white plate restores the legibility halo behind σ-line labels while
    # letting bars show through (Issue 1/7 swung opaque → fully transparent;
    # this is the compromise: text reads cleanly, bars are not occluded).
    _label_bbox = dict(facecolor="white", alpha=0.5, edgecolor="none", pad=1.0)
    ax.text(
        xmin + 0.05,
        SIGMA_MEAN + 0.014,
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
        SIGMA_P90 + 0.014,
        f"p90 σ = {SIGMA_P90:.2f}",
        ha="left",
        va="bottom",
        color="#444444",
        fontsize=9,
        bbox=_label_bbox,
        zorder=6,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [GEX_CT_DISPLAY.get(ct, ct) for ct in ct_order],
        rotation=35,
        ha="right",
    )
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel("RMSE (log1p GEX)", labelpad=12)
    ax.set_title("GEX RMSE", pad=4)
    place_legend(
        ax,
        position="bottom",
        ncol=3,
        fontsize=7,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.8,
        borderpad=0.2,
    )
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max(1.0, SIGMA_P90 + 0.10))


# ── Panel F: gene coverage + same-gene-set accuracy ──────────────────────────


def _compute_cg_restricted_r() -> float:
    """
    Compute CITEgeist SACE flat_r when restricted to Cell2Location's gene set.

    This is the one metric not available in GEX_RERUN_CSV (which only evaluates
    each method on its own gene set). All other values come from load_gex_all_summary().

    Reads from:
        GEX_CG_LAYERS_DIR   — CITEgeist SACE layer CSVs  (spots × genes)
        GEX_C2L_LAYERS_DIR  — Cell2Location layer CSVs    (spots × genes)
        GEX_GT_LAYERS_DIR   — Xenium single-cell GT CSVs  (spots × genes)

    Raises FileNotFoundError if any required layer directory is missing.
    """
    for label, path in [
        ("CITEgeist SACE layers", GEX_CG_LAYERS_DIR),
        ("Cell2Location layers", GEX_C2L_LAYERS_DIR),
        ("GT GEX layers", GEX_GT_LAYERS_DIR),
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"{label} directory not found: {path}\n"
                "Re-run the SACE / Cell2Location benchmark to regenerate layer CSVs."
            )

    # Collect C2L gene set (upper-cased for case-insensitive matching)
    c2l_genes: set = set()
    for reg in range(5):
        layer_dir = GEX_C2L_LAYERS_DIR / f"Xenium_region_{reg}" / "layers"
        if not layer_dir.exists():
            continue
        for f in layer_dir.glob("*.csv"):
            df = pd.read_csv(f, index_col=0)
            c2l_genes |= set(df.columns.str.upper())

    if not c2l_genes:
        raise FileNotFoundError(
            f"No Cell2Location layer CSVs found under {GEX_C2L_LAYERS_DIR}.\n"
            "Re-run Cell2Location benchmark to regenerate layer outputs."
        )

    type_names = [
        "B_cells",
        "CD4pos_T_cells",
        "CD8pos_T_cells",
        "Endothelial",
        "Epithelial",
        "Fibroblasts",
        "Macrophages",
    ]
    restr_pred: List[np.ndarray] = []
    restr_gt: List[np.ndarray] = []

    for reg in range(5):
        cg_reg = GEX_CG_LAYERS_DIR / f"Xenium_region_{reg}" / "layers"
        gt_reg = GEX_GT_LAYERS_DIR / f"Xenium_region_{reg}"
        for tname in type_names:
            cg_fp = cg_reg / f"{tname}_layer.csv"
            gt_fp = gt_reg / f"{tname}_GT.csv"
            if not cg_fp.exists() or not gt_fp.exists():
                continue
            gt_df = pd.read_csv(gt_fp, index_col=0)
            gt_df.columns = gt_df.columns.str.upper()
            gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
            cg_df = pd.read_csv(cg_fp, index_col=0)
            cg_df.columns = cg_df.columns.str.upper()
            cg_df = cg_df.loc[:, ~cg_df.columns.duplicated()]
            common_spots = sorted(set(cg_df.index) & set(gt_df.index))
            if len(common_spots) < 3:
                continue
            rg = sorted(set(cg_df.columns) & set(gt_df.columns) & c2l_genes)
            if not rg:
                continue
            pv = np.log1p(cg_df.loc[common_spots, rg].values.astype(float))
            gv = np.log1p(gt_df.loc[common_spots, rg].values.astype(float))
            m = (pv.ravel() > 0) | (gv.ravel() > 0)
            restr_pred.append(pv.ravel()[m])
            restr_gt.append(gv.ravel()[m])

    if not restr_pred:
        raise RuntimeError(
            "No overlapping spots/genes found between CITEgeist layers and C2L gene set. "
            "Check that layer CSVs exist and share spot IDs."
        )

    r, _ = pearsonr(np.concatenate(restr_pred), np.concatenate(restr_gt))
    return float(r)


def draw_panel_F(ax: plt.Axes) -> None:
    """
    Panel F: CITEgeist vs Cell2Location flat_r on own gene set vs same 180-gene set.

    Illustrates that SACE's accuracy advantage holds — and widens — when both
    methods are evaluated on exactly the same genes (Cell2Location's 180 training
    genes), ruling out gene-set selection as a confound.
    """
    cg_color = METHOD_COLORS["CITEgeist"]
    c2l_color = METHOD_COLORS["Cell2Location"]

    # Values from CSV (auto-updated when compare_xenium_gex_methods.py is re-run)
    gex_summary = load_gex_all_summary()
    cg_full_r = gex_summary["CITEgeist_SACE"]["flat_r"]
    c2l_full_r = gex_summary["Cell2Location"]["flat_r"]
    cg_n = gex_summary["CITEgeist_SACE"]["n_genes"]
    c2l_n = gex_summary["Cell2Location"]["n_genes"]

    # Restricted r: live-computed from layer CSVs (not in GEX_RERUN_CSV)
    cg_restr_r = _compute_cg_restricted_r()

    x = np.array([0.0, 2.5])
    width = 0.32

    cg_vals = [cg_full_r, cg_restr_r]
    c2l_vals = [c2l_full_r, c2l_full_r]

    ax.bar(
        x - width / 2,
        cg_vals,
        width=width,
        color=cg_color,
        label="CITEgeist",
        zorder=3,
        edgecolor="black",
        linewidth=0.6,
    )
    ax.bar(x + width / 2, c2l_vals, width=width, color=c2l_color, label="C2L", zorder=3)

    # Bracket + delta over "same genes" group — placed in the clear upper
    # region (~0.92) well above the bar tops (~0.73) so it doesn't compete
    # with any other annotation for vertical space.
    y_top = 0.92
    delta = cg_restr_r - c2l_full_r
    ax.annotate(
        "",
        xy=(2.5 - width / 2, y_top),
        xytext=(2.5 + width / 2, y_top),
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.0),
    )
    ax.text(2.5, y_top + 0.018, f"Δ = +{delta:.3f}", ha="center", va="bottom", fontweight="bold", fontsize=8)

    for xpos, cg_v, c2l_v in zip(x, cg_vals, c2l_vals):
        ax.text(
            xpos - width / 2,
            cg_v + 0.03,
            f"{cg_v:.3f}",
            ha="center",
            va="bottom",
            color=cg_color,
            fontweight="bold",
            fontsize=7,
            zorder=5,
        )
        ax.text(
            xpos + width / 2,
            c2l_v + 0.03,
            f"{c2l_v:.3f}",
            ha="center",
            va="bottom",
            color=c2l_color,
            fontweight="bold",
            fontsize=7,
            zorder=5,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        ["Own genes", f"Same {c2l_n}"],
        fontsize=8,
    )
    ax.set_ylabel("Flat Pearson r", labelpad=16)
    ax.set_title("Matched genes", pad=4, fontsize=10)
    # Headroom for bar-top value labels and the bracket annotation.
    ax.set_ylim(0, 1.15)
    # External-bottom legend (below x-tick labels).
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.34),
        borderaxespad=0,
        frameon=False,
        fontsize=7,
        ncol=1,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.8,
        borderpad=0.2,
    )
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)


# ── Save helpers ──────────────────────────────────────────────────────────────


def save_panels(
    results_6: Dict,
    results_7: Dict,
    gex_df: Optional[pd.DataFrame],
    panel_sizes: Optional[Dict] = None,
) -> None:
    # Convert list format [w_mm, h_mm] → dict format {w_mm: ..., h_mm: ...} for PanelContext
    sizes_dict: Optional[Dict] = None
    if panel_sizes:
        sizes_dict = {
            label: {"w_mm": v[0], "h_mm": v[1]} if isinstance(v, (list, tuple)) else v
            for label, v in panel_sizes.items()
        }

    # A
    with PanelContext("A", sizes_dict, PANELS_DIR, default_figsize=(6, 4)) as (fig, ax):
        draw_panel_A(ax, results_6, results_7)

    # B
    with PanelContext("B", sizes_dict, PANELS_DIR, default_figsize=(8, 5)) as (fig, ax):
        draw_panel_B(ax, results_6)

    # Panel C (concordance scatter) removed — redundant with heatmap B.

    # Current panel layout (v34 E<->F relabel). Draw-function names are
    # historical and do NOT track the letter they are emitted as — the
    # (label, draw fn) pairing below is the source of truth:
    #   C = per-cell-type GEX flat Pearson r        (draw_panel_E)
    #   D = per-cell-type GEX RMSE w/ σ ref lines   (draw_panel_gex_rmse)
    #   E = SACE vs Cell2Location matched gene-set  (draw_panel_F)
    #   F = full-width spatial pie maps             (draw_panel_D_into_axes)

    # C — per-type GEX flat Pearson r
    with PanelContext("C", sizes_dict, PANELS_DIR, default_figsize=(9.5, 4)) as (fig, ax):
        draw_panel_E(ax, gex_df)

    # D — per-type GEX RMSE w/ σ reference lines
    with PanelContext("D", sizes_dict, PANELS_DIR, default_figsize=(9.5, 4)) as (fig, ax):
        draw_panel_gex_rmse(ax, gex_df)

    # F — full-width spatial pie maps; standalone figure (relabelled E->F in v34)
    with PanelContext(
        "F", sizes_dict, PANELS_DIR, default_figsize=(16, 4.5), figonly=True, constrained_layout=False
    ) as fig:
        gs_f_outer = gridspec.GridSpec(1, 1, figure=fig, left=0.04, right=0.97, top=0.88, bottom=0.04)
        draw_panel_D_into_axes(fig, gs_f_outer[0, 0])

    # E — SACE vs Cell2Location matched gene-set comparison (promoted from supp S15 in v33; relabelled F->E in v34)
    with PanelContext("E", sizes_dict, PANELS_DIR, default_figsize=(6, 4)) as (fig, ax):
        draw_panel_F(ax)


# ── CLI ───────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = panel_argparser(description=__doc__)
    args = parser.parse_args()
    panel_sizes = load_panel_sizes(args)

    print("Loading 6-type benchmark results …")
    results_6 = load_method_results(DIR_6TYPE, METHOD_FILES_6TYPE)
    print(f"  {len(results_6)} methods: {list(results_6)}")

    print("Loading 7-type benchmark results …")
    results_7 = load_method_results(DIR_7TYPE, METHOD_FILES_7TYPE)
    print(f"  {len(results_7)} methods: {list(results_7)}")

    print("Loading GEX rerun data …")
    gex_df = load_gex_rerun_data()
    if gex_df is not None:
        print(f"  Methods: {sorted(gex_df['method'].unique())}")
    else:
        print(f"  Not found at {GEX_RERUN_CSV}; panel E will be empty.")

    print("Saving individual panels …")
    save_panels(results_6, results_7, gex_df, panel_sizes=panel_sizes)
    print(f"  Written to {PANELS_DIR}")

    print("Done.")


if __name__ == "__main__":
    main()
