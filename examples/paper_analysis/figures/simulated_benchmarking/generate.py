#!/usr/bin/env python3
"""
Simulated Benchmarking Figure — CITEgeist Patterns v10.

Panels:
  A — Proportion benchmark (high_seg + mixed), all methods, Pearson r ± std
  B — GEX RMSE bar chart, methods with GEX output
  C — Spatial pie maps: high_seg condition (GT + 4 methods incl. Seurat)
  D — Spatial pie maps: mixed condition (GT + 4 methods incl. Seurat)
  (v42: the tissue robustness scatter, draw_panel_C, is a supplementary figure: generate_supp_tissue_robustness.py)

v41 (env FIG_LABELSET=v41): reference-based methods only (CITEgeist is not evaluated in simulation), on the original
sim (GT replicates/, method props regen/canonical/, summary results/sim_summary.json).
  A — 9-type per-type mean r +- s.e.m. over replicates (bars; primary), with the 6-type merged per-type mean r
        +- s.e.m. overlaid as markers (secondary). Source: the per-type JSON written by
        benchmarks/simulation_benchmarking/score_calsim_per_type.py --reference-only (CALSIM_PER_TYPE_JSON,
        default results/calsim_per_type_canonical.json).
  B — GEX RMSE per replicate (Cell2Location, Tangram) from results/sim_summary.json.
  C — high-seg vs mixed 9-type per-type mean r (same JSON as A).
  D/E — pie maps (GT + the five reference-based methods) from regen/canonical/<method>/<condition>/Wu_rep_0/props.csv.

Usage:
    python generate.py               # panels + composite
    python generate.py --panels-only # panels only (fast)
"""

from __future__ import annotations

import json
import os
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

# Every sim summary this generator reads: the stamped sim_summary.json (single writer
# regen/src/build_sim_summary.py).
BENCH = PROJECT / "benchmarks" / "simulation_benchmarking" / "results"
# Env FIG_LABELSET: unset/"" = v40 (sim_summary.json, regen/canonical, replicates);
# "v41" = the same sim and summary, reference-based methods only (CITEgeist rows are not read).
FIG_LABELSET = os.environ.get("FIG_LABELSET", "")
SIM_SUMMARY_FILE = BENCH / "sim_summary.json"
SIM_SUMMARY = {
    "props_high_seg": SIM_SUMMARY_FILE,
    "props_mixed": SIM_SUMMARY_FILE,
    "gex_high_seg": SIM_SUMMARY_FILE,
    "gex_mixed": SIM_SUMMARY_FILE,
}
if FIG_LABELSET not in ("", "v41", "v42"):
    raise ValueError(f"FIG_LABELSET must be '', 'v41' or 'v42', got {FIG_LABELSET!r}")
# v42 = v41 behaviour: the simulation sources (sim_summary, calsim per-type) are the same reused, unchanged files
# (CITEgeist is not evaluated in simulation, so the v2 rollout does not touch this figure's numbers).
_V41 = FIG_LABELSET in ("v41", "v42")
CANONICAL_DIRNAME = "canonical"
REPLICATES_DIRNAME = "replicates"
# v41 per-type proportion scores (9-type primary, 6-type merged secondary), one file for all panels that plot r.
CALSIM_PER_TYPE_JSON = Path(
    os.environ.get("CALSIM_PER_TYPE_JSON", str(BENCH / "calsim_per_type_canonical_v42_blocksigma.json"))
)
# sim_summary.json names our method "CITEgeist-SACE"; figures label it "CITEgeist".
_JSON_METHOD = {"CITEgeist-SACE": "CITEgeist"}
PANELS_DIR = SCRIPT_DIR / "panels"
COMPOSITE_DIR = SCRIPT_DIR / "composite"
PANELS_DIR.mkdir(parents=True, exist_ok=True)
COMPOSITE_DIR.mkdir(parents=True, exist_ok=True)

# Methods to show. v41: reference-based methods only (CITEgeist is not evaluated in simulation).
PROP_METHODS = (["CITEgeist"] if not _V41 else []) + ["Cell2Location", "RCTD", "CARD", "Seurat", "Tangram"]
GEX_METHODS = (["CITEgeist"] if not _V41 else []) + ["Cell2Location", "Tangram"]

CITEGEIST_COLOR = METHOD_COLORS["CITEgeist"]


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def _sim_rows(json_path: Path, kind: str, condition: str) -> pd.DataFrame:
    """Per-replicate rows of one kind ("prop"/"gex") and condition from sim_summary.json."""
    if not json_path.exists():
        raise FileNotFoundError(f"Missing benchmark data: {json_path}")
    with open(json_path) as f:
        rows = json.load(f)["rows"]
    df = pd.DataFrame([r for r in rows if r["kind"] == kind and r["condition"] == condition])
    if df.empty:
        raise ValueError(f"{json_path}: no {kind} rows for condition {condition!r}")
    df["Method"] = df["method"].replace(_JSON_METHOD)
    return df.rename(columns={"replicate": "Replicate"})


def load_prop(json_path: Path, condition: str) -> pd.DataFrame:
    """Per-replicate mean Pearson r (one row per Method x Replicate)."""
    df = _sim_rows(json_path, "prop", condition)
    df = df[df["Method"].isin(PROP_METHODS)]
    return df[["Method", "Replicate", "corr"]].reset_index(drop=True)


def load_gex(json_path: Path, condition: str) -> pd.DataFrame:
    """Per-replicate Average RMSE for each GEX method (proportional-split controls excluded)."""
    df = _sim_rows(json_path, "gex", condition)
    df = df[~df["control"].astype(bool) & df["Method"].isin(GEX_METHODS)]
    df = df.rename(columns={"average_rmse": "Average RMSE"})
    return df[["Method", "Replicate", "Average RMSE"]].reset_index(drop=True)


# per-type JSON method key aliases -> figure method name
_PER_TYPE_ALIASES = {
    "CITEgeist": ("CITEgeist", "CITEgeist-SACE"),
    "Cell2Location": ("Cell2Location", "C2L"),
    "RCTD": ("RCTD",),
    "CARD": ("CARD",),
    "Seurat": ("Seurat",),
    "Tangram": ("Tangram",),
}
SIM_CONDITIONS = ("high_seg", "mixed")
PER_TYPE_FIELDS = ("per_type_r_9", "per_type_r_9_sem", "per_type_r_6", "per_type_r_6_sem", "n_reps")


def _flat_row(c: dict) -> dict[str, float]:
    """One (method, condition) cell in the flat field set, from either input shape."""
    row: dict[str, float] = {"n_reps": float(c["n_reps"])}
    for g in ("9", "6"):
        v = c[f"per_type_r_{g}"]
        if isinstance(v, dict):  # scorer payload: {mean, sd, sem, per_rep}
            row[f"per_type_r_{g}"], row[f"per_type_r_{g}_sem"] = float(v["mean"]), float(v["sem"])
        else:
            row[f"per_type_r_{g}"], row[f"per_type_r_{g}_sem"] = float(v), float(c[f"per_type_r_{g}_sem"])
    return row


def parse_per_type(d: dict, where: str = "per-type JSON") -> dict[str, dict[str, dict[str, float]]]:
    """{method: {condition: {per_type_r_9, per_type_r_9_sem, per_type_r_6, per_type_r_6_sem, n_reps}}} (pure).

    Reads either shape: the scorer's payload (score_calsim_per_type.py: methods.<Method>.<cond>.per_type_r_{9,6} =
    {mean, sd, sem, per_rep}, n_reps) or the flat manuscript-data shape (<Method>.<cond>.{per_type_r_9,
    per_type_r_9_sem, per_type_r_6, per_type_r_6_sem, n_reps}). Method keys may be CITEgeist or CITEgeist-SACE,
    Cell2Location or C2L. Every PROP_METHODS entry must be present for both conditions with finite fields and
    n_reps >= 2; nothing is filled in. No SPOTS / QP keys are read even if present.
    """
    if isinstance(d.get("methods"), dict):
        d = d["methods"]
    out: dict[str, dict[str, dict[str, float]]] = {}
    for method in PROP_METHODS:
        keys = [k for k in _PER_TYPE_ALIASES[method] if k in d]
        if len(keys) != 1:
            raise KeyError(f"{where}: expected exactly one key of {_PER_TYPE_ALIASES[method]} for {method}, got {keys}")
        block = d[keys[0]]
        out[method] = {}
        for cond in SIM_CONDITIONS:
            if cond not in block:
                raise KeyError(f"{where}: {keys[0]} has no {cond!r} block")
            row = _flat_row(block[cond])
            if not all(np.isfinite(v) for v in row.values()):
                raise ValueError(f"{where}: non-finite field in {keys[0]}.{cond}: {row}")
            if row["n_reps"] < 2:
                raise ValueError(f"{where}: {keys[0]}.{cond} n_reps={row['n_reps']} < 2")
            out[method][cond] = row
    return out


def check_per_type_against_data(per_type: dict, data: dict, ndigits: int = 3) -> None:
    """Plotted per-type r / s.e.m. must equal the sim_reference.json summary (pure)."""
    assert "SPOTS" not in data, "sim_reference.json carries a SPOTS row"
    assert "CITEgeist" not in data, "sim_reference.json carries a CITEgeist row (v41: reference-based methods only)"
    parsed = parse_per_type(data, "sim_reference.json")
    for method, conds in per_type.items():
        for cond, row in conds.items():
            for f in ("per_type_r_9", "per_type_r_9_sem", "per_type_r_6", "per_type_r_6_sem"):
                got, want = round(row[f], ndigits), round(parsed[method][cond][f], ndigits)
                assert got == want, f"{method}.{cond}.{f}: plotted {got} != sim_reference.json {want}"


_SIGMA_BLOCKS = ("gt_gex_sigma", "sigma", "gex_sigma")
DERIVED_METRICS_JSON = Path("/path/to/derived_metrics/derived_metrics.json")  # optional; used only if present


def sigma_from(per_type_raw: dict | None, derived: dict | None) -> tuple[float, float]:
    """(mean, p90) per-gene GT sigma of the simulation (pure).

    Prefers the sigma block of the per-type JSON (score_calsim_per_type.py); otherwise derived_metrics.json
    sigma_sim_{mean,p90} (the values the text reports). When both exist they must agree to 3 decimals.
    """
    block = None
    for k in _SIGMA_BLOCKS:
        b = per_type_raw.get(k) if per_type_raw else None
        if not isinstance(b, dict):
            continue
        # scorer payload: mean_expressed / p90_expressed (GT-expressed genes, the v40 SIGMA_*_EXPRESSED rule)
        for mk, pk in (("mean_expressed", "p90_expressed"), ("mean", "p90")):
            if mk in b and pk in b:
                block = (float(b[mk]), float(b[pk]))
                break
        if block is not None:
            break
    dm = None
    if derived and "sigma_sim_mean" in derived and "sigma_sim_p90" in derived:
        dm = (float(derived["sigma_sim_mean"]), float(derived["sigma_sim_p90"]))
    if block is None and dm is None:
        raise KeyError(
            f"no sim GT sigma: none of {_SIGMA_BLOCKS} in the per-type JSON and no derived_metrics sigma_sim_*"
        )
    if block is not None and dm is not None:
        assert (round(block[0], 3), round(block[1], 3)) == dm, f"sigma per-type JSON {block} != derived_metrics {dm}"
    return block if block is not None else dm


def load_sigma_v41() -> tuple[float, float]:
    raw = json.loads(CALSIM_PER_TYPE_JSON.read_text()) if CALSIM_PER_TYPE_JSON.exists() else None
    derived = json.loads(DERIVED_METRICS_JSON.read_text()) if DERIVED_METRICS_JSON.exists() else None
    return sigma_from(raw, derived)


def load_per_type_v41() -> dict[str, dict[str, dict[str, float]]]:
    """Per-type scores of the reference-based methods in simulation (gated against the manuscript-data mirror when it exists)."""
    if not CALSIM_PER_TYPE_JSON.exists():
        raise FileNotFoundError(
            f"calibrated-sim per-type scores not found: {CALSIM_PER_TYPE_JSON}\n"
            "Run benchmarks/simulation_benchmarking/score_calsim_per_type.py (or set CALSIM_PER_TYPE_JSON)."
        )
    with open(CALSIM_PER_TYPE_JSON) as fh:
        per_type = parse_per_type(json.load(fh), str(CALSIM_PER_TYPE_JSON))
    # The data-freshness gate against the text-data summary is maintainer-only and not shipped.
    return per_type


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
    gt_csv = SIM_ROOT / REPLICATES_DIRNAME / condition / "ST_sim" / "Wu_ST_0_prop.csv"
    if not gt_csv.exists():
        raise FileNotFoundError(f"GT prop CSV not found: {gt_csv}")
    gt = pd.read_csv(gt_csv, index_col="spot")
    coords = gt[["spot_x", "spot_y"]].copy()
    gt_props = _threshold_renorm(gt)

    # ── Method predictions: canonical per-replicate props (replicate 0) ──────
    # Same files the sim_summary.json stats were computed from
    # (sim_summary.json canonical_manifest; regen/canonical/<method>/<condition>/Wu_rep_<r>/props.csv).
    def _canon(method: str) -> pd.DataFrame:
        p = SIM_ROOT / "regen" / CANONICAL_DIRNAME / method / condition / "Wu_rep_0" / "props.csv"
        if not p.exists():
            raise FileNotFoundError(f"{method} canonical props CSV not found: {p}")
        return pd.read_csv(p, index_col=0)

    canon_ids = {"CITEgeist": "CITEgeist-SACE"}  # figure name -> canonical dir name (others identical)
    preds = {m: _canon(canon_ids.get(m, m)) for m in PROP_METHODS}

    # ── Align to common spot index ─────────────────────────────────────────────
    common = gt_props.index
    for df in preds.values():
        common = common.intersection(df.index)

    def _fmt(df: pd.DataFrame) -> pd.DataFrame:
        props = _threshold_renorm(df.loc[common])
        props = props.reindex(columns=CELL_TYPES_ORDERED, fill_value=0.0)
        return pd.concat([coords.loc[common], props], axis=1)

    return {"Ground Truth": _fmt(gt.loc[common]), **{m: _fmt(df) for m, df in preds.items()}}


def _draw_spatial_row(
    fig: plt.Figure,
    subplot_spec,
    condition: str,
) -> None:
    """Draw pie-map subplots (GT | each of PROP_METHODS) into subplot_spec (v41: GT + five reference-based methods)."""
    data = _load_sim_spatial(condition)
    method_keys = ["Ground Truth"] + PROP_METHODS
    method_titles = list(method_keys)

    gs_sub = gridspec.GridSpecFromSubplotSpec(
        1,
        len(method_keys),
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


def draw_panel_A_v41(ax: plt.Axes, per_type: dict | None = None) -> None:
    """v41: 9-type per-type mean r +- s.e.m. (bars), 6-type merged per-type mean r +- s.e.m. (markers)."""
    per_type = per_type if per_type is not None else load_per_type_v41()
    bar_w = 0.35
    ekw = dict(ecolor="#444", elinewidth=0.8, capsize=3, capthick=0.8)
    lows = []
    for i, method in enumerate(PROP_METHODS):
        col = METHOD_COLORS.get(method, PALETTE["neutral"])
        for off, cond, alpha, hatch in ((-bar_w / 2, "high_seg", 0.92, None), (bar_w / 2, "mixed", 0.75, "///")):
            row = per_type[method][cond]
            ax.bar(
                i + off,
                row["per_type_r_9"],
                bar_w,
                color=col,
                alpha=alpha,
                edgecolor="black" if method == "CITEgeist" else "#555",
                linewidth=1.2 if method == "CITEgeist" else 0.8,
                yerr=row["per_type_r_9_sem"],
                error_kw=ekw,
                zorder=3,
                hatch=hatch,
            )
            ax.errorbar(
                i + off,
                row["per_type_r_6"],
                yerr=row["per_type_r_6_sem"],
                fmt="D",
                markersize=3.5,
                markerfacecolor="white",
                markeredgecolor="black",
                markeredgewidth=0.8,
                ecolor="black",
                elinewidth=0.7,
                capsize=1.5,
                zorder=5,
            )
            lows += [row["per_type_r_9"] - row["per_type_r_9_sem"], row["per_type_r_6"] - row["per_type_r_6_sem"]]

    ax.set_xticks(np.arange(len(PROP_METHODS)))
    _abbr = {"Cell2Location": "C2L"}
    ax.set_xticklabels(
        [_abbr.get(METHOD_DISPLAY.get(m, m), METHOD_DISPLAY.get(m, m)) for m in PROP_METHODS],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=7.5,
    )
    ax.set_ylabel("Per-type mean r", labelpad=12)
    ymin = min(0.0, min(lows) - 0.04)
    ax.set_ylim(ymin, 1.05)
    if ymin < 0:
        ax.axhline(0, color="#444", linewidth=0.8, zorder=2)
    ax.set_title("Proportion benchmark", fontsize=9)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    place_legend(
        ax,
        handles=[
            Patch(facecolor="#888", alpha=0.92, label="High seg"),
            Patch(facecolor="#888", alpha=0.45, hatch="///", label="Mixed"),
            Line2D(
                [], [], marker="D", linestyle="none", markerfacecolor="white", markeredgecolor="black", label="6-type"
            ),
        ],
        position="bottom",
        fontsize=8,
        ncol=3,
        # Panel A is ~53.5 mm wide; the default columnspacing=2.0 leaves big gaps
        # between the three entries and clips the last ("6-type") off the panel edge.
        # Tighten spacing so the whole legend fits within the panel (v41 R1).
        columnspacing=0.6,
        handlelength=1.2,
        handletextpad=0.3,
        borderpad=0.25,
    )
    ax.yaxis.grid(True, linestyle="--", alpha=0.35, zorder=0)
    ax.set_axisbelow(True)


def draw_panel_A(ax: plt.Axes) -> None:
    """Proportion benchmark: grouped bar chart, high_seg (solid) vs mixed (hatched)."""
    if _V41:
        draw_panel_A_v41(ax)
        return
    hs = load_prop(SIM_SUMMARY["props_high_seg"], "high_seg")
    mx = load_prop(SIM_SUMMARY["props_mixed"], "mixed")

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
    hs = load_gex(SIM_SUMMARY["gex_high_seg"], "high_seg")
    mx = load_gex(SIM_SUMMARY["gex_mixed"], "mixed")

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
    if _V41:  # v41: per-gene sigma of the sim GT GEX from the per-type JSON (not the v40 literals above)
        SIGMA_MEAN, SIGMA_P90 = load_sigma_v41()
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
    # Tangram's native-projection GEX RMSE (~0.8-0.9) sits far above the other methods, so the
    # y-limit follows the tallest bar + error bar; otherwise the bar is cut off at 0.5.
    bar_top = max((hs_s["mean"] + hs_s["std"].fillna(0)).max(), (mx_s["mean"] + mx_s["std"].fillna(0)).max())
    ax.set_ylim(0, max(0.50, SIGMA_P90 + 0.07, bar_top + 0.05))
    ax.set_yticks(np.arange(0, ax.get_ylim()[1], 0.25))


def draw_panel_C(ax: plt.Axes) -> None:
    """Tissue robustness scatter: high_seg r vs mixed r per method."""
    if _V41:  # 9-type per-type mean r (the primary metric of panel A)
        per_type = load_per_type_v41()
        hs_mean = pd.Series({m: per_type[m]["high_seg"]["per_type_r_9"] for m in PROP_METHODS})
        mx_mean = pd.Series({m: per_type[m]["mixed"]["per_type_r_9"] for m in PROP_METHODS})
    else:
        hs = load_prop(SIM_SUMMARY["props_high_seg"], "high_seg")
        mx = load_prop(SIM_SUMMARY["props_mixed"], "mixed")
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
        # per-type r is bounded by 1: cap the axes there; labels are kept inside by adjust_text
        lo = min(min(xvals), min(yvals)) - 0.08
        hi = 1.0
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)

    # y = x reference line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "--", color="#aaaaaa", linewidth=0.8, zorder=1)
    # labelpad=8 clears tick labels (framework was flagging axis_title_tick=7)
    _metric = "per-type r" if _V41 else "Pearson r"
    ax.set_xlabel(f"High-seg {_metric}", labelpad=8)
    ax.set_ylabel(f"Mixed tissue {_metric}", labelpad=16)
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

    # v42: panel C (tissue robustness scatter, draw_panel_C) moved to the supplement (supplementary/generate_supp_tissue_robustness.py)
    for label, fn in [("A", draw_panel_A), ("B", draw_panel_B)]:
        with PanelContext(label, sizes_dict, PANELS_DIR) as (fig, ax):
            fn(ax)

    # C and D: spatial pie map rows (standalone — need figonly + constrained_layout=False)
    for label, condition in [("C", "high_seg"), ("D", "mixed")]:
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
