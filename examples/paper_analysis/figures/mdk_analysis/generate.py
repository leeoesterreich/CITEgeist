#!/usr/bin/env python3
"""
Midkine (MDK) Figure — CITEgeist Patterns v21.

Panels (8 flat, A-H):
  A — Spatial grid (MDK map + proportion maps side-by-side)
  B — COMMOT MDK-NCL compartment-level hurdle forest (A12.9)
  C — Bivariate Moran's I forest (12 specimens)
  D — ELISA MDK secretion (MCF7 vs T47D)
  E — IF micrographs (MCF7 only, top 44% crop)
  F — D538G fold-change barplot
  G — ChIP-seq heatmap

Usage:
    python generate.py               # panels only (no composite)
    python generate.py --panels-only  # same (composite removed)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.image import imread

SCRIPT_DIR = Path(__file__).resolve().parent
# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(SCRIPT_DIR.parent))
from _shared.generate_utils import PanelContext, load_panel_sizes, panel_argparser, place_annotation, place_legend
from _shared.spatial_utils import manuscript_protocol_sha12
from _shared.style import apply_style
from _shared import v2_runs

apply_style()
plt.rcParams.update({"svg.fonttype": "none"})

PANELS_DIR = SCRIPT_DIR / "panels"
PANELS_DIR.mkdir(parents=True, exist_ok=True)

MDK_FIGS = PROJECT / "mdk_analysis" / "outputs" / "figures"
MDK_IF = PROJECT / "mdk_analysis" / "outputs" / "if_analysis"

PANEL_PATHS: dict[str, Path] = {
    # v42: written by generate_panel_A_rasters.py from the stamped B5 compartments (never the old if_analysis/panel_A_grid.png)
    "A_grid": SCRIPT_DIR / "outputs" / "panel_A_grid_v42.png",
    "F": PANELS_DIR / "panel_F.jpg",
}

# CSV data sources for vector panels
HURDLE_RESULTS_JSON = Path("/path/to/midkine_data/hurdle_results.json")
MORANS_CSV = PROJECT / "mdk_analysis/outputs/if_analysis/bivariate_morans_results.csv"
# Fig 6D: recomputed from the raw plate 2 readings by scripts/elisa_table.py (4PL standard curve, blank-subtracted, per ug/uL cell protein, Welch).
# (scripts/elisa_table.py in the authoring tree; not shipped)
ELISA_SCRIPTS = Path("/path/to/elisa_scripts")
ELISA_CSV = Path("/path/to/midkine_data/elisa_plate2_wells.csv")
ELISA_YLABEL = "MDK (ng/mL per \u03bcg/\u03bcL cell protein)"
# Produced by build_er_chipseq_table.py (max MACS peak score within hg38 TSS +/- 10 kb, capped at 3100).
CHIPSEQ_CSV = SCRIPT_DIR / "er_chipseq" / "er_chipseq_peak_scores.csv"

# Competitor-GEX comparison (promoted from supplementary biology_validation → Fig 5 C/D)
# Def C (A12) MDK bivariate Moran's I sources. The def C mean (0.0606) SUPERSEDES
# the retired 0.254 from bivariate_morans_results.csv (A12.10/A12.14).
# Do NOT revert to the legacy CSV.
MORANS_DEF_C_JSON = Path("/path/to/midkine_data/morans_results.json")
# v42: stage B9 (run id from v42_run_ids.yaml; the v40 leaf name m1 does not exist under the v2 protocol).
MORANS_COMPARATORS_DEF_C_JSON = v2_runs.results_dir("B9") / "morans_results.json"
# Fig 5D source: stamped A12 pr2 program-recovery registry (all three methods, one rule).
PROGRAM_RECOVERY_JSON = v2_runs.results_dir("B11") / "program_recovery.json"
TISSUE_JSON = Path("/path/to/midkine_data/tissue.json")
METHOD_DISPLAY_CMP = {"cell2location": "Cell2Loc.", "tangram": "Tangram"}
METHOD_COLORS_CMP = {"CITEgeist": "#c51b8a", "cell2location": "#2171b5", "tangram": "#fe9929"}
PROGRAM_COLORS_CMP = {"full": "#4CAF50", "partial": "#FF9800", "absent": "#F44336"}

# Panel D: patient sample labels and colors
SAMPLE_SHORT_D = {
    "HCC22-088-P1-S1": "P1-S1 (Bx)",
    "HCC22-088-P1-S2": "P1-S2 (Sx)",
    "HCC22-088-P2-S1": "P2-S1 (Bx)",
    "HCC22-088-P2-S2": "P2-S2 (Sx)",
    "HCC22-088-P3-S1_A": "P3-S1 (Bx)",
    "HCC22-088-P3-S2": "P3-S2 (Sx)",
    "HCC22-088-P4-S2_1i_rep": "P4-S2 (Sx)",
    "HCC22-088-P5-S1": "P5-S1 (Bx)",
    "HCC22-088-P5-S2_F_rep": "P5-S2 (Sx)",
    "HCC22-088-P6-S1": "P6-S1 (Bx)",
    "HCC22-088-P6-S2_D": "P6-S2 (Sx)",
}
SAMPLE_ORDER_D = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]
PATIENT_COLORS_D = {
    "P1": "#d62728",
    "P2": "#1f77b4",
    "P3": "#2ca02c",
    "P4": "#ff7f0e",
    "P5": "#9467bd",
    "P6": "#17becf",
}


def _trim_image_whitespace(img: np.ndarray, tol: float = 0.995, pad: int = 15) -> np.ndarray:
    if img.ndim == 2:
        rgb = np.stack([img, img, img], axis=-1)
    else:
        rgb = img[..., :3]
    bright = np.all(rgb >= tol, axis=-1)
    keep = ~bright
    if not keep.any():
        return img
    ys, xs = np.where(keep)
    y0 = max(int(ys.min()) - pad, 0)
    y1 = min(int(ys.max()) + pad + 1, img.shape[0])
    x0 = max(int(xs.min()) - pad, 0)
    x1 = min(int(xs.max()) + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def _placeholder(ax: plt.Axes, title: str, msg: str = "") -> None:
    ax.set_facecolor("#f0f0f0")
    ax.text(
        0.5,
        0.5,
        f"{title}\n{msg}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        color="#555",
        style="italic",
        wrap=True,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#cccccc")


def _show_img(ax: plt.Axes, path: Path) -> None:
    if path.exists():
        img = imread(str(path))
        img = _trim_image_whitespace(img)
        h, w = img.shape[:2]
        ax.imshow(img, aspect="equal", rasterized=True)
        ax.set_xlim(0, w)
        ax.set_ylim(h, 0)  # origin='upper' so y increases downward
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    else:
        _placeholder(ax, path.name, f"(missing: {path.name})")


def draw_panel_A(ax: plt.Axes) -> None:
    """MDK expression (row 1) and secretory score (row 2) of the eligible cancer compartments, v2 (B5), three specimens."""
    if not PANEL_PATHS["A_grid"].exists():  # no placeholder: a missing raster must stop the build, not ship a grey box
        raise FileNotFoundError(f"{PANEL_PATHS['A_grid']}: run generate_panel_A_rasters.py first")
    _show_img(ax, PANEL_PATHS["A_grid"])


HURDLE_KEY = "a12_compartment_hurdle"
# (sensitivity run id, display label) in top-to-bottom order; R00 is the primary analysis.
HURDLE_FOREST_RUNS = [
    ("R00", "Primary"),
    ("R10", "Between-spot only"),
    ("R01", "Depth-only norm."),
    ("R02", "Total-count norm."),
    ("R11", "Nuclei-free"),
]
HURDLE_OR_RUNS = ("R00", "R11")


def load_hurdle_a12(path: Path | None = None) -> dict:
    """Return the ``a12_compartment_hurdle`` block of hurdle_results.json."""
    import json

    with open(path or HURDLE_RESULTS_JSON, encoding="utf-8") as f:
        data = json.load(f)
    if HURDLE_KEY not in data:
        raise KeyError(f"{HURDLE_KEY!r} missing from hurdle results JSON")
    return data[HURDLE_KEY]


def _require(d: dict, key: str, where: str):
    if not isinstance(d, dict) or key not in d or d[key] is None:
        raise KeyError(f"hurdle results: missing {where}[{key!r}]")
    return d[key]


def hurdle_panel_rows(a12: dict) -> list[dict]:
    """Forest-plot rows for Figure 5B, top to bottom (pure function; raises on any missing key).

    Each row: ``key, label, kind`` ('primary' | 'sensitivity' | 'control'), ``gmr``,
    ``lo``/``hi`` (nominal 95% interval; primary only, else None), and
    ``part_or`` (participation OR; R00 and R11 only, else None).
    The last row is the 3x-MDK positive control (``injected_gmr``).
    """
    head = _require(a12, "headline", "a12")
    sens = _require(a12, "sensitivity_gmr", "a12")
    sens_or = _require(a12, "sensitivity_participation_or", "a12")
    rows = []
    for run, label in HURDLE_FOREST_RUNS:
        row = {
            "key": run,
            "label": label,
            "kind": "primary" if run == "R00" else "sensitivity",
            "gmr": float(_require(sens, run, "sensitivity_gmr")),
            "lo": None,
            "hi": None,
            "part_or": float(_require(sens_or, run, "sensitivity_participation_or")) if run in HURDLE_OR_RUNS else None,
        }
        if run == "R00":
            row["lo"] = float(_require(head, "ci_lo", "headline"))
            row["hi"] = float(_require(head, "ci_hi", "headline"))
        rows.append(row)
    pc = _require(a12, "positive_control", "a12")
    if pc.get("status") != "ok":
        raise ValueError(f"hurdle results: positive_control status is {pc.get('status')!r}, not 'ok'")
    fold = float(_require(pc, "fold_applied", "positive_control"))
    rows.append(
        {
            "key": "PC",
            "label": f"Positive control\n({fold:g}× MDK)",
            "kind": "control",
            "gmr": float(_require(pc, "injected_gmr", "positive_control")),
            "lo": None,
            "hi": None,
            "part_or": None,
        }
    )
    return rows


def hurdle_null_band(a12: dict):
    """Shifted-label null 95% band (q025, q975) of the primary GMR; None when the negative control is not_estimable."""
    nc = _require(a12, "negative_control", "a12")
    if nc.get("status") == "not_estimable":
        return None
    s = _require(nc, "summary", "negative_control")
    return float(_require(s, "q025", "negative_control.summary")), float(
        _require(s, "q975", "negative_control.summary")
    )


def draw_panel_B(ax: plt.Axes) -> None:
    """COMMOT MDK-NCL compartment-level hurdle (A12.9), P4-S2.

    Left (main) axis: forest plot of the magnitude geometric-mean ratio (GMR,
    D538G-proxy-high / low) among active senders on a log axis, dashed line at 1.
    Primary analysis (R00) carries its nominal 95% block-bootstrap interval; the
    four sensitivity runs (between-spot only, depth-only, total-count, nuclei-free)
    are points. A gray band marks the 2.5-97.5% range of the shifted-label null;
    the open diamond is the injected 3x-MDK positive control.
    Right axis (shares the row positions): participation odds ratio for the
    primary and nuclei-free runs, dashed line at 1.
    No p-values or significance marks. Data: hurdle_results.json
    (``a12_compartment_hurdle``).
    """
    a12 = load_hurdle_a12()
    rows = hurdle_panel_rows(a12)
    band = hurdle_null_band(a12)
    head = a12["headline"]
    ink = "#333333"

    ctx = getattr(PanelContext, "_active", None)
    fonts = getattr(ctx, "_fonts", None) or {}
    ann_pt = fonts.get("annotation_pt", 9)

    fig = ax.get_figure()
    ax.set_visible(False)
    y_pos = [0, 1, 2, 3, 4, 5.6]
    ax_m = fig.add_axes([0.40, 0.25, 0.34, 0.54])
    ax_p = fig.add_axes([0.83, 0.25, 0.16, 0.54], sharey=ax_m)

    if band is not None:
        ax_m.axvspan(band[0], band[1], color="#E3E3E3", zorder=0, linewidth=0)
    ax_m.axvline(1.0, color="gray", linestyle="--", linewidth=1.0, zorder=1)
    for y, r in zip(y_pos, rows):
        if r["kind"] == "control":
            ax_m.plot(
                r["gmr"],
                y,
                marker="D",
                markersize=6,
                markerfacecolor="white",
                markeredgecolor=ink,
                markeredgewidth=1.2,
                linestyle="none",
                zorder=3,
            )
            continue
        if r["lo"] is not None:
            ax_m.plot([r["lo"], r["hi"]], [y, y], color=ink, linewidth=1.4, solid_capstyle="butt", zorder=2)
        ax_m.plot(
            r["gmr"],
            y,
            marker="o",
            markersize=6 if r["kind"] == "primary" else 5,
            color=ink,
            linestyle="none",
            zorder=3,
        )
    ax_m.set_xscale("log")
    ax_m.set_xlim(0.5, 4.0)
    ax_m.set_xticks([0.5, 1, 2, 4])
    ax_m.set_xticklabels(["0.5", "1", "2", "4"])
    ax_m.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax_m.set_yticks(y_pos)
    ax_m.set_yticklabels([r["label"] for r in rows])
    ax_m.set_ylim(6.2, -0.6)
    ax_m.set_xlabel("Sender GMR\n(high/low)", labelpad=3)
    ax_m.spines["top"].set_visible(False)
    ax_m.spines["right"].set_visible(False)
    if band is not None:
        ax_m.text(
            band[0],
            1.02,
            "Shifted-label 95%",
            transform=ax_m.get_xaxis_transform(),
            ha="left",
            va="bottom",
            fontsize=ann_pt,
            color="#666666",
        )
    fig.text(
        0.40,
        0.985,
        f"P4-S2, MDK-NCL\nn = {head['n_active_high']} vs {head['n_active_low']} active senders",
        ha="left",
        va="top",
        fontsize=ann_pt,
    )

    for y, r in zip(y_pos, rows):
        if r["part_or"] is not None:
            ax_p.plot(
                r["part_or"],
                y,
                marker="o",
                markersize=6 if r["kind"] == "primary" else 5,
                color=ink,
                linestyle="none",
                zorder=3,
            )
    ax_p.axvline(1.0, color="gray", linestyle="--", linewidth=1.0, zorder=1)
    ax_p.set_xscale("log")
    ax_p.set_xlim(0.4, 2.5)
    ax_p.set_xticks([0.5, 1, 2])
    ax_p.set_xticklabels(["0.5", "1", "2"])
    ax_p.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax_p.tick_params(axis="y", left=False, labelleft=False)
    ax_p.set_xlabel("Particip.\nOR", labelpad=3)
    ax_p.spines["top"].set_visible(False)
    ax_p.spines["right"].set_visible(False)
    ax_p.spines["left"].set_visible(False)


def _get_patient_d(sample: str) -> str:
    """Extract patient ID (P1-P6) from sample name."""
    for p in PATIENT_COLORS_D:
        if p in sample:
            return p
    return "?"


def draw_panel_C(ax: plt.Axes) -> None:
    """Bivariate Moran's I (MDK vs secretory) across 11 specimens (cohort)."""
    df = pd.read_csv(MORANS_CSV).set_index("sample").reindex(SAMPLE_ORDER_D).reset_index()
    n = len(df)
    y_pos = np.arange(n)
    colors_list = [PATIENT_COLORS_D[_get_patient_d(s)] for s in df["sample"]]

    bars = ax.barh(y_pos, df["morans_I"].fillna(0), height=0.7, color=colors_list, alpha=0.85)

    for i in range(n):
        I_val = df["morans_I"].iloc[i]
        I_val = 0 if pd.isna(I_val) else I_val
        p_val = df["p_value"].iloc[i]
        if pd.isna(p_val):
            marker, c = "nd", "gray"
        elif p_val <= 0.001:
            marker, c = "***", "black"
        elif p_val <= 0.01:
            marker, c = "**", "black"
        elif p_val <= 0.05:
            marker, c = "*", "black"
        else:
            marker, c = "ns", "gray"
        ax.text(I_val + 0.006, y_pos[i], marker, va="center", ha="left", color=c)

    valid_I = df["morans_I"].dropna()
    mean_I = valid_I.mean()
    ax.axvline(mean_I, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.text(mean_I + 0.006, n - 0.4, f"mean = {mean_I:.3f}", va="top", color="black", alpha=0.8)

    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    ax.set_yticks(y_pos)
    ax.set_yticklabels([SAMPLE_SHORT_D.get(s, s) for s in df["sample"]])
    ax.set_xlabel("Moran's I", labelpad=8)
    # Push y-tick labels away from spine to prevent left-edge clipping.
    ax.tick_params(axis="y", pad=15)
    # Right-margin headroom so significance markers and "mean = ..." text are not clipped.
    ax.set_xlim(0, valid_I.max() * 1.40)

    for i, s in enumerate(df["sample"]):
        if "P4" in s:
            bars[i].set_edgecolor("#ff7f0e")
            bars[i].set_linewidth(2)


def draw_panel_E(ax: plt.Axes) -> None:
    """ELISA MDK secretion across 8 conditions with per-replicate points."""
    from scipy.stats import ttest_ind

    if str(ELISA_SCRIPTS) not in sys.path:
        sys.path.insert(0, str(ELISA_SCRIPTS))
    import elisa_table

    df = pd.DataFrame(elisa_table.sample_values(elisa_table.load_table(ELISA_CSV)))
    df["group"] = df["line"] + " " + df["genotype"] + " " + df["e2"] + "E2"
    group_order = [
        f"{ln} {gt} {e}E2"
        for ln in elisa_table.LINES
        for gt, e in (("WT", "-"), ("WT", "+"), ("D538G", "-"), ("D538G", "+"))
    ]
    n = len(group_order)
    colors_list = ["#66C2A5", "#FC8D62", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#8DA0CB", "#E78AC3"]

    means, sems = [], []
    group_vals: list[np.ndarray] = []
    for i, grp in enumerate(group_order):
        vals = df.loc[df["group"] == grp, "mdk_norm"].values
        group_vals.append(vals)
        means.append(float(vals.mean()))
        sems.append(float(vals.std(ddof=1) / np.sqrt(len(vals))))
        rng = np.random.default_rng(seed=42 + i)
        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            color=colors_list[i],
            edgecolors="black",
            linewidths=0.5,
            s=20,
            alpha=0.7,
            zorder=3,
        )

    ax.bar(
        np.arange(n),
        means,
        yerr=sems,
        width=0.6,
        color=[colors_list[i] for i in range(n)],
        alpha=0.5,
        capsize=3,
        error_kw={"linewidth": 1.2},
    )

    ax.set_xticks(np.arange(n))
    # Per-bar tick = genotype x E2 only (short -> the leftmost 45deg label no
    # longer overhangs/clips the fixed-canvas edge); cell line moves to spanning
    # group headers below (panel D is too narrow for full per-bar labels).
    short_labels = [g.replace("MCF7 ", "").replace("T47D ", "") for g in group_order]
    ax.set_xticklabels(
        short_labels,
        rotation=45,
        ha="right",
        fontsize=7,
    )
    ax.set_ylabel(ELISA_YLABEL, labelpad=8)

    # Cell-line group headers spanning bars 0-3 (MCF7) and 4-7 (T47D), placed
    # ABOVE the frame (column-group convention) -- the short tick labels consume
    # the bottom margin, so headers go on top. x in data coords, y in axes frac.
    xax_tr = ax.get_xaxis_transform()
    for x_center, line_name in ((1.5, "MCF7"), (5.5, "T47D")):
        ax.text(
            x_center,
            1.10,
            line_name,
            transform=xax_tr,
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
            clip_on=False,
        )
        ax.plot(
            [x_center - 1.4, x_center + 1.4],
            [1.07, 1.07],
            transform=xax_tr,
            color="black",
            linewidth=0.8,
            clip_on=False,
        )

    # Significance bars: WT vs D538G within each cell line × E2 condition.
    # group_order = [MCF7 WT -E2, MCF7 WT +E2, MCF7 D538G -E2, MCF7 D538G +E2,
    #                T47D WT -E2, T47D WT +E2, T47D D538G -E2, T47D D538G +E2].
    # WT vs D538G holds ±E2 fixed → pair indices step by 2 within each cell line.
    pairs = [(0, 2), (1, 3), (4, 6), (5, 7)]
    bar_tops = [m + s for m, s in zip(means, sems)]
    data_max = max(bar_tops)
    tick = data_max * 0.04
    # Stagger bracket heights: within each cell line, pair (0,2)/(4,6) covers a
    # wider x-span than (1,3)/(5,7) and the two cross — stack the inner pair
    # higher so the brackets and bars stay visually separated. Both heights
    # measured from the global max bar top so they don't dip near taller bars.
    global_top = max(bar_tops)
    base_offset = tick * 3.5
    stagger_extra = tick * 5.5
    bar_y = [global_top + base_offset + (stagger_extra if i % 2 == 1 else 0) for i, _ in enumerate(pairs)]

    def _stars(p: float) -> str:
        if p < 1e-4:
            return "****"
        if p < 1e-3:
            return "***"
        if p < 1e-2:
            return "**"
        if p < 5e-2:
            return "*"
        return "ns"

    for (a, b), y in zip(pairs, bar_y):
        p = float(ttest_ind(group_vals[a], group_vals[b], equal_var=False).pvalue)
        ax.plot([a, a, b, b], [y, y + tick * 0.4, y + tick * 0.4, y], lw=0.8, color="black")
        label = _stars(p)
        ax.text(
            (a + b) / 2.0,
            y + tick * 0.85,
            label,
            ha="center",
            va="bottom",
            fontsize=7 if label == "ns" else 9,
            fontfamily="Arial",
        )

    ax.set_ylim(0, max(bar_y) + tick * 6.0)
    ax.set_xlim(-0.55, n - 0.45)
    # Let constrained_layout (PanelContext default) handle the bottom margin.


def draw_panel_F(ax: plt.Axes) -> None:
    """MCF7-only crop of IF montage (top 44%)."""
    path = PANEL_PATHS["F"]
    if not path.exists():
        _placeholder(ax, "IF Montage (MCF7)", f"(missing: {path.name})")
        return
    img = plt.imread(str(path))
    if img.ndim == 3 and img.shape[2] == 4:
        img = img[..., :3]
    h = img.shape[0]
    img = img[: int(h * 0.44), :, :]
    ax.imshow(img, aspect="equal", rasterized=True)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def draw_panel_G(ax: plt.Axes) -> None:
    """Secretory gene log2 fold-change (D538G/WT) in MCF7 vs T47D (GSE89888)."""
    import json as _json

    tpm_path = PROJECT / "mdk_analysis/data/GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    mapping_path = PROJECT / "mdk_analysis/data/gene_id_mapping.json"

    with open(mapping_path, encoding="utf-8") as f:
        sym_to_id = _json.load(f)
    id_to_sym = {str(v): k for k, v in sym_to_id.items()}

    tpm_df = pd.read_csv(tpm_path, sep="\t", index_col=0)
    tpm_df.index = tpm_df.index.astype(str)

    groups = {
        "MCF7_WT": ["GSM2392606", "GSM2392607", "GSM2392608", "GSM2392609"],
        "MCF7_D538G": ["GSM2392614", "GSM2392615", "GSM2392616", "GSM2392617"],
        "T47D_WT": ["GSM2392582", "GSM2392583", "GSM2392584", "GSM2392585"],
        "T47D_D538G": ["GSM2392590", "GSM2392591", "GSM2392592", "GSM2392593"],
    }

    # Secretory machinery genes + MDK
    genes = ["HSP90B1", "HSPA5", "PDIA4", "CALR", "CANX", "ATF6", "MAN1A1", "MDK"]
    gene_ids = {g: str(sym_to_id[g]) for g in genes if g in sym_to_id}

    mcf7_fc, t47d_fc = [], []
    valid_genes = []
    for g in genes:
        gid = gene_ids.get(g)
        if gid is None or gid not in tpm_df.index:
            continue
        mcf7_wt_mean = tpm_df.loc[gid, groups["MCF7_WT"]].mean()
        mcf7_d_mean = tpm_df.loc[gid, groups["MCF7_D538G"]].mean()
        t47d_wt_mean = tpm_df.loc[gid, groups["T47D_WT"]].mean()
        t47d_d_mean = tpm_df.loc[gid, groups["T47D_D538G"]].mean()
        mcf7_fc.append(np.log2(mcf7_d_mean / mcf7_wt_mean) if mcf7_wt_mean > 0 else 0)
        t47d_fc.append(np.log2(t47d_d_mean / t47d_wt_mean) if t47d_wt_mean > 0 else 0)
        valid_genes.append(g)

    n = len(valid_genes)
    y_pos = np.arange(n)
    bar_h = 0.35

    ax.barh(y_pos - bar_h / 2, mcf7_fc, height=bar_h, color="#E07B39", alpha=0.85, label="MCF7")
    ax.barh(y_pos + bar_h / 2, t47d_fc, height=bar_h, color="#C87DC8", alpha=0.85, label="T47D")

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(valid_genes, fontsize=7)
    ax.tick_params(axis="y", pad=4)
    # crowded by x-ticks + xlabel; even at bbox_to_anchor=(0.5, -0.32) the
    # legend overlapped y-axis gene tick labels (perceptual framework caught
    # legend_tick=10). Move legend ABOVE the axes — there's ~15mm of gutter
    # above panel G with no other content.
    ax.set_xlabel("log₂(D538G / WT)", labelpad=5)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=2,
        frameon=False,
        fontsize=8,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.8,
        borderpad=0.2,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()


def draw_panel_H(ax: plt.Axes) -> None:
    """ChIP-seq heatmap: secretory program genes across 8 conditions."""
    df = pd.read_csv(CHIPSEQ_CSV)

    # Sort genes: those with any ER binding (row sum > 0) first, then unbound
    cond_cols = [c for c in df.columns if c != "gene"]
    df["_row_sum"] = df[cond_cols].sum(axis=1)
    df["_has_binding"] = df["_row_sum"] > 0
    df = df.sort_values(["_has_binding", "_row_sum"], ascending=[False, False]).reset_index(drop=True)

    genes = df["gene"].values
    data = df[cond_cols].values.astype(float)
    has_binding = df["_has_binding"].values

    row_means = data.mean(axis=1, keepdims=True)
    row_stds = data.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    normalized = (data - row_means) / row_stds

    # Insert a blank visual spacer row between bound and unbound gene groups
    # so the divider is unambiguously *between* rows (no risk of overlapping
    # a y-tick label as a strikethrough). We splice a NaN row into the
    # displayed matrix only — the underlying data are unchanged.
    n_bound = int(has_binding.sum())
    if 0 < n_bound < len(genes):
        spacer = np.full((1, normalized.shape[1]), np.nan)
        display_norm = np.vstack([normalized[:n_bound], spacer, normalized[n_bound:]])
        display_genes = list(genes[:n_bound]) + [""] + list(genes[n_bound:])
        # Row-index map back into the original data array (None for spacer).
        orig_row = list(range(n_bound)) + [None] + list(range(n_bound, len(genes)))
    else:
        display_norm = normalized
        display_genes = list(genes)
        orig_row = list(range(len(genes)))

    # aspect="auto": cells widen to the panel so 3-4 digit peak scores fit inside them.
    _ = ax.imshow(display_norm, aspect="auto", cmap="RdYlBu_r", vmin=-2, vmax=2, rasterized=True)
    ax.set_xticks(range(len(cond_cols)))
    ax.set_xticklabels([c.replace("_", " ") for c in cond_cols], rotation=90, ha="center", fontsize=7)
    ax.set_yticks(range(len(display_genes)))
    ax.set_yticklabels(display_genes, fontsize=8)
    ax.tick_params(axis="y", pad=2)

    for i_disp, gi in enumerate(orig_row):
        if gi is None:
            continue
        for j in range(len(cond_cols)):
            val = int(data[gi, j])
            tc = "white" if abs(normalized[gi, j]) > 1 else "black"
            ax.text(j, i_disp, str(val), ha="center", va="center", color=tc, fontsize=8)

    # Italic group annotations placed in the left margin, vertically centered
    # within each group (no line — the blank row IS the separator).
    # clip_on=False allows the rotated text to render outside the axes boundary
    # without being clipped by the panel edge (fixes 4G y-axis clipping).
    if 0 < n_bound < len(genes):
        bound_mid = (n_bound - 1) / 2.0
        unbound_mid = n_bound + 1 + (len(genes) - n_bound - 1) / 2.0
        # x in points left of the axes (cells are no longer square, so a data-unit offset would drift).
        for y_mid, label in ((bound_mid, "ER binding"), (unbound_mid, "No binding")):
            ax.annotate(
                label,
                xy=(0, y_mid),
                xycoords=("axes fraction", "data"),
                xytext=(-46, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=7.4,
                fontstyle="italic",
                rotation=90,
                annotation_clip=False,
            )


def _citegeist_mdk_morans():
    """CITEgeist MDK bivariate Moran's I: (mean, sem, per-sample valid df).

    Def C (A12) source (mean I = 0.0606). The retired 0.254 mean from
    mdk_analysis/outputs/if_analysis/bivariate_morans_results.csv is superseded
    and deliberately NOT used.
    """
    import json

    data = json.loads(MORANS_DEF_C_JSON.read_text())
    valid = pd.DataFrame([s for s in data["specimens"] if s["status"] == "ok"])[["sample", "I"]].rename(
        columns={"I": "morans_I"}
    )
    vals = valid["morans_I"].values
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(len(vals))), valid


def draw_panel_morans_methods(ax: plt.Axes) -> None:
    """Panel C — MDK bivariate Moran's I by method (CITEgeist / C2L / Tangram, shared n=10 estimable set; Def C)."""
    import json

    comp_data = json.loads(MORANS_COMPARATORS_DEF_C_JSON.read_text())
    shared = comp_data["comparators"]["shared"]
    morans = pd.DataFrame(
        [
            {"gene": "MDK", "method": m, "sample": sid, "I": sp[m]["I"]}
            for sid, sp in shared.items()
            for m in sp
            if isinstance(sp[m], dict) and "I" in sp[m]
        ]
    )
    cg_I, cg_sem, _ = _citegeist_mdk_morans()
    mdk = morans[morans["gene"] == "MDK"]
    stats = mdk.groupby("method")["I"].agg(["mean", "sem"]).reset_index()

    labels, means, sems, colors = [], [], [], []
    for m in ["CITEgeist", "cell2location", "tangram"]:
        if m == "CITEgeist":
            labels.append("CITEgeist")
            means.append(cg_I)
            sems.append(cg_sem)
            colors.append(METHOD_COLORS_CMP["CITEgeist"])
            continue
        row = stats[stats["method"] == m]
        if row.empty:
            continue
        labels.append(METHOD_DISPLAY_CMP.get(m, m))
        means.append(float(row["mean"].iloc[0]))
        sems.append(float(row["sem"].iloc[0]))
        colors.append(METHOD_COLORS_CMP[m])

    x = np.arange(len(labels))
    bars = ax.bar(x, means, width=0.55, color=colors, edgecolor="black", linewidth=0.6)
    for i, (mv, sv) in enumerate(zip(means, sems)):
        if sv > 0:
            ax.errorbar(x[i], mv, yerr=sv, fmt="none", color="black", capsize=4, linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("MDK bivariate Moran's I", labelpad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, mv, sv in zip(bars, means, sems):
        ax.text(bar.get_x() + bar.get_width() / 2, mv + sv + 0.004, f"{mv:.3f}", ha="center", va="bottom")

    _top = max(
        m + s for m, s in zip(means, sems)
    )  # clear the tallest error bar (and its value label) with the first bracket
    y_max = _top * 1.14
    step = _top * 0.16
    off = 0.0
    inference = comp_data["comparators"][
        "inference"
    ]  # stamped B9: patient-level sign-flip, Holm over the 2 comparators (two-sided)
    for comp in ["cell2location", "tangram"]:
        disp = METHOD_DISPLAY_CMP.get(comp, comp)
        if disp not in labels:
            continue
        ci = labels.index(disp)
        inf = inference[comp]
        sig = f"Holm p = {inf['holm_p']:.2f} ({inf['n_patients']} patients)"
        by = y_max + off
        ax.plot([x[0], x[0], x[ci], x[ci]], [by - 0.005, by, by, by - 0.005], color="black", linewidth=0.8)
        ax.text((x[0] + x[ci]) / 2, by + 0.005, sig, ha="center", va="bottom")
        off += step
    ax.set_ylim(0, y_max + off + step * 0.5)


# Moran's I robustness (Fig 5C): mean real I beside the two matched-control means, per adjustment form (B22, descriptive).
# (plain-English label, key of tissue.json::moran_robustness.arms). Forms: D1 = intercept + spot depth; D2 = depth deciles;
# D4 = depth + tumor content (qp_prop); D5 = depth + RNA per cell (hurdle.covariate_matrix columns 0,2 / 0,1).
ROBUSTNESS_FORMS = [
    ("No adjustment", "check11|full_counts_none"),
    ("Depth", "check11|full_counts_D1"),
    ("Depth\n(deciles)", "check11|full_counts_D2"),
    ("Depth +\ntumor content", "check13|full_counts_D4"),
    ("Depth +\nRNA per cell", "check13|full_counts_D5"),
    ("15 neighbors", "check15|none_k15"),
]
ROBUSTNESS_COLORS = {"real": "#c51b8a", "ctrl_a": "#9e9e9e", "ctrl_b": "#d9d9d9"}


def moran_robustness_rows(tissue: dict, forms=None) -> list:
    """[(label, mean_real, mean_ctrl_a, mean_ctrl_b, n_specimens)] from tissue.json (pure; raises on any gap).

    Checks the B22 provenance (stage, descriptive_only, k_primary == 6) and that every plotted arm has a mean for the real
    statistic and both matched-control families on all of its specimens. Plots means, never ratios.
    """
    mr = tissue["moran_robustness"]
    if mr.get("stage") != "B22" or mr.get("descriptive_only") is not True or mr.get("k_primary") != 6:
        raise ValueError(
            f"moran_robustness provenance is not B22 / descriptive / k=6: { {k: mr.get(k) for k in ('stage', 'descriptive_only', 'k_primary')} }"
        )
    rows = []
    for label, key in forms or ROBUSTNESS_FORMS:
        arm = mr["arms"][key]
        vals = [arm.get("mean_real"), arm.get("mean_ctrl_a"), arm.get("mean_ctrl_b")]
        if any(v is None for v in vals) or arm.get("n_with_controls") != arm.get("n_specimens"):
            raise ValueError(
                f"moran_robustness arm {key!r}: missing mean or controls on {arm.get('n_with_controls')}/{arm.get('n_specimens')} specimens"
            )
        rows.append((label, float(vals[0]), float(vals[1]), float(vals[2]), int(arm["n_specimens"])))
    if not rows:
        raise ValueError("no robustness forms to plot")
    return rows


def draw_panel_moran_robustness(ax: plt.Axes) -> None:
    """Fig 5C: mean MDK-secretory bivariate Moran's I (real) next to the two matched random-gene control means, by form."""
    rows = moran_robustness_rows(json.loads(TISSUE_JSON.read_text()))
    n = {r[4] for r in rows}
    x = np.arange(len(rows))
    w = 0.26
    series = [
        ("Real", 1, ROBUSTNESS_COLORS["real"]),
        ("Random gene for MDK", 2, ROBUSTNESS_COLORS["ctrl_a"]),
        ("Random genes for secretory", 3, ROBUSTNESS_COLORS["ctrl_b"]),
    ]
    for i, (lab, col, color) in enumerate(series):
        ax.bar(
            x + (i - 1) * w,
            [r[col] for r in rows],
            width=w,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            label=lab,
            zorder=3,
        )
    ax.axhline(0, color="black", linewidth=0.6, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([r[0] for r in rows])
    ax.set_ylabel(f"Mean bivariate Moran's I\n({'/'.join(str(v) for v in sorted(n))} specimens)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    place_legend(ax, position="bottom", frameon=False, ncol=3)


NCL_TYPE_LABELS = {
    "B_Cells": "B cells",
    "CD4_T_Cells": "CD4\nT cells",
    "CD8_T_Cells": "CD8\nT cells",
    "Cancer_Cells": "Cancer\ncells",
    "Endothelial": "Endothelial",
    "Fibroblasts": "Fibroblasts",
    "Macrophages": "Macrophages",
}


def mdk_ncl_rows(tissue: dict) -> tuple:
    """(rows, headline specimen, min_active, n_specimens) from tissue.json::mdk_ncl_celltype_map (pure; raises on any gap).

    Each row: label, sender / receiver median share across specimens, the headline specimen's share, and the number of
    specimens with at least ``min_active`` active compartments of that type (sender, receiver). Sorted by sender median.
    """
    m = tissue["mdk_ncl_celltype_map"]
    if m.get("stage") != "B18" or m.get("descriptive_only") is not True:
        raise ValueError("mdk_ncl_celltype_map provenance is not B18 / descriptive")
    rows = []
    for e in m["cohort_mdk_ncl"]:
        if e["cell_type"] not in NCL_TYPE_LABELS:
            raise KeyError(f"no plain-English label for cell type {e['cell_type']!r}")
        s, r = e["s"], e["r"]
        rows.append(
            {
                "label": NCL_TYPE_LABELS[e["cell_type"]],
                "s_med": float(s["median_share_across_specimens"]),
                "r_med": float(r["median_share_across_specimens"]),
                "s_p4": float(s["P4-S2_share"]),
                "r_p4": float(r["P4-S2_share"]),
                "s_n": int(s["n_specimens_ge30_active"]),
                "r_n": int(r["n_specimens_ge30_active"]),
                "n_with": int(s["n_specimens_with_type"]),
            }
        )
    if not rows:
        raise ValueError("no cell types in mdk_ncl_celltype_map")
    rows.sort(key=lambda r: -r["s_med"])
    return rows, m["headline_specimen"], int(m["min_active"]), int(m["n_specimens"])


def draw_panel_mdk_ncl_celltypes(ax: plt.Axes) -> None:
    """Fig 5B: MDK-NCL sender and receiver score share by cell type (median across specimens; P4-S2 marked; B18, descriptive)."""
    rows, head, min_active, n_spec = mdk_ncl_rows(json.loads(TISSUE_JSON.read_text()))
    n = len(rows)
    y = np.arange(n)[::-1]
    h = 0.36
    c_s, c_r = "#E07B39", "#5B8DBE"
    ax.barh(
        y + h / 2,
        [r["s_med"] for r in rows],
        height=h,
        color=c_s,
        edgecolor="black",
        linewidth=0.5,
        label="Sender",
        zorder=3,
    )
    ax.barh(
        y - h / 2,
        [r["r_med"] for r in rows],
        height=h,
        color=c_r,
        edgecolor="black",
        linewidth=0.5,
        label="Receiver",
        zorder=3,
    )
    ax.scatter(
        [r["s_p4"] for r in rows],
        y + h / 2,
        marker="D",
        s=22,
        color="white",
        edgecolor="black",
        linewidth=0.8,
        zorder=5,
    )
    ax.scatter(
        [r["r_p4"] for r in rows],
        y - h / 2,
        marker="D",
        s=22,
        color="white",
        edgecolor="black",
        linewidth=0.8,
        zorder=5,
        label="P4-S2",
    )
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=9)
    ax.set_xlim(0, 0.95)
    ax.set_xticks([0, 0.2, 0.4, 0.6])
    ax.set_xlabel("Share of MDK-NCL score")
    # specimens with enough active compartments of each type (small-type shares need caution)
    ax.text(
        0.9,
        1.02,
        f"Specimens with >={min_active}\nactive (of {n_spec}):\nsender / receiver",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )
    for yi, r in zip(y, rows):
        ax.text(
            0.9, yi, f"{r['s_n']} / {r['r_n']}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.6, zorder=0)
    place_legend(ax, position="bottom", frameon=False, ncol=3)


SHORT = (
    lambda sid: sid.replace("HCC22-088-", "")
    .replace("_1i_rep", "")
    .replace("_F_rep", "")
    .replace("_D", "")
    .replace("_A", "")
)


def moran_control_rows(cp: dict) -> list:
    """Per-specimen rows from the stamped B17 file (pure): sample, real_I, and p5 / median / p95 of both matched controls.

    Every specimen with status ok is returned (none dropped); the stage must be B17 and descriptive_only.
    """
    if cp.get("stage") != "B17" or cp.get("descriptive_only") is not True or cp.get("k") != 6:
        raise ValueError("control_pairs provenance is not B17 / descriptive / k=6")
    rows = []
    for sid, e in cp["per_specimen"].items():
        if e.get("status") != "ok":
            continue
        a, b = e["control_a_random_mdk"], e["control_b_random_secretory"]
        rows.append(
            {
                "sample": sid,
                "real": float(e["real_I"]),
                "a": (a["p5"], a["median"], a["p95"], bool(a["real_gt_p95"])),
                "b": (b["p5"], b["median"], b["p95"], bool(b["real_gt_p95"])),
            }
        )
    if not rows:
        raise ValueError("no estimable specimens in control_pairs")
    return rows


def draw_panel_moran_controls(ax: plt.Axes) -> None:
    """Fig 5C: real MDK-secretory Moran's I per specimen (dot) against the 5th-95th percentile range of its 200 matched
    random-gene controls, for both control families. Every estimable specimen is shown; the counts above p95 are printed.
    """
    cp = json.loads((v2_runs.results_dir("B17") / "control_pairs_results.json").read_text())
    rows = moran_control_rows(cp)
    n = len(rows)
    y = np.arange(n)[::-1]
    off = 0.17
    for i, r in enumerate(rows):
        for sign, key, color in ((+1, "a", "#9e9e9e"), (-1, "b", "#c9c9c9")):
            lo, med, hi, _ = r[key]
            ax.plot([lo, hi], [y[i] + sign * off] * 2, color=color, linewidth=3.2, solid_capstyle="butt", zorder=2)
            ax.plot([med], [y[i] + sign * off], marker="|", color="black", markersize=5, zorder=3)
        ax.scatter(
            [r["real"]], [y[i]], s=34, color=ROBUSTNESS_COLORS["real"], edgecolor="black", linewidth=0.5, zorder=4
        )
    ax.axvline(0, color="black", linewidth=0.5, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT(r["sample"]) for r in rows])
    ax.set_xlabel("Moran's I (bars: 5th-95th pct)")
    na, nb = sum(r["a"][3] for r in rows), sum(r["b"][3] for r in rows)
    mr, ma, mb = (
        float(np.mean(v)) for v in ([r["real"] for r in rows], [r["a"][1] for r in rows], [r["b"][1] for r in rows])
    )
    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=ROBUSTNESS_COLORS["real"],
            markeredgecolor="black",
            markersize=6,
            label=f"Real (mean {mr:.3f})",
        ),
        Line2D([0], [0], color="#9e9e9e", linewidth=3.2, label=f"Random MDK gene: above {na}/{n}"),
        Line2D([0], [0], color="#c9c9c9", linewidth=3.2, label=f"Random secretory genes: above {nb}/{n}"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def program_recovery_rows(pr: dict) -> tuple:
    """Return (n_specimens, [(method, full, partial, absent), ...]) from the pr2 registry.

    Raises on a missing key, a non-ok (not-estimable/failed) specimen, or counts
    that do not sum to n_specimens. Order: CITEgeist, cell2location, tangram.
    """
    n = pr["n_specimens"]
    rows = []
    for m in ("CITEgeist", "cell2location", "tangram"):
        t = pr["totals"][m]
        if t["n_not_estimable"] or t["n_failed"] or t["n_specimens_ok"] != n:
            raise ValueError(f"program recovery: {m} has non-ok specimens")
        f, p, a = t["full"], t["partial"], t["absent"]
        if f + p + a != n:
            raise ValueError(f"program recovery: {m} full+partial+absent != n_specimens={n}")
        rows.append((m, f, p, a))
    return n, rows


def draw_panel_program_recovery(ax: plt.Axes) -> None:
    """MDK NMF program recovery by method (stacked specimen counts, of N).

    Rendered as the single panel of Supplementary Figure program_recovery_by_method
    (generate_supp_program_recovery_by_method.py); no longer a main-figure panel.
    """
    pr = json.loads(PROGRAM_RECOVERY_JSON.read_text())
    if pr.get("_code_dirty") is not False:
        raise ValueError(f"{PROGRAM_RECOVERY_JSON}: _code_dirty is {pr.get('_code_dirty')!r}")
    n, rows = program_recovery_rows(pr)
    labels = [METHOD_DISPLAY_CMP.get(r[0], r[0]) for r in rows]
    x = np.arange(len(labels))
    bottom = np.zeros(len(rows))
    for i, (cl, key) in enumerate([("Full", "full"), ("Partial", "partial"), ("Absent", "absent")]):
        vals = np.array([r[i + 1] for r in rows], dtype=float)
        ax.bar(
            x,
            vals,
            bottom=bottom,
            width=0.65,
            color=PROGRAM_COLORS_CMP[key],
            edgecolor="black",
            linewidth=0.6,
            label=cl,
        )
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel(f"Specimens (of {n})")
    ax.set_ylim(0, n)
    ax.set_yticks(range(0, n + 1, 2))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    place_legend(ax, position="right", title="Recovery", frameon=False)


def panel_specs_for(figure):
    """Return the [(output_label, draw_fn), ...] for a midkine sub-figure.

    Pure — no rendering, no data access — so it is unit-testable on a login node.
    Draw-fn names are historical and do not track output letters.
    """
    if figure == "spatial":
        return [
            ("A", draw_panel_A),  # MDK / secretory spatial maps (v2 B5 compartments)
            ("B", draw_panel_mdk_ncl_celltypes),  # MDK-NCL sender / receiver share by cell type (B18, descriptive)
            ("C", draw_panel_moran_controls),  # real I per specimen vs matched random-gene controls (B17)
            # v42: hurdle (draw_panel_B), the 3-method Moran's comparison and the across-forms robustness panel moved to the supplement
        ]
    if figure == "mechanism":
        return [
            ("A", draw_panel_H),  # ChIP-seq
            ("B", draw_panel_G),  # secretory fold-change (RNA)
            # ("C", ...) IF montage — owned by generate_if_panel.py
            ("D", draw_panel_E),  # ELISA
        ]
    raise ValueError(f"unknown figure {figure!r}")


def save_panels(panel_sizes: dict | None = None, figure: str = "spatial", panels_dir=None) -> None:
    """Save panels for one midkine sub-figure as JPG + SVG.

    The draw-function names are historical and do NOT track the output letter —
    the (label, fn) tuples returned by panel_specs_for() are the source of truth
    for content↔letter mapping.

    panels_dir overrides the module PANELS_DIR so the spatial (Fig 5) and
    mechanism (Fig 6) sub-figures don't clobber each other's panel_*.svg.

    Panel A uses imshow (spatial). The IF montage panel is produced by
    generate_if_panel.py directly from raw z-stack TIFFs and is intentionally
    NOT regenerated here.
    """
    out_dir = panels_dir if panels_dir is not None else PANELS_DIR
    panel_specs = panel_specs_for(figure)
    # Convert list format [w_mm, h_mm] → dict format {w_mm: ..., h_mm: ...} for PanelContext
    sizes_dict: dict | None = None
    if panel_sizes:
        sizes_dict = {
            label: {"w_mm": v[0], "h_mm": v[1]} if isinstance(v, (list, tuple)) else v
            for label, v in panel_sizes.items()
        }
    for label, draw_fn in panel_specs:
        # ELISA (mechanism label D) renders without constrained_layout (creates sub-axes).
        kwargs = {"constrained_layout": False} if (figure == "mechanism" and label == "D") else {}
        with PanelContext(label, sizes_dict, out_dir, **kwargs) as (fig, ax):
            draw_fn(ax)


def main() -> None:
    parser = panel_argparser("Generate midkine (MDK) figure panels.")
    parser.add_argument("--figure", choices=["spatial", "mechanism"], default="spatial")
    parser.add_argument(
        "--panels-dir", default=None, help="Override output dir (keeps spatial/mechanism panels separate)."
    )
    args = parser.parse_args()
    panel_sizes = load_panel_sizes(args)

    from pathlib import Path

    pdir = Path(args.panels_dir) if args.panels_dir else None
    print(f"Saving {args.figure} panels...")
    save_panels(panel_sizes=panel_sizes, figure=args.figure, panels_dir=pdir)
    print("Done.")


if __name__ == "__main__":
    main()
