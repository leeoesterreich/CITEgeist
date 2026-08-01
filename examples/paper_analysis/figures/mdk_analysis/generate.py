#!/usr/bin/env python3
"""
Midkine (MDK) Figure — CITEgeist Patterns v21.

Panels (8 flat, A-H):
  A — Spatial grid (MDK map + proportion maps side-by-side)
  B — COMMOT MDK-NCL interaction boxplot
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

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from matplotlib.image import imread

SCRIPT_DIR = Path(__file__).resolve().parent
# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
sys.path.insert(0, str(SCRIPT_DIR.parent))
from _shared.generate_utils import PanelContext, load_panel_sizes, panel_argparser, place_annotation, place_legend
from _shared.style import apply_style

apply_style()
plt.rcParams.update({"svg.fonttype": "none"})

PANELS_DIR = SCRIPT_DIR / "panels"
PANELS_DIR.mkdir(parents=True, exist_ok=True)

MDK_FIGS = PROJECT / "mdk_analysis" / "outputs" / "figures"
MDK_IF = PROJECT / "mdk_analysis" / "outputs" / "if_analysis"

PANEL_PATHS: dict[str, Path] = {
    "A_grid": MDK_IF / "panel_A_grid.png",
    "F": PANELS_DIR / "panel_F.jpg",
}

# CSV data sources for vector panels
HURDLE_CSV = PROJECT / "mdk_analysis/outputs/commot_p4s2/hurdle_results_summary.csv"
HURDLE_JSON = PROJECT / "mdk_analysis/outputs/commot_p4s2/hurdle_results.json"
MORANS_CSV = PROJECT / "mdk_analysis/outputs/if_analysis/bivariate_morans_results.csv"
ELISA_CSV = PROJECT / "mdk_analysis/outputs/if_analysis/elisa_data.csv"
CHIPSEQ_CSV = PROJECT / "mdk_analysis/_archive/mdk_saturation_pipeline/outputs/tables/er_chipseq_all_conditions.csv"
TRAJECTORY_CSV = PROJECT / "mdk_analysis/outputs/spatial/trajectory_deltas.csv"

# Competitor-GEX comparison (promoted from supplementary biology_validation → Fig 5 C/D)
COMPETITOR_DIR = PROJECT / "output" / "competitor_gex" / "comparison"
MORANS_COMPARE_CSV = COMPETITOR_DIR / "morans_comparison.csv"
FAILURE_CSV = COMPETITOR_DIR / "failure_mode_summary.csv"
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
    """Spatial deconvolution overview: 2x4 grid of Cancer Luminal proportion
    (row 1) and MDK expression (row 2) across 4 patient samples.
    """
    _show_img(ax, PANEL_PATHS["A_grid"])


def draw_panel_B(ax: plt.Axes) -> None:
    """COMMOT hurdle model (P4, single patient): two significant axes.

    Left  = overall signaling participation (pooled across pathways),
            OR=4.01, p=8.7e-11 (Fisher's exact).
    Right = MDK sender magnitude among active senders (NCL 2.08x, SDC2 2.01x;
            Mann-Whitney FDR<0.05).
    The n.s. MDK-NCL *participation* bar was removed: it made the panel read
    as null even though the finding (a ~2x magnitude increase) is significant.
    MDK-NCL participation (39% vs 36%, p=0.47) is reported in the caption.
    """
    import json

    # Load data
    hurdle_df = pd.read_csv(HURDLE_CSV)
    with open(HURDLE_JSON, encoding="utf-8") as f:
        hurdle_json = json.load(f)

    # --- Set up two subplots ---
    fig = ax.get_figure()
    pos = ax.get_position()
    ax.set_visible(False)

    # Reserve a bottom strip for the shared legend so it does not overlap bars.
    ax_left = fig.add_axes([pos.x0, pos.y0 + pos.height * 0.18, pos.width * 0.45, pos.height * 0.70])
    ax_right = fig.add_axes(
        [pos.x0 + pos.width * 0.55, pos.y0 + pos.height * 0.18, pos.width * 0.45, pos.height * 0.70]
    )

    color_mut = "#E07B39"  # D538G-high (orange)
    color_wt = "#5B8DBE"  # D538G-low (blue)

    # --- Left: overall (pooled) signaling participation -----------------
    # Only the pooled "total" participation is shown: it is the significant
    # participation result. MDK is significant in magnitude, not engagement
    # frequency, so its effect is shown on the right instead.
    total_mut = hurdle_json["total_sender_participation"]["mut_active_pct"]
    total_wt = hurdle_json["total_sender_participation"]["wt_active_pct"]
    total_p = hurdle_json["total_sender_participation"]["fisher_p"]
    total_or = hurdle_json["total_sender_participation"]["fisher_OR"]

    bar_w = 0.4
    ax_left.bar([-bar_w / 2 - 0.02], [total_mut], width=bar_w, color=color_mut, label="D538G-high")
    ax_left.bar([bar_w / 2 + 0.02], [total_wt], width=bar_w, color=color_wt, label="D538G-low")

    ax_left.set_xticks([0])
    ax_left.set_xticklabels(["All pathways"])
    ax_left.set_xlim(-0.7, 0.7)
    ax_left.set_ylabel("Sender participation (%)", labelpad=8)
    ax_left.set_ylim(0, 118)
    ax_left.set_title("Overall signaling")
    # Legend placed below the left sub-axes so the bars stay visible.
    ax_left.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=1,
        frameon=False,
    )
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)

    # Significance annotation for the pooled participation comparison.
    def _annot_sig(ax_t, x_pos, y_top, p_val):
        stars = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
        ax_t.text(x_pos, y_top + 2, stars, ha="center", va="bottom", fontweight="bold" if stars != "ns" else "normal")

    y_top = max(total_mut, total_wt)
    _annot_sig(ax_left, 0, y_top, total_p)
    ax_left.text(0, y_top + 11, f"OR = {total_or:.1f}", ha="center", va="bottom")

    # --- Right: MDK sender magnitude among active senders ----------------
    mdk_rows = hurdle_df[hurdle_df["pathway"].str.startswith("s-MDK-")].copy()
    # Exclude ITGA6_ITGB1 (no active senders)
    mdk_rows = mdk_rows[~mdk_rows["pathway"].str.contains("ITGA6_ITGB1")]
    mdk_rows = mdk_rows.dropna(subset=["sender_fc_nonzero"])

    receptors = [p.replace("s-MDK-", "") for p in mdk_rows["pathway"]]
    fc_vals = mdk_rows["sender_fc_nonzero"].values
    fdr_vals = mdk_rows["mw_fdr_nonzero"].values

    colors = ["#d62728" if fdr < 0.05 else "#999999" for fdr in fdr_vals]
    y_pos = np.arange(len(receptors))

    ax_right.barh(y_pos, fc_vals, color=colors, height=0.6, alpha=0.85)
    ax_right.axvline(1.0, color="gray", linestyle="--", linewidth=1.2)
    ax_right.set_yticks(y_pos)
    ax_right.set_yticklabels(receptors)
    ax_right.set_xlabel("Fold change (D538G/WT)", labelpad=8)
    ax_right.set_title("MDK sender strength")
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)

    # Annotate fold-change with FDR significance stars (Mann-Whitney on
    # nonzero senders). Significant bars are red, n.s. grey.
    def _stars(fdr):
        return "***" if fdr < 0.001 else ("**" if fdr < 0.01 else ("*" if fdr < 0.05 else ""))

    for i, (fc, fdr) in enumerate(zip(fc_vals, fdr_vals)):
        star = _stars(fdr)
        label = f"{fc:.2f}x {star}".rstrip() if star else f"{fc:.2f}x (ns)"
        ax_right.text(fc + 0.05, i, label, va="center", ha="left")
    if len(fc_vals):
        ax_right.set_xlim(0, max(fc_vals) * 1.40)


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

    df = pd.read_csv(ELISA_CSV)
    group_order = df["group"].drop_duplicates().tolist()
    n = len(group_order)
    colors_list = ["#66C2A5", "#FC8D62", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#8DA0CB", "#E78AC3"]

    means, sems = [], []
    group_vals: list[np.ndarray] = []
    for i, grp in enumerate(group_order):
        vals = df.loc[df["group"] == grp, "mdk_pg_ml"].values
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
    ax.set_ylabel("MDK (pg/mL)", labelpad=8)

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
    stagger_extra = tick * 3.0
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

    _ = ax.imshow(display_norm, aspect="equal", cmap="RdYlBu_r", vmin=-2, vmax=2, rasterized=True)
    ax.set_xticks(range(len(cond_cols)))
    ax.set_xticklabels([c.replace("_", " ") for c in cond_cols], rotation=90, ha="right", fontsize=7)
    ax.set_yticks(range(len(display_genes)))
    ax.set_yticklabels(display_genes, fontsize=8)
    ax.tick_params(axis="y", pad=2)

    for i_disp, gi in enumerate(orig_row):
        if gi is None:
            continue
        for j in range(len(cond_cols)):
            val = int(data[gi, j])
            tc = "white" if abs(normalized[gi, j]) > 1 else "black"
            ax.text(j, i_disp, str(val), ha="center", va="center", color=tc)

    # Italic group annotations placed in the left margin, vertically centered
    # within each group (no line — the blank row IS the separator).
    # clip_on=False allows the rotated text to render outside the axes boundary
    # without being clipped by the panel edge (fixes 4G y-axis clipping).
    if 0 < n_bound < len(genes):
        bound_mid = (n_bound - 1) / 2.0
        unbound_mid = n_bound + 1 + (len(genes) - n_bound - 1) / 2.0
        ax.text(
            -3.8,
            bound_mid,
            "ER binding",
            ha="right",
            va="center",
            fontsize=7.4,
            fontstyle="italic",
            rotation=90,
            clip_on=False,
        )
        ax.text(
            -3.8,
            unbound_mid,
            "No binding",
            ha="right",
            va="center",
            fontsize=7.4,
            fontstyle="italic",
            rotation=90,
            clip_on=False,
        )


def draw_panel_I(ax: plt.Axes) -> None:
    """Patient MDK-secretory trajectory: biopsy to surgery."""
    df = pd.read_csv(TRAJECTORY_CSV)
    patient_meta = {
        "P1": {"color": "#d62728", "marker": "s", "is_prog": True},
        "P2": {"color": "#1f77b4", "marker": "o", "is_prog": False},
        "P3": {"color": "#2ca02c", "marker": "o", "is_prog": False},
        "P5": {"color": "#9467bd", "marker": "o", "is_prog": False},
        "P6": {"color": "#17becf", "marker": "o", "is_prog": False},
    }

    endpoints: list[tuple[str, float, dict]] = []
    for _, row in df.iterrows():
        p = row["patient"]
        if p not in patient_meta:
            continue
        m = patient_meta[p]
        prog = m["is_prog"]
        lw = 2.5 if prog else 1.5
        ms = 9 if prog else 7
        alpha = 1.0 if prog else 0.6
        z = 10 if prog else 5

        ax.plot([0, 1], [row["biopsy_I"], row["surgery_I"]], color=m["color"], linewidth=lw, alpha=alpha, zorder=z)
        ax.scatter(
            0,
            row["biopsy_I"],
            s=ms**2,
            color=m["color"],
            marker=m["marker"],
            edgecolors=m["color"],
            facecolors="white",
            linewidths=1.5,
            alpha=alpha,
            zorder=z + 1,
        )
        ax.scatter(1, row["surgery_I"], s=ms**2, color=m["color"], marker=m["marker"], alpha=alpha, zorder=z + 1)
        endpoints.append((p, float(row["surgery_I"]), m))

    endpoints.sort(key=lambda e: e[1])
    all_vals = list(df["biopsy_I"]) + list(df["surgery_I"])
    y_range = max(all_vals) - min(all_vals)
    # 0.12 (was 0.06): the smaller gap left near-coincident endpoints (P2/P5)
    # visually tight without triggering a nudge; 0.12×range forces a legible gap
    # between adjacent endpoint labels (with leader lines back to true points).
    min_sep = max(y_range * 0.12, 0.025)
    placed_ys: list[float] = []
    for p, y_orig, meta in endpoints:
        y = y_orig
        if placed_ys and y - placed_ys[-1] < min_sep:
            y = placed_ys[-1] + min_sep
        placed_ys.append(y)
        # Draw a thin leader line when the label has been nudged away from
        # the true endpoint so readers can trace which point the label belongs to.
        if abs(y - y_orig) > 1e-9:
            ax.annotate(
                p,
                xy=(1.05, y_orig),
                xytext=(1.05, y),
                va="center",
                color=meta["color"],
                fontweight="bold" if meta["is_prog"] else "normal",
                alpha=1.0 if meta["is_prog"] else 0.6,
                fontsize=7,
                arrowprops=dict(arrowstyle="-", lw=0.4, color="0.6"),
            )
        else:
            ax.text(
                1.05,
                y,
                p,
                va="center",
                color=meta["color"],
                fontweight="bold" if meta["is_prog"] else "normal",
                alpha=1.0 if meta["is_prog"] else 0.6,
                fontsize=7,
            )

    mean_delta = float(df["delta_I"].mean())
    # Paired one-sided Wilcoxon signed-rank: H1 is biopsy < surgery (delta > 0).
    # n=5 all-positive deltas saturate at p = 1/2^5 = 0.03125 (formally significant).
    from scipy.stats import wilcoxon

    deltas = df["delta_I"].dropna().values
    if len(deltas) >= 5:
        _, wilcoxon_p = wilcoxon(deltas, alternative="greater", zero_method="wilcox")
        # Single-line annotation: the prior two-line form was misread by the
        # composite validator as two overlapping x-tick labels. Keep it compact
        # (full methodological detail lives in the caption) so it does not crowd
        # the data area or overflow the right edge.
        annot = f"Δ̄ = +{mean_delta:.3f}; Wilcoxon p = {wilcoxon_p:.3f}"
    else:
        annot = f"Δ̄ = +{mean_delta:.3f}"

    max_i = max(all_vals) if all_vals else 1.0
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Biopsy", "Surgery"])
    ax.set_ylabel("Moran's I (MDK–secretory)", labelpad=8)
    # Extra right margin hosts de-collided endpoint labels; extra top headroom
    # gives the annotation room above the data.
    ax.set_xlim(-0.18, 1.55)
    ax.set_ylim(-0.02, max_i * 1.55)

    # Annotation in upper-LEFT corner so the long "Wilcoxon 1-sided p = ..."
    # string never extends past the panel's right edge (T17 fix).
    place_annotation(
        ax,
        0.02,
        0.97,
        annot,
        ha="left",
        va="top",
        fontsize=7,
    )

    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], marker="s", color="#d62728", markersize=8, label="Progressor", linewidth=2),
        Line2D([0], [0], marker="o", color="gray", markersize=7, label="Responder", linewidth=1.5, alpha=0.6),
    ]
    # External-bottom legend matches Fig3 C/E/F treatment so the
    # Progressor/Responder key never overlaps the trajectory lines.
    ax.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.32),
        ncol=2,
        frameon=False,
        fontsize=8,
    )


def _sig_label_cmp(p: float) -> str:
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return "n.s."


def _citegeist_mdk_morans():
    """CITEgeist MDK bivariate Moran's I: (mean, sem, per-sample valid df)."""
    df = pd.read_csv(MORANS_CSV)
    valid = df.dropna(subset=["morans_I"])
    vals = valid["morans_I"].values
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(len(vals))), valid


def _paired_wilcoxon_vs_citegeist(competitor: str, morans: pd.DataFrame, cg_valid: pd.DataFrame):
    """One-sided paired Wilcoxon (CITEgeist > competitor) on matched specimens."""
    cg = cg_valid[["sample", "morans_I"]].rename(columns={"morans_I": "I_cg"})
    comp = morans[(morans["gene"] == "MDK") & (morans["method"] == competitor)][["sample", "I"]].rename(
        columns={"I": "I_comp"}
    )
    merged = cg.merge(comp, on="sample", how="inner").dropna(subset=["I_cg", "I_comp"])
    n = len(merged)
    if n < 5:
        return np.nan, n
    try:
        return wilcoxon(merged["I_cg"].values, merged["I_comp"].values, alternative="greater").pvalue, n
    except ValueError:
        return np.nan, n


def draw_panel_morans_methods(ax: plt.Axes) -> None:
    """Panel C — MDK bivariate Moran's I by method (CITEgeist 11 / C2L 12 / Tangram 12)."""
    morans = pd.read_csv(MORANS_COMPARE_CSV)
    cg_I, cg_sem, cg_valid = _citegeist_mdk_morans()
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
    ax.set_ylabel("MDK bivariate Moran's I")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, mv in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.008, f"{mv:.3f}", ha="center", va="bottom")

    y_max = max(means) * 1.20
    step = max(means) * 0.14
    off = 0.0
    for comp in ["cell2location", "tangram"]:
        disp = METHOD_DISPLAY_CMP.get(comp, comp)
        if disp not in labels:
            continue
        ci = labels.index(disp)
        p, n = _paired_wilcoxon_vs_citegeist(comp, morans, cg_valid)
        sig = f"n.s. (n={n})" if np.isnan(p) else f"{_sig_label_cmp(p)} (p={p:.4f}, n={n})"
        by = y_max + off
        ax.plot([x[0], x[0], x[ci], x[ci]], [by - 0.005, by, by, by - 0.005], color="black", linewidth=0.8)
        ax.text((x[0] + x[ci]) / 2, by + 0.005, sig, ha="center", va="bottom")
        off += step
    ax.set_ylim(0, y_max + off + step * 0.5)


def draw_panel_program_recovery(ax: plt.Axes) -> None:
    """Panel D — MDK NMF program recovery by method (out of 12)."""
    failure = pd.read_csv(FAILURE_CSV)
    prog = failure[failure["analysis"] == "programs_mdk"].copy().sort_values("method")
    prog = prog[prog["method"].isin(["cell2location", "tangram"])]
    labels = [METHOD_DISPLAY_CMP.get(m, m) for m in prog["method"]]
    x = np.arange(len(labels))
    bar_w = 0.65 / 3
    cats = [
        ("Full", prog["full"].tolist(), PROGRAM_COLORS_CMP["full"]),
        ("Partial", prog["partial"].tolist(), PROGRAM_COLORS_CMP["partial"]),
        ("Absent", prog["absent"].tolist(), PROGRAM_COLORS_CMP["absent"]),
    ]
    for (cl, vals, c), boff in zip(cats, [-bar_w, 0, bar_w]):
        ax.bar(x + boff, vals, width=bar_w * 0.9, color=c, edgecolor="black", linewidth=0.6, label=cl)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("MDK programs (of 12)")
    ax.set_yticks(range(0, 13, 2))
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
            ("A", draw_panel_A),  # spatial grid
            ("B", draw_panel_B),  # COMMOT hurdle
            ("C", draw_panel_morans_methods),  # 3-method Moran's (promoted S9-A; forest detail → Supp S6)
            ("D", draw_panel_program_recovery),  # NMF program recovery (promoted S9-B)
            ("E", draw_panel_I),  # trajectory (was D)
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
