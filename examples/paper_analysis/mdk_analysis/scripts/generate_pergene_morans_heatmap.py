#!/usr/bin/env python
"""
Generate Figure 4D: Per-gene bivariate Moran's I heatmap (MDK vs secretory genes).

Reads pergene_bivariate_morans.csv, produces a gene x sample heatmap with
FDR significance annotations and a side panel comparing secretory vs control.
"""

import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

PROJECT_DIR = "/path/to/CITEgeist_analysis"
INPUT_CSV = os.path.join(PROJECT_DIR, "mdk_analysis/outputs/if_analysis/pergene_bivariate_morans.csv")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "figures/export/Figure4")

SECRETORY_GENES = ["HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4", "PDIA6", "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1"]

# Short sample labels
SAMPLE_SHORT = {
    "HCC22-088-P1-S1": "P1-S1\n(Bx, Prog)",
    "HCC22-088-P1-S2": "P1-S2\n(Sx, Prog)",
    "HCC22-088-P2-S1": "P2-S1\n(Bx, Resp)",
    "HCC22-088-P2-S2": "P2-S2\n(Sx, Resp)",
    "HCC22-088-P3-S1_A": "P3-S1\n(Bx, Resp)",
    "HCC22-088-P3-S2": "P3-S2\n(Sx, Resp)",
    "HCC22-088-P4-S1": "P4-S1\n(Bx, Prog)",
    "HCC22-088-P4-S2_1i_rep": "P4-S2\n(Sx, Prog)",
    "HCC22-088-P5-S1": "P5-S1\n(Bx, Resp)",
    "HCC22-088-P5-S2_F_rep": "P5-S2\n(Sx, Resp)",
    "HCC22-088-P6-S1": "P6-S1\n(Bx, Resp)",
    "HCC22-088-P6-S2_D": "P6-S2\n(Sx, Resp)",
}

SAMPLE_ORDER = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)
    sec_df = df[df["gene_type"] == "secretory"].copy()
    ctrl_df = df[df["gene_type"] == "control"].copy()

    # Pivot to gene x sample matrix
    pivot_I = sec_df.pivot(index="gene", columns="sample", values="morans_I")
    pivot_q = sec_df.pivot(index="gene", columns="sample", values="fdr_q")

    # Reorder
    pivot_I = pivot_I.reindex(index=SECRETORY_GENES, columns=SAMPLE_ORDER)
    pivot_q = pivot_q.reindex(index=SECRETORY_GENES, columns=SAMPLE_ORDER)

    # --- Figure: 2-panel layout ---
    fig = plt.figure(figsize=(10, 5.5))

    # Left: heatmap (wider)
    ax_heat = fig.add_axes([0.08, 0.18, 0.58, 0.72])
    # Right: secretory vs control comparison
    ax_comp = fig.add_axes([0.74, 0.18, 0.22, 0.72])

    # === Heatmap ===
    cmap = plt.cm.YlOrRd
    vmin, vmax = 0, 0.65
    im = ax_heat.imshow(pivot_I.values, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax, interpolation="nearest")

    # FDR significance markers
    for i in range(pivot_I.shape[0]):
        for j in range(pivot_I.shape[1]):
            q = pivot_q.iloc[i, j]
            val = pivot_I.iloc[i, j]
            if pd.isna(val):
                ax_heat.text(j, i, "—", ha="center", va="center", fontsize=7, color="gray")
            elif pd.notna(q) and q < 0.001:
                ax_heat.text(
                    j,
                    i,
                    f"{val:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    fontweight="bold",
                    color="white" if val > 0.35 else "black",
                )
            elif pd.notna(q) and q < 0.05:
                ax_heat.text(
                    j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="white" if val > 0.35 else "black"
                )
            else:
                # Not significant
                ax_heat.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="gray", fontstyle="italic")

    # Axis labels
    ax_heat.set_xticks(range(len(SAMPLE_ORDER)))
    ax_heat.set_xticklabels([SAMPLE_SHORT[s] for s in SAMPLE_ORDER], fontsize=7, rotation=0, ha="center")
    ax_heat.set_yticks(range(len(SECRETORY_GENES)))
    ax_heat.set_yticklabels(SECRETORY_GENES, fontsize=9)

    # Highlight P4-S1 column (weakest)
    p4s1_idx = SAMPLE_ORDER.index("HCC22-088-P4-S1")
    ax_heat.add_patch(
        Rectangle(
            (p4s1_idx - 0.5, -0.5),
            1,
            len(SECRETORY_GENES),
            linewidth=2,
            edgecolor="red",
            facecolor="none",
            linestyle="--",
        )
    )

    # Colorbar
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.02, pad=0.02)
    cbar.set_label("Bivariate Moran's I", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax_heat.set_title(
        "MDK co-localization with secretory genes\n(bold = FDR < 0.001, italic gray = ns)", fontsize=10, pad=8
    )

    # === Comparison panel: secretory vs control ===
    # Per-sample mean I for secretory and control
    sec_means = sec_df.groupby("sample")["morans_I"].mean().reindex(SAMPLE_ORDER)
    ctrl_means = ctrl_df.groupby("sample")["morans_I"].mean().reindex(SAMPLE_ORDER)

    x = np.arange(len(SAMPLE_ORDER))
    width = 0.35
    bars1 = ax_comp.barh(x - width / 2, sec_means.values, width, color="#d62728", alpha=0.8, label="Secretory")
    bars2 = ax_comp.barh(x + width / 2, ctrl_means.values, width, color="#7f7f7f", alpha=0.6, label="Control")

    ax_comp.set_yticks(x)
    ax_comp.set_yticklabels([])  # shared with heatmap
    ax_comp.set_xlabel("Mean Moran's I", fontsize=9)
    ax_comp.set_title("Secretory vs\ncontrol genes", fontsize=10, pad=8)
    ax_comp.legend(fontsize=7, loc="lower right")
    ax_comp.set_xlim(0, 0.55)
    ax_comp.tick_params(labelsize=8)
    ax_comp.invert_yaxis()

    # Add overall means as vertical lines
    ax_comp.axvline(sec_means.mean(), color="#d62728", linestyle="--", linewidth=1, alpha=0.7)
    ax_comp.axvline(ctrl_means.mean(), color="#7f7f7f", linestyle="--", linewidth=1, alpha=0.7)

    # Summary text
    n_sig = (sec_df["fdr_q"] < 0.05).sum()
    n_total = sec_df["fdr_q"].notna().sum()
    fig.text(
        0.5,
        0.02,
        f"Secretory: {n_sig}/{n_total} pairs FDR < 0.05, mean I = {sec_df['morans_I'].mean():.3f}  |  "
        f"Control: mean I = {ctrl_df['morans_I'].mean():.3f}",
        ha="center",
        fontsize=9,
        style="italic",
    )

    # Save
    for ext in ["png", "svg"]:
        out_path = os.path.join(OUTPUT_DIR, f"Fig4D_pergene_morans_heatmap.{ext}")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {out_path}")

    plt.close(fig)
    print("\nDone!")


if __name__ == "__main__":
    main()
