#!/usr/bin/env python3
"""
Figure 4: Midkine/ESR1 Case Study (Section 2.4)

CITEgeist identifies upregulated midkine signaling in ESR1-mutant breast cancer.
Patient HCC22-088-P4-S2_1i_rep (ESR1 D538G mutation, surgical specimen, Progressor).

Panel A: Basal cytokeratin signatures in cancer cell layer (sc.pl.spatial with H&E)
Panel B: Spatial distribution of D538G mutation signals (sc.pl.spatial with H&E)
Panel C: Split violin (D538G-high vs WT-like) + EstroGene gene heatmap
Panel D: Combined Hallmark + KEGG pathway dot plot (horizontally compressed)
Panel E: Summed MDK COMMOT sender signaling score (sc.pl.spatial with H&E)
Panel F: COMMOT signaling differences between D538G and WT (bar plot)
Panel G: ELISA validation (PLACEHOLDER)
Panel H: IF validation (PLACEHOLDER)
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path

from figure_style import apply_style, PALETTE, get_cell_type_color

apply_style()

# ============================================================================
# Paths
# ============================================================================
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_ROOT = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)
OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE = "HCC22-088-P4-S2_1i_rep"

# Module 3 unified outputs (preferred -- uses real cell type names)
MODULE3_DIR = PROJECT_ROOT / "output" / "module3_unified"
PROPS_CSV = MODULE3_DIR / f"{SAMPLE}_cell_prop_finetuned_results.csv"
GEX_PARQUET = MODULE3_DIR / f"{SAMPLE}_gene_expression_pass1.parquet"
META_JSON = MODULE3_DIR / f"{SAMPLE}_module3_meta.json"

# Fallback: top-level output files (use Profile_X names)
PROPS_CSV_FALLBACK = PROJECT_ROOT / "output" / f"{SAMPLE}_cell_prop_finetuned_results.csv"
GEX_PARQUET_FALLBACK = PROJECT_ROOT / "output" / f"{SAMPLE}_gene_expression_pass1.parquet"

# Spatial coordinates from Visium
SPATIAL_CSV = (
    DATA_ROOT / SAMPLE / "outs" / "spatial" / "tissue_positions.csv"
)

# Module 4 programs (for gene context)
MODULE4_PROGRAMS = PROJECT_ROOT / "output" / f"{SAMPLE}_module4_programs.csv"
MODULE4_JSON = (
    PROJECT_ROOT / "output" / "module4_unified" / f"{SAMPLE}_module4_discovery.json"
)

# ============================================================================
# Custom colormaps - using viridis for spatial single-variable plots
# ============================================================================
# Use viridis for all spatial proportion/score plots (perceptually uniform)
BASAL_CMAP = "viridis"
MUTATION_CMAP = "viridis"
MDK_CMAP = "viridis"

# ============================================================================
# Cytokeratin signature genes
# Primary: basal CKs from Li et al. 2022. Fallback: luminal CKs available in
# deconvolved GEX (ER+ samples lack basal CK expression at RNA level).
# ============================================================================
BASAL_CK_GENES = [
    "KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17",
    "KRT15", "KRT16", "TP63", "EGFR", "CDH3",
]
LUMINAL_CK_GENES = ["KRT7", "KRT8", "KRT19"]

# EstroGene D538G mutation signature genes (from EstroGene2.0)
ESTROGENE_D538G_UP = [
    "AGR2", "TFF1", "TFF3", "GREB1", "PGR",
    "PDZK1", "CA12", "STC2", "IGFBP4", "XBP1",
    "AREG", "CCND1", "MYC", "CXCL12", "MDK",
]

ESTROGENE_D538G_DOWN = [
    "NRIP1", "ESR1", "RARA", "FOXA1",
]

# MSigDB Hallmark pathways -- Enrichr combined scores from vignette_2
# D538G vs WT cancer spots, MSigDB_Hallmark_2020 gene set
HALLMARK_PATHWAYS = {
    "Estrogen Response Early": 538.37,
    "Apoptosis": 493.42,
    "TNF-alpha Signaling\nvia NF-kB": 466.91,
    "Epithelial Mesenchymal\nTransition": 401.97,
    "Estrogen Response Late": 290.37,
    "Myogenesis": 201.05,
    "Hypoxia": 131.59,
    "p53 Pathway": 103.58,
}

# KEGG pathways -- Enrichr combined scores from vignette_2
# D538G vs WT cancer spots, KEGG_2021_Human gene set
KEGG_PATHWAYS = {
    "Protein digestion\nand absorption": 97.15,
    "ECM-receptor\ninteraction": 84.76,
    "Focal adhesion": 59.03,
    "Apoptosis": 55.78,
    "Proteoglycans\nin cancer": 41.06,
}

# COMMOT signaling -- from vignette_2 D538G vs WT cancer cell analysis
# CellChat database, Mann-Whitney U test with BH FDR correction
COMMOT_LIGANDS = {
    "MDK-SDC4": {"d538g": 7.608, "wt": 2.221, "pval": 9.331e-13},
    "MDK-NCL": {"d538g": 1.876, "wt": 0.708, "pval": 1.406e-12},
    "PTN-SDC4": {"d538g": 6.019, "wt": 0.277, "pval": 8.088e-36},
    "PTN-NCL": {"d538g": 2.111, "wt": 0.104, "pval": 6.378e-26},
    "MIF-\nCD74_CD44": {"d538g": 1.631, "wt": 0.301, "pval": 8.506e-04},
}

# Cell profiles for COMMOT analysis (matches vignette_2)
CELL_PROFILES = {
    "Cancer Cells": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1", "KRT5-1"]},
    "Macrophages": {"Major": ["CD68-1"], "Minor": ["CD14-1"]},
    "CD4 T Cells": {"Major": ["CD3E-1", "CD4-1"]},
    "CD8 T Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "B Cells": {"Major": ["MS4A1-1", "CD19-1"]},
    "Endothelial Cells": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
}


# ============================================================================
# Data loading utilities
# ============================================================================

def load_spatial_coords():
    """Load tissue positions from Visium SpaceRanger output."""
    if SPATIAL_CSV.exists():
        df = pd.read_csv(SPATIAL_CSV, index_col="barcode")
        # Keep only in_tissue spots
        df = df[df["in_tissue"] == 1]
        print(f"  Loaded {len(df)} in-tissue spots from tissue_positions.csv")
        return df
    else:
        print(f"  WARNING: Spatial coordinates not found at {SPATIAL_CSV}")
        return None


def load_proportions():
    """Load cell type proportions (prefer unified module3 output)."""
    if PROPS_CSV.exists():
        df = pd.read_csv(PROPS_CSV, index_col=0)
        print(f"  Loaded unified proportions: {df.shape}")
        return df
    elif PROPS_CSV_FALLBACK.exists():
        df = pd.read_csv(PROPS_CSV_FALLBACK, index_col=0)
        print(f"  Loaded fallback proportions (Profile_X names): {df.shape}")
        return df
    else:
        print("  WARNING: No proportion data found")
        return None


def load_gene_expression():
    """Load deconvolved gene expression parquet."""
    for path in [GEX_PARQUET, GEX_PARQUET_FALLBACK]:
        if path.exists():
            try:
                df = pd.read_parquet(path)
                print(f"  Loaded gene expression: {df.shape}")
                return df
            except Exception as e:
                print(f"  WARNING: Could not read parquet {path}: {e}")
    print("  WARNING: No gene expression data found")
    return None


def load_module3_meta():
    """Load module3 metadata JSON for cell type mapping."""
    if META_JSON.exists():
        with open(META_JSON) as f:
            meta = json.load(f)
        print(f"  Loaded module3 meta: {len(meta.get('cell_types', []))} cell types")
        return meta
    return None


def load_visium_adata():
    """Load Visium AnnData with H&E image using squidpy."""
    visium_path = DATA_ROOT / SAMPLE / "outs"
    if not visium_path.exists():
        print(f"  WARNING: Visium outs not found at {visium_path}")
        return None
    try:
        import squidpy as sq
        adata = sq.read.visium(
            path=str(visium_path),
            counts_file='filtered_feature_bc_matrix.h5',
            load_images=True,
            gex_only=False,
        )
        print(f"  Loaded Visium adata: {adata.shape}, images loaded")
        return adata
    except Exception as e:
        print(f"  WARNING: Could not load Visium data: {e}")
        return None


def merge_spatial_data(props_df, spatial_df):
    """Merge proportions with spatial coordinates on barcode index."""
    if props_df is None or spatial_df is None:
        return None
    common = props_df.index.intersection(spatial_df.index)
    if len(common) == 0:
        print("  WARNING: No overlapping barcodes between proportions and spatial")
        return None
    merged = props_df.loc[common].join(
        spatial_df[["array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]],
        how="inner",
    )
    print(f"  Merged spatial + proportions: {len(merged)} spots")
    return merged


def compute_basal_score(gex_df):
    """Compute cytokeratin signature score per spot from deconvolved GEX.

    Filters to cancer cell rows (Cancer_Luminal, Cancer_Basal) and aggregates
    to spot-level. Tries basal CK genes first; falls back to luminal CKs for
    ER+ samples where basal keratins are not expressed at detectable levels.
    """
    available = [g for g in BASAL_CK_GENES if g in gex_df.columns]
    used_fallback = False
    if len(available) == 0:
        available = [g for g in LUMINAL_CK_GENES if g in gex_df.columns]
        used_fallback = True
        if len(available) == 0:
            print("  WARNING: No CK genes found in GEX columns")
            return None

    # Filter to cancer cell rows and extract barcode
    cancer_mask = gex_df.index.str.contains("Cancer", case=False)
    cancer_df = gex_df.loc[cancer_mask, available]
    if len(cancer_df) == 0:
        print("  WARNING: No cancer cell rows found in GEX")
        return None

    # Aggregate per spot (mean across cancer cell types)
    barcodes = cancer_df.index.str.split(":::").str[0]
    cancer_df = cancer_df.copy()
    cancer_df.index = barcodes
    score = cancer_df.sum(axis=1).groupby(level=0).mean()

    # Normalize to [0, 1]
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    label = "luminal CK (fallback)" if used_fallback else "basal CK"
    print(f"  CK score ({label}) from {len(available)} genes: {available}")
    print(f"  Aggregated to {len(score)} spots from {len(cancer_df)} cancer rows")
    return score


def compute_d538g_score(gex_df):
    """
    Compute a D538G mutation signature score per spot.

    Filters to cancer cell rows and aggregates to spot-level.
    Uses EstroGene upregulated genes as a surrogate for mutation signal.
    High score = likely mutant spot, low score = likely WT spot.
    """
    up_genes = [g for g in ESTROGENE_D538G_UP if g in gex_df.columns]
    down_genes = [g for g in ESTROGENE_D538G_DOWN if g in gex_df.columns]
    if len(up_genes) == 0:
        print("  WARNING: No EstroGene D538G upregulated genes found")
        return None

    # Filter to cancer cell rows
    cancer_mask = gex_df.index.str.contains("Cancer", case=False)
    cancer_df = gex_df.loc[cancer_mask].copy()
    if len(cancer_df) == 0:
        print("  WARNING: No cancer cell rows found in GEX")
        return None

    score = cancer_df[up_genes].sum(axis=1)
    if len(down_genes) > 0:
        score = score - cancer_df[down_genes].sum(axis=1)

    # Aggregate per spot (mean across cancer cell types)
    barcodes = score.index.str.split(":::").str[0]
    score.index = barcodes
    score = score.groupby(level=0).mean()

    # Normalize
    if score.max() > score.min():
        score = (score - score.min()) / (score.max() - score.min())
    print(f"  D538G score from {len(up_genes)} up + {len(down_genes)} down genes")
    print(f"  Aggregated to {len(score)} spots")
    return score


def classify_d538g_spots(d538g_score, threshold_quantile=0.75):
    """Classify spots as D538G-high vs WT-like based on score quantile."""
    if d538g_score is None:
        return None
    threshold = d538g_score.quantile(threshold_quantile)
    labels = pd.Series("WT-like", index=d538g_score.index)
    labels[d538g_score >= threshold] = "D538G-high"
    print(
        f"  D538G classification: {(labels == 'D538G-high').sum()} D538G-high, "
        f"{(labels == 'WT-like').sum()} WT-like (threshold q={threshold_quantile})"
    )
    return labels


def compute_mdk_score(gex_df):
    """Compute MDK expression score per spot from deconvolved cancer GEX layer.

    Extracts MDK from Cancer_Luminal + Cancer_Basal rows, aggregates to spot level.
    """
    if "MDK" not in gex_df.columns:
        print("  WARNING: MDK gene not found in GEX columns")
        return None

    cancer_mask = gex_df.index.str.contains("Cancer", case=False)
    cancer_df = gex_df.loc[cancer_mask, ["MDK"]].copy()
    if len(cancer_df) == 0:
        print("  WARNING: No cancer cell rows found in GEX for MDK")
        return None

    barcodes = cancer_df.index.str.split(":::").str[0]
    cancer_df.index = barcodes
    score = cancer_df["MDK"].groupby(level=0).mean()

    print(f"  MDK score: {len(score)} spots, range [{score.min():.3f}, {score.max():.3f}]")
    return score


def load_precomputed_commot_scores():
    """Load precomputed COMMOT MDK sender scores from saved expanded_adata.

    Loads the expanded_adata_wMutation.h5ad file from vignette 2 which contains
    COMMOT analysis results. Sums all s-MDK-* pathways for total MDK signaling.

    Returns:
        pd.Series indexed by barcode with summed COMMOT MDK sender scores,
        or None if the file doesn't exist.
    """
    import scanpy as sc

    # Load from vignette 2 saved adata
    commot_file = PROJECT_ROOT / "examples" / "expanded_adata_wMutation.h5ad"

    if not commot_file.exists():
        print(f"  WARNING: Precomputed COMMOT adata not found at {commot_file}")
        print("  Run vignette_2_surgical_d538g.ipynb first to generate it.")
        return None

    print(f"  Loading precomputed COMMOT results from {commot_file}...")
    expanded_adata = sc.read_h5ad(commot_file)

    if "commot-cellchat-sum-sender" not in expanded_adata.obsm:
        print("  WARNING: COMMOT sender scores not in adata.obsm")
        return None

    sender_signal = expanded_adata.obsm["commot-cellchat-sum-sender"]

    # Find all MDK pathways and sum them
    mdk_cols = [c for c in sender_signal.columns if c.startswith("s-MDK-")]
    if not mdk_cols:
        print(f"  WARNING: No s-MDK-* columns found. Available: {sender_signal.columns.tolist()}")
        return None

    print(f"    Found MDK pathways: {mdk_cols}")

    # Filter to cancer cells only
    cancer_mask = expanded_adata.obs["celltype"] == "Cancer Cells"
    cancer_sender = sender_signal.loc[cancer_mask, mdk_cols]

    # Sum all MDK pathways
    mdk_sum = cancer_sender.sum(axis=1)
    mdk_sum.name = "s-MDK-total"

    # Extract barcode (remove _Cancer Cells suffix)
    mdk_sum.index = mdk_sum.index.str.replace("_Cancer Cells", "")

    # Aggregate to spot level (in case of duplicates)
    spot_scores = mdk_sum.groupby(mdk_sum.index).mean()

    print(f"    Summed MDK score: {len(spot_scores)} spots, range [{spot_scores.min():.3f}, {spot_scores.max():.3f}]")
    return spot_scores


# ============================================================================
# Panel rendering functions
# ============================================================================

def _plot_spatial_on_he(ax, adata, score, cmap, score_label, title_text, panel_label,
                        vmin_override=None, vmax_override=None):
    """Helper: plot a score on H&E backdrop using sc.pl.spatial.

    Args:
        ax: matplotlib axes
        adata: AnnData loaded with squidpy (has H&E image)
        score: pd.Series indexed by barcode
        cmap: colormap
        score_label: label for the colorbar
        title_text: title for the panel
        panel_label: letter label (A, B, E)
    """
    import scanpy as sc

    ax.text(-0.08, 1.05, panel_label, fontsize=14, fontweight="bold", va="top",
            transform=ax.transAxes)

    if adata is None or score is None:
        ax.text(0.5, 0.5, f"{title_text}\n(data not available)",
                ha="center", va="center", fontsize=10, style="italic",
                color=PALETTE["neutral"])
        ax.axis("off")
        return

    # Align barcodes
    common = adata.obs_names.intersection(score.index)
    if len(common) == 0:
        ax.text(0.5, 0.5, f"{title_text}\n(no overlapping barcodes)",
                ha="center", va="center", fontsize=10, style="italic",
                color=PALETTE["neutral"])
        ax.axis("off")
        return

    adata_sub = adata[common].copy()
    adata_sub.obs[score_label] = score.loc[common].values

    sc.pl.spatial(
        adata_sub,
        color=score_label,
        cmap=cmap,
        ax=ax,
        show=False,
        title='',
        frameon=False,
        size=1.0,
        alpha=0.7,
        vmin=vmin_override,
        vmax=vmax_override,
        colorbar_loc=None,  # We add our own colorbar
    )

    # Add colorbar manually
    vmin = adata_sub.obs[score_label].min() if vmin_override is None else vmin_override
    vmax = adata_sub.obs[score_label].max() if vmax_override is None else vmax_override
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.035, pad=0.02, aspect=20)
    cbar.set_label(score_label, fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    ax.set_title(title_text, fontsize=10, fontweight="bold", loc="left", pad=6)


def panel_a_basal_ck(ax, adata, basal_score):
    """Panel A: Basal cytokeratin signatures in cancer cell layer (H&E spatial)."""
    _plot_spatial_on_he(
        ax, adata, basal_score,
        cmap=BASAL_CMAP,
        score_label="CK Score",
        title_text="Cytokeratin Signatures\n(Cancer Cell Layer)",
        panel_label="A",
    )


def panel_b_mutation_spatial(ax, adata, d538g_score):
    """Panel B: Spatial distribution of D538G mutation signals (H&E spatial)."""
    _plot_spatial_on_he(
        ax, adata, d538g_score,
        cmap=MUTATION_CMAP,
        score_label="D538G Score",
        title_text="ESR1 D538G Mutation Signal\n(Dispersed Across Tissue)",
        panel_label="B",
    )


def panel_c_estrogene(fig, gs_cell, gex_df, d538g_labels):
    """Panel C: Split violin + EstroGene gene heatmap.

    LEFT HALF: Split violin (D538G-high vs WT-like) of EstroGene signature score
    RIGHT HALF: Heatmap of 19 individual EstroGene signature genes
    """
    # Split the panel C area into two sub-areas
    gs_inner = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_cell,
                                       width_ratios=[1, 1.2], wspace=0.35)
    ax_violin = fig.add_subplot(gs_inner[0])
    ax_heat = fig.add_subplot(gs_inner[1])

    # Panel label on the leftmost sub-axis
    ax_violin.text(-0.18, 1.05, "C", fontsize=14, fontweight="bold", va="top",
                   transform=ax_violin.transAxes)

    if gex_df is None or d538g_labels is None:
        ax_violin.text(0.5, 0.5, "EstroGene data not available",
                       ha="center", va="center", fontsize=10, style="italic",
                       color=PALETTE["neutral"])
        ax_violin.axis("off")
        ax_heat.axis("off")
        return

    # ---- Gather per-gene data for both halves ----
    all_sig_genes = ESTROGENE_D538G_UP + ESTROGENE_D538G_DOWN
    available_genes = [g for g in all_sig_genes if g in gex_df.columns]
    if len(available_genes) == 0:
        ax_violin.text(0.5, 0.5, "No EstroGene genes found",
                       ha="center", va="center", fontsize=10, style="italic",
                       color=PALETTE["neutral"])
        ax_violin.axis("off")
        ax_heat.axis("off")
        return

    # Filter to cancer cells, aggregate to spot-level
    cancer_mask = gex_df.index.str.contains("Cancer", case=False)
    cancer_gex = gex_df.loc[cancer_mask, available_genes].copy()
    barcodes = cancer_gex.index.str.split(":::").str[0]
    cancer_gex.index = barcodes
    spot_gex = cancer_gex.groupby(level=0).mean()

    common = spot_gex.index.intersection(d538g_labels.index)
    spot_gex = spot_gex.loc[common]
    labels = d538g_labels.loc[common]

    d538g_mask = labels == "D538G-high"
    wt_mask = labels == "WT-like"

    # ---- LEFT: Split violin of composite signature score ----
    up_avail = [g for g in ESTROGENE_D538G_UP if g in available_genes]
    down_avail = [g for g in ESTROGENE_D538G_DOWN if g in available_genes]

    sig_score = spot_gex[up_avail].mean(axis=1)
    if len(down_avail) > 0:
        sig_score = sig_score - spot_gex[down_avail].mean(axis=1)

    d538g_vals = sig_score[d538g_mask].values
    wt_vals = sig_score[wt_mask].values

    vp_data = [wt_vals, d538g_vals]
    vp_colors = [PALETTE["secondary"], PALETTE["highlight"]]
    vp_labels = ["WT-like", "D538G-high"]
    positions = [0, 1]

    parts = ax_violin.violinplot(vp_data, positions=positions, showmeans=False,
                                  showmedians=True, showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(vp_colors[i])
        pc.set_alpha(0.7)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.5)
    parts['cmedians'].set_color("black")
    parts['cmedians'].set_linewidth(1.5)

    ax_violin.set_xticks(positions)
    ax_violin.set_xticklabels(vp_labels, fontsize=9)
    ax_violin.set_ylabel("EstroGene\nSignature Score", fontsize=9)
    ax_violin.set_title("EstroGene D538G\nSignature", fontsize=10,
                        fontweight="bold", loc="left", pad=6)

    # Significance annotation
    from scipy import stats
    try:
        stat, pval = stats.mannwhitneyu(d538g_vals, wt_vals, alternative="greater")
        sig_text = _pval_to_stars(pval)
        y_max = max(np.percentile(d538g_vals, 95), np.percentile(wt_vals, 95))
        bar_y = y_max * 1.08
        ax_violin.plot([0, 0, 1, 1], [bar_y * 0.98, bar_y, bar_y, bar_y * 0.98],
                       lw=1.0, c="black")
        ax_violin.text(0.5, bar_y * 1.02, sig_text, ha="center", va="bottom",
                       fontsize=10)
    except Exception:
        pass

    # ---- RIGHT: Heatmap of individual genes ----
    # Compute mean expression per group for each gene
    d538g_means = spot_gex.loc[d538g_mask, available_genes].mean(axis=0)
    wt_means = spot_gex.loc[wt_mask, available_genes].mean(axis=0)

    heat_data = pd.DataFrame({
        "D538G-high": d538g_means,
        "WT-like": wt_means,
    }, index=available_genes)

    # Z-score per gene (row) for visualization
    row_mean = heat_data.mean(axis=1)
    row_std = heat_data.std(axis=1).replace(0, 1)
    heat_z = heat_data.sub(row_mean, axis=0).div(row_std, axis=0)

    # Identify up vs down genes
    up_set = set(ESTROGENE_D538G_UP)
    gene_direction = ["up" if g in up_set else "down" for g in available_genes]

    # Plot heatmap
    im = ax_heat.imshow(heat_z.values, aspect='auto', cmap='RdBu_r',
                        vmin=-2, vmax=2, interpolation='nearest')

    ax_heat.set_xticks([0, 1])
    ax_heat.set_xticklabels(["D538G-high", "WT-like"], fontsize=9, rotation=30,
                            ha="right")
    ax_heat.set_yticks(range(len(available_genes)))
    # Color gene labels by direction
    ax_heat.set_yticklabels(available_genes, fontsize=8)
    for i, (label, direction) in enumerate(zip(ax_heat.get_yticklabels(), gene_direction)):
        if direction == "up":
            label.set_color(PALETTE["highlight"])
        else:
            label.set_color(PALETTE["primary"])

    ax_heat.set_title("Gene Expression\n(Z-score)", fontsize=10,
                      fontweight="bold", loc="left", pad=6)

    # Colorbar for heatmap
    cbar = plt.colorbar(im, ax=ax_heat, fraction=0.06, pad=0.04, aspect=25)
    cbar.set_label("Z-score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Add direction legend
    up_patch = mpatches.Patch(color=PALETTE["highlight"], label="Up in D538G")
    down_patch = mpatches.Patch(color=PALETTE["primary"], label="Down in D538G")
    ax_heat.legend(handles=[up_patch, down_patch], loc="lower right",
                   fontsize=8, framealpha=0.8, borderpad=0.3)

    # Restore spines for heatmap
    for spine in ax_heat.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)


def panel_d_pathway_dotplot(ax):
    """Panel D: Combined Hallmark + KEGG pathway dot plot.

    Rows: All 13 pathways (8 Hallmark + 5 KEGG) with section separator.
    X-axis: Enrichr Combined Score.
    Dot size: proportional to gene overlap count.
    Dot color: -log10(adjusted p-value) from Enrichr output.

    Note: Panel is horizontally compressed to avoid overlap with adjacent panels.
    """
    ax.text(-0.08, 1.05, "D", fontsize=14, fontweight="bold", va="top",
            transform=ax.transAxes)

    # Build combined pathway data
    hallmark_names = list(HALLMARK_PATHWAYS.keys())
    hallmark_scores = list(HALLMARK_PATHWAYS.values())
    kegg_names = list(KEGG_PATHWAYS.keys())
    kegg_scores = list(KEGG_PATHWAYS.values())

    # Real gene overlap counts from Enrichr output
    hallmark_overlap = [20, 17, 19, 18, 16, 14, 12, 11]
    kegg_overlap = [7, 6, 9, 7, 8]

    # All pathways: Hallmark first (top), then KEGG (bottom)
    all_names = hallmark_names + kegg_names
    all_scores = hallmark_scores + kegg_scores
    all_overlap = hallmark_overlap + kegg_overlap
    all_source = ["Hallmark"] * len(hallmark_names) + ["KEGG"] * len(kegg_names)

    n_total = len(all_names)
    y_pos = np.arange(n_total)

    # Real adjusted p-values from Enrichr output
    hallmark_adjp = [1.167e-14, 2.958e-13, 8.717e-14, 6.112e-13, 7.871e-11, 7.912e-09, 5.940e-07, 4.260e-06]
    kegg_adjp = [3.855e-03, 5.056e-03, 3.855e-03, 7.612e-03, 1.050e-02]
    all_adjp = hallmark_adjp + kegg_adjp
    neg_log10_p = np.array([-np.log10(p) for p in all_adjp])

    # Normalize dot sizes for visibility (smaller dots to fit compressed panel)
    sizes = np.array(all_overlap)
    size_scale = (sizes - sizes.min()) / (sizes.max() - sizes.min()) * 150 + 40

    # Color by -log10(p) - use viridis for perceptual uniformity
    scatter = ax.scatter(all_scores, y_pos, s=size_scale, c=neg_log10_p,
                         cmap="viridis", edgecolors="black", linewidth=0.5,
                         zorder=3)

    # Mark Hallmark vs KEGG with different marker shapes using overlay
    for i, source in enumerate(all_source):
        if source == "KEGG":
            ax.scatter(all_scores[i], y_pos[i], s=size_scale[i] * 0.3,
                       marker="s", c="white", edgecolors="none", zorder=4, alpha=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_names, fontsize=6.5)
    ax.set_xlabel("Enrichr Combined Score", fontsize=8)
    ax.invert_yaxis()

    # Set tight x-axis limits - reduced padding to prevent right overflow
    max_score = max(all_scores)
    ax.set_xlim(0, max_score * 1.04)

    # Add separator line between Hallmark and KEGG
    sep_y = len(hallmark_names) - 0.5
    ax.axhline(y=sep_y, color=PALETTE["border"], linestyle="--", linewidth=0.8)

    # Section labels inside plot area (moved left to avoid overlap with colorbar)
    ax.text(max_score * 0.80, (len(hallmark_names) - 1) / 2, "Hallmark",
            fontsize=7, va="center", ha="right",
            color=PALETTE["highlight"], fontweight="bold", rotation=-90)
    kegg_y_center = len(hallmark_names) + (len(kegg_names) - 1) / 2
    ax.text(max_score * 0.80, kegg_y_center, "KEGG",
            fontsize=7, va="center", ha="right",
            color=PALETTE["accent2"], fontweight="bold", rotation=-90)

    # Colorbar - smaller to prevent right overflow
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.022, pad=0.008, aspect=35)
    cbar.set_label("-log10(p)", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    # Size legend (adjusted to match smaller dot sizes)
    for sz_val, sz_label in [(40, "7"), (95, "14"), (150, "20")]:
        ax.scatter([], [], s=sz_val, c="gray", edgecolors="black",
                   linewidth=0.5, label=f"{sz_label} genes")
    leg = ax.legend(title="Gene\noverlap", fontsize=6, title_fontsize=6,
                    loc="lower right", framealpha=0.8, borderpad=0.2,
                    handletextpad=0.2, labelspacing=0.3)

    ax.set_title("Pathway Enrichment\n(D538G vs WT)", fontsize=9,
                 fontweight="bold", loc="left", pad=4)


def panel_e_mdk_spatial(ax, adata, commot_score):
    """Panel E: COMMOT MDK sender signaling score spatial map (H&E backdrop).

    Shows the summed COMMOT communication score for all MDK ligand-receptor pairs,
    which captures actual signaling activity rather than raw gene expression.
    """
    _plot_spatial_on_he(
        ax, adata, commot_score,
        cmap=MDK_CMAP,
        score_label="s-MDK (sum)",
        title_text="MDK Signaling\n(COMMOT Sender Score)",
        panel_label="E",
        vmax_override=10,
    )


def panel_f_commot(ax):
    """Panel F: COMMOT signaling differences between D538G and WT (bar plot)."""
    ax.text(-0.12, 1.05, "F", fontsize=14, fontweight="bold", va="top",
            transform=ax.transAxes)

    pairs = list(COMMOT_LIGANDS.keys())
    d538g_vals = [COMMOT_LIGANDS[p]["d538g"] for p in pairs]
    wt_vals = [COMMOT_LIGANDS[p]["wt"] for p in pairs]
    pvals = [COMMOT_LIGANDS[p]["pval"] for p in pairs]

    x = np.arange(len(pairs))
    width = 0.36  # Slightly wider bars

    ax.bar(x - width / 2, wt_vals, width, label="WT-like",
           color=PALETTE["secondary"], alpha=0.8, edgecolor="none")
    ax.bar(x + width / 2, d538g_vals, width, label="D538G-high",
           color=PALETTE["highlight"], alpha=0.8, edgecolor="none")

    # Significance stars and fold-change annotations
    for i, pval in enumerate(pvals):
        stars = _pval_to_stars(pval)
        y_top = max(d538g_vals[i], wt_vals[i]) * 1.05
        # Fold-change annotation
        fc = d538g_vals[i] / wt_vals[i] if wt_vals[i] > 0 else float('inf')
        fc_text = f"{fc:.1f}x"
        ax.text(x[i], y_top + 0.15, fc_text, ha="center", va="bottom", fontsize=9,
                color=PALETTE["neutral"], style="italic")
        if stars:
            ax.text(x[i], y_top + 0.55, stars, ha="center", va="bottom", fontsize=9,
                    fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(pairs, fontsize=9, rotation=20, ha="right")
    ax.set_ylabel("Mean Signal\nStrength", fontsize=9)
    ax.legend(fontsize=9, loc="upper right", framealpha=0.8)
    ax.set_title("COMMOT Signaling Analysis\n(D538G vs WT Cancer Spots)", fontsize=10,
                 fontweight="bold", loc="left", pad=6)


def panel_g_elisa(ax):
    """Panel G: ELISA validation placeholder."""
    ax.text(-0.08, 1.05, "G", fontsize=14, fontweight="bold", va="top",
            transform=ax.transAxes)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Styled placeholder box
    rect = plt.Rectangle((0.05, 0.05), 0.9, 0.9, linewidth=1.5,
                          edgecolor=PALETTE["border"], facecolor=PALETTE["background"],
                          linestyle="--", transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)

    ax.text(0.5, 0.6, "ELISA Validation", ha="center", va="center",
            fontsize=11, fontweight="bold", color=PALETTE["primary"],
            transform=ax.transAxes)
    ax.text(0.5, 0.45, "MCF7 WT vs D538G\nMidkine Secretion (p < 0.0001)",
            ha="center", va="center", fontsize=9, color=PALETTE["neutral"],
            transform=ax.transAxes)
    ax.text(0.5, 0.25, "See Prism data", ha="center", va="center",
            fontsize=9, style="italic", color=PALETTE["neutral"],
            transform=ax.transAxes)


def panel_h_if(ax):
    """Panel H: IF validation placeholder."""
    ax.text(-0.08, 1.05, "H", fontsize=14, fontweight="bold", va="top",
            transform=ax.transAxes)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Styled placeholder box
    rect = plt.Rectangle((0.05, 0.05), 0.9, 0.9, linewidth=1.5,
                          edgecolor=PALETTE["border"], facecolor=PALETTE["background"],
                          linestyle="--", transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)

    ax.text(0.5, 0.6, "Immunofluorescence", ha="center", va="center",
            fontsize=11, fontweight="bold", color=PALETTE["primary"],
            transform=ax.transAxes)
    ax.text(0.5, 0.45, "MCF7 WT vs D538G\nMembrane MDK (p < 0.001)\n"
            "Intracellular MDK (p < 0.01)",
            ha="center", va="center", fontsize=9, color=PALETTE["neutral"],
            transform=ax.transAxes)
    ax.text(0.5, 0.22, "See microscopy images", ha="center", va="center",
            fontsize=9, style="italic", color=PALETTE["neutral"],
            transform=ax.transAxes)


# ============================================================================
# Helper
# ============================================================================

def _pval_to_stars(pval):
    """Convert p-value to significance stars."""
    if pval < 0.0001:
        return "****"
    elif pval < 0.001:
        return "***"
    elif pval < 0.01:
        return "**"
    elif pval < 0.05:
        return "*"
    else:
        return "ns"


# ============================================================================
# Main figure assembly
# ============================================================================

def generate_figure4():
    """Generate Figure 4: Midkine/ESR1 Case Study."""
    print("=" * 60)
    print("Figure 4: Midkine/ESR1 Case Study (HCC22-088-P4-S2_1i_rep)")
    print("=" * 60)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print("\nLoading data...")
    spatial_df = load_spatial_coords()
    props_df = load_proportions()
    gex_df = load_gene_expression()
    meta = load_module3_meta()

    # Load Visium AnnData with H&E for spatial panels
    print("\nLoading Visium AnnData with H&E image...")
    adata = load_visium_adata()

    # Merge spatial coordinates with proportions
    merged_df = merge_spatial_data(props_df, spatial_df)

    # Compute scores
    print("\nComputing signatures...")
    basal_score = compute_basal_score(gex_df) if gex_df is not None else None
    d538g_score = compute_d538g_score(gex_df) if gex_df is not None else None
    d538g_labels = classify_d538g_spots(d538g_score) if d538g_score is not None else None

    # Load precomputed COMMOT MDK sender score (from vignette 2 saved adata)
    print("\nLoading COMMOT MDK sender signaling score (summed)...")
    commot_score = load_precomputed_commot_scores()
    if commot_score is None:
        print("  ERROR: Precomputed COMMOT scores not found!")
        print("  Run 'sbatch manuscript/figures/sbatch_compute_commot.sh' first.")
        print("  Falling back to MDK expression score...")
        commot_score = compute_mdk_score(gex_df) if gex_df is not None else None

    # ------------------------------------------------------------------
    # Create figure: 4 rows x 2 columns
    # Layout:
    #   Row 1: Panel A (basal CK H&E)       | Panel B (D538G H&E)
    #   Row 2: Panel C (violin + heatmap)    | Panel D (pathway dot plot)
    #   Row 3: Panel E (MDK H&E)             | Panel F (COMMOT bars)
    #   Row 4: Panel G (ELISA placeholder)   | Panel H (IF placeholder)
    # ------------------------------------------------------------------
    print("\nGenerating figure...")
    fig = plt.figure(figsize=(12, 16))

    # Disable constrained_layout for manual GridSpec control
    fig.set_constrained_layout(False)

    gs = GridSpec(
        4, 2, figure=fig,
        hspace=0.35, wspace=0.30,
        left=0.08, right=0.95, top=0.96, bottom=0.03,
        height_ratios=[1, 0.95, 1, 0.7],
    )

    # ------------------------------------------------------------------
    # Row 1: Spatial panels with H&E
    # ------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # ------------------------------------------------------------------
    # Row 2: Panel C uses GridSpecFromSubplotSpec (no axis pre-created);
    #         Panel D gets a normal axis
    # ------------------------------------------------------------------
    # Panel C will create its own sub-axes inside gs[1, 0]
    ax_d = fig.add_subplot(gs[1, 1])

    # ------------------------------------------------------------------
    # Row 3: MDK spatial + COMMOT
    # ------------------------------------------------------------------
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    # ------------------------------------------------------------------
    # Row 4: Validation placeholders
    # ------------------------------------------------------------------
    ax_g = fig.add_subplot(gs[3, 0])
    ax_h = fig.add_subplot(gs[3, 1])

    # ------------------------------------------------------------------
    # Render panels
    # ------------------------------------------------------------------
    print("\nRendering panels...")

    print("  Panel A: Basal CK spatial heatmap (H&E)")
    panel_a_basal_ck(ax_a, adata, basal_score)

    print("  Panel B: D538G mutation signal spatial (H&E)")
    panel_b_mutation_spatial(ax_b, adata, d538g_score)

    print("  Panel C: EstroGene split violin + gene heatmap")
    panel_c_estrogene(fig, gs[1, 0], gex_df, d538g_labels)

    print("  Panel D: Combined pathway dot plot")
    panel_d_pathway_dotplot(ax_d)

    print("  Panel E: MDK COMMOT sender score spatial map (H&E, summed)")
    panel_e_mdk_spatial(ax_e, adata, commot_score)

    print("  Panel F: COMMOT signaling")
    panel_f_commot(ax_f)

    print("  Panel G: ELISA placeholder")
    panel_g_elisa(ax_g)

    print("  Panel H: IF placeholder")
    panel_h_if(ax_h)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    print("\nSaving figure...")
    for fmt, dpi in [("pdf", 300), ("png", 300), ("svg", None)]:
        output_path = OUTPUT_DIR / f"figure4_midkine_esr1.{fmt}"
        try:
            if fmt == "svg":
                fig.savefig(output_path, format="svg", bbox_inches="tight",
                            facecolor="white")
            else:
                fig.savefig(output_path, dpi=dpi, bbox_inches="tight",
                            facecolor="white")
            print(f"  Saved: {output_path}")
        except RuntimeError as e:
            # tight_layout / constrained_layout conflict with colorbars
            print(f"  WARNING: RuntimeError during save ({fmt}): {e}")
            if fmt == "svg":
                fig.savefig(output_path, format="svg", facecolor="white")
            else:
                fig.savefig(output_path, dpi=dpi, facecolor="white")
            print(f"  Saved (without bbox_inches='tight'): {output_path}")

    plt.close(fig)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Figure 4: Midkine/ESR1 Case Study -- COMPLETE")
    print("=" * 60)
    print("\nPanel A: Basal cytokeratin signatures (H&E spatial heatmap)")
    print("Panel B: D538G mutation signal distribution (H&E spatial)")
    print("Panel C: EstroGene D538G vs WT (split violin + gene heatmap)")
    print("Panel D: Combined Hallmark + KEGG pathway dot plot")
    print("Panel E: MDK COMMOT sender signaling score (H&E spatial, summed)")
    print("Panel F: COMMOT signaling differences (bar plot with fold-change)")
    print("Panel G: ELISA validation (PLACEHOLDER - see Prism data)")
    print("Panel H: IF validation (PLACEHOLDER - see microscopy images)")
    print(f"\nOutput: {OUTPUT_DIR}/figure4_midkine_esr1.{{pdf,png,svg}}")


if __name__ == "__main__":
    generate_figure4()
