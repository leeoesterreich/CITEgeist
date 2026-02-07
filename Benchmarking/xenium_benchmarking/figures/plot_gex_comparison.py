"""
Fair GEX deconvolution comparison: CITEgeist vs scResolve (Region 0).

Approach:
  1. Hungarian-match scResolve cells to GT cells (spatial + expression similarity)
  2. Transfer GT cell type labels to matched scResolve cells
  3. Aggregate scResolve cell GEX by GT-matched cell type into compartments
  4. If scResolve under/over-reports cells, unmatched cells count as wrong
     (missing = 0 expression in compartment, extra = ignored)
  5. Evaluate BOTH methods on ALL spots (zeros for missing compartments)
  6. Generate figure: spatial pie charts (top) + GEX metric bars (bottom)
"""

import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from scipy import stats
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csc_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking")
PSEUDO_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
OUTPUT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))

REGION = 0

# ── Constants ──────────────────────────────────────────────────────────────────
COLLAPSE_MAPPING = {
    "B cells": "B cells", "B_cells": "B cells",
    "CD8+ T cells": "CD8+ T cells", "CD8pos_T_cells": "CD8+ T cells",
    "Proliferating T": "CD8+ T cells", "Proliferating_T": "CD8+ T cells",
    "CD4+ T cells": "CD4+ T cells", "CD4pos_T_cells": "CD4+ T cells",
    "Mixed Immune": "CD4+ T cells", "Mixed_Immune": "CD4+ T cells",
    "Macrophages": "Macrophages",
    "Epithelial": "Epithelial",
    "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial", "Vascular_Stromal": "Endothelial",
    "Myofibroblasts": "Fibroblasts", "Stromal": "Fibroblasts",
    "Fibroblasts": "Fibroblasts",
}

ACHIEVABLE_7 = ["B cells", "CD4+ T cells", "CD8+ T cells",
                "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]

CELL_COLORS = {
    "B cells": "#1b9e77", "CD4+ T cells": "#d95f02", "CD8+ T cells": "#7570b3",
    "Macrophages": "#e7298a", "Endothelial": "#66a61e",
    "Epithelial": "#e6ab02", "Fibroblasts": "#a6761d",
}

GT_TO_ACHIEVABLE_7 = {
    "B cells": "B cells", "Mixed Immune": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells", "Proliferating T": "CD8+ T cells",
    "Macrophages": "Macrophages", "Endothelial": "Endothelial",
    "Vascular Stromal": "Endothelial", "Epithelial": "Epithelial",
    "Myofibroblasts": "Fibroblasts", "Stromal": "Fibroblasts",
}


# ── Step 1: Load ground truth GEX (genes x spots) ─────────────────────────────

def load_gt_gex():
    """Load and collapse ground truth GEX into achievable-7 compartments."""
    gt_dir = PSEUDO_DIR / "data" / "ground_truth_gex" / f"Xenium_region_{REGION}"
    gt_gex = {}
    for f in gt_dir.glob("*_GT.csv"):
        ct_name = f.stem.replace("_GT", "")
        ct_collapsed = COLLAPSE_MAPPING.get(ct_name, ct_name)
        if ct_collapsed not in ACHIEVABLE_7:
            continue
        df = pd.read_csv(f, index_col=0)  # genes x spots
        if ct_collapsed in gt_gex:
            gt_gex[ct_collapsed] = gt_gex[ct_collapsed].add(df, fill_value=0)
        else:
            gt_gex[ct_collapsed] = df.copy()
    logger.info(f"GT GEX: {len(gt_gex)} cell types loaded")
    return gt_gex


# ── Step 2: Load CITEgeist GEX layers ─────────────────────────────────────────

def load_citegeist_gex():
    """Load CITEgeist per-cell-type GEX layers (spots x genes -> genes x spots)."""
    layers_dir = (BASE_DIR / "CITEgeist" / "output_achievable_7_no_unknown" /
                  f"Xenium_region_{REGION}_pass1" / "layers" / "pass1")
    gex = {}
    for f in layers_dir.glob("*_layer_pass1.csv"):
        ct_name = f.stem.replace("_layer_pass1", "").replace("_", " ")
        if ct_name == "CD4  T cells":
            ct_name = "CD4+ T cells"
        if ct_name == "CD8  T cells":
            ct_name = "CD8+ T cells"
        if ct_name not in ACHIEVABLE_7:
            continue
        df = pd.read_csv(f, index_col=0)  # spots x genes
        gex[ct_name] = df.T  # genes x spots
    logger.info(f"CITEgeist GEX: {len(gex)} cell types loaded")
    return gex


# ── Step 3: Region-wide Hungarian-matched scResolve GEX ───────────────────────

MATCH_RADIUS = 100.0  # microns — spatial radius for candidate pairs


def compute_scresolve_hungarian_gex():
    """
    Region-wide Hungarian matching of scResolve cells to GT cells.

    Instead of matching within each spot independently, matches across the
    entire region using spatial radius filtering + connected component
    decomposition for tractable Hungarian assignment.

    Steps:
      1. Load GT and scResolve cell positions + expression
      2. Build bipartite graph: edges between GT-scResolve pairs within radius
      3. Find connected components in the bipartite graph
      4. Run Hungarian matching per component (cosine distance on expression)
      5. Transfer GT cell types to matched scResolve cells
      6. Aggregate into spot-level compartments (0 for missing)
    """
    import networkx as nx

    region_name = f"Xenium_region_{REGION}"

    # Load achievable-7 cell type annotations
    cell_types_path = PSEUDO_DIR / "data_rna_gt" / "cell_types_achievable7.csv"
    cell_types_df = pd.read_csv(cell_types_path, index_col=0)
    logger.info(f"Cell type annotations: {len(cell_types_df)}")

    # Load original Xenium cells (expression)
    logger.info("Loading original Xenium cell matrix...")
    h5_path = XENIUM_DIR / "cell_feature_matrix.h5"
    with h5py.File(h5_path, 'r') as f:
        data = f['matrix/data'][:]
        indices = f['matrix/indices'][:]
        indptr = f['matrix/indptr'][:]
        shape = f['matrix/shape'][:]
        gene_names = f['matrix/features/name'][:].astype(str)
        feature_types = f['matrix/features/feature_type'][:].astype(str)
        barcodes = f['matrix/barcodes'][:].astype(str)

    X_sparse = csc_matrix((data, indices, indptr), shape=shape)
    orig_X = X_sparse.T.toarray()
    gex_mask = feature_types == 'Gene Expression'
    orig_genes = gene_names[gex_mask]
    orig_X = orig_X[:, gex_mask]
    barcode_to_idx = {b: i for i, b in enumerate(barcodes)}
    logger.info(f"Original Xenium: {orig_X.shape[0]} cells x {orig_X.shape[1]} genes")

    # Load original cell spatial coordinates
    logger.info("Loading original cell coordinates...")
    cells_df = pd.read_csv(XENIUM_DIR / "cells.csv.gz",
                           usecols=['cell_id', 'x_centroid', 'y_centroid'])
    cells_df = cells_df.set_index('cell_id')

    # Load cell-to-spot mapping (to identify region cells and their spot assignments)
    mapping_df = pd.read_csv(PSEUDO_DIR / "data" / "cell_to_spot_mapping.csv", index_col=0)
    mapping_df = mapping_df[mapping_df['region_id'] == REGION]
    mapping_df = mapping_df[mapping_df['spot_idx'] >= 0]
    logger.info(f"Region {REGION}: {len(mapping_df)} mapped GT cells")

    # Get GT cell indices, coordinates, and barcodes
    gt_barcodes = mapping_df.index.values
    gt_global_idx = np.array([barcode_to_idx[bc] for bc in gt_barcodes
                              if bc in barcode_to_idx])
    gt_barcodes_valid = np.array([bc for bc in gt_barcodes if bc in barcode_to_idx])
    gt_coords = cells_df.loc[gt_barcodes_valid, ['x_centroid', 'y_centroid']].values
    gt_spot_ids = mapping_df.loc[gt_barcodes_valid, 'spot_id'].values
    logger.info(f"GT cells with coords: {len(gt_barcodes_valid)}")

    # Load scResolve cells
    cells_path = (BASE_DIR / "scResolve" / "output_morphology" / region_name /
                  "segmentation" / "cells_exp_loc.h5ad")
    scres_adata = sc.read_h5ad(cells_path)
    scres_X = scres_adata.X.toarray() if hasattr(scres_adata.X, 'toarray') else scres_adata.X
    scres_coords = scres_adata.obs[['x', 'y']].values
    scres_genes = scres_adata.var_names.values
    logger.info(f"scResolve: {len(scres_adata)} cells x {len(scres_genes)} genes")

    # Load spot info for output
    spots_adata = sc.read_h5ad(
        PSEUDO_DIR / "data" / "h5ad_objects" / f"{region_name}_GEX.h5ad")
    spot_coords = spots_adata.obsm['spatial']
    spot_ids = spots_adata.obs_names.values
    all_spots = list(spot_ids)

    # Find common genes
    common_genes = sorted(list(set(orig_genes) & set(scres_genes)))
    orig_gene_idx = np.array([np.where(orig_genes == g)[0][0] for g in common_genes])
    scres_gene_idx = np.array([np.where(scres_genes == g)[0][0] for g in common_genes])
    logger.info(f"Common genes: {len(common_genes)}")

    # ── Region-wide spatial radius matching ──
    logger.info(f"Building spatial candidate pairs (radius={MATCH_RADIUS}μm)...")

    # KDTree on GT cell coordinates
    gt_tree = cKDTree(gt_coords)
    # For each scResolve cell, find GT cells within radius
    candidate_pairs = gt_tree.query_ball_point(scres_coords, r=MATCH_RADIUS)

    # Build bipartite graph for connected component decomposition
    # Nodes: "gt_i" for GT cells, "sr_j" for scResolve cells
    G = nx.Graph()
    n_edges = 0
    for sr_idx, gt_neighbors in enumerate(candidate_pairs):
        if len(gt_neighbors) > 0:
            for gt_idx in gt_neighbors:
                G.add_edge(f"gt_{gt_idx}", f"sr_{sr_idx}")
                n_edges += 1

    logger.info(f"Bipartite graph: {G.number_of_nodes()} nodes, {n_edges} edges")

    # Find connected components
    components = list(nx.connected_components(G))
    logger.info(f"Connected components: {len(components)}")

    # Log component size distribution
    comp_sizes = [len(c) for c in components]
    logger.info(f"  Component sizes: min={min(comp_sizes)}, max={max(comp_sizes)}, "
                f"median={np.median(comp_sizes):.0f}, mean={np.mean(comp_sizes):.0f}")

    # Initialize result
    gex_result = {}
    for ct in ACHIEVABLE_7:
        gex_result[ct] = pd.DataFrame(0.0, index=common_genes, columns=all_spots)

    # Assign scResolve cells to nearest spot (for aggregation)
    spot_tree = cKDTree(spot_coords)
    _, scres_spot_idx = spot_tree.query(scres_coords)
    scres_spot_ids = spot_ids[scres_spot_idx]

    # Pre-compute normalized expression for cosine distance
    gt_X_common = orig_X[gt_global_idx][:, orig_gene_idx]
    gt_norms = np.linalg.norm(gt_X_common, axis=1, keepdims=True)
    gt_norms[gt_norms == 0] = 1
    gt_X_normed = gt_X_common / gt_norms

    scres_X_common = scres_X[:, scres_gene_idx]
    scres_norms = np.linalg.norm(scres_X_common, axis=1, keepdims=True)
    scres_norms[scres_norms == 0] = 1
    scres_X_normed = scres_X_common / scres_norms

    n_matched = 0
    n_unmatched_gt = 0
    n_extra_scres = 0

    logger.info("Running Hungarian matching per connected component...")

    for ci, component in enumerate(components):
        # Separate GT and scResolve nodes
        gt_nodes = sorted([int(n.split('_')[1]) for n in component if n.startswith('gt_')])
        sr_nodes = sorted([int(n.split('_')[1]) for n in component if n.startswith('sr_')])

        if len(gt_nodes) == 0 or len(sr_nodes) == 0:
            if len(gt_nodes) > 0:
                n_unmatched_gt += len(gt_nodes)
            if len(sr_nodes) > 0:
                n_extra_scres += len(sr_nodes)
            continue

        # Build dense cost matrix for this component
        gt_expr = gt_X_normed[gt_nodes]
        sr_expr = scres_X_normed[sr_nodes]
        cost = cdist(gt_expr, sr_expr, metric='cosine')
        cost = np.nan_to_num(cost, nan=1.0)

        n_gt, n_sr = len(gt_nodes), len(sr_nodes)

        if n_gt <= n_sr:
            row_ind, col_ind = linear_sum_assignment(cost)
            n_extra_scres += n_sr - n_gt
        else:
            col_ind_t, row_ind_t = linear_sum_assignment(cost.T)
            row_ind, col_ind = row_ind_t, col_ind_t
            n_unmatched_gt += n_gt - n_sr

        # Transfer cell types and aggregate
        for k in range(len(row_ind)):
            gt_local = gt_nodes[row_ind[k]]
            sr_local = sr_nodes[col_ind[k]]

            gt_bc = gt_barcodes_valid[gt_local]
            if gt_bc not in cell_types_df.index:
                continue
            cell_type = cell_types_df.loc[gt_bc, 'cell_type']
            if cell_type not in ACHIEVABLE_7:
                continue

            # Assign to the spot of the matched scResolve cell
            spot_id = scres_spot_ids[sr_local]
            scres_cell_expr = scres_X_common[sr_local]
            gex_result[cell_type][spot_id] += scres_cell_expr
            n_matched += 1

        if (ci + 1) % 500 == 0:
            logger.info(f"  {ci+1}/{len(components)} components processed...")

    logger.info(f"Region-wide matching complete: {n_matched} matched, "
                f"{n_unmatched_gt} unmatched GT cells, "
                f"{n_extra_scres} extra scResolve cells (ignored)")

    return gex_result


# ── Step 4: Compute fair GEX metrics on ALL spots ─────────────────────────────

def compute_fair_gex_metrics(gt_gex, pred_gex, method_name):
    """Compute per-spot gene profile correlation on jointly-active spots.

    For each cell type:
      1. Identify spots where BOTH GT and prediction have nonzero expression
      2. At each such spot, compute Pearson r across genes (log1p values)
         — this measures gene profile similarity independent of library size
      3. Report mean per-spot r, number of active spots, and coverage
         (fraction of GT-active spots that the method also predicts active)

    This avoids the zero-inflation problem of flattening sparse matrices
    and directly measures gene expression profile quality.
    """
    results = {}
    for ct in ACHIEVABLE_7:
        if ct not in gt_gex:
            continue
        gt_df = gt_gex[ct]
        if ct not in pred_gex:
            results[ct] = {
                'spot_pearson_r': 0.0, 'coverage': 0.0,
                'n_genes': 0, 'n_active_spots': 0, 'n_gt_active': 0,
            }
            continue

        pred_df = pred_gex[ct]

        # Align to common genes and ALL GT spots
        common_genes = gt_df.index.intersection(pred_df.index)
        all_spots = gt_df.columns

        if len(common_genes) == 0:
            results[ct] = {
                'spot_pearson_r': 0.0, 'coverage': 0.0,
                'n_genes': 0, 'n_active_spots': 0, 'n_gt_active': 0,
            }
            continue

        gt_vals = gt_df.loc[common_genes, all_spots].values       # genes x spots
        pred_aligned = pred_df.reindex(index=common_genes, columns=all_spots, fill_value=0.0)
        pred_vals = pred_aligned.values                            # genes x spots

        # Identify active spots
        gt_active = gt_vals.sum(axis=0) > 0
        pred_active = pred_vals.sum(axis=0) > 0
        both_active = gt_active & pred_active
        n_gt_active = int(gt_active.sum())
        n_both_active = int(both_active.sum())
        coverage = n_both_active / n_gt_active if n_gt_active > 0 else 0.0

        # Per-spot gene profile Pearson r on jointly-active spots
        spot_rs = []
        if n_both_active > 0:
            gt_active_log = np.log1p(gt_vals[:, both_active])
            pred_active_log = np.log1p(pred_vals[:, both_active])
            for j in range(n_both_active):
                gvec = gt_active_log[:, j]
                pvec = pred_active_log[:, j]
                if np.std(gvec) > 0 and np.std(pvec) > 0:
                    r, _ = stats.pearsonr(gvec, pvec)
                    spot_rs.append(r)

        mean_spot_r = float(np.mean(spot_rs)) if spot_rs else 0.0

        results[ct] = {
            'spot_pearson_r': mean_spot_r,
            'coverage': float(coverage),
            'n_genes': int(len(common_genes)),
            'n_active_spots': n_both_active,
            'n_gt_active': n_gt_active,
        }
        logger.info(f"  {method_name} - {ct}: spot r={mean_spot_r:.4f}, "
                     f"coverage={coverage:.2%} ({n_both_active}/{n_gt_active} spots), "
                     f"{len(common_genes)} genes")

    # Overall
    valid = [v for k, v in results.items() if k != '_overall']
    all_r = [v['spot_pearson_r'] for v in valid]
    all_cov = [v['coverage'] for v in valid]
    results['_overall'] = {
        'mean_spot_pearson_r': float(np.mean(all_r)) if all_r else 0.0,
        'mean_coverage': float(np.mean(all_cov)) if all_cov else 0.0,
    }
    return results


# ── Step 5: Spatial pie chart drawing (shared with other script) ───────────────

def load_spatial_coords(h5ad_path):
    with h5py.File(h5ad_path, 'r') as f:
        coords = f['obsm']['spatial'][:]
        names = [x.decode() if isinstance(x, bytes) else x for x in f['obs']['_index'][:]]
    return pd.DataFrame(coords, index=names, columns=['x', 'y'])


def collapse_gt_to_achievable_7(gt_df):
    collapsed = pd.DataFrame(0.0, index=gt_df.index, columns=ACHIEVABLE_7)
    for src, dst in GT_TO_ACHIEVABLE_7.items():
        if src in gt_df.columns:
            collapsed[dst] += gt_df[src]
    row_sums = collapsed.sum(axis=1).replace(0, 1)
    return collapsed.div(row_sums, axis=0)


def load_proportions():
    """Load proportion data for spatial pie charts."""
    h5ad = PSEUDO_DIR / "data" / "h5ad_objects" / f"Xenium_region_{REGION}_GEX.h5ad"
    coords = load_spatial_coords(h5ad)

    gt_raw = pd.read_csv(
        PSEUDO_DIR / "data_granular_gt" / "ground_truth" / f"Xenium_region_{REGION}_prop.csv",
        index_col=0)
    gt = collapse_gt_to_achievable_7(gt_raw)

    pred = pd.read_csv(
        BASE_DIR / "CITEgeist" / "output_achievable_7_no_unknown" /
        f"Xenium_region_{REGION}_cell_prop_global_results.csv", index_col=0)
    if "Unknown" in pred.columns:
        pred = pred.drop(columns=["Unknown"])
    pred = pred[ACHIEVABLE_7]
    pred = pred.div(pred.sum(axis=1).replace(0, 1), axis=0)

    sr = pd.read_csv(
        BASE_DIR / "scResolve" / "output_morphology" / f"Xenium_region_{REGION}" /
        f"Xenium_region_{REGION}_fair_proportions.csv", index_col=0)
    sr = sr[ACHIEVABLE_7]
    sr = sr.div(sr.sum(axis=1).replace(0, 1), axis=0)

    common = coords.index.intersection(gt.index).intersection(
        pred.index).intersection(sr.index)
    return coords.loc[common], gt.loc[common], pred.loc[common], sr.loc[common]


def draw_pie_at(ax, x, y, proportions, colors, radius=30):
    mask = proportions > 0.01
    if mask.sum() == 0:
        mask[proportions.argmax()] = True
    props = proportions[mask]
    cols = [colors[i] for i in range(len(proportions)) if mask[i]]
    angles = props / props.sum() * 360.0
    start_angle = 90
    for angle, color in zip(angles, cols):
        theta1 = start_angle
        theta2 = start_angle + angle
        wedge = mpatches.Wedge((x, y), radius, theta1, theta2,
                               facecolor=color, edgecolor='none', linewidth=0)
        ax.add_patch(wedge)
        start_angle = theta2


def compute_pie_radius(x, y):
    from scipy.spatial import distance
    sample_idx = np.random.choice(len(x), min(200, len(x)), replace=False)
    sample_coords = np.column_stack([x[sample_idx], y[sample_idx]])
    nn_dists = []
    for i in range(len(sample_coords)):
        d = distance.cdist([sample_coords[i]], sample_coords)[0]
        d[i] = np.inf
        nn_dists.append(d.min())
    return np.median(nn_dists) * 0.40


def draw_spatial_panel(ax, coords, proportions_df, title, radius, xlim, ylim):
    x = coords['x'].values
    y = coords['y'].values
    color_list = [CELL_COLORS[ct] for ct in ACHIEVABLE_7]
    prop_values = proportions_df[ACHIEVABLE_7].values
    for i in range(len(x)):
        draw_pie_at(ax, x[i], y[i], prop_values[i], color_list, radius=radius)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=14, fontweight='bold', fontfamily='sans-serif', pad=10)
    ax.axis('off')


# ── Step 6: Draw combined figure ──────────────────────────────────────────────

def draw_combined_figure(coords, gt_prop, cg_prop, sr_prop,
                         cg_metrics, sr_metrics, output_path):
    """Top: 3 spatial pie panels. Bottom: per-cell-type GEX Pearson r bars."""
    fig = plt.figure(figsize=(22, 18))
    gs = fig.add_gridspec(2, 3, height_ratios=[2.2, 1], hspace=0.2, wspace=0.15)

    # ── Top: Spatial pie charts ──
    ax_gt = fig.add_subplot(gs[0, 0])
    ax_cg = fig.add_subplot(gs[0, 1])
    ax_sr = fig.add_subplot(gs[0, 2])

    x = coords['x'].values
    y = coords['y'].values
    radius = compute_pie_radius(x, y)
    pad = radius * 3
    xlim = (x.min() - pad, x.max() + pad)
    ylim = (y.min() - pad, y.max() + pad)

    logger.info("Drawing spatial panels...")
    draw_spatial_panel(ax_gt, coords, gt_prop, "Ground Truth", radius, xlim, ylim)
    draw_spatial_panel(ax_cg, coords, cg_prop, "CITEgeist", radius, xlim, ylim)
    draw_spatial_panel(ax_sr, coords, sr_prop, "scResolve", radius, xlim, ylim)

    # Legend for spatial panels
    legend_handles = [
        mpatches.Patch(facecolor=CELL_COLORS[ct], edgecolor='#333333',
                       linewidth=0.5, label=ct)
        for ct in ACHIEVABLE_7
    ]
    fig.legend(handles=legend_handles, loc='upper center', ncol=7,
               fontsize=11, frameon=True, edgecolor='#cccccc',
               bbox_to_anchor=(0.5, 0.505))

    # ── Bottom: GEX metric bar charts ──
    gs_bottom = gs[1, :].subgridspec(1, 2, wspace=0.35)

    cell_types_with_data = [ct for ct in ACHIEVABLE_7
                            if ct in cg_metrics or ct in sr_metrics]
    x_pos = np.arange(len(cell_types_with_data))
    width = 0.35

    # Left panel: Per-spot gene profile Pearson r (on jointly-active spots)
    ax_pct = fig.add_subplot(gs_bottom[0])
    cg_vals = [cg_metrics.get(ct, {}).get('spot_pearson_r', 0) for ct in cell_types_with_data]
    sr_vals = [sr_metrics.get(ct, {}).get('spot_pearson_r', 0) for ct in cell_types_with_data]

    bars1 = ax_pct.bar(x_pos - width/2, cg_vals, width, label='CITEgeist',
                       color='#4A90D9', edgecolor='#333333', linewidth=0.5)
    bars2 = ax_pct.bar(x_pos + width/2, sr_vals, width, label='scResolve',
                       color='#E85D75', edgecolor='#333333', linewidth=0.5)
    ax_pct.set_xticks(x_pos)
    ax_pct.set_xticklabels([ct.replace(' ', '\n') for ct in cell_types_with_data],
                           fontsize=9, fontfamily='sans-serif')
    ax_pct.set_ylabel('Mean Pearson r', fontsize=11, fontfamily='sans-serif')
    ax_pct.set_title('Gene Profile Correlation\n(per-spot, active spots only)',
                     fontsize=13, fontweight='bold', fontfamily='sans-serif')
    ax_pct.legend(fontsize=10, frameon=True)
    ax_pct.axhline(y=0, color='#999999', linewidth=0.5, linestyle='--')
    ax_pct.spines['top'].set_visible(False)
    ax_pct.spines['right'].set_visible(False)

    # Right panel: Coverage (fraction of GT-active spots that method predicts active)
    ax_cov = fig.add_subplot(gs_bottom[1])
    cg_cov = [cg_metrics.get(ct, {}).get('coverage', 0) for ct in cell_types_with_data]
    sr_cov = [sr_metrics.get(ct, {}).get('coverage', 0) for ct in cell_types_with_data]

    bars3 = ax_cov.bar(x_pos - width/2, cg_cov, width, label='CITEgeist',
                       color='#4A90D9', edgecolor='#333333', linewidth=0.5)
    bars4 = ax_cov.bar(x_pos + width/2, sr_cov, width, label='scResolve',
                       color='#E85D75', edgecolor='#333333', linewidth=0.5)
    ax_cov.set_xticks(x_pos)
    ax_cov.set_xticklabels([ct.replace(' ', '\n') for ct in cell_types_with_data],
                           fontsize=9, fontfamily='sans-serif')
    ax_cov.set_ylabel('Coverage', fontsize=11, fontfamily='sans-serif')
    ax_cov.set_title('Cell Type Coverage\n(fraction of GT spots with prediction)',
                     fontsize=13, fontweight='bold', fontfamily='sans-serif')
    ax_cov.legend(fontsize=10, frameon=True)
    ax_cov.set_ylim(0, 1.05)
    ax_cov.spines['top'].set_visible(False)
    ax_cov.spines['right'].set_visible(False)

    fig.suptitle("Region 0: CITEgeist vs scResolve — GEX Deconvolution (Fair Comparison)",
                 fontsize=18, fontweight='bold', fontfamily='sans-serif', y=0.98)

    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    logger.info(f"Saved: {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    logger.info("=" * 60)
    logger.info("Fair GEX Comparison: CITEgeist vs scResolve (Region 0)")
    logger.info("=" * 60)

    # Load GEX data
    logger.info("\n--- Loading GT GEX ---")
    gt_gex = load_gt_gex()

    logger.info("\n--- Loading CITEgeist GEX ---")
    citegeist_gex = load_citegeist_gex()

    logger.info("\n--- Computing scResolve Hungarian-matched GEX ---")
    scresolve_gex = compute_scresolve_hungarian_gex()

    # Compute fair metrics (all spots)
    logger.info("\n--- Computing fair metrics ---")
    logger.info("CITEgeist:")
    cg_metrics = compute_fair_gex_metrics(gt_gex, citegeist_gex, "CITEgeist")
    logger.info("scResolve:")
    sr_metrics = compute_fair_gex_metrics(gt_gex, scresolve_gex, "scResolve")

    logger.info(f"\nCITEgeist  mean spot Pearson r: {cg_metrics['_overall']['mean_spot_pearson_r']:.4f}")
    logger.info(f"scResolve  mean spot Pearson r: {sr_metrics['_overall']['mean_spot_pearson_r']:.4f}")
    logger.info(f"CITEgeist  mean coverage: {cg_metrics['_overall']['mean_coverage']:.2%}")
    logger.info(f"scResolve  mean coverage: {sr_metrics['_overall']['mean_coverage']:.2%}")

    # Save metrics
    results_path = BASE_DIR / "evaluation" / "results_gex_comparison_fair.json"
    results_clean = {}
    for name, metrics in [("CITEgeist", cg_metrics), ("scResolve", sr_metrics)]:
        results_clean[name] = {}
        for ct, m in metrics.items():
            results_clean[name][ct] = {
                k: float(v) if isinstance(v, (np.floating, float)) else v
                for k, v in m.items()
            }
    with open(results_path, 'w') as f:
        json.dump(results_clean, f, indent=2)
    logger.info(f"Saved metrics: {results_path}")

    # Load proportion data for spatial panels
    logger.info("\n--- Loading proportion data for spatial panels ---")
    coords, gt_prop, cg_prop, sr_prop = load_proportions()

    # Draw figure
    logger.info("\n--- Drawing figure ---")
    fig_path = OUTPUT_DIR / "benchmark_gex_comparison.png"
    draw_combined_figure(coords, gt_prop, cg_prop, sr_prop,
                         cg_metrics, sr_metrics, str(fig_path))

    logger.info("\nDone!")
