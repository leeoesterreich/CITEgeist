#!/usr/bin/env python
"""
Analyze what ground truth looks like when requiring BOTH protein AND RNA agreement.
"""

import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import confusion_matrix

XENIUM_DIR = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')
PSEUDOVISIUM_DIR = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium')

# Simplified 6 cell types (T cells combined)
CELL_TYPE_ORDER = ["B cells", "T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]

RNA_CLUSTER_MAP = {
    1: "CD8+ T cells", 2: "Macrophages", 3: "Mixed Immune", 4: "Epithelial",
    5: "Fibroblasts", 6: "Stromal", 7: "Endothelial", 8: "B cells",
    9: "Proliferating T", 10: "Vascular Stromal",
}

RNA_SIMPLIFIED_MAP = {
    "CD8+ T cells": "T cells", "Macrophages": "Macrophages", "Mixed Immune": "Macrophages",
    "Epithelial": "Epithelial", "Fibroblasts": "Fibroblasts", "Stromal": "Fibroblasts",
    "Endothelial": "Endothelial", "B cells": "B cells", "Proliferating T": "T cells",
    "Vascular Stromal": "Endothelial",
}


def get_threshold(expr_df, marker, percentile):
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0


def classify_cells_protein(expr_df):
    """Protein-based classification with current thresholds."""
    cell_types = pd.Series('Unknown', index=expr_df.index)

    # Current thresholds (p50/p25/p75/p90)
    CD20_thresh = get_threshold(expr_df, 'CD20', 25)
    CD3E_thresh = get_threshold(expr_df, 'CD3E', 50)
    CD68_thresh = get_threshold(expr_df, 'CD68', 50)
    CD31_thresh = get_threshold(expr_df, 'CD31', 50)
    PanCK_thresh = get_threshold(expr_df, 'PanCK', 25)
    ECad_thresh = get_threshold(expr_df, 'E-Cadherin', 90)
    alphaSMA_thresh = get_threshold(expr_df, 'alphaSMA', 75)

    # B cells
    b_cells = expr_df['CD20'] > CD20_thresh
    cell_types[b_cells] = 'B cells'

    # T cells (combined CD4+ and CD8+)
    t_cells = (expr_df['CD3E'] > CD3E_thresh) & (cell_types == 'Unknown')
    cell_types[t_cells] = 'T cells'

    # Macrophages
    cd3e_neg = expr_df['CD3E'] < CD3E_thresh
    macrophages = (expr_df['CD68'] > CD68_thresh) & cd3e_neg & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'

    # Endothelial
    cd68_neg = expr_df['CD68'] < CD68_thresh
    cd31_pos = expr_df['CD31'] > CD31_thresh
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'

    # Epithelial
    panck_pos = expr_df['PanCK'] > PanCK_thresh
    ecad_high = expr_df['E-Cadherin'] > ECad_thresh
    epithelial = (panck_pos | ecad_high) & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'

    # Fibroblasts
    asma_high = expr_df['alphaSMA'] > alphaSMA_thresh
    fibroblasts = asma_high & ~cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[fibroblasts] = 'Fibroblasts'

    return cell_types


def main():
    print("=" * 70)
    print("Consensus Ground Truth Analysis (Protein + RNA Agreement)")
    print("=" * 70)

    # Load protein data
    print("\n1. Loading protein expression data...")
    adata = sc.read_10x_h5(XENIUM_DIR / 'cell_feature_matrix.h5', gex_only=False)
    protein_mask = adata.var['feature_types'] == 'Protein Expression'
    adata_protein = adata[:, protein_mask].copy()

    X = adata_protein.X.toarray()
    proteins = list(adata_protein.var_names)
    cell_ids = [str(x) for x in adata_protein.obs_names]
    prot_df = pd.DataFrame(X, index=cell_ids, columns=proteins)
    print(f"   Loaded {len(prot_df)} cells")

    # Load RNA clusters
    print("\n2. Loading RNA clusters...")
    analysis_tar = XENIUM_DIR / 'analysis.tar.gz'
    with tarfile.open(analysis_tar, 'r:gz') as tar:
        f = tar.extractfile('analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv')
        clusters_df = pd.read_csv(f)
    clusters_df['cell_id'] = clusters_df['Barcode'].astype(str)
    clusters_df = clusters_df.set_index('cell_id')
    rna_types = clusters_df['Cluster'].map(RNA_CLUSTER_MAP).map(RNA_SIMPLIFIED_MAP)

    # Get protein-based types
    print("\n3. Classifying cells by protein...")
    protein_types = classify_cells_protein(prot_df)

    # Find consensus cells
    print("\n4. Finding consensus cells (both modalities agree)...")
    common_cells = protein_types.index.intersection(rna_types.index)
    protein_aligned = protein_types.loc[common_cells]
    rna_aligned = rna_types.loc[common_cells]

    # Consensus = both agree AND not Unknown
    consensus_mask = (protein_aligned == rna_aligned) & (protein_aligned != 'Unknown')
    consensus_cells = common_cells[consensus_mask]
    consensus_types = protein_aligned.loc[consensus_cells]

    print(f"\n   Total cells: {len(common_cells):,}")
    print(f"   Protein Unknown: {(protein_aligned == 'Unknown').sum():,} ({100*(protein_aligned == 'Unknown').sum()/len(common_cells):.1f}%)")
    print(f"   Consensus cells: {len(consensus_cells):,} ({100*len(consensus_cells)/len(common_cells):.1f}%)")

    # Per-type breakdown
    print("\n5. Per-cell-type consensus counts:")
    print("-" * 60)
    print(f"{'Cell Type':<15} {'Protein':<10} {'RNA':<10} {'Consensus':<10} {'% Kept':<10}")
    print("-" * 60)

    for ct in CELL_TYPE_ORDER:
        n_protein = (protein_aligned == ct).sum()
        n_rna = (rna_aligned == ct).sum()
        n_consensus = (consensus_types == ct).sum()
        pct_kept = 100 * n_consensus / n_protein if n_protein > 0 else 0
        print(f"{ct:<15} {n_protein:<10,} {n_rna:<10,} {n_consensus:<10,} {pct_kept:<10.1f}%")

    # Load spatial coordinates to compute spot-level proportions
    print("\n6. Computing spot-level consensus GT proportions...")
    cells_tar = XENIUM_DIR / 'cells.csv.gz'
    coords_df = pd.read_csv(cells_tar, compression='gzip')
    coords_df['cell_id'] = coords_df['cell_id'].astype(str)
    coords_df = coords_df.set_index('cell_id')

    # Load pseudo-Visium spots
    spots_df = pd.read_csv(PSEUDOVISIUM_DIR / 'data' / 'pseudo_visium_spots.csv')
    spot_barcodes = spots_df['barcode'].values
    spot_coords = spots_df[['x', 'y']].values

    # Assign consensus cells to spots (within 27.5um radius)
    SPOT_RADIUS = 27.5
    from scipy.spatial import KDTree

    # Get coordinates for consensus cells only
    consensus_coords = coords_df.loc[coords_df.index.intersection(consensus_cells)]
    cell_xy = consensus_coords[['x_centroid', 'y_centroid']].values
    cell_ids_with_coords = consensus_coords.index.values

    spot_tree = KDTree(spot_coords)
    cell_tree = KDTree(cell_xy)

    # Build spot-to-cells mapping
    spot_cells = {barcode: [] for barcode in spot_barcodes}
    indices = spot_tree.query_ball_tree(cell_tree, r=SPOT_RADIUS)

    for spot_idx, cell_indices in enumerate(indices):
        barcode = spot_barcodes[spot_idx]
        for cell_idx in cell_indices:
            cell_id = cell_ids_with_coords[cell_idx]
            spot_cells[barcode].append(cell_id)

    # Compute consensus proportions
    consensus_gt = pd.DataFrame(0.0, index=spot_barcodes, columns=CELL_TYPE_ORDER)

    for barcode, cell_list in spot_cells.items():
        if len(cell_list) > 0:
            types = consensus_types.loc[cell_list]
            type_counts = types.value_counts()
            total = len(cell_list)
            for ct in CELL_TYPE_ORDER:
                if ct in type_counts.index:
                    consensus_gt.loc[barcode, ct] = type_counts[ct] / total

    # Summary statistics
    print("\n7. Consensus GT summary:")
    print("-" * 60)

    # Cells per spot
    cells_per_spot = [len(v) for v in spot_cells.values()]
    print(f"   Consensus cells per spot: mean={np.mean(cells_per_spot):.1f}, median={np.median(cells_per_spot):.0f}")

    # Proportion sums
    row_sums = consensus_gt.sum(axis=1)
    print(f"   Proportion row sums: mean={row_sums.mean():.3f}, std={row_sums.std():.3f}")
    print(f"   Spots with 0 cells: {(row_sums == 0).sum()}")

    # Per-type mean proportions
    print("\n   Mean proportions per cell type:")
    for ct in CELL_TYPE_ORDER:
        print(f"      {ct}: {consensus_gt[ct].mean():.3f}")

    # Compare with protein-only GT
    print("\n8. Loading protein-only GT for comparison...")
    protein_gt = pd.read_csv(PSEUDOVISIUM_DIR / 'data_protein_gt' / 'ground_truth_proportions.csv', index_col=0)

    # Align columns
    protein_gt_aligned = protein_gt.reindex(columns=CELL_TYPE_ORDER, fill_value=0)

    print("\n   Comparison: Protein-only vs Consensus GT")
    print("-" * 60)
    print(f"{'Cell Type':<15} {'Protein GT':<12} {'Consensus GT':<12} {'Diff':<10}")
    print("-" * 60)
    for ct in CELL_TYPE_ORDER:
        prot_mean = protein_gt_aligned[ct].mean()
        cons_mean = consensus_gt[ct].mean()
        diff = cons_mean - prot_mean
        print(f"{ct:<15} {prot_mean:<12.3f} {cons_mean:<12.3f} {diff:>+.3f}")

    # Save consensus GT
    output_dir = PSEUDOVISIUM_DIR / 'data_consensus_gt'
    output_dir.mkdir(exist_ok=True)
    consensus_gt.to_csv(output_dir / 'ground_truth_proportions.csv')
    print(f"\n9. Saved consensus GT to {output_dir / 'ground_truth_proportions.csv'}")

    # Key insight
    print("\n" + "=" * 70)
    print("KEY INSIGHTS")
    print("=" * 70)
    print(f"""
1. Consensus keeps {100*len(consensus_cells)/len(common_cells):.1f}% of cells (only those where both agree)
2. This eliminates the 'Unknown' problem entirely
3. Consensus GT will have higher confidence but fewer cells per spot

Implications for benchmarking:
- Consensus GT is the 'gold standard' (both modalities agree)
- Methods that predict well on consensus GT are truly accurate
- But: fewer cells per spot → noisier spot-level proportions
""")


if __name__ == '__main__':
    main()
