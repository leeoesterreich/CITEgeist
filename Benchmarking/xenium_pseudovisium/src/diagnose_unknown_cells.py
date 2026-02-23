#!/usr/bin/env python
"""
Diagnose protein expression in "Unknown" cells to assess gating stringency.

Questions:
1. Are Unknown cells truly negative or borderline?
2. What would happen with looser thresholds?
3. Which markers are these cells closest to?
"""

import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

XENIUM_DIR = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')
OUTPUT_DIR = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/analysis')


def get_threshold(expr_df, marker, percentile=50):
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0


def main():
    print("=" * 70)
    print("Unknown Cell Protein Expression Diagnosis")
    print("=" * 70)

    # Load protein data
    print("\n1. Loading data...")
    adata = sc.read_10x_h5(XENIUM_DIR / 'cell_feature_matrix.h5', gex_only=False)
    protein_mask = adata.var['feature_types'] == 'Protein Expression'
    adata_protein = adata[:, protein_mask].copy()

    X = adata_protein.X.toarray()
    proteins = list(adata_protein.var_names)
    cell_ids = [str(x) for x in adata_protein.obs_names]
    prot_df = pd.DataFrame(X, index=cell_ids, columns=proteins)
    print(f"   Loaded {len(prot_df)} cells x {len(proteins)} proteins")

    # Current thresholds
    print("\n2. Current gating thresholds (percentile of non-zero):")
    thresholds = {
        'CD20': get_threshold(prot_df, 'CD20', 25),
        'CD3E': get_threshold(prot_df, 'CD3E', 50),
        'CD68': get_threshold(prot_df, 'CD68', 50),
        'CD31': get_threshold(prot_df, 'CD31', 50),
        'PanCK': get_threshold(prot_df, 'PanCK', 25),
        'E-Cadherin': get_threshold(prot_df, 'E-Cadherin', 90),
        'alphaSMA': get_threshold(prot_df, 'alphaSMA', 75),
    }
    for marker, thresh in thresholds.items():
        print(f"   {marker}: {thresh:.1f}")

    # Apply gating to identify Unknown cells
    print("\n3. Identifying Unknown cells with current gating...")
    cell_types = pd.Series('Unknown', index=prot_df.index)

    # B cells
    b_cells = prot_df['CD20'] > thresholds['CD20']
    cell_types[b_cells] = 'B cells'

    # T cells
    t_cells = (prot_df['CD3E'] > thresholds['CD3E']) & (cell_types == 'Unknown')
    cell_types[t_cells] = 'T cells'

    # Macrophages
    macrophages = (prot_df['CD68'] > thresholds['CD68']) & (prot_df['CD3E'] < thresholds['CD3E']) & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'

    # Endothelial
    endothelial = (prot_df['CD31'] > thresholds['CD31']) & (prot_df['CD68'] < thresholds['CD68']) & (prot_df['CD3E'] < thresholds['CD3E']) & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'

    # Epithelial
    epithelial = ((prot_df['PanCK'] > thresholds['PanCK']) | (prot_df['E-Cadherin'] > thresholds['E-Cadherin'])) & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'

    # Fibroblasts
    cd31_pos = prot_df['CD31'] > thresholds['CD31']
    cd68_neg = prot_df['CD68'] < thresholds['CD68']
    cd3e_neg = prot_df['CD3E'] < thresholds['CD3E']
    fibroblasts = (prot_df['alphaSMA'] > thresholds['alphaSMA']) & ~cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[fibroblasts] = 'Fibroblasts'

    unknown_mask = cell_types == 'Unknown'
    n_unknown = unknown_mask.sum()
    print(f"   Unknown cells: {n_unknown} ({100*n_unknown/len(cell_types):.1f}%)")

    # Analyze Unknown cells
    unknown_df = prot_df[unknown_mask]

    print("\n4. Protein expression in Unknown cells:")
    print("   " + "-" * 60)
    print(f"   {'Marker':<15} {'Mean':>10} {'Median':>10} {'% > 0':>10} {'Threshold':>10}")
    print("   " + "-" * 60)

    key_markers = ['CD20', 'CD3E', 'CD68', 'CD31', 'PanCK', 'E-Cadherin', 'alphaSMA', 'Vimentin', 'HLA-DR']
    for marker in key_markers:
        if marker in unknown_df.columns:
            vals = unknown_df[marker]
            mean = vals.mean()
            median = vals.median()
            pct_pos = 100 * (vals > 0).mean()
            thresh = thresholds.get(marker, 'N/A')
            thresh_str = f"{thresh:.1f}" if isinstance(thresh, float) else thresh
            print(f"   {marker:<15} {mean:>10.1f} {median:>10.1f} {pct_pos:>9.1f}% {thresh_str:>10}")

    # Check how close Unknown cells are to each threshold
    print("\n5. How close are Unknown cells to thresholds?")
    print("   (% of Unknown cells above X% of threshold)")
    print("   " + "-" * 50)

    for marker in ['CD20', 'CD3E', 'CD68', 'CD31', 'PanCK', 'alphaSMA']:
        thresh = thresholds[marker]
        vals = unknown_df[marker]
        pct_above_50 = 100 * (vals > 0.5 * thresh).mean()
        pct_above_75 = 100 * (vals > 0.75 * thresh).mean()
        pct_above_90 = 100 * (vals > 0.90 * thresh).mean()
        print(f"   {marker:<12}: >50%={pct_above_50:>5.1f}%  >75%={pct_above_75:>5.1f}%  >90%={pct_above_90:>5.1f}%")

    # What if we used lower thresholds?
    print("\n6. Impact of lowering thresholds (to p25 instead of p50):")
    lower_thresholds = {
        'CD20': get_threshold(prot_df, 'CD20', 10),  # Very low for sparse
        'CD3E': get_threshold(prot_df, 'CD3E', 25),
        'CD68': get_threshold(prot_df, 'CD68', 25),
        'CD31': get_threshold(prot_df, 'CD31', 25),
        'PanCK': get_threshold(prot_df, 'PanCK', 10),
        'E-Cadherin': get_threshold(prot_df, 'E-Cadherin', 75),
        'alphaSMA': get_threshold(prot_df, 'alphaSMA', 50),
    }

    print(f"   {'Marker':<15} {'Current':>10} {'Lower':>10} {'Unknowns rescued':>20}")
    print("   " + "-" * 55)

    for marker in ['CD3E', 'CD68', 'CD31', 'PanCK', 'alphaSMA']:
        current = thresholds[marker]
        lower = lower_thresholds[marker]
        # Count unknowns that would be rescued
        rescued = ((unknown_df[marker] > lower) & (unknown_df[marker] <= current)).sum()
        print(f"   {marker:<15} {current:>10.1f} {lower:>10.1f} {rescued:>15} ({100*rescued/n_unknown:.1f}%)")

    # Load RNA clusters to see what RNA calls these Unknown cells
    print("\n7. What RNA clustering calls Unknown cells:")
    analysis_tar = XENIUM_DIR / 'analysis.tar.gz'
    with tarfile.open(analysis_tar, 'r:gz') as tar:
        f = tar.extractfile('analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv')
        clusters_df = pd.read_csv(f)
    clusters_df['cell_id'] = clusters_df['Barcode'].astype(str)
    clusters_df = clusters_df.set_index('cell_id')

    RNA_CLUSTER_MAP = {
        1: "CD8+ T cells", 2: "Macrophages", 3: "Mixed Immune", 4: "Epithelial",
        5: "Fibroblasts", 6: "Stromal", 7: "Endothelial", 8: "B cells",
        9: "Proliferating T", 10: "Vascular Stromal",
    }

    unknown_ids = unknown_df.index.intersection(clusters_df.index)
    unknown_clusters = clusters_df.loc[unknown_ids, 'Cluster'].map(RNA_CLUSTER_MAP)
    cluster_dist = unknown_clusters.value_counts()

    print("   " + "-" * 40)
    for ct, count in cluster_dist.items():
        print(f"   {ct:<20}: {count:>8} ({100*count/len(unknown_clusters):.1f}%)")

    # For each RNA-classified Unknown cell type, what's their protein expression?
    print("\n8. Protein expression in Unknown cells BY RNA cluster:")
    print("   (Are they borderline for the expected marker?)")
    print("   " + "-" * 70)

    rna_to_expected_marker = {
        'Epithelial': 'PanCK',
        'Macrophages': 'CD68',
        'Fibroblasts': 'alphaSMA',
        'Stromal': 'Vimentin',
        'Mixed Immune': 'CD68',
        'CD8+ T cells': 'CD3E',
        'Endothelial': 'CD31',
    }

    for rna_type in ['Epithelial', 'Macrophages', 'Fibroblasts', 'CD8+ T cells', 'Endothelial']:
        mask = unknown_clusters == rna_type
        if mask.sum() == 0:
            continue
        subset_ids = unknown_ids[mask]
        subset = prot_df.loc[subset_ids]

        expected_marker = rna_to_expected_marker.get(rna_type)
        if expected_marker and expected_marker in subset.columns:
            vals = subset[expected_marker]
            thresh = thresholds.get(expected_marker, 0)
            mean = vals.mean()
            pct_above_thresh = 100 * (vals > thresh).mean()
            pct_above_half = 100 * (vals > 0.5 * thresh).mean()
            print(f"   RNA={rna_type:<15} | {expected_marker}: mean={mean:.1f}, {pct_above_half:.0f}% above 50% thresh, {pct_above_thresh:.0f}% above thresh")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
Key findings:
1. {n_unknown} cells ({100*n_unknown/len(cell_types):.0f}%) are Unknown by protein gating
2. RNA classifies them mostly as Epithelial/Macrophages/Fibroblasts
3. Check above for whether they're borderline or truly negative
""")


if __name__ == '__main__':
    main()
