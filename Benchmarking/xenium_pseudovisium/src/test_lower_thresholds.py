#!/usr/bin/env python
"""
Test concordance with lower protein gating thresholds.
Compare current (p50) vs lower (p25) thresholds.
"""

import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import confusion_matrix

XENIUM_DIR = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')
PSEUDOVISIUM_DIR = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium')

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


def classify_cells(expr_df, threshold_percentiles):
    """Classify cells with given threshold percentiles."""
    cell_types = pd.Series('Unknown', index=expr_df.index)

    # Get thresholds
    CD20_thresh = get_threshold(expr_df, 'CD20', threshold_percentiles.get('CD20', 25))
    CD3E_thresh = get_threshold(expr_df, 'CD3E', threshold_percentiles.get('CD3E', 50))
    CD68_thresh = get_threshold(expr_df, 'CD68', threshold_percentiles.get('CD68', 50))
    CD31_thresh = get_threshold(expr_df, 'CD31', threshold_percentiles.get('CD31', 50))
    PanCK_thresh = get_threshold(expr_df, 'PanCK', threshold_percentiles.get('PanCK', 25))
    ECad_thresh = get_threshold(expr_df, 'E-Cadherin', threshold_percentiles.get('E-Cadherin', 90))
    alphaSMA_thresh = get_threshold(expr_df, 'alphaSMA', threshold_percentiles.get('alphaSMA', 75))

    # B cells
    b_cells = expr_df['CD20'] > CD20_thresh
    cell_types[b_cells] = 'B cells'

    # T cells
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


def calculate_concordance(protein_types, rna_types):
    """Calculate concordance metrics."""
    common_cells = protein_types.index.intersection(rna_types.index)
    protein_aligned = protein_types.loc[common_cells]
    rna_aligned = rna_types.loc[common_cells]

    all_types = CELL_TYPE_ORDER + ['Unknown']
    cm = confusion_matrix(protein_aligned, rna_aligned, labels=all_types)

    total = len(common_cells)
    diagonal = np.diag(cm).sum()
    overall = diagonal / total * 100

    # Per-type
    per_type = {}
    for i, ct in enumerate(all_types):
        if cm[i, :].sum() > 0:
            per_type[ct] = cm[i, i] / cm[i, :].sum() * 100

    # Unknown rate
    unknown_rate = (protein_aligned == 'Unknown').sum() / len(protein_aligned) * 100

    return overall, per_type, unknown_rate


def main():
    print("=" * 70)
    print("Testing Lower Protein Gating Thresholds")
    print("=" * 70)

    # Load data
    print("\n1. Loading data...")
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

    # Test different threshold settings
    threshold_configs = {
        'Current (p50/p25/p75/p90)': {
            'CD20': 25, 'CD3E': 50, 'CD68': 50, 'CD31': 50,
            'PanCK': 25, 'E-Cadherin': 90, 'alphaSMA': 75
        },
        'Lower (p25 for main markers)': {
            'CD20': 10, 'CD3E': 25, 'CD68': 25, 'CD31': 25,
            'PanCK': 10, 'E-Cadherin': 75, 'alphaSMA': 50
        },
        'Very Low (p10 for main markers)': {
            'CD20': 5, 'CD3E': 10, 'CD68': 10, 'CD31': 10,
            'PanCK': 5, 'E-Cadherin': 50, 'alphaSMA': 25
        },
    }

    print("\n3. Testing threshold configurations...")
    print("=" * 70)

    results = []
    for config_name, thresholds in threshold_configs.items():
        print(f"\n--- {config_name} ---")

        # Show actual threshold values
        print("   Thresholds:")
        for marker in ['CD3E', 'CD68', 'CD31', 'PanCK', 'alphaSMA']:
            thresh_val = get_threshold(prot_df, marker, thresholds[marker])
            print(f"      {marker} (p{thresholds[marker]}): {thresh_val:.1f}")

        # Classify
        protein_types = classify_cells(prot_df, thresholds)

        # Calculate concordance
        overall, per_type, unknown_rate = calculate_concordance(protein_types, rna_types)

        print(f"\n   Results:")
        print(f"      Unknown rate: {unknown_rate:.1f}%")
        print(f"      Overall concordance: {overall:.1f}%")
        print(f"      Per-type concordance:")
        for ct in CELL_TYPE_ORDER:
            if ct in per_type:
                print(f"         {ct}: {per_type[ct]:.1f}%")

        results.append({
            'config': config_name,
            'unknown_rate': unknown_rate,
            'overall_concordance': overall,
            'per_type': per_type
        })

    # Summary comparison
    print("\n" + "=" * 70)
    print("SUMMARY COMPARISON")
    print("=" * 70)
    print(f"\n{'Config':<35} {'Unknown%':>10} {'Concordance%':>12}")
    print("-" * 60)
    for r in results:
        print(f"{r['config']:<35} {r['unknown_rate']:>10.1f} {r['overall_concordance']:>12.1f}")

    # Trade-off analysis
    print("\n" + "=" * 70)
    print("TRADE-OFF ANALYSIS")
    print("=" * 70)

    current = results[0]
    lower = results[1]

    print(f"\nLowering thresholds from p50 to p25:")
    print(f"   Unknown rate: {current['unknown_rate']:.1f}% → {lower['unknown_rate']:.1f}% (Δ = {lower['unknown_rate'] - current['unknown_rate']:.1f}%)")
    print(f"   Concordance:  {current['overall_concordance']:.1f}% → {lower['overall_concordance']:.1f}% (Δ = {lower['overall_concordance'] - current['overall_concordance']:.1f}%)")

    print("\n   Per-type changes:")
    for ct in CELL_TYPE_ORDER:
        curr = current['per_type'].get(ct, 0)
        low = lower['per_type'].get(ct, 0)
        delta = low - curr
        direction = "↑" if delta > 0 else "↓" if delta < 0 else "="
        print(f"      {ct}: {curr:.1f}% → {low:.1f}% ({direction} {abs(delta):.1f}%)")


if __name__ == '__main__':
    main()
