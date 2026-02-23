#!/usr/bin/env python
"""
Analyze what RNA cell types the protein-Unknown cells belong to.

Hypothesis: Unknown cells (low protein) still have RNA signal that
RNA-based methods can predict, but CITEgeist can't.
"""

import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

XENIUM_DIR = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')
PSEUDOVISIUM_DIR = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium')

# RNA cluster mapping (same as rna_cell_types.py)
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

CELL_TYPE_ORDER = ["B cells", "T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]


def get_threshold(expr_df, marker, percentile):
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0


def classify_cells_protein(expr_df):
    """Protein-based classification."""
    cell_types = pd.Series('Unknown', index=expr_df.index)

    CD20_thresh = get_threshold(expr_df, 'CD20', 25)
    CD3E_thresh = get_threshold(expr_df, 'CD3E', 50)
    CD68_thresh = get_threshold(expr_df, 'CD68', 50)
    CD31_thresh = get_threshold(expr_df, 'CD31', 50)
    PanCK_thresh = get_threshold(expr_df, 'PanCK', 25)
    ECad_thresh = get_threshold(expr_df, 'E-Cadherin', 90)
    alphaSMA_thresh = get_threshold(expr_df, 'alphaSMA', 75)

    # Hierarchical gating
    b_cells = expr_df['CD20'] > CD20_thresh
    cell_types[b_cells] = 'B cells'

    t_cells = (expr_df['CD3E'] > CD3E_thresh) & (cell_types == 'Unknown')
    cell_types[t_cells] = 'T cells'

    cd3e_neg = expr_df['CD3E'] < CD3E_thresh
    macrophages = (expr_df['CD68'] > CD68_thresh) & cd3e_neg & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'

    cd68_neg = expr_df['CD68'] < CD68_thresh
    cd31_pos = expr_df['CD31'] > CD31_thresh
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'

    panck_pos = expr_df['PanCK'] > PanCK_thresh
    ecad_high = expr_df['E-Cadherin'] > ECad_thresh
    epithelial = (panck_pos | ecad_high) & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'

    asma_high = expr_df['alphaSMA'] > alphaSMA_thresh
    fibroblasts = asma_high & ~cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[fibroblasts] = 'Fibroblasts'

    return cell_types


def main():
    print("=" * 70)
    print("Analyzing RNA Distribution of Protein-Unknown Cells")
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

    # Get protein classifications
    print("\n3. Classifying cells by protein...")
    protein_types = classify_cells_protein(prot_df)

    # Align
    common_cells = protein_types.index.intersection(rna_types.index)
    protein_aligned = protein_types.loc[common_cells]
    rna_aligned = rna_types.loc[common_cells]

    # Separate Unknown vs classified
    unknown_mask = protein_aligned == 'Unknown'
    n_unknown = unknown_mask.sum()
    n_classified = (~unknown_mask).sum()

    print(f"\n4. Cell counts:")
    print(f"   Total cells: {len(common_cells):,}")
    print(f"   Protein-classified: {n_classified:,} ({100*n_classified/len(common_cells):.1f}%)")
    print(f"   Protein-Unknown: {n_unknown:,} ({100*n_unknown/len(common_cells):.1f}%)")

    # What RNA types are the Unknown cells?
    print("\n" + "=" * 70)
    print("5. RNA cell type distribution of PROTEIN-UNKNOWN cells")
    print("=" * 70)

    unknown_rna = rna_aligned[unknown_mask]
    unknown_rna_counts = unknown_rna.value_counts()

    print(f"\n   {'RNA Cell Type':<20} {'Count':>10} {'% of Unknown':>15} {'% of RNA Type':>15}")
    print("   " + "-" * 60)

    for ct in CELL_TYPE_ORDER:
        if ct in unknown_rna_counts.index:
            count = unknown_rna_counts[ct]
            pct_unknown = 100 * count / n_unknown
            # What fraction of this RNA type is Unknown in protein?
            total_rna_type = (rna_aligned == ct).sum()
            pct_of_type = 100 * count / total_rna_type
            print(f"   {ct:<20} {count:>10,} {pct_unknown:>14.1f}% {pct_of_type:>14.1f}%")

    # Key question: which RNA cell types are MOSTLY Unknown in protein?
    print("\n" + "=" * 70)
    print("6. What fraction of each RNA type is protein-Unknown?")
    print("   (High % = RNA method can see them, CITEgeist can't)")
    print("=" * 70)

    print(f"\n   {'RNA Cell Type':<20} {'Total':>10} {'Unknown':>10} {'Classified':>10} {'% Unknown':>10}")
    print("   " + "-" * 60)

    type_unknown_pct = {}
    for ct in CELL_TYPE_ORDER:
        ct_mask = rna_aligned == ct
        total = ct_mask.sum()
        unknown_in_type = (unknown_mask & ct_mask).sum()
        classified = total - unknown_in_type
        pct_unknown = 100 * unknown_in_type / total if total > 0 else 0
        type_unknown_pct[ct] = pct_unknown
        print(f"   {ct:<20} {total:>10,} {unknown_in_type:>10,} {classified:>10,} {pct_unknown:>9.1f}%")

    # Compare with CITEgeist performance gap
    print("\n" + "=" * 70)
    print("7. Correlation with CITEgeist performance gap")
    print("=" * 70)

    # From the benchmark results (protein GT r - RNA GT r for CITEgeist_Continuous)
    # These are approximate from the output
    protein_gt_r = {
        "B cells": 0.70,  # CITEgeist does well on protein GT
        "T cells": 0.65,
        "Macrophages": 0.75,
        "Endothelial": 0.65,
        "Epithelial": 0.60,
        "Fibroblasts": 0.50,
    }
    rna_gt_r = {
        "B cells": 0.64,
        "T cells": 0.69,
        "Macrophages": 0.66,
        "Endothelial": 0.50,
        "Epithelial": 0.54,
        "Fibroblasts": 0.43,
    }

    print(f"\n   {'Cell Type':<20} {'% Unknown':>10} {'Prot GT r':>10} {'RNA GT r':>10} {'Δ r':>10}")
    print("   " + "-" * 60)

    for ct in CELL_TYPE_ORDER:
        pct_unk = type_unknown_pct.get(ct, 0)
        pr = protein_gt_r.get(ct, 0)
        rr = rna_gt_r.get(ct, 0)
        delta = rr - pr
        print(f"   {ct:<20} {pct_unk:>9.1f}% {pr:>10.2f} {rr:>10.2f} {delta:>+10.2f}")

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("""
If high '% Unknown' correlates with negative Δr (CITEgeist worse on RNA GT),
this confirms the hypothesis: cells that are protein-low but RNA-classifiable
create a disadvantage for CITEgeist when evaluated against RNA GT.
""")


if __name__ == '__main__':
    main()
