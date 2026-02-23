#!/usr/bin/env python
"""
Diagnose CD4+ T cell discrepancy between protein gating and RNA clustering.

Questions to answer:
1. Are protein-gated CD4+ T cells actually CD4+ by protein expression?
2. What does CD4 GENE expression look like in these cells?
3. Is our gating threshold wrong?
4. What does the RNA cluster 9 (Proliferating T) actually look like?
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
    print("=" * 60)
    print("CD4+ T Cell Discrepancy Diagnosis")
    print("=" * 60)

    # Load protein data
    print("\n1. Loading protein expression...")
    adata = sc.read_10x_h5(XENIUM_DIR / 'cell_feature_matrix.h5', gex_only=False)

    protein_mask = adata.var['feature_types'] == 'Protein Expression'
    adata_protein = adata[:, protein_mask].copy()

    X_prot = adata_protein.X.toarray()
    proteins = list(adata_protein.var_names)
    cell_ids = [str(x) for x in adata_protein.obs_names]

    prot_df = pd.DataFrame(X_prot, index=cell_ids, columns=proteins)
    print(f"   Loaded {len(prot_df)} cells x {len(proteins)} proteins")

    # Load RNA data
    print("\n2. Loading RNA expression...")
    rna_mask = adata.var['feature_types'] == 'Gene Expression'
    adata_rna = adata[:, rna_mask].copy()

    X_rna = adata_rna.X.toarray()
    genes = list(adata_rna.var_names)

    rna_df = pd.DataFrame(X_rna, index=cell_ids, columns=genes)
    print(f"   Loaded {len(rna_df)} cells x {len(genes)} genes")

    # Load RNA clusters
    print("\n3. Loading RNA k-means clusters...")
    analysis_tar = XENIUM_DIR / 'analysis.tar.gz'
    with tarfile.open(analysis_tar, 'r:gz') as tar:
        f = tar.extractfile('analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv')
        clusters_df = pd.read_csv(f)

    clusters_df['cell_id'] = clusters_df['Barcode'].astype(str)
    clusters_df = clusters_df.set_index('cell_id')

    # Check thresholds
    print("\n4. Checking protein thresholds...")
    CD3E_thresh = get_threshold(prot_df, 'CD3E', 50)
    CD4_thresh = get_threshold(prot_df, 'CD4', 50)
    CD8A_thresh = get_threshold(prot_df, 'CD8A', 50)

    print(f"   CD3E threshold (p50): {CD3E_thresh:.1f}")
    print(f"   CD4 threshold (p50): {CD4_thresh:.1f}")
    print(f"   CD8A threshold (p50): {CD8A_thresh:.1f}")

    # Apply protein gating for CD4+ T cells
    print("\n5. Applying protein gating for T cells...")
    cd3e_pos = prot_df['CD3E'] > CD3E_thresh
    cd4_pos = prot_df['CD4'] > CD4_thresh
    cd8a_pos = prot_df['CD8A'] > CD8A_thresh
    cd8a_neg = prot_df['CD8A'] < CD8A_thresh

    # Our gating: CD4+ T cells = CD3E+ CD4+ CD8A-
    protein_cd4_tcells = cd3e_pos & cd4_pos & cd8a_neg
    # CD8+ T cells = CD3E+ CD8A+
    protein_cd8_tcells = cd3e_pos & cd8a_pos

    print(f"   Protein-gated CD4+ T cells: {protein_cd4_tcells.sum()}")
    print(f"   Protein-gated CD8+ T cells: {protein_cd8_tcells.sum()}")

    # Check actual expression in protein-gated CD4+ T cells
    print("\n6. Protein expression in protein-gated CD4+ T cells:")
    cd4_cells = prot_df[protein_cd4_tcells]
    print(f"   CD3E: mean={cd4_cells['CD3E'].mean():.1f}, median={cd4_cells['CD3E'].median():.1f}")
    print(f"   CD4:  mean={cd4_cells['CD4'].mean():.1f}, median={cd4_cells['CD4'].median():.1f}")
    print(f"   CD8A: mean={cd4_cells['CD8A'].mean():.1f}, median={cd4_cells['CD8A'].median():.1f}")

    # Check RNA cluster distribution of protein-gated CD4+ T cells
    print("\n7. RNA cluster distribution of protein-gated CD4+ T cells:")
    cd4_cell_ids = prot_df.index[protein_cd4_tcells]
    cd4_clusters = clusters_df.loc[cd4_cell_ids.intersection(clusters_df.index), 'Cluster']
    cluster_dist = cd4_clusters.value_counts().sort_index()

    RNA_CLUSTER_MAP = {
        1: "CD8+ T cells",
        2: "Macrophages",
        3: "Mixed Immune",
        4: "Epithelial",
        5: "Fibroblasts",
        6: "Stromal",
        7: "Endothelial",
        8: "B cells",
        9: "Proliferating T",
        10: "Vascular Stromal",
    }

    for cluster, count in cluster_dist.items():
        pct = 100 * count / len(cd4_clusters)
        print(f"   Cluster {cluster} ({RNA_CLUSTER_MAP[cluster]}): {count} ({pct:.1f}%)")

    # Check CD4 GENE expression in these cells
    print("\n8. CD4 GENE expression in protein-gated CD4+ T cells:")
    if 'CD4' in rna_df.columns:
        cd4_gene = rna_df.loc[cd4_cell_ids.intersection(rna_df.index), 'CD4']
        print(f"   CD4 gene: mean={cd4_gene.mean():.2f}, median={cd4_gene.median():.2f}")
        print(f"   CD4 gene > 0: {(cd4_gene > 0).sum()} ({100*(cd4_gene > 0).mean():.1f}%)")
    else:
        print("   CD4 gene not in panel!")
        # Check what T cell genes are available
        t_genes = [g for g in genes if 'CD3' in g or 'CD4' in g or 'CD8' in g]
        print(f"   Available T cell genes: {t_genes}")

    # Check CD8A GENE expression in protein-gated CD4+ T cells
    print("\n9. CD8A GENE expression in protein-gated CD4+ T cells:")
    if 'CD8A' in rna_df.columns:
        cd8a_gene = rna_df.loc[cd4_cell_ids.intersection(rna_df.index), 'CD8A']
        print(f"   CD8A gene: mean={cd8a_gene.mean():.2f}, median={cd8a_gene.median():.2f}")
        print(f"   CD8A gene > 0: {(cd8a_gene > 0).sum()} ({100*(cd8a_gene > 0).mean():.1f}%)")

    # Check CD3E GENE expression
    print("\n10. CD3E GENE expression in protein-gated CD4+ T cells:")
    if 'CD3E' in rna_df.columns:
        cd3e_gene = rna_df.loc[cd4_cell_ids.intersection(rna_df.index), 'CD3E']
        print(f"   CD3E gene: mean={cd3e_gene.mean():.2f}, median={cd3e_gene.median():.2f}")
        print(f"   CD3E gene > 0: {(cd3e_gene > 0).sum()} ({100*(cd3e_gene > 0).mean():.1f}%)")

    # Compare to RNA cluster 9 (Proliferating T) - our "CD4+ T cells" RNA mapping
    print("\n11. RNA Cluster 9 (Proliferating T) - what we map to CD4+ T cells:")
    cluster9_ids = clusters_df[clusters_df['Cluster'] == 9].index
    cluster9_prot = prot_df.loc[cluster9_ids.intersection(prot_df.index)]

    print(f"   Total cells: {len(cluster9_prot)}")
    print(f"   CD3E protein: mean={cluster9_prot['CD3E'].mean():.1f}")
    print(f"   CD4 protein: mean={cluster9_prot['CD4'].mean():.1f}")
    print(f"   CD8A protein: mean={cluster9_prot['CD8A'].mean():.1f}")

    # How many cluster 9 cells pass protein CD4+ T cell gating?
    cluster9_cd4_gated = (
        (cluster9_prot['CD3E'] > CD3E_thresh) &
        (cluster9_prot['CD4'] > CD4_thresh) &
        (cluster9_prot['CD8A'] < CD8A_thresh)
    )
    print(f"   Pass CD4+ T cell gating: {cluster9_cd4_gated.sum()} ({100*cluster9_cd4_gated.mean():.1f}%)")

    # Compare to RNA cluster 1 (CD8+ T cells)
    print("\n12. RNA Cluster 1 (CD8+ T cells) - where protein CD4+ T cells actually go:")
    cluster1_ids = clusters_df[clusters_df['Cluster'] == 1].index
    cluster1_prot = prot_df.loc[cluster1_ids.intersection(prot_df.index)]

    print(f"   Total cells: {len(cluster1_prot)}")
    print(f"   CD3E protein: mean={cluster1_prot['CD3E'].mean():.1f}")
    print(f"   CD4 protein: mean={cluster1_prot['CD4'].mean():.1f}")
    print(f"   CD8A protein: mean={cluster1_prot['CD8A'].mean():.1f}")

    # How many cluster 1 cells pass protein CD4+ T cell gating?
    cluster1_cd4_gated = (
        (cluster1_prot['CD3E'] > CD3E_thresh) &
        (cluster1_prot['CD4'] > CD4_thresh) &
        (cluster1_prot['CD8A'] < CD8A_thresh)
    )
    print(f"   Pass CD4+ T cell gating: {cluster1_cd4_gated.sum()} ({100*cluster1_cd4_gated.mean():.1f}%)")

    # CRITICAL CHECK: Is the CD4 threshold too low?
    print("\n13. CRITICAL: CD4 protein distribution analysis")
    print(f"   All cells CD4 > 0: {(prot_df['CD4'] > 0).sum()} ({100*(prot_df['CD4'] > 0).mean():.1f}%)")
    print(f"   CD4 percentiles (non-zero):")
    cd4_nonzero = prot_df['CD4'][prot_df['CD4'] > 0]
    for p in [10, 25, 50, 75, 90, 95, 99]:
        print(f"      p{p}: {np.percentile(cd4_nonzero, p):.1f}")

    # CRITICAL CHECK: Are we using the wrong marker?
    print("\n14. Available markers in protein panel:")
    for p in proteins:
        if 'CD' in p or 'T cell' in p.lower():
            print(f"   {p}")

    print("\n" + "=" * 60)
    print("DIAGNOSIS SUMMARY")
    print("=" * 60)


if __name__ == '__main__':
    main()
