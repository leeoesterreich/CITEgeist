#!/usr/bin/env python
"""
08_spatial_comparison.py

Question: Can the secretory program be demonstrated in spatial data,
          and does the patient resemble MCF7 or T47D?
Inputs:   Raw spatial data (h5), GSE89888 bulk RNA-seq
Outputs:  Comparison metrics and conclusions
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import stats
from scipy.spatial.distance import euclidean, cosine

# Paths - Note: squidpy.read.visium expects the parent directory (it adds /outs internally)
SPATIAL_DATA = Path("/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/HCC22-088-P4-S2_1i_rep")
CITEGEIST_OUTPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk/citegeist_output")
PIPELINE_DIR = Path(__file__).parent.parent
OUTPUT_DIR = PIPELINE_DIR / "outputs"
TABLES_DIR = OUTPUT_DIR / "tables"

# Cancer cell proportion threshold for filtering
CANCER_CELL_THRESHOLD = 0.3  # Spots with >30% cancer cells

# Secretory/chaperone genes of interest
SECRETORY_GENES = [
    'HSP90B1',  # Key ER-resident chaperone (the main target)
    'HSPA5',    # BiP/GRP78
    'CALR',     # Calreticulin
    'CANX',     # Calnexin
    'PDIA4',    # Protein disulfide isomerase
    'PDIA6',    # Protein disulfide isomerase
    'SEC61A1',  # Translocon
    'SEC61B',   # Translocon
    'XBP1',     # UPR master regulator
    'ATF6',     # UPR sensor
    'ERO1A',    # Oxidoreductin
    'MDK',      # Midkine (the secreted output)
]

def load_cell_proportions():
    """Load CITEgeist cell proportion results"""
    prop_path = CITEGEIST_OUTPUT / "HCC22-088-P4-S2_1i_rep_cell_prop_finetuned_results.csv"
    prop_df = pd.read_csv(prop_path, index_col=0)
    return prop_df


def load_spatial_data(filter_cancer_cells=True):
    """Load raw spatial transcriptomics data, optionally filtering to cancer cell spots"""
    print("=" * 80)
    print("LOADING SPATIAL DATA")
    print("=" * 80)

    # Load directly using scanpy (h5 file is in outs subdirectory)
    h5_path = SPATIAL_DATA / "outs" / "filtered_feature_bc_matrix.h5"
    print(f"Loading from: {h5_path}")
    adata = sc.read_10x_h5(str(h5_path), gex_only=False)
    adata.var_names_make_unique()

    print(f"Loaded: {adata.shape[0]} spots, {adata.shape[1]} features")

    # Check feature types
    if 'feature_types' in adata.var.columns:
        print(f"\nFeature types: {adata.var['feature_types'].value_counts().to_dict()}")
        # Filter to gene expression only
        gex_mask = adata.var['feature_types'] == 'Gene Expression'
        adata_gex = adata[:, gex_mask].copy()
        print(f"After GEX filter: {adata_gex.shape[0]} spots, {adata_gex.shape[1]} genes")
        adata = adata_gex

    # Filter to Cancer Cells if requested
    if filter_cancer_cells:
        print(f"\n--- Filtering to Cancer Cell spots (>{CANCER_CELL_THRESHOLD*100:.0f}% Cancer Cells) ---")
        prop_df = load_cell_proportions()

        # Get spots with high cancer cell proportion
        cancer_spots = prop_df[prop_df['Cancer Cells'] > CANCER_CELL_THRESHOLD].index.tolist()

        # Filter adata to only cancer cell spots
        common_spots = [s for s in cancer_spots if s in adata.obs_names]
        adata_cancer = adata[common_spots, :].copy()

        print(f"Cancer cell spots (>{CANCER_CELL_THRESHOLD*100:.0f}%): {len(cancer_spots)}")
        print(f"Spots in both datasets: {len(common_spots)}")
        print(f"After filter: {adata_cancer.shape[0]} spots, {adata_cancer.shape[1]} genes")

        # Report cancer cell proportion stats for filtered spots
        filtered_props = prop_df.loc[common_spots, 'Cancer Cells']
        print(f"Cancer cell proportion in filtered spots: mean={filtered_props.mean():.2f}, median={filtered_props.median():.2f}")

        adata = adata_cancer

    # Check which secretory genes are present
    present_genes = [g for g in SECRETORY_GENES if g in adata.var_names]
    missing_genes = [g for g in SECRETORY_GENES if g not in adata.var_names]

    print(f"\nSecretory genes present: {len(present_genes)}/{len(SECRETORY_GENES)}")
    print(f"  Present: {present_genes}")
    if missing_genes:
        print(f"  Missing: {missing_genes}")

    return adata, present_genes


def load_bulk_rna_data():
    """Load GSE89888 bulk RNA-seq data for MCF7 and T47D"""
    print("\n" + "=" * 80)
    print("LOADING BULK RNA-SEQ DATA (GSE89888)")
    print("=" * 80)

    # Load chaperone expression data from pipeline outputs
    expr_df = pd.read_csv(TABLES_DIR / "chaperone_expression.csv")

    # Create profiles for each condition
    profiles = {
        'MCF7_WT': expr_df.set_index('gene')['MCF7_WT_TPM'].to_dict(),
        'MCF7_D538G': expr_df.set_index('gene')['MCF7_D538G_TPM'].to_dict(),
        'T47D_WT': expr_df.set_index('gene')['T47D_WT_TPM'].to_dict(),
        'T47D_D538G': expr_df.set_index('gene')['T47D_D538G_TPM'].to_dict(),
    }

    print("\nBulk RNA profiles loaded:")
    for name, profile in profiles.items():
        print(f"  {name}: {len(profile)} genes")

    return profiles


def compare_spatial_to_bulk(adata, present_genes, bulk_profiles):
    """Compare spatial expression to bulk RNA profiles"""
    print("\n" + "=" * 80)
    print("COMPARING SPATIAL DATA TO BULK PROFILES")
    print("=" * 80)

    # Get common genes between spatial and bulk
    bulk_genes = list(bulk_profiles['MCF7_WT'].keys())
    common_genes = [g for g in present_genes if g in bulk_genes]

    print(f"\nCommon genes for comparison: {len(common_genes)}")
    print(f"  {common_genes}")

    if len(common_genes) < 3:
        print("\nWARNING: Too few genes for robust comparison")
        return None

    # Normalize the data first to handle zero counts
    adata_norm = adata.copy()
    sc.pp.filter_cells(adata_norm, min_counts=1)  # Remove empty spots
    sc.pp.normalize_total(adata_norm, target_sum=1e4)

    # Extract spatial expression (mean across all spots)
    spatial_expr = {}
    for gene in common_genes:
        if gene in adata_norm.var_names:
            gene_idx = adata_norm.var_names.get_loc(gene)
            # Get normalized counts
            counts = adata_norm.X[:, gene_idx].toarray().flatten() if hasattr(adata_norm.X, 'toarray') else adata_norm.X[:, gene_idx]
            spatial_expr[gene] = np.mean(counts)

    print("\nSpatial expression (mean CPM across spots):")
    for gene, expr in spatial_expr.items():
        print(f"  {gene}: {expr:.2f}")

    # Create vectors for comparison
    spatial_vec = np.array([spatial_expr.get(g, 0) for g in common_genes])

    # Normalize for comparison (z-score)
    spatial_vec_z = (spatial_vec - np.mean(spatial_vec)) / (np.std(spatial_vec) + 1e-10)

    results = []
    for name, profile in bulk_profiles.items():
        bulk_vec = np.array([profile.get(g, 0) for g in common_genes])
        bulk_vec_z = (bulk_vec - np.mean(bulk_vec)) / (np.std(bulk_vec) + 1e-10)

        # Calculate similarity metrics
        eucl_dist = euclidean(spatial_vec_z, bulk_vec_z)
        cos_sim = 1 - cosine(spatial_vec_z, bulk_vec_z)
        pearson_r, pearson_p = stats.pearsonr(spatial_vec, bulk_vec)
        spearman_r, spearman_p = stats.spearmanr(spatial_vec, bulk_vec)

        results.append({
            'condition': name,
            'euclidean_distance': eucl_dist,
            'cosine_similarity': cos_sim,
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p,
        })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('euclidean_distance')

    print("\n" + "-" * 60)
    print("SIMILARITY TO BULK PROFILES (lower distance = more similar)")
    print("-" * 60)
    print(results_df.to_string(index=False))

    return results_df


def analyze_gene_expression_stats(adata, present_genes):
    """Analyze expression statistics for secretory genes"""
    print("\n" + "=" * 80)
    print("EXPRESSION STATISTICS IN SPATIAL DATA")
    print("=" * 80)

    # Normalize data
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)

    results = []
    for gene in present_genes:
        if gene in adata_norm.var_names:
            gene_idx = adata_norm.var_names.get_loc(gene)
            expr = adata_norm.X[:, gene_idx].toarray().flatten() if hasattr(adata_norm.X, 'toarray') else adata_norm.X[:, gene_idx]
            results.append({
                'gene': gene,
                'mean_expr': np.mean(expr),
                'median_expr': np.median(expr),
                'std_expr': np.std(expr),
                'pct_expressing': np.sum(expr > 0) / len(expr) * 100
            })

    stats_df = pd.DataFrame(results)
    stats_df = stats_df.sort_values('mean_expr', ascending=False)

    print("\nExpression statistics (CPM):")
    print(stats_df.to_string(index=False))

    return stats_df


def main():
    print("=" * 80)
    print("SCRIPT 08: SPATIAL COMPARISON TO CELL LINES")
    print("=" * 80)

    # Load data
    adata, present_genes = load_spatial_data()
    bulk_profiles = load_bulk_rna_data()

    # Compare spatial to bulk
    comparison_df = compare_spatial_to_bulk(adata, present_genes, bulk_profiles)

    # Analyze expression statistics
    stats_df = analyze_gene_expression_stats(adata, present_genes)

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    if comparison_df is not None:
        # Find most similar condition
        most_similar = comparison_df.iloc[0]['condition']

        # Check if patient is more similar to MCF7 or T47D
        mcf7_dist = comparison_df[comparison_df['condition'].str.contains('MCF7')]['euclidean_distance'].min()
        t47d_dist = comparison_df[comparison_df['condition'].str.contains('T47D')]['euclidean_distance'].min()

        print(f"\nMost similar condition: {most_similar}")
        print(f"\nMCF7 vs T47D comparison:")
        print(f"  Closest MCF7 distance: {mcf7_dist:.3f}")
        print(f"  Closest T47D distance: {t47d_dist:.3f}")

        if mcf7_dist < t47d_dist:
            ratio = t47d_dist / mcf7_dist
            print(f"\n  --> Patient is {ratio:.1f}x MORE SIMILAR to MCF7 than T47D")
            print(f"\n  This supports the hypothesis that the patient's secretory")
            print(f"  context resembles MCF7 (where D538G causes MDK UP)")
        else:
            ratio = mcf7_dist / t47d_dist
            print(f"\n  --> Patient is {ratio:.1f}x MORE SIMILAR to T47D than MCF7")

        # Save results
        comparison_df.to_csv(TABLES_DIR / "spatial_bulk_comparison.csv", index=False)
        stats_df.to_csv(TABLES_DIR / "spatial_expression_stats.csv", index=False)

        print(f"\nSaved: {TABLES_DIR / 'spatial_bulk_comparison.csv'}")
        print(f"Saved: {TABLES_DIR / 'spatial_expression_stats.csv'}")

    print("\n" + "=" * 80)
    print("SCRIPT 08 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
