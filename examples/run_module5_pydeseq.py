#!/usr/bin/env python3
"""
Module 5 PyDESeq2 Analysis
==========================

Performs differential expression analysis comparing responder vs progressor samples.

Strategy:
1. Load gene expression from all patient samples
2. Assign response status based on patient ID (P1, P5 = responder; P2-P4, P6 = progressor)
3. Run PyDESeq2 comparing responder vs progressor spots
4. Pathway enrichment with GSEApy

Based on Module 5 response analysis:
- 3 responder-enriched programs
- 4 progressor-enriched programs

Author: CITEgeist manuscript analysis
Date: 2026-02-06
"""

import os
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Bioinformatics
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import gseapy as gp

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_response_status(sample_name: str) -> str:
    """
    Determine response status for a sample based on patient ID.

    Patient response mapping (from clinical data):
    - P1, P5: Responders
    - P2, P3, P4, P6: Progressors
    """
    if 'P1-' in sample_name or 'P5-' in sample_name:
        return 'responder'
    else:
        return 'progressor'


def load_sample_expression(sample_name: str, output_dir: Path) -> Optional[pd.DataFrame]:
    """Load gene expression for a sample."""
    gex_file = output_dir / f"{sample_name}_gene_expression_pass1.parquet"
    if not gex_file.exists():
        logger.warning(f"Gene expression file not found: {gex_file}")
        return None
    return pd.read_parquet(gex_file)


def subsample_spots(df: pd.DataFrame, n_spots: int = 1000, random_state: int = 42) -> pd.DataFrame:
    """Subsample spots if dataset is too large."""
    if len(df) > n_spots:
        return df.sample(n=n_spots, random_state=random_state)
    return df


def run_pydeseq2(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    contrast_column: str = 'condition',
    n_cpus: int = 16,
    min_gene_counts: int = 10,
    min_nonzero_spots: int = 5,
    max_genes: int = 2000
) -> pd.DataFrame:
    """
    Run PyDESeq2 differential expression analysis.

    Filtering strategy:
    1. Remove genes with very low total counts
    2. Remove genes expressed in too few spots
    3. Keep top variable genes for computational efficiency
    """
    logger.info("Preparing data for PyDESeq2...")

    # Ensure integer counts (fill NaN, clip negatives, round)
    counts_int = counts_df.fillna(0).clip(lower=0).round().astype(int)
    logger.info(f"  Initial: {counts_int.shape[0]} spots, {counts_int.shape[1]} genes")

    # Filter genes with very low total counts
    gene_counts = counts_int.sum(axis=0)
    valid_genes = gene_counts[gene_counts >= min_gene_counts].index
    counts_int = counts_int[valid_genes]
    logger.info(f"  After min count filter (>={min_gene_counts}): {counts_int.shape[1]} genes")

    # Filter genes expressed in too few spots (must be non-zero in at least N spots)
    nonzero_spots = (counts_int > 0).sum(axis=0)
    expressed_genes = nonzero_spots[nonzero_spots >= min_nonzero_spots].index
    counts_int = counts_int[expressed_genes]
    logger.info(f"  After expression filter (>={min_nonzero_spots} spots): {counts_int.shape[1]} genes")

    # Keep top variable genes for computational efficiency
    if counts_int.shape[1] > max_genes:
        # Use coefficient of variation for better variable gene selection
        gene_mean = counts_int.mean(axis=0)
        gene_std = counts_int.std(axis=0)
        gene_cv = gene_std / (gene_mean + 1)  # Add 1 to avoid division by zero
        top_var_genes = gene_cv.nlargest(max_genes).index
        counts_int = counts_int[top_var_genes]
        logger.info(f"  After variance filter (top {max_genes}): {counts_int.shape[1]} genes")

    # Remove any spots with zero counts across remaining genes
    spot_counts = counts_int.sum(axis=1)
    valid_spots = spot_counts[spot_counts > 0].index
    counts_int = counts_int.loc[valid_spots]
    metadata_df = metadata_df.loc[valid_spots]
    logger.info(f"  After spot filter (>0 counts): {counts_int.shape[0]} spots")

    logger.info(f"Running PyDESeq2 with {counts_int.shape[0]} spots, {counts_int.shape[1]} genes")

    if counts_int.shape[0] < 10:
        logger.error("Too few spots for analysis")
        return pd.DataFrame()

    # Initialize inference
    inference = DefaultInference(n_cpus=n_cpus)

    # Create DESeq2 dataset
    logger.info("Creating DESeq2 dataset...")
    dds = DeseqDataSet(
        counts=counts_int,
        metadata=metadata_df,
        design=f"~{contrast_column}",
        inference=inference,
        refit_cooks=False  # Disable Cooks outlier refitting for speed
    )

    # Fit model with deseq2() which handles all steps properly
    logger.info("Fitting DESeq2 model (this may take a while)...")
    dds.deseq2()
    logger.info("  DESeq2 fitting complete!")

    # Extract statistics
    stat_res = DeseqStats(
        dds,
        contrast=[contrast_column, 'responder', 'progressor'],
        alpha=0.05,
        inference=inference
    )
    stat_res.summary()

    # Get results
    results_df = stat_res.results_df
    results_df = results_df.sort_values('padj')

    return results_df


def run_pathway_analysis(
    de_results: pd.DataFrame,
    output_dir: Path,
    prefix: str = 'all_spots'
) -> Dict:
    """Run pathway enrichment analysis on DE genes."""
    # Get significant genes
    sig_genes = de_results[de_results['padj'] < 0.05]

    if len(sig_genes) == 0:
        logger.warning("No significant genes found")
        return {}

    # Split into up and down regulated
    up_genes = sig_genes[sig_genes['log2FoldChange'] > 0].index.tolist()
    down_genes = sig_genes[sig_genes['log2FoldChange'] < 0].index.tolist()

    results = {}
    gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021', 'MSigDB_Hallmark_2020']

    for gene_set in gene_sets:
        for direction, gene_list in [('responder_up', up_genes), ('progressor_up', down_genes)]:
            if len(gene_list) < 5:
                continue

            try:
                enr = gp.enrichr(
                    gene_list=gene_list[:500],  # Limit to top 500
                    gene_sets=gene_set,
                    organism='Human',
                    outdir=None
                )

                if not enr.results.empty:
                    enr_df = enr.results[enr.results['Adjusted P-value'] < 0.05]
                    if not enr_df.empty:
                        out_file = output_dir / f"{prefix}_{direction}_{gene_set}.csv"
                        enr_df.to_csv(out_file, index=False)
                        results[f"{direction}_{gene_set}"] = len(enr_df)
                        logger.info(f"  {direction} {gene_set}: {len(enr_df)} significant pathways")

            except Exception as e:
                logger.warning(f"Enrichment failed for {direction} {gene_set}: {e}")

    return results


def main():
    """Main analysis pipeline."""

    # Paths
    base_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")
    output_dir = base_dir / "output"
    results_dir = base_dir / "examples" / "output_module5_pydeseq"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Sample list
    samples = [
        "HCC22-088-P1-S1", "HCC22-088-P1-S2",
        "HCC22-088-P2-S1", "HCC22-088-P2-S2",
        "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
        "HCC22-088-P4-S1", "HCC22-088-P4-S2", "HCC22-088-P4-S2_1i_rep",
        "HCC22-088-P5-S1", "HCC22-088-P5-S2", "HCC22-088-P5-S2_F_rep",
        "HCC22-088-P6-S1", "HCC22-088-P6-S2_D"
    ]

    logger.info("="*60)
    logger.info("Module 5 PyDESeq2 Analysis")
    logger.info("Comparing Responder (P1, P5) vs Progressor (P2-P4, P6)")
    logger.info("="*60)

    # Collect expression data from all samples
    responder_dfs = []
    progressor_dfs = []

    for sample in samples:
        gex_df = load_sample_expression(sample, output_dir)
        if gex_df is None:
            continue

        response = get_response_status(sample)

        # Subsample to manage memory (use more spots for sparse data)
        gex_subsampled = subsample_spots(gex_df, n_spots=1000)
        gex_subsampled['sample'] = sample
        gex_subsampled['condition'] = response

        logger.info(f"  {sample}: {len(gex_subsampled)} spots ({response})")

        if response == 'responder':
            responder_dfs.append(gex_subsampled)
        else:
            progressor_dfs.append(gex_subsampled)

    # Check data availability
    if not responder_dfs or not progressor_dfs:
        logger.error("Insufficient data for comparison")
        return

    # Combine
    responder_all = pd.concat(responder_dfs, ignore_index=True)
    progressor_all = pd.concat(progressor_dfs, ignore_index=True)
    all_spots = pd.concat([responder_all, progressor_all], ignore_index=True)

    logger.info(f"\nTotal: {len(responder_all)} responder spots, {len(progressor_all)} progressor spots")

    # Prepare for PyDESeq2
    meta_cols = ['sample', 'condition']
    gene_cols = [c for c in all_spots.columns if c not in meta_cols]

    counts_df = all_spots[gene_cols].copy()
    counts_df.index = [f"spot_{i}" for i in range(len(counts_df))]

    metadata_df = all_spots[meta_cols].copy()
    metadata_df.index = counts_df.index

    # Diagnostic: check data characteristics
    logger.info(f"\nData diagnostics:")
    logger.info(f"  Total genes in data: {len(gene_cols)}")
    logger.info(f"  Data value range: {counts_df.values.min():.2f} - {counts_df.values.max():.2f}")
    logger.info(f"  Mean value: {counts_df.values.mean():.4f}")
    logger.info(f"  Fraction of zeros: {(counts_df.values == 0).mean():.2%}")

    # Check how many genes have reasonable expression
    gene_sums = counts_df.sum(axis=0)
    logger.info(f"  Genes with sum > 10: {(gene_sums > 10).sum()}")
    logger.info(f"  Genes with sum > 100: {(gene_sums > 100).sum()}")
    logger.info(f"  Genes with sum > 1000: {(gene_sums > 1000).sum()}")

    # Run differential expression
    logger.info("\nRunning PyDESeq2...")
    de_results = run_pydeseq2(counts_df, metadata_df, n_cpus=16)

    if de_results.empty:
        logger.error("PyDESeq2 analysis failed")
        return

    # Save DE results
    de_file = results_dir / "responder_vs_progressor_de_results.csv"
    de_results.to_csv(de_file)
    logger.info(f"Saved DE results to {de_file}")

    # Summary statistics
    n_sig = (de_results['padj'] < 0.05).sum()
    n_up = ((de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 0)).sum()
    n_down = ((de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < 0)).sum()

    logger.info(f"\nSignificant genes (padj < 0.05): {n_sig}")
    logger.info(f"  Up in responder: {n_up}")
    logger.info(f"  Up in progressor: {n_down}")

    # Top genes
    if n_up > 0:
        logger.info("\nTop 15 responder-enriched genes:")
        top_up = de_results[(de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 0)].head(15)
        for gene, row in top_up.iterrows():
            logger.info(f"  {gene}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")

    if n_down > 0:
        logger.info("\nTop 15 progressor-enriched genes:")
        top_down = de_results[(de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < 0)].head(15)
        for gene, row in top_down.iterrows():
            logger.info(f"  {gene}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")

    # Pathway analysis
    logger.info("\nRunning pathway enrichment...")
    pathway_results = run_pathway_analysis(de_results, results_dir)

    # Save summary
    summary = {
        'n_responder_spots': len(responder_all),
        'n_progressor_spots': len(progressor_all),
        'n_genes_tested': len(de_results),
        'n_significant': int(n_sig),
        'n_responder_up': int(n_up),
        'n_progressor_up': int(n_down),
        'pathway_results': pathway_results
    }

    if n_up > 0:
        summary['top_responder_genes'] = de_results[
            (de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 0)
        ].head(20).index.tolist()

    if n_down > 0:
        summary['top_progressor_genes'] = de_results[
            (de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < 0)
        ].head(20).index.tolist()

    summary_file = results_dir / "module5_pydeseq_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"\nSaved summary to {summary_file}")

    logger.info("\n" + "="*60)
    logger.info("Module 5 PyDESeq2 Analysis Complete")
    logger.info("="*60)


if __name__ == "__main__":
    main()
