#!/usr/bin/env python3
"""
Module 5 PyDESeq2 Analysis - Pseudo-bulk Approach
==================================================

Performs differential expression analysis comparing responder vs progressor samples
using pseudo-bulk aggregation (sum expression within each sample).

This is the standard approach for multi-sample spatial/single-cell DE analysis.

Author: CITEgeist manuscript analysis
Date: 2026-02-06
"""

import os
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

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
    """Determine response status based on patient ID.

    Patient response mapping (from clinical data):
    - P1, P4: Progressors (2 patients)
    - P2, P3, P5, P6: Responders (4 patients)
    """
    if 'P1-' in sample_name or 'P4-' in sample_name:
        return 'progressor'
    return 'responder'


def get_timepoint(sample_name: str) -> str:
    """Determine timepoint from sample name.

    S1 = biopsy (pre-treatment), S2 = surgical (post-treatment).
    """
    if '-S1' in sample_name:
        return 'biopsy'
    return 'surgical'


def load_sample_expression(sample_name: str, output_dir: Path) -> Optional[pd.DataFrame]:
    """Load gene expression for a sample."""
    # Try NB output path first (per-sample subdirectory)
    gex_file = output_dir / "module3_nb" / sample_name / f"{sample_name}_gene_expression_pass1.parquet"
    if not gex_file.exists():
        # Fall back to flat output directory
        gex_file = output_dir / f"{sample_name}_gene_expression_pass1.parquet"
    if not gex_file.exists():
        logger.warning(f"Gene expression file not found: {gex_file}")
        return None
    return pd.read_parquet(gex_file)


def create_pseudobulk(samples: List[str], output_dir: Path) -> tuple:
    """
    Create pseudo-bulk expression matrix by summing counts within each sample.

    Returns:
        counts_df: genes x samples matrix
        metadata_df: sample metadata with condition and timepoint
    """
    sample_sums = {}
    conditions = {}
    timepoints = {}

    for sample in samples:
        gex_df = load_sample_expression(sample, output_dir)
        if gex_df is None:
            continue

        # Sum expression across all spots in this sample
        sample_sum = gex_df.sum(axis=0)
        sample_sums[sample] = sample_sum
        conditions[sample] = get_response_status(sample)
        timepoints[sample] = get_timepoint(sample)

        logger.info(f"  {sample}: {len(gex_df)} spots -> pseudobulk ({conditions[sample]}, {timepoints[sample]})")

    # Create counts matrix (genes x samples)
    counts_df = pd.DataFrame(sample_sums).T

    # Create metadata
    metadata_df = pd.DataFrame({
        'sample': list(conditions.keys()),
        'condition': list(conditions.values()),
        'timepoint': [timepoints[s] for s in conditions.keys()]
    })
    metadata_df = metadata_df.set_index('sample')

    return counts_df, metadata_df


def run_pydeseq2_pseudobulk(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    min_count: int = 10,
    n_cpus: int = 8
) -> pd.DataFrame:
    """
    Run PyDESeq2 on pseudo-bulk data.
    """
    logger.info(f"Input: {counts_df.shape[0]} samples, {counts_df.shape[1]} genes")

    # Ensure integer counts
    counts_int = counts_df.fillna(0).clip(lower=0).round().astype(int)

    # Filter low-count genes (must have >= min_count total across samples)
    gene_sums = counts_int.sum(axis=0)
    valid_genes = gene_sums[gene_sums >= min_count].index
    counts_int = counts_int[valid_genes]
    logger.info(f"After count filter (>={min_count}): {counts_int.shape[1]} genes")

    # Filter genes with zero variance
    gene_var = counts_int.var(axis=0)
    var_genes = gene_var[gene_var > 0].index
    counts_int = counts_int[var_genes]
    logger.info(f"After variance filter: {counts_int.shape[1]} genes")

    logger.info(f"Running PyDESeq2 with {counts_int.shape[0]} samples, {counts_int.shape[1]} genes")

    # Initialize
    inference = DefaultInference(n_cpus=n_cpus)

    # Create DESeq2 dataset
    dds = DeseqDataSet(
        counts=counts_int,
        metadata=metadata_df,
        design="~condition + timepoint",
        inference=inference,
        refit_cooks=False
    )

    # Fit model
    logger.info("Fitting DESeq2 model...")
    dds.deseq2()
    logger.info("DESeq2 fitting complete!")

    # Extract statistics
    stat_res = DeseqStats(
        dds,
        contrast=['condition', 'responder', 'progressor'],
        alpha=0.05,
        inference=inference
    )
    stat_res.summary()

    # Get results
    results_df = stat_res.results_df.copy()
    results_df = results_df.sort_values('padj')

    return results_df


def run_pathway_analysis(de_results: pd.DataFrame, output_dir: Path) -> Dict:
    """Run pathway enrichment on DE genes."""
    sig_genes = de_results[de_results['padj'] < 0.05]

    if len(sig_genes) == 0:
        logger.warning("No significant genes for pathway analysis")
        return {}

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
                    gene_list=gene_list[:500],
                    gene_sets=gene_set,
                    organism='Human',
                    outdir=None
                )

                if not enr.results.empty:
                    enr_df = enr.results[enr.results['Adjusted P-value'] < 0.05]
                    if not enr_df.empty:
                        out_file = output_dir / f"pseudobulk_{direction}_{gene_set}.csv"
                        enr_df.to_csv(out_file, index=False)
                        results[f"{direction}_{gene_set}"] = len(enr_df)
                        logger.info(f"  {direction} {gene_set}: {len(enr_df)} pathways")

            except Exception as e:
                logger.warning(f"Enrichment failed for {direction} {gene_set}: {e}")

    return results


def main():
    """Main analysis pipeline using pseudo-bulk approach."""

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
        "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
        "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
        "HCC22-088-P6-S1", "HCC22-088-P6-S2_D"
    ]

    logger.info("="*60)
    logger.info("Module 5 PyDESeq2 Analysis (Pseudo-bulk)")
    logger.info("Comparing Responder (P2, P3, P5, P6) vs Progressor (P1, P4)")
    logger.info("Design: ~condition + timepoint (biopsy vs surgical)")
    logger.info("="*60)

    # Create pseudo-bulk matrix
    logger.info("\nCreating pseudo-bulk expression matrix...")
    counts_df, metadata_df = create_pseudobulk(samples, output_dir)

    n_responder = (metadata_df['condition'] == 'responder').sum()
    n_progressor = (metadata_df['condition'] == 'progressor').sum()
    logger.info(f"\nPseudo-bulk samples: {n_responder} responders, {n_progressor} progressors")

    # Run DE analysis
    logger.info("\nRunning differential expression...")
    de_results = run_pydeseq2_pseudobulk(counts_df, metadata_df)

    if de_results.empty:
        logger.error("DE analysis failed")
        return

    # Save results
    de_file = results_dir / "pseudobulk_de_results.csv"
    de_results.to_csv(de_file)
    logger.info(f"Saved DE results to {de_file}")

    # Summary
    n_sig = (de_results['padj'] < 0.05).sum()
    n_up = ((de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 0)).sum()
    n_down = ((de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < 0)).sum()

    logger.info(f"\n{'='*40}")
    logger.info(f"Significant genes (padj < 0.05): {n_sig}")
    logger.info(f"  Up in responders: {n_up}")
    logger.info(f"  Up in progressors: {n_down}")
    logger.info(f"{'='*40}")

    # Top genes
    if n_up > 0:
        logger.info("\nTop 20 responder-enriched genes:")
        top_up = de_results[(de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 0)].head(20)
        for gene, row in top_up.iterrows():
            logger.info(f"  {gene}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")

    if n_down > 0:
        logger.info("\nTop 20 progressor-enriched genes:")
        top_down = de_results[(de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < 0)].head(20)
        for gene, row in top_down.iterrows():
            logger.info(f"  {gene}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")

    # Pathway analysis
    logger.info("\nRunning pathway enrichment...")
    pathway_results = run_pathway_analysis(de_results, results_dir)

    # Save summary
    summary = {
        'method': 'pseudo-bulk',
        'n_responder_samples': int(n_responder),
        'n_progressor_samples': int(n_progressor),
        'n_genes_tested': len(de_results),
        'n_significant': int(n_sig),
        'n_responder_up': int(n_up),
        'n_progressor_up': int(n_down),
        'pathway_results': pathway_results
    }

    if n_up > 0:
        summary['top_responder_genes'] = de_results[
            (de_results['padj'] < 0.05) & (de_results['log2FoldChange'] > 0)
        ].head(50).index.tolist()

    if n_down > 0:
        summary['top_progressor_genes'] = de_results[
            (de_results['padj'] < 0.05) & (de_results['log2FoldChange'] < 0)
        ].head(50).index.tolist()

    summary_file = results_dir / "module5_pydeseq_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"\nSaved summary to {summary_file}")

    logger.info("\n" + "="*60)
    logger.info("Module 5 PyDESeq2 Analysis Complete")
    logger.info("="*60)


if __name__ == "__main__":
    main()
