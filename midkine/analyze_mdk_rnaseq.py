#!/usr/bin/env python
"""
Analyze MDK contextual genes in GSE89888 RNA-seq data.

This script validates the spatial MDK contextual genes by checking
for interaction effects (D538G response differs between MCF7 and T47D).
"""

import os
import gzip
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind
import GEOparse
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"
TPM_FILE = os.path.join(DATA_DIR, "GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
OUTPUT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk"

# MDK contextual genes from spatial analysis
MDK_CONTEXT_GENES = [
    "MUC1", "XBP1", "AZGP1", "TFF3", "FOS", "ANKRD30A", "TPM1", "MCL1",
    "BTG2", "MLPH", "CRIP1", "WFDC2", "B4GALT1", "OAZ1", "CLU", "TRPS1",
    "EVL", "JUN", "RHOB", "MDK", "GATA3", "TFF1", "FOXA1", "ESR1"
]


def load_tpm_data():
    """Load TPM expression data."""
    print("Loading TPM data...")
    df = pd.read_csv(TPM_FILE, sep='\t', compression='gzip', index_col=0)
    print(f"Loaded: {df.shape[0]} genes x {df.shape[1]} samples")
    return df


def get_sample_metadata():
    """Fetch sample metadata from GEO."""
    print("Fetching sample metadata from GEO...")

    gse = GEOparse.get_GEO(geo="GSE89888", destdir=DATA_DIR, silent=True)

    samples = []
    for gsm_name, gsm in gse.gsms.items():
        meta = gsm.metadata
        title = meta.get('title', [''])[0]

        # Parse title: format is "CellLine-treatment-genotype-repN"
        # Example: "T47D-veh-WT-rep1", "MCF7-E2-D538G-rep2"
        parts = title.split('-')

        if len(parts) >= 4:
            cell_line = parts[0]  # T47D or MCF7
            treatment = parts[1]  # veh or E2
            genotype = parts[2]   # WT, D538G, Y537S
        else:
            cell_line = ""
            treatment = ""
            genotype = ""

        samples.append({
            'sample_id': gsm_name,
            'title': title,
            'cell_line': cell_line,
            'genotype': genotype,
            'treatment': treatment,
        })

    meta_df = pd.DataFrame(samples)
    print(f"Got metadata for {len(meta_df)} samples")
    print("\nSample breakdown:")
    print(meta_df.groupby(['cell_line', 'genotype', 'treatment']).size())

    return meta_df, gse


def map_gene_ids_to_symbols(tpm_df, gse):
    """Map NCBI gene IDs to gene symbols using the GPL platform annotation."""
    print("\nMapping gene IDs to symbols...")

    # Try to get gene symbol mapping from GPL
    gpl_name = list(gse.gpls.keys())[0] if gse.gpls else None

    if gpl_name:
        gpl = gse.gpls[gpl_name]
        # Check if there's annotation table
        if hasattr(gpl, 'table') and gpl.table is not None:
            annot = gpl.table
            print(f"GPL annotation columns: {annot.columns.tolist()[:10]}")

    # Alternative: use mygene for ID mapping
    try:
        import mygene
        mg = mygene.MyGeneInfo()

        # Get unique gene IDs
        gene_ids = tpm_df.index.astype(str).tolist()

        # Query in batches
        print(f"Querying {len(gene_ids)} gene IDs...")
        results = mg.querymany(gene_ids, scopes='entrezgene',
                              fields='symbol', species='human',
                              returnall=True)

        # Build mapping
        id_to_symbol = {}
        for hit in results['out']:
            if 'symbol' in hit:
                id_to_symbol[str(hit['query'])] = hit['symbol']

        print(f"Mapped {len(id_to_symbol)} gene IDs to symbols")

        # Rename index
        tpm_df.index = tpm_df.index.astype(str)
        tpm_df_mapped = tpm_df.copy()
        tpm_df_mapped.index = tpm_df_mapped.index.map(lambda x: id_to_symbol.get(x, x))

        return tpm_df_mapped

    except ImportError:
        print("mygene not installed. Using gene IDs directly.")
        return tpm_df


def analyze_interaction_effects(tpm_df, meta_df, genes_of_interest):
    """
    Analyze interaction effects for genes of interest.

    Interaction = (D538G - WT effect in MCF7) - (D538G - WT effect in T47D)
    Positive interaction = gene responds more to D538G in MCF7 than T47D
    """
    print("\n" + "="*60)
    print("INTERACTION EFFECT ANALYSIS")
    print("="*60)

    # Filter to vehicle treatment (baseline comparison)
    veh_mask = meta_df['treatment'].str.lower().str.contains('veh|vehicle|dmso', na=False)
    if veh_mask.sum() == 0:
        print("No vehicle samples found, using all samples")
        meta_subset = meta_df
    else:
        meta_subset = meta_df[veh_mask]
        print(f"Filtered to {len(meta_subset)} vehicle-treated samples")

    results = []

    for gene in genes_of_interest:
        if gene not in tpm_df.index:
            continue

        expr = tpm_df.loc[gene]

        gene_results = {'gene': gene}

        for cell_line in ['MCF7', 'T47D']:
            # Get samples for this cell line
            cl_mask = meta_subset['cell_line'] == cell_line

            # Split by genotype (WT vs D538G)
            wt_mask = cl_mask & (meta_subset['genotype'] == 'WT')
            mut_mask = cl_mask & (meta_subset['genotype'] == 'D538G')

            wt_samples = meta_subset[wt_mask]['sample_id'].tolist()
            mut_samples = meta_subset[mut_mask]['sample_id'].tolist()

            # Get expression values
            wt_vals = expr[[s for s in wt_samples if s in expr.index]].values
            mut_vals = expr[[s for s in mut_samples if s in expr.index]].values

            if len(wt_vals) > 0 and len(mut_vals) > 0:
                gene_results[f'{cell_line}_wt_mean'] = np.mean(wt_vals)
                gene_results[f'{cell_line}_mut_mean'] = np.mean(mut_vals)
                gene_results[f'{cell_line}_fc'] = (np.mean(mut_vals) + 0.1) / (np.mean(wt_vals) + 0.1)
                gene_results[f'{cell_line}_diff'] = np.mean(mut_vals) - np.mean(wt_vals)

                # Statistical test
                if len(wt_vals) >= 2 and len(mut_vals) >= 2:
                    try:
                        _, pval = ttest_ind(mut_vals, wt_vals)
                        gene_results[f'{cell_line}_pval'] = pval
                    except:
                        gene_results[f'{cell_line}_pval'] = np.nan

        # Compute interaction
        if 'MCF7_diff' in gene_results and 'T47D_diff' in gene_results:
            gene_results['interaction'] = gene_results['MCF7_diff'] - gene_results['T47D_diff']
            gene_results['interaction_fc'] = gene_results['MCF7_fc'] / gene_results['T47D_fc']

        results.append(gene_results)

    results_df = pd.DataFrame(results)

    # Sort by interaction effect
    if 'interaction' in results_df.columns:
        results_df = results_df.sort_values('interaction', ascending=False)

    return results_df


def main():
    print("="*60)
    print("MDK CONTEXTUAL GENE VALIDATION (GSE89888)")
    print("="*60)

    # Load data
    tpm_df = load_tpm_data()
    meta_df, gse = get_sample_metadata()

    # Map gene IDs to symbols
    tpm_df = map_gene_ids_to_symbols(tpm_df, gse)

    # Save metadata
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    meta_df.to_csv(os.path.join(OUTPUT_DIR, "GSE89888_metadata.csv"), index=False)

    # Check MDK expression first
    print("\n" + "="*60)
    print("MDK EXPRESSION CHECK")
    print("="*60)

    if 'MDK' in tpm_df.index:
        mdk_expr = tpm_df.loc['MDK']
        print(f"MDK found in data. Mean expression: {mdk_expr.mean():.2f}")

        # Quick check by condition
        for gsm in meta_df['sample_id'][:5]:
            title = meta_df[meta_df['sample_id'] == gsm]['title'].values[0]
            if gsm in mdk_expr.index:
                print(f"  {gsm}: {mdk_expr[gsm]:.2f} ({title[:50]})")
    else:
        print("MDK not found by symbol. Checking if numeric ID works...")
        # MDK NCBI gene ID is 4192
        if '4192' in tpm_df.index.astype(str):
            print("Found MDK by ID (4192)")

    # Analyze interaction effects
    results_df = analyze_interaction_effects(tpm_df, meta_df, MDK_CONTEXT_GENES)

    # Display results
    print("\n" + "="*60)
    print("INTERACTION RESULTS (MCF7-specific D538G response)")
    print("="*60)

    cols_to_show = ['gene', 'MCF7_fc', 'T47D_fc', 'interaction', 'MCF7_pval', 'T47D_pval']
    cols_available = [c for c in cols_to_show if c in results_df.columns]

    print("\nTop genes with POSITIVE interaction (respond MORE in MCF7):")
    print(results_df[cols_available].head(15).to_string(index=False))

    print("\nTop genes with NEGATIVE interaction (respond MORE in T47D):")
    print(results_df[cols_available].tail(10).to_string(index=False))

    # Save full results
    results_df.to_csv(os.path.join(OUTPUT_DIR, "mdk_context_interaction_analysis.csv"), index=False)
    print(f"\nResults saved to {OUTPUT_DIR}/mdk_context_interaction_analysis.csv")

    # Specific check for MDK
    print("\n" + "="*60)
    print("MDK INTERACTION EFFECT")
    print("="*60)

    if 'MDK' in results_df['gene'].values:
        mdk_row = results_df[results_df['gene'] == 'MDK'].iloc[0].to_dict()
        mcf7_wt = mdk_row.get('MCF7_wt_mean', np.nan)
        mcf7_mut = mdk_row.get('MCF7_mut_mean', np.nan)
        mcf7_fc = mdk_row.get('MCF7_fc', np.nan)
        t47d_wt = mdk_row.get('T47D_wt_mean', np.nan)
        t47d_mut = mdk_row.get('T47D_mut_mean', np.nan)
        t47d_fc = mdk_row.get('T47D_fc', np.nan)
        interaction = mdk_row.get('interaction', np.nan)

        print(f"MDK in MCF7: WT={mcf7_wt:.2f}, D538G={mcf7_mut:.2f}, FC={mcf7_fc:.2f}")
        print(f"MDK in T47D: WT={t47d_wt:.2f}, D538G={t47d_mut:.2f}, FC={t47d_fc:.2f}")
        print(f"Interaction: {interaction:.2f}")

        if mdk_row.get('interaction', 0) > 0:
            print("\n✓ CONFIRMED: MDK shows MCF7-specific D538G response")
        else:
            print("\n✗ MDK does NOT show MCF7-specific response in this analysis")

    # Biological interpretation
    print("\n" + "="*60)
    print("BIOLOGICAL INTERPRETATION")
    print("="*60)

    # Genes with positive interaction (MCF7-specific)
    if 'interaction' in results_df.columns:
        mcf7_specific = results_df[results_df['interaction'] > 0]['gene'].tolist()
        t47d_specific = results_df[results_df['interaction'] < 0]['gene'].tolist()

        print(f"\nGenes with MCF7-specific D538G response: {', '.join(mcf7_specific[:10])}")
        print(f"Genes with T47D-specific D538G response: {', '.join(t47d_specific[:10])}")

        # Check key contextual factors
        key_genes = ['XBP1', 'GATA3', 'FOS', 'JUN', 'MUC1', 'TFF1', 'FOXA1']
        print("\nKey contextual factor analysis:")
        for gene in key_genes:
            if gene in results_df['gene'].values:
                row = results_df[results_df['gene'] == gene].iloc[0]
                interaction = row.get('interaction', 0)
                direction = "MCF7-specific" if interaction > 0 else "T47D-specific"
                print(f"  {gene}: interaction={interaction:.2f} ({direction})")


if __name__ == "__main__":
    main()
