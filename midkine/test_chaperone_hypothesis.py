#!/usr/bin/env python
"""
Critically test the chaperone hypothesis for MCF7-specific MDK secretion.

Tests:
1. Do HSP90B1/HSPA5 actually show MCF7-specific upregulation?
2. What are the baseline levels? (Same issue as XBP1?)
3. Are there other confounding factors?
4. Does the ChIP-seq support ER regulation of these genes?
"""

import os
import gzip
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"
TPM_FILE = os.path.join(DATA_DIR, "GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
OUTPUT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk"

# Chaperone and secretory pathway genes to test
CHAPERONE_GENES = [
    'HSP90B1',   # GRP94 - ER chaperone
    'HSPA5',     # BiP/GRP78 - ER chaperone
    'CALR',      # Calreticulin - ER chaperone
    'CANX',      # Calnexin - ER chaperone
    'PDIA4',     # Protein disulfide isomerase
    'PDIA6',     # Protein disulfide isomerase
    'ERO1A',     # ER oxidoreductin
    'ERO1B',     # ER oxidoreductin
    'SEC61A1',   # Translocon component
    'SEC61B',    # Translocon component
    'SRP54',     # Signal recognition particle
    'SSR1',      # Signal sequence receptor
    'XBP1',      # UPR transcription factor
    'ATF6',      # UPR transcription factor
    'ATF4',      # UPR transcription factor
    'DDIT3',     # CHOP - UPR stress marker
    'MDK',       # Midkine itself
]


def load_and_map_data():
    """Load TPM data and map gene IDs to symbols."""
    print("Loading TPM data...")
    df = pd.read_csv(TPM_FILE, sep='\t', compression='gzip', index_col=0)
    print(f"Loaded: {df.shape[0]} genes x {df.shape[1]} samples")

    # Map gene IDs to symbols
    try:
        import mygene
        mg = mygene.MyGeneInfo()

        gene_ids = df.index.astype(str).tolist()
        print(f"Querying {len(gene_ids)} gene IDs...")
        results = mg.querymany(gene_ids, scopes='entrezgene',
                              fields='symbol', species='human',
                              returnall=True)

        id_to_symbol = {}
        for hit in results['out']:
            if 'symbol' in hit:
                id_to_symbol[str(hit['query'])] = hit['symbol']

        print(f"Mapped {len(id_to_symbol)} gene IDs to symbols")

        df.index = df.index.astype(str)
        df.index = df.index.map(lambda x: id_to_symbol.get(x, x))

    except ImportError:
        print("mygene not available")

    return df


def get_sample_groups(df):
    """Parse sample names to get experimental groups."""
    # Sample naming: GSM2392606 = MCF7-veh-WT-rep1
    # Based on metadata file

    groups = {
        'MCF7_WT_veh': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
        'MCF7_D538G_veh': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
        'T47D_WT_veh': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
        'T47D_D538G_veh': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
    }

    return groups


def analyze_gene(df, gene, groups):
    """Analyze a single gene across all conditions."""
    if gene not in df.index:
        return None

    expr = df.loc[gene]

    results = {'gene': gene}

    for group_name, samples in groups.items():
        valid_samples = [s for s in samples if s in expr.index]
        if valid_samples:
            vals = expr[valid_samples].values.astype(float)
            results[f'{group_name}_mean'] = np.mean(vals)
            results[f'{group_name}_std'] = np.std(vals)
            results[f'{group_name}_vals'] = vals

    # Calculate fold changes and statistics
    for cell_line in ['MCF7', 'T47D']:
        wt_key = f'{cell_line}_WT_veh_vals'
        mut_key = f'{cell_line}_D538G_veh_vals'

        if wt_key in results and mut_key in results:
            wt_vals = results[wt_key]
            mut_vals = results[mut_key]

            wt_mean = np.mean(wt_vals)
            mut_mean = np.mean(mut_vals)

            # Fold change (with pseudocount for low expressors)
            results[f'{cell_line}_fc'] = (mut_mean + 0.1) / (wt_mean + 0.1)
            results[f'{cell_line}_diff'] = mut_mean - wt_mean

            # T-test
            if len(wt_vals) >= 2 and len(mut_vals) >= 2:
                _, pval = ttest_ind(mut_vals, wt_vals)
                results[f'{cell_line}_pval'] = pval

    # Interaction effect
    if 'MCF7_diff' in results and 'T47D_diff' in results:
        results['interaction'] = results['MCF7_diff'] - results['T47D_diff']
        results['interaction_fc'] = results.get('MCF7_fc', 1) / results.get('T47D_fc', 1)

    return results


def main():
    print("="*70)
    print("CRITICAL TEST: CHAPERONE HYPOTHESIS")
    print("="*70)

    df = load_and_map_data()
    groups = get_sample_groups(df)

    # Analyze all chaperone genes
    print("\n" + "="*70)
    print("CHAPERONE & SECRETORY PATHWAY GENE ANALYSIS")
    print("="*70)

    all_results = []

    for gene in CHAPERONE_GENES:
        result = analyze_gene(df, gene, groups)
        if result:
            all_results.append(result)

    results_df = pd.DataFrame(all_results)

    # Display results
    print("\n" + "-"*70)
    print("EXPRESSION VALUES (TPM, vehicle condition)")
    print("-"*70)
    print(f"{'Gene':<12} {'MCF7-WT':>10} {'MCF7-D538G':>12} {'MCF7-FC':>8} {'T47D-WT':>10} {'T47D-D538G':>12} {'T47D-FC':>8}")
    print("-"*70)

    for _, row in results_df.iterrows():
        gene = row['gene']
        mcf7_wt = row.get('MCF7_WT_veh_mean', np.nan)
        mcf7_mut = row.get('MCF7_D538G_veh_mean', np.nan)
        mcf7_fc = row.get('MCF7_fc', np.nan)
        t47d_wt = row.get('T47D_WT_veh_mean', np.nan)
        t47d_mut = row.get('T47D_D538G_veh_mean', np.nan)
        t47d_fc = row.get('T47D_fc', np.nan)

        print(f"{gene:<12} {mcf7_wt:>10.1f} {mcf7_mut:>12.1f} {mcf7_fc:>8.2f} {t47d_wt:>10.1f} {t47d_mut:>12.1f} {t47d_fc:>8.2f}")

    # Critical test 1: Do HSP90B1/HSPA5 show MCF7-specific UPREGULATION?
    print("\n" + "="*70)
    print("CRITICAL TEST 1: MCF7-Specific UPREGULATION")
    print("="*70)
    print("For the hypothesis to hold, we need:")
    print("  - MCF7 FC > 1 (upregulated in D538G)")
    print("  - T47D FC < 1 (downregulated in D538G)")
    print("  - Interaction > 0 (MCF7 effect > T47D effect)")
    print("-"*70)

    for gene in ['HSP90B1', 'HSPA5', 'CALR', 'CANX']:
        if gene in results_df['gene'].values:
            row = results_df[results_df['gene'] == gene].iloc[0]
            mcf7_fc = row.get('MCF7_fc', np.nan)
            t47d_fc = row.get('T47D_fc', np.nan)
            interaction = row.get('interaction', np.nan)
            mcf7_pval = row.get('MCF7_pval', np.nan)
            t47d_pval = row.get('T47D_pval', np.nan)

            mcf7_up = "YES" if mcf7_fc > 1 else "NO"
            t47d_down = "YES" if t47d_fc < 1 else "NO"
            mcf7_specific = "YES" if interaction > 0 else "NO"

            print(f"\n{gene}:")
            print(f"  MCF7 FC = {mcf7_fc:.3f} (p={mcf7_pval:.4f}) - Upregulated? {mcf7_up}")
            print(f"  T47D FC = {t47d_fc:.3f} (p={t47d_pval:.4f}) - Downregulated? {t47d_down}")
            print(f"  Interaction = {interaction:.2f} - MCF7-specific? {mcf7_specific}")

    # Critical test 2: Baseline expression comparison (XBP1 problem)
    print("\n" + "="*70)
    print("CRITICAL TEST 2: BASELINE EXPRESSION (XBP1 Problem Check)")
    print("="*70)
    print("If T47D has higher baseline, the gene cannot explain MCF7 specificity")
    print("-"*70)

    for gene in ['HSP90B1', 'HSPA5', 'XBP1', 'MDK']:
        if gene in results_df['gene'].values:
            row = results_df[results_df['gene'] == gene].iloc[0]
            mcf7_wt = row.get('MCF7_WT_veh_mean', np.nan)
            t47d_wt = row.get('T47D_WT_veh_mean', np.nan)
            ratio = mcf7_wt / t47d_wt if t47d_wt > 0 else np.nan

            higher = "MCF7" if mcf7_wt > t47d_wt else "T47D"

            print(f"{gene}: MCF7-WT={mcf7_wt:.1f}, T47D-WT={t47d_wt:.1f}, Ratio={ratio:.2f} ({higher} higher)")

    # Critical test 3: Are the chaperones actually changing in the right direction?
    print("\n" + "="*70)
    print("CRITICAL TEST 3: DIRECTION OF CHANGE")
    print("="*70)
    print("Secretion UP in MCF7-D538G requires chaperones to go UP")
    print("Secretion DOWN in T47D-D538G requires chaperones to go DOWN")
    print("-"*70)

    # For MCF7: need FC > 1
    # For T47D: need FC < 1

    mcf7_consistent = []
    t47d_consistent = []
    both_consistent = []

    for gene in CHAPERONE_GENES:
        if gene in results_df['gene'].values:
            row = results_df[results_df['gene'] == gene].iloc[0]
            mcf7_fc = row.get('MCF7_fc', 1)
            t47d_fc = row.get('T47D_fc', 1)

            mcf7_ok = mcf7_fc > 1  # Should go UP in MCF7
            t47d_ok = t47d_fc < 1  # Should go DOWN in T47D

            if mcf7_ok:
                mcf7_consistent.append(gene)
            if t47d_ok:
                t47d_consistent.append(gene)
            if mcf7_ok and t47d_ok:
                both_consistent.append(gene)

    print(f"\nGenes UP in MCF7-D538G (consistent): {mcf7_consistent}")
    print(f"Genes DOWN in T47D-D538G (consistent): {t47d_consistent}")
    print(f"Genes with BOTH patterns (fully consistent): {both_consistent}")

    # MDK itself
    print("\n" + "="*70)
    print("MDK mRNA vs SECRETION PARADOX")
    print("="*70)

    if 'MDK' in results_df['gene'].values:
        row = results_df[results_df['gene'] == 'MDK'].iloc[0]
        mcf7_wt = row.get('MCF7_WT_veh_mean', np.nan)
        mcf7_mut = row.get('MCF7_D538G_veh_mean', np.nan)
        t47d_wt = row.get('T47D_WT_veh_mean', np.nan)
        t47d_mut = row.get('T47D_D538G_veh_mean', np.nan)

        print(f"MDK mRNA:")
        print(f"  MCF7: WT={mcf7_wt:.1f} -> D538G={mcf7_mut:.1f} (FC={mcf7_mut/mcf7_wt:.2f})")
        print(f"  T47D: WT={t47d_wt:.1f} -> D538G={t47d_mut:.1f} (FC={t47d_mut/t47d_wt:.2f})")
        print(f"\nMDK secretion (from spatial data):")
        print(f"  MCF7: WT=low -> D538G=high (UP)")
        print(f"  T47D: WT=high -> D538G=low (DOWN)")
        print(f"\nParadox: mRNA and secretion go in OPPOSITE directions!")

    # Final verdict
    print("\n" + "="*70)
    print("FINAL VERDICT ON CHAPERONE HYPOTHESIS")
    print("="*70)

    if len(both_consistent) > 0:
        print(f"\nGenes supporting the hypothesis: {both_consistent}")
        print("These show: UP in MCF7-D538G AND DOWN in T47D-D538G")
    else:
        print("\nNO genes show the expected pattern (UP in MCF7, DOWN in T47D)")
        print("The chaperone hypothesis is NOT SUPPORTED by the data.")

    # Save results
    # Clean up vals columns before saving
    save_cols = [c for c in results_df.columns if not c.endswith('_vals')]
    results_df[save_cols].to_csv(os.path.join(OUTPUT_DIR, "chaperone_hypothesis_test.csv"), index=False)
    print(f"\nResults saved to {OUTPUT_DIR}/chaperone_hypothesis_test.csv")


if __name__ == "__main__":
    main()
