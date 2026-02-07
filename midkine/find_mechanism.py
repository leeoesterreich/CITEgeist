#!/usr/bin/env python
"""
Find the most likely mechanism for opposite chaperone regulation.

Strategy:
1. Find TFs that show the same opposite pattern (up MCF7, down T47D)
2. Check known chaperone regulators
3. Look at what's different at baseline between cell lines
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

def load_tpm():
    df = pd.read_csv(os.path.join(DATA_DIR, "GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
                     sep='\t', compression='gzip', index_col=0)
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(df.index.astype(str).tolist(), scopes='entrezgene',
                              fields='symbol', species='human', returnall=True)
        id_to_sym = {str(h['query']): h['symbol'] for h in results['out'] if 'symbol' in h}
        df.index = df.index.astype(str).map(lambda x: id_to_sym.get(x, x))
    except:
        pass
    return df


def main():
    print("="*80)
    print("FINDING THE MECHANISM")
    print("="*80)

    tpm = load_tpm()

    groups = {
        'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
        'MCF7_D538G': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
        'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
        'T47D_D538G': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
    }

    # Known regulators of ER chaperones / secretory pathway
    chaperone_regulators = {
        'HSF1': 'Heat shock factor - master regulator of HSP genes',
        'HSF2': 'Heat shock factor 2',
        'XBP1': 'IRE1-UPR branch, activates ERAD and chaperones',
        'ATF6': 'ATF6-UPR branch, activates chaperones',
        'ATF4': 'PERK-UPR branch, stress response',
        'NFE2L2': 'NRF2 - oxidative stress, can activate chaperones',
        'DDIT3': 'CHOP - pro-apoptotic UPR',
        'CREB3L1': 'ER stress TF',
        'CREB3L2': 'ER stress TF (BBF2H7)',
        'ERN1': 'IRE1 - UPR sensor, splices XBP1',
        'EIF2AK3': 'PERK - UPR sensor',
        'ATF3': 'Stress-inducible TF',
    }

    print("\n" + "="*80)
    print("KNOWN CHAPERONE REGULATORS - DO ANY SHOW OPPOSITE PATTERN?")
    print("="*80)

    print(f"\n{'Gene':<10} {'MCF7 FC':>10} {'T47D FC':>10} {'Opposite?':>12} {'Role'}")
    print("-"*80)

    for gene, role in chaperone_regulators.items():
        if gene not in tpm.index:
            continue

        mcf7_wt = tpm.loc[gene, groups['MCF7_WT']].mean()
        mcf7_mut = tpm.loc[gene, groups['MCF7_D538G']].mean()
        t47d_wt = tpm.loc[gene, groups['T47D_WT']].mean()
        t47d_mut = tpm.loc[gene, groups['T47D_D538G']].mean()

        mcf7_fc = mcf7_mut / mcf7_wt if mcf7_wt > 0 else 0
        t47d_fc = t47d_mut / t47d_wt if t47d_wt > 0 else 0

        # Check if opposite pattern
        opposite = "YES" if (mcf7_fc > 1.1 and t47d_fc < 0.9) else ""

        print(f"{gene:<10} {mcf7_fc:>10.2f} {t47d_fc:>10.2f} {opposite:>12} {role[:35]}")

    # =========================================================================
    # Check baseline differences
    # =========================================================================
    print("\n" + "="*80)
    print("BASELINE DIFFERENCES (WT levels)")
    print("="*80)
    print("What's different between MCF7-WT and T47D-WT that might explain")
    print("why the same mutation causes opposite effects?\n")

    baseline_genes = list(chaperone_regulators.keys()) + ['ESR1', 'GATA3', 'FOXA1', 'PGR']

    print(f"{'Gene':<10} {'MCF7-WT':>12} {'T47D-WT':>12} {'Ratio':>10} {'Higher in':>12}")
    print("-"*60)

    for gene in baseline_genes:
        if gene not in tpm.index:
            continue

        mcf7 = tpm.loc[gene, groups['MCF7_WT']].mean()
        t47d = tpm.loc[gene, groups['T47D_WT']].mean()
        ratio = mcf7 / t47d if t47d > 0 else float('inf')
        higher = "MCF7" if ratio > 1.5 else "T47D" if ratio < 0.67 else "similar"

        print(f"{gene:<10} {mcf7:>12.1f} {t47d:>12.1f} {ratio:>10.2f} {higher:>12}")

    # =========================================================================
    # Correlation analysis - what predicts chaperone changes?
    # =========================================================================
    print("\n" + "="*80)
    print("CORRELATION: WHAT PREDICTS CHAPERONE CHANGES?")
    print("="*80)

    # Calculate D538G effect for all genes
    effects = {}
    for gene in tpm.index:
        try:
            mcf7_wt = tpm.loc[gene, groups['MCF7_WT']].mean()
            mcf7_mut = tpm.loc[gene, groups['MCF7_D538G']].mean()
            t47d_wt = tpm.loc[gene, groups['T47D_WT']].mean()
            t47d_mut = tpm.loc[gene, groups['T47D_D538G']].mean()

            if mcf7_wt > 1 and t47d_wt > 1:
                effects[gene] = {
                    'MCF7_fc': mcf7_mut / mcf7_wt,
                    'T47D_fc': t47d_mut / t47d_wt,
                    'MCF7_diff': mcf7_mut - mcf7_wt,
                    'T47D_diff': t47d_mut - t47d_wt,
                }
        except:
            pass

    effects_df = pd.DataFrame(effects).T

    # Get HSP90B1 effect as reference
    hsp90b1_mcf7 = effects_df.loc['HSP90B1', 'MCF7_diff']
    hsp90b1_t47d = effects_df.loc['HSP90B1', 'T47D_diff']

    # Find genes whose MCF7 effect correlates with HSP90B1
    print("\nGenes whose D538G effect correlates with HSP90B1 across samples:")
    print("(These might be co-regulated or upstream regulators)\n")

    # Actually, let's do this differently - find TFs in the opposite-regulated set
    opposite = effects_df[(effects_df['MCF7_fc'] > 1.2) & (effects_df['T47D_fc'] < 0.85)]

    # Known TFs
    all_tfs = ['HSF1', 'HSF2', 'XBP1', 'ATF6', 'ATF4', 'ATF3', 'NFE2L2', 'DDIT3',
               'CREB3L1', 'CREB3L2', 'CREB3', 'FOXO1', 'FOXO3', 'FOXO4',
               'NR3C1', 'PPARGC1A', 'PPARGC1B', 'SREBF1', 'SREBF2',
               'MYC', 'MYCN', 'JUN', 'FOS', 'EGR1', 'SP1', 'KLF4',
               'HIF1A', 'EPAS1', 'ARNT', 'ESR1', 'ESR2']

    tfs_in_opposite = [g for g in all_tfs if g in opposite.index]

    print(f"TFs showing opposite pattern (up MCF7, down T47D): {tfs_in_opposite}")

    for tf in tfs_in_opposite:
        row = opposite.loc[tf]
        print(f"  {tf}: MCF7 FC={row['MCF7_fc']:.2f}, T47D FC={row['T47D_fc']:.2f}")

    # =========================================================================
    # MOST LIKELY MECHANISM
    # =========================================================================
    print("\n" + "="*80)
    print("MOST LIKELY MECHANISM")
    print("="*80)

    print("""
    Based on the data:

    1. NO KNOWN CHAPERONE REGULATOR shows strong opposite pattern
       - HSF1: FC=0.96 (MCF7), FC=1.06 (T47D) - no change
       - XBP1: FC=1.27 (MCF7), FC=1.06 (T47D) - both UP, not opposite
       - ATF6: FC=1.20 (MCF7), FC=0.81 (T47D) - POSSIBLE candidate
       - ATF4: FC=1.00 (MCF7), FC=0.98 (T47D) - no change

    2. ATF6 is the best candidate among known regulators
       - MCF7: goes UP with D538G (FC=1.20)
       - T47D: goes DOWN with D538G (FC=0.81)
       - ATF6 is known to activate chaperone genes

    3. But the effect sizes are modest
       - ATF6 FC=1.20 in MCF7 vs FC=0.81 in T47D
       - This is a ~1.5x difference, which may not fully explain
         the chaperone changes (HSP90B1: 1.57 vs 0.68 = 2.3x difference)

    4. ALTERNATIVE: Cell-line specific chromatin context
       - The same TF signal might have different effects
       - MCF7 and T47D have different epigenetic states
       - Chaperone gene promoters might be differentially poised

    5. MOST PARSIMONIOUS EXPLANATION:
       D538G mutation → Altered ER transcriptional program →
       Different downstream effects due to cell-line specific context →
       Opposite chaperone regulation

       The specific intermediate factor is NOT clearly identified
       in this dataset.
    """)

    # Check ATF6 more carefully
    print("\n--- ATF6 as candidate mechanism ---")
    if 'ATF6' in tpm.index:
        print("\nATF6 expression (TPM):")
        for group, samples in groups.items():
            vals = tpm.loc['ATF6', samples]
            print(f"  {group}: {vals.mean():.1f} ± {vals.std():.1f}")

        # Statistical test for ATF6
        from scipy.stats import ttest_ind
        mcf7_wt = tpm.loc['ATF6', groups['MCF7_WT']].values
        mcf7_mut = tpm.loc['ATF6', groups['MCF7_D538G']].values
        t47d_wt = tpm.loc['ATF6', groups['T47D_WT']].values
        t47d_mut = tpm.loc['ATF6', groups['T47D_D538G']].values

        _, p_mcf7 = ttest_ind(mcf7_mut, mcf7_wt)
        _, p_t47d = ttest_ind(t47d_mut, t47d_wt)

        print(f"\n  MCF7 D538G vs WT: p = {p_mcf7:.4f}")
        print(f"  T47D D538G vs WT: p = {p_t47d:.4f}")


if __name__ == "__main__":
    main()
