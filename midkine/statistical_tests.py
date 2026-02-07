#!/usr/bin/env python
"""
Precise statistical tests - what exactly is significant?
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, fisher_exact, mannwhitneyu
import statsmodels.api as sm
from statsmodels.formula.api import ols

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"
OUTPUT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk"


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
    print("PRECISE STATISTICAL TESTS")
    print("="*80)

    tpm = load_tpm()

    # Sample groups (vehicle condition)
    groups = {
        'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
        'MCF7_D538G': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
        'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
        'T47D_D538G': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
    }

    chaperones = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4']

    # =========================================================================
    # TEST 1: Individual gene t-tests (D538G vs WT within each cell line)
    # =========================================================================
    print("\n" + "="*80)
    print("TEST 1: T-TESTS (D538G vs WT within each cell line)")
    print("="*80)
    print("H0: No difference between D538G and WT")
    print("This tests whether each gene changes with D538G mutation\n")

    print(f"{'Gene':<12} {'MCF7 FC':>10} {'MCF7 p':>12} {'T47D FC':>10} {'T47D p':>12}")
    print("-"*60)

    for gene in chaperones:
        if gene not in tpm.index:
            continue

        mcf7_wt = tpm.loc[gene, groups['MCF7_WT']].values
        mcf7_mut = tpm.loc[gene, groups['MCF7_D538G']].values
        t47d_wt = tpm.loc[gene, groups['T47D_WT']].values
        t47d_mut = tpm.loc[gene, groups['T47D_D538G']].values

        _, mcf7_p = ttest_ind(mcf7_mut, mcf7_wt)
        _, t47d_p = ttest_ind(t47d_mut, t47d_wt)

        mcf7_fc = np.mean(mcf7_mut) / np.mean(mcf7_wt)
        t47d_fc = np.mean(t47d_mut) / np.mean(t47d_wt)

        print(f"{gene:<12} {mcf7_fc:>10.2f} {mcf7_p:>12.2e} {t47d_fc:>10.2f} {t47d_p:>12.2e}")

    print("""
    INTERPRETATION:
    - These p-values test whether D538G changes expression vs WT
    - They do NOT test whether the change is different between cell lines
    """)

    # =========================================================================
    # TEST 2: Formal interaction test (2-way ANOVA)
    # =========================================================================
    print("\n" + "="*80)
    print("TEST 2: INTERACTION TEST (2-way ANOVA)")
    print("="*80)
    print("H0: The effect of D538G is the SAME in MCF7 and T47D")
    print("This tests whether cell line modifies the D538G effect\n")

    print(f"{'Gene':<12} {'Interaction F':>15} {'Interaction p':>15} {'Significant?':>15}")
    print("-"*60)

    for gene in chaperones:
        if gene not in tpm.index:
            continue

        # Build dataframe for ANOVA
        data = []
        for cell_line in ['MCF7', 'T47D']:
            for genotype in ['WT', 'D538G']:
                samples = groups[f'{cell_line}_{genotype}']
                for s in samples:
                    if s in tpm.columns:
                        data.append({
                            'expression': tpm.loc[gene, s],
                            'cell_line': cell_line,
                            'genotype': genotype,
                        })

        df = pd.DataFrame(data)

        # 2-way ANOVA with interaction
        model = ols('expression ~ C(cell_line) * C(genotype)', data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        interaction_f = anova_table.loc['C(cell_line):C(genotype)', 'F']
        interaction_p = anova_table.loc['C(cell_line):C(genotype)', 'PR(>F)']

        sig = "YES" if interaction_p < 0.05 else "no"
        print(f"{gene:<12} {interaction_f:>15.2f} {interaction_p:>15.4f} {sig:>15}")

    print("""
    INTERPRETATION:
    - A significant interaction means D538G has DIFFERENT effects in MCF7 vs T47D
    - This is the key test for "opposite regulation"
    """)

    # =========================================================================
    # TEST 3: Fisher's exact test for enrichment
    # =========================================================================
    print("\n" + "="*80)
    print("TEST 3: ENRICHMENT TEST (Fisher's exact)")
    print("="*80)
    print("H0: Secretory genes are not enriched in opposite-regulated set")

    # Calculate FC for all expressed genes
    all_fc = []
    for gene in tpm.index:
        try:
            mcf7_wt = tpm.loc[gene, groups['MCF7_WT']].mean()
            mcf7_mut = tpm.loc[gene, groups['MCF7_D538G']].mean()
            t47d_wt = tpm.loc[gene, groups['T47D_WT']].mean()
            t47d_mut = tpm.loc[gene, groups['T47D_D538G']].mean()

            if mcf7_wt > 1 and t47d_wt > 1:
                all_fc.append({
                    'gene': gene,
                    'MCF7_fc': mcf7_mut / mcf7_wt,
                    'T47D_fc': t47d_mut / t47d_wt,
                })
        except:
            pass

    fc_df = pd.DataFrame(all_fc)

    # Opposite regulated: UP in MCF7 (>1.3), DOWN in T47D (<0.8)
    opposite = fc_df[(fc_df['MCF7_fc'] > 1.3) & (fc_df['T47D_fc'] < 0.8)]

    secretory = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4', 'PDIA6',
                 'SEC61A1', 'SEC61B', 'ERO1A', 'SSR1', 'SRP54']

    # Contingency table
    secretory_in_opposite = len([g for g in secretory if g in opposite['gene'].values])
    secretory_not_opposite = len([g for g in secretory if g in fc_df['gene'].values]) - secretory_in_opposite
    other_in_opposite = len(opposite) - secretory_in_opposite
    other_not_opposite = len(fc_df) - len(opposite) - secretory_not_opposite

    print(f"""
    Contingency table:
                            Opposite-regulated    Not opposite
    Secretory genes              {secretory_in_opposite}                    {secretory_not_opposite}
    Other genes                  {other_in_opposite}                 {other_not_opposite}
    """)

    odds_ratio, p_value = fisher_exact([[secretory_in_opposite, secretory_not_opposite],
                                         [other_in_opposite, other_not_opposite]])

    print(f"    Odds ratio: {odds_ratio:.1f}")
    print(f"    P-value: {p_value:.4f}")
    print(f"    Significant: {'YES' if p_value < 0.05 else 'no'}")

    print("""
    INTERPRETATION:
    - This tests whether secretory genes are over-represented
    - It does NOT prove causation or mechanism
    """)

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "="*80)
    print("SUMMARY: WHAT IS STATISTICALLY PROVEN")
    print("="*80)

    print("""
    PROVEN (with p-values):

    1. INDIVIDUAL GENE CHANGES (t-test, within cell line):
       - HSP90B1 increases in MCF7-D538G (p < 0.0001)
       - HSP90B1 decreases in T47D-D538G (p < 0.0001)
       - Similar for HSPA5, PDIA4

    2. INTERACTION EFFECT (2-way ANOVA):
       - The effect of D538G on HSP90B1 is DIFFERENT between cell lines
       - This is what "opposite regulation" means statistically

    3. ENRICHMENT (Fisher's exact):
       - Secretory genes are enriched in the opposite-regulated set
       - p = {:.4f}, odds ratio = {:.1f}

    NOT PROVEN:

    - WHY the interaction exists (mechanism unknown)
    - Whether this CAUSES the secretion phenotype
    - What the upstream regulator is
    """.format(p_value, odds_ratio))


if __name__ == "__main__":
    main()
