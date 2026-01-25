#!/usr/bin/env python
"""
02_analyze_chaperone_expression.py

Datasets: GSE89888 (RNA-seq, D538G vs WT in MCF7 and T47D)
Question: Do chaperones show opposite regulation in MCF7 vs T47D?
Inputs:   data/geo/GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
Outputs:  outputs/tables/chaperone_expression.csv
          outputs/tables/interaction_stats.csv
          outputs/figures/fig2_chaperone_heatmap.pdf
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
from scipy.stats import ttest_ind
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

DATA_DIR = PIPELINE_DIR / config['paths']['data_dir']
OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']
CHAPERONES = config['parameters']['chaperone_genes']


def load_tpm():
    """Load GSE89888 TPM data with gene symbol mapping."""
    filepath = DATA_DIR / config['datasets']['GSE89888']
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)

    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(df.index.astype(str).tolist(), scopes='entrezgene',
                              fields='symbol', species='human', returnall=True)
        id_to_sym = {str(h['query']): h['symbol'] for h in results['out'] if 'symbol' in h}
        df.index = df.index.astype(str).map(lambda x: id_to_sym.get(x, x))
    except Exception as e:
        print(f"Warning: Gene mapping failed: {e}")

    return df


def validate_inputs():
    """Check required files exist."""
    filepath = DATA_DIR / config['datasets']['GSE89888']
    if not filepath.exists():
        sys.exit(f"Missing required file: {filepath}")


def main():
    print("=" * 80)
    print("SCRIPT 02: CHAPERONE EXPRESSION ANALYSIS (GSE89888)")
    print("=" * 80)

    validate_inputs()

    # Sample groups (vehicle condition)
    groups = {
        'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
        'MCF7_D538G': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
        'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
        'T47D_D538G': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
    }

    print("\nLoading GSE89888 RNA-seq data...")
    tpm = load_tpm()

    # Calculate expression and fold changes
    results = []
    for gene in CHAPERONES:
        if gene not in tpm.index:
            print(f"Warning: {gene} not found")
            continue

        mcf7_wt = tpm.loc[gene, groups['MCF7_WT']].mean()
        mcf7_d538g = tpm.loc[gene, groups['MCF7_D538G']].mean()
        t47d_wt = tpm.loc[gene, groups['T47D_WT']].mean()
        t47d_d538g = tpm.loc[gene, groups['T47D_D538G']].mean()

        mcf7_fc = mcf7_d538g / mcf7_wt if mcf7_wt > 0 else np.nan
        t47d_fc = t47d_d538g / t47d_wt if t47d_wt > 0 else np.nan

        results.append({
            'gene': gene,
            'MCF7_WT_TPM': mcf7_wt,
            'MCF7_D538G_TPM': mcf7_d538g,
            'MCF7_FC': mcf7_fc,
            'MCF7_direction': 'UP' if mcf7_fc > 1.1 else 'DOWN' if mcf7_fc < 0.9 else 'NC',
            'T47D_WT_TPM': t47d_wt,
            'T47D_D538G_TPM': t47d_d538g,
            'T47D_FC': t47d_fc,
            'T47D_direction': 'UP' if t47d_fc > 1.1 else 'DOWN' if t47d_fc < 0.9 else 'NC',
            'opposite_regulation': (mcf7_fc > 1.1 and t47d_fc < 0.9) or (mcf7_fc < 0.9 and t47d_fc > 1.1)
        })

    expr_df = pd.DataFrame(results)
    print(f"\n{expr_df.to_string()}")

    # 2-way ANOVA for interaction effect
    print("\n" + "-" * 40)
    print("2-WAY ANOVA: Cell line x Genotype interaction")
    print("-" * 40)

    interaction_results = []
    for gene in CHAPERONES:
        if gene not in tpm.index:
            continue

        # Build data for ANOVA
        data = []
        for cell_line in ['MCF7', 'T47D']:
            for genotype in ['WT', 'D538G']:
                samples = groups[f'{cell_line}_{genotype}']
                for s in samples:
                    if s in tpm.columns:
                        data.append({
                            'expression': tpm.loc[gene, s],
                            'cell_line': cell_line,
                            'genotype': genotype
                        })

        df_anova = pd.DataFrame(data)
        model = ols('expression ~ C(cell_line) * C(genotype)', data=df_anova).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        interaction_p = anova_table.loc['C(cell_line):C(genotype)', 'PR(>F)']
        interaction_results.append({
            'gene': gene,
            'interaction_pvalue': interaction_p,
            'significant': interaction_p < 0.05
        })
        print(f"{gene}: interaction p = {interaction_p:.4f} {'***' if interaction_p < 0.005 else '**' if interaction_p < 0.01 else '*' if interaction_p < 0.05 else ''}")

    interaction_df = pd.DataFrame(interaction_results)

    # Save outputs
    expr_df.to_csv(OUTPUT_DIR / "tables" / "chaperone_expression.csv", index=False)
    interaction_df.to_csv(OUTPUT_DIR / "tables" / "interaction_stats.csv", index=False)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(8, 6))
    heatmap_data = expr_df[['gene', 'MCF7_FC', 'T47D_FC']].set_index('gene')
    heatmap_data.columns = ['MCF7 D538G/WT', 'T47D D538G/WT']

    # Log2 transform for visualization
    heatmap_log = np.log2(heatmap_data)

    sns.heatmap(heatmap_log, annot=heatmap_data.round(2), fmt='', cmap='RdBu_r',
                center=0, ax=ax, cbar_kws={'label': 'log2(Fold Change)'})
    ax.set_title('Chaperone Expression: D538G vs WT\n(GSE89888)')
    plt.tight_layout()

    fig.savefig(OUTPUT_DIR / "figures" / "fig2_chaperone_heatmap.pdf")
    fig.savefig(OUTPUT_DIR / "figures" / "fig2_chaperone_heatmap.png", dpi=300)
    plt.close()

    print(f"\nSaved: outputs/tables/chaperone_expression.csv")
    print(f"Saved: outputs/tables/interaction_stats.csv")
    print(f"Saved: outputs/figures/fig2_chaperone_heatmap.pdf")

    # Summary
    opposite_count = expr_df['opposite_regulation'].sum()
    sig_count = interaction_df['significant'].sum()
    print(f"\nSummary: {opposite_count}/{len(CHAPERONES)} show opposite regulation")
    print(f"         {sig_count}/{len(CHAPERONES)} significant interaction (p<0.05)")

    print("\n" + "=" * 80)
    print("SCRIPT 02 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
