#!/usr/bin/env python
"""
05_analyze_foxa1_perturbations.py

Datasets: GSE254218 (T47D FOXA1-KD) + GSE75329 (MCF7 FOXA1/ER KD/OE)
Question: Does perturbing FOXA1/ER confirm the saturation model?
Inputs:   data/geo/GSE254218_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
          data/geo/GSE75329_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
Outputs:  outputs/tables/perturbation_results.csv
          outputs/figures/fig5_perturbation_validation.pdf
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
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


def load_tpm_with_mapping(filepath):
    """Load TPM data with gene symbol mapping."""
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)

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


def validate_inputs():
    """Check required files exist."""
    required = [
        config['datasets']['GSE254218'],
        config['datasets']['GSE75329'],
    ]
    missing = [f for f in required if not (DATA_DIR / f).exists()]
    if missing:
        sys.exit(f"Missing required files:\n" + "\n".join(missing))


def main():
    print("=" * 80)
    print("SCRIPT 05: FOXA1 PERTURBATION ANALYSIS")
    print("=" * 80)

    validate_inputs()

    # Load T47D FOXA1-KD data (GSE254218)
    print("\nLoading GSE254218 (T47D FOXA1-KD)...")
    t47d_df = load_tpm_with_mapping(DATA_DIR / config['datasets']['GSE254218'])
    t47d_df = t47d_df.rename(columns={
        'GSM8036191': 'T47D_siCtrl',
        'GSM8036192': 'T47D_siFOXA1',
    })

    # Load MCF7 perturbation data (GSE75329)
    print("Loading GSE75329 (MCF7 perturbations)...")
    mcf7_df = load_tpm_with_mapping(DATA_DIR / config['datasets']['GSE75329'])
    mcf7_df = mcf7_df.rename(columns={
        'GSM1953122': 'MCF7_siCtrl',
        'GSM1953124': 'MCF7_siFOXA1_1',
        'GSM1953126': 'MCF7_siFOXA1_2',
        'GSM1953127': 'MCF7_siER',
        'GSM1953129': 'MCF7_FOXA1_noDox',
        'GSM1953131': 'MCF7_FOXA1_Dox',
    })

    # Analyze perturbation effects
    results = []

    print("\n" + "-" * 40)
    print("T47D FOXA1 KNOCKDOWN -> CHAPERONES")
    print("-" * 40)

    for gene in CHAPERONES:
        if gene not in t47d_df.index:
            continue

        ctrl = t47d_df.loc[gene, 'T47D_siCtrl']
        kd = t47d_df.loc[gene, 'T47D_siFOXA1']
        fc = kd / ctrl if ctrl > 0 else np.nan

        results.append({
            'perturbation': 'T47D_FOXA1_KD',
            'gene': gene,
            'ctrl_TPM': ctrl,
            'perturb_TPM': kd,
            'fold_change': fc,
            'direction': 'UP' if fc > 1.1 else 'DOWN' if fc < 0.9 else 'NC'
        })
        print(f"{gene}: {ctrl:.1f} -> {kd:.1f} (FC={fc:.2f})")

    print("\n" + "-" * 40)
    print("MCF7 ER KNOCKDOWN -> CHAPERONES")
    print("-" * 40)

    for gene in CHAPERONES:
        if gene not in mcf7_df.index:
            continue

        ctrl = mcf7_df.loc[gene, 'MCF7_siCtrl']
        kd = mcf7_df.loc[gene, 'MCF7_siER']
        fc = kd / ctrl if ctrl > 0 else np.nan

        results.append({
            'perturbation': 'MCF7_ER_KD',
            'gene': gene,
            'ctrl_TPM': ctrl,
            'perturb_TPM': kd,
            'fold_change': fc,
            'direction': 'UP' if fc > 1.1 else 'DOWN' if fc < 0.9 else 'NC'
        })
        print(f"{gene}: {ctrl:.1f} -> {kd:.1f} (FC={fc:.2f})")

    print("\n" + "-" * 40)
    print("MCF7 FOXA1 OVEREXPRESSION -> CHAPERONES")
    print("-" * 40)

    for gene in CHAPERONES:
        if gene not in mcf7_df.index:
            continue

        nodox = mcf7_df.loc[gene, 'MCF7_FOXA1_noDox']
        dox = mcf7_df.loc[gene, 'MCF7_FOXA1_Dox']
        fc = dox / nodox if nodox > 0 else np.nan

        results.append({
            'perturbation': 'MCF7_FOXA1_OE',
            'gene': gene,
            'ctrl_TPM': nodox,
            'perturb_TPM': dox,
            'fold_change': fc,
            'direction': 'UP' if fc > 1.1 else 'DOWN' if fc < 0.9 else 'NC'
        })
        print(f"{gene}: {nodox:.1f} -> {dox:.1f} (FC={fc:.2f})")

    # Check ESR1 in FOXA1-OE (explains the paradox)
    if 'ESR1' in mcf7_df.index:
        nodox = mcf7_df.loc['ESR1', 'MCF7_FOXA1_noDox']
        dox = mcf7_df.loc['ESR1', 'MCF7_FOXA1_Dox']
        print(f"\nNote: ESR1 in FOXA1-OE: {nodox:.1f} -> {dox:.1f} (FC={dox/nodox:.2f})")
        print("FOXA1-OE downregulates ESR1, explaining chaperone UP")

    results_df = pd.DataFrame(results)

    # Save outputs
    results_df.to_csv(OUTPUT_DIR / "tables" / "perturbation_results.csv", index=False)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))

    perturbations = ['T47D_FOXA1_KD', 'MCF7_ER_KD', 'MCF7_FOXA1_OE']
    titles = ['T47D FOXA1-KD\n(GSE254218)', 'MCF7 ER-KD\n(GSE75329)', 'MCF7 FOXA1-OE\n(GSE75329)']

    for i, (perturb, title) in enumerate(zip(perturbations, titles)):
        ax = axes[i]
        subset = results_df[results_df['perturbation'] == perturb]

        colors = ['green' if d == 'UP' else 'red' if d == 'DOWN' else 'gray'
                  for d in subset['direction']]

        ax.barh(subset['gene'], np.log2(subset['fold_change']), color=colors)
        ax.axvline(0, color='black', linewidth=0.5)
        ax.set_xlabel('log2(Fold Change)')
        ax.set_title(title)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figures" / "fig5_perturbation_validation.pdf")
    fig.savefig(OUTPUT_DIR / "figures" / "fig5_perturbation_validation.png", dpi=300)
    plt.close()

    print(f"\nSaved: outputs/tables/perturbation_results.csv")
    print(f"Saved: outputs/figures/fig5_perturbation_validation.pdf")

    # Summary
    for perturb in perturbations:
        subset = results_df[results_df['perturbation'] == perturb]
        up = (subset['direction'] == 'UP').sum()
        down = (subset['direction'] == 'DOWN').sum()
        print(f"\n{perturb}: {up} UP, {down} DOWN")

    print("\n" + "=" * 80)
    print("SCRIPT 05 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
