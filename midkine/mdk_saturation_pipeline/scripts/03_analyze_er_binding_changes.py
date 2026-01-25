#!/usr/bin/env python
"""
03_analyze_er_binding_changes.py

Datasets: GSE125117 (ER ChIP-seq, D538G vs WT)
Question: Does ER binding change oppositely at chaperone loci?
Inputs:   data/geo/GSM3563751_MCF7_WT_E2_peaks.bed.gz (and 3 others)
Outputs:  outputs/tables/binding_changes.csv
          outputs/figures/fig3_binding_changes.pdf
"""

import os
import sys
from pathlib import Path
import gzip
import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

DATA_DIR = PIPELINE_DIR / config['paths']['data_dir']
OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']
GENE_COORDS = {k: tuple(v) for k, v in config['gene_coordinates_hg38'].items()}


def load_peaks(filepath):
    """Load BED peak file."""
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                peaks.append({
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2])
                })
    return pd.DataFrame(peaks)


def count_peaks_in_region(peaks_df, chrom, start, end):
    """Count peaks overlapping a genomic region."""
    if peaks_df is None or len(peaks_df) == 0:
        return 0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[(chr_peaks['start'] < end) & (chr_peaks['end'] > start)]
    return len(overlapping)


def validate_inputs():
    """Check required files exist."""
    required = [
        config['datasets']['GSE125117_MCF7_WT'],
        config['datasets']['GSE125117_MCF7_D538G'],
        config['datasets']['GSE125117_T47D_WT'],
        config['datasets']['GSE125117_T47D_D538G'],
    ]
    missing = [f for f in required if not (DATA_DIR / f).exists()]
    if missing:
        sys.exit(f"Missing required files:\n" + "\n".join(missing))


def main():
    print("=" * 80)
    print("SCRIPT 03: ER BINDING CHANGES (GSE125117)")
    print("=" * 80)

    validate_inputs()

    # Load peak files
    print("\nLoading ER ChIP-seq peaks...")
    peaks = {
        'MCF7_WT': load_peaks(DATA_DIR / config['datasets']['GSE125117_MCF7_WT']),
        'MCF7_D538G': load_peaks(DATA_DIR / config['datasets']['GSE125117_MCF7_D538G']),
        'T47D_WT': load_peaks(DATA_DIR / config['datasets']['GSE125117_T47D_WT']),
        'T47D_D538G': load_peaks(DATA_DIR / config['datasets']['GSE125117_T47D_D538G']),
    }

    for name, df in peaks.items():
        print(f"  {name}: {len(df)} peaks")

    # Global binding changes
    print("\n" + "-" * 40)
    print("GLOBAL ER BINDING CHANGES")
    print("-" * 40)

    mcf7_global_delta = len(peaks['MCF7_D538G']) - len(peaks['MCF7_WT'])
    t47d_global_delta = len(peaks['T47D_D538G']) - len(peaks['T47D_WT'])
    print(f"MCF7: {len(peaks['MCF7_WT'])} -> {len(peaks['MCF7_D538G'])} ({mcf7_global_delta:+d} peaks)")
    print(f"T47D: {len(peaks['T47D_WT'])} -> {len(peaks['T47D_D538G'])} ({t47d_global_delta:+d} peaks)")

    # Binding at chaperone loci
    print("\n" + "-" * 40)
    print("BINDING AT CHAPERONE LOCI")
    print("-" * 40)

    results = []
    chaperones = config['parameters']['chaperone_genes'][:5]  # Top 5

    for gene in chaperones:
        if gene not in GENE_COORDS:
            continue

        chrom, start, end = GENE_COORDS[gene]

        mcf7_wt = count_peaks_in_region(peaks['MCF7_WT'], chrom, start, end)
        mcf7_d538g = count_peaks_in_region(peaks['MCF7_D538G'], chrom, start, end)
        t47d_wt = count_peaks_in_region(peaks['T47D_WT'], chrom, start, end)
        t47d_d538g = count_peaks_in_region(peaks['T47D_D538G'], chrom, start, end)

        mcf7_delta = mcf7_d538g - mcf7_wt
        t47d_delta = t47d_d538g - t47d_wt

        results.append({
            'gene': gene,
            'MCF7_WT_peaks': mcf7_wt,
            'MCF7_D538G_peaks': mcf7_d538g,
            'MCF7_delta': mcf7_delta,
            'T47D_WT_peaks': t47d_wt,
            'T47D_D538G_peaks': t47d_d538g,
            'T47D_delta': t47d_delta,
            'opposite_binding': (mcf7_delta < 0 and t47d_delta > 0) or (mcf7_delta > 0 and t47d_delta < 0)
        })

        print(f"{gene}: MCF7 {mcf7_delta:+d}, T47D {t47d_delta:+d}")

    binding_df = pd.DataFrame(results)

    # Save outputs
    binding_df.to_csv(OUTPUT_DIR / "tables" / "binding_changes.csv", index=False)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Global binding bar chart
    ax1 = axes[0]
    x = ['MCF7', 'T47D']
    wt_vals = [len(peaks['MCF7_WT']), len(peaks['T47D_WT'])]
    d538g_vals = [len(peaks['MCF7_D538G']), len(peaks['T47D_D538G'])]

    width = 0.35
    ax1.bar([0 - width/2, 1 - width/2], wt_vals, width, label='WT', color='steelblue')
    ax1.bar([0 + width/2, 1 + width/2], d538g_vals, width, label='D538G', color='coral')
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(x)
    ax1.set_ylabel('Number of ER peaks')
    ax1.set_title('Global ER Binding')
    ax1.legend()

    # Chaperone locus binding
    ax2 = axes[1]
    genes = binding_df['gene'].tolist()
    mcf7_deltas = binding_df['MCF7_delta'].tolist()
    t47d_deltas = binding_df['T47D_delta'].tolist()

    x_pos = np.arange(len(genes))
    ax2.bar(x_pos - width/2, mcf7_deltas, width, label='MCF7', color='steelblue')
    ax2.bar(x_pos + width/2, t47d_deltas, width, label='T47D', color='coral')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(genes, rotation=45, ha='right')
    ax2.set_ylabel('Delta ER peaks (D538G - WT)')
    ax2.set_title('ER Binding at Chaperone Loci')
    ax2.axhline(0, color='black', linewidth=0.5)
    ax2.legend()

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figures" / "fig3_binding_changes.pdf")
    fig.savefig(OUTPUT_DIR / "figures" / "fig3_binding_changes.png", dpi=300)
    plt.close()

    print(f"\nSaved: outputs/tables/binding_changes.csv")
    print(f"Saved: outputs/figures/fig3_binding_changes.pdf")

    opposite_count = binding_df['opposite_binding'].sum()
    print(f"\nSummary: {opposite_count}/{len(results)} genes show opposite binding changes")

    print("\n" + "=" * 80)
    print("SCRIPT 03 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
