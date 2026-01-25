#!/usr/bin/env python
"""
04_quantify_saturation.py

Datasets: GSE254216 (ATAC-seq) + GSE72249 (FOXA1 ChIP-seq)
Question: Why do MCF7 and T47D respond differently to D538G?
Inputs:   data/geo/GSM8036144_MCF7_summits.bed.gz (ATAC)
          data/geo/GSM1858624_GH917_FoxA1_unt_MCF7_ChIP_hg19.bedgraph.gz
Outputs:  outputs/tables/saturation_metrics.csv
          outputs/figures/fig4_saturation_model.pdf
"""

import os
import sys
from pathlib import Path
import gzip
import subprocess
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


def load_peaks(filepath):
    """Load BED peak file."""
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and not line.startswith('track'):
                try:
                    peaks.append({
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2])
                    })
                except ValueError:
                    continue
    return pd.DataFrame(peaks)


def get_bedgraph_signal(bedgraph_file, chrom, start, end):
    """Get mean signal from bedgraph using awk."""
    cmd = f"zcat {bedgraph_file} | awk '$1==\"{chrom}\" && $2 >= {start} && $3 <= {end} {{sum+=$4*($3-$2); len+=$3-$2}} END {{if(len>0) print sum/len; else print 0}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    try:
        return float(result.stdout.strip())
    except:
        return 0.0


def validate_inputs():
    """Check required files exist."""
    required = [
        config['datasets']['GSE254216_MCF7'],
        config['datasets']['GSE254216_T47D'],
        config['datasets']['GSE125117_MCF7_WT'],
        config['datasets']['GSE125117_T47D_WT'],
    ]
    missing = [f for f in required if not (DATA_DIR / f).exists()]
    if missing:
        sys.exit(f"Missing required files:\n" + "\n".join(missing))


def main():
    print("=" * 80)
    print("SCRIPT 04: SATURATION QUANTIFICATION")
    print("=" * 80)

    validate_inputs()

    # Load ATAC peaks (chromatin accessibility)
    print("\nLoading ATAC-seq data (GSE254216)...")
    atac = {
        'MCF7': load_peaks(DATA_DIR / config['datasets']['GSE254216_MCF7']),
        'T47D': load_peaks(DATA_DIR / config['datasets']['GSE254216_T47D']),
    }

    # Load ER ChIP peaks
    print("Loading ER ChIP-seq data (GSE125117)...")
    er = {
        'MCF7': load_peaks(DATA_DIR / config['datasets']['GSE125117_MCF7_WT']),
        'T47D': load_peaks(DATA_DIR / config['datasets']['GSE125117_T47D_WT']),
    }

    # Calculate saturation metrics
    print("\n" + "-" * 40)
    print("SATURATION METRICS")
    print("-" * 40)

    metrics = []
    for cell_line in ['MCF7', 'T47D']:
        atac_peaks = len(atac[cell_line])
        er_peaks = len(er[cell_line])
        occupancy = (er_peaks / atac_peaks * 100) if atac_peaks > 0 else 0

        metrics.append({
            'cell_line': cell_line,
            'ATAC_peaks': atac_peaks,
            'ER_peaks': er_peaks,
            'ER_occupancy_pct': occupancy,
            'available_sites_pct': 100 - occupancy
        })

        print(f"{cell_line}:")
        print(f"  ATAC peaks (open chromatin): {atac_peaks:,}")
        print(f"  ER peaks: {er_peaks:,}")
        print(f"  ER occupancy: {occupancy:.1f}%")
        print(f"  Available sites: {100 - occupancy:.1f}%")

    metrics_df = pd.DataFrame(metrics)

    # FOXA1 binding at chaperones (if files exist)
    print("\n" + "-" * 40)
    print("FOXA1 BINDING AT CHAPERONES (GSE72249)")
    print("-" * 40)

    foxa1_mcf7_file = DATA_DIR / config['datasets'].get('GSE72249_MCF7_FOXA1', '')
    foxa1_t47d_file = DATA_DIR / config['datasets'].get('GSE72249_T47D_FOXA1', '')

    foxa1_results = []
    if foxa1_mcf7_file.exists() and foxa1_t47d_file.exists():
        gene_coords_hg19 = {k: tuple(v) for k, v in config['gene_coordinates_hg19'].items()}

        for gene in ['HSP90B1', 'HSPA5', 'CALR']:
            if gene not in gene_coords_hg19:
                continue
            chrom, start, end = gene_coords_hg19[gene]

            mcf7_signal = get_bedgraph_signal(foxa1_mcf7_file, chrom, start, end)
            t47d_signal = get_bedgraph_signal(foxa1_t47d_file, chrom, start, end)
            ratio = t47d_signal / mcf7_signal if mcf7_signal > 0 else 0

            foxa1_results.append({
                'gene': gene,
                'MCF7_FOXA1': mcf7_signal,
                'T47D_FOXA1': t47d_signal,
                'T47D_MCF7_ratio': ratio
            })
            print(f"{gene}: MCF7={mcf7_signal:.2f}, T47D={t47d_signal:.2f}, ratio={ratio:.2f}")
    else:
        print("FOXA1 ChIP-seq files not found, skipping")

    foxa1_df = pd.DataFrame(foxa1_results) if foxa1_results else pd.DataFrame()

    # Save outputs
    metrics_df.to_csv(OUTPUT_DIR / "tables" / "saturation_metrics.csv", index=False)
    if not foxa1_df.empty:
        foxa1_df.to_csv(OUTPUT_DIR / "tables" / "foxa1_chaperone_binding.csv", index=False)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Saturation comparison
    ax1 = axes[0]
    x = ['MCF7\n(Saturated)', 'T47D\n(Unsaturated)']
    occupancy = metrics_df['ER_occupancy_pct'].tolist()
    available = metrics_df['available_sites_pct'].tolist()

    ax1.bar(x, occupancy, label='ER occupied', color='coral')
    ax1.bar(x, available, bottom=occupancy, label='Available', color='lightgray')
    ax1.set_ylabel('% of open chromatin')
    ax1.set_title('Chromatin Occupancy\n(ER peaks / ATAC peaks)')
    ax1.legend()

    # Global peaks comparison
    ax2 = axes[1]
    x = ['MCF7', 'T47D']
    atac_vals = [len(atac['MCF7']), len(atac['T47D'])]
    er_vals = [len(er['MCF7']), len(er['T47D'])]

    width = 0.35
    ax2.bar([0 - width/2, 1 - width/2], atac_vals, width, label='ATAC', color='steelblue')
    ax2.bar([0 + width/2, 1 + width/2], er_vals, width, label='ER', color='coral')
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(x)
    ax2.set_ylabel('Number of peaks')
    ax2.set_title('Global Peaks')
    ax2.legend()

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figures" / "fig4_saturation_model.pdf")
    fig.savefig(OUTPUT_DIR / "figures" / "fig4_saturation_model.png", dpi=300)
    plt.close()

    print(f"\nSaved: outputs/tables/saturation_metrics.csv")
    print(f"Saved: outputs/figures/fig4_saturation_model.pdf")

    print("\n" + "=" * 80)
    print("SCRIPT 04 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
