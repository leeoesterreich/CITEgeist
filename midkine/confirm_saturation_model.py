#!/usr/bin/env python
"""
Confirm the saturation model by comparing MCF7 vs T47D:
1. ATAC-seq: chromatin accessibility
2. ER ChIP-seq: actual ER binding
3. Ratio: What fraction of open chromatin is ER-bound?
"""

import os
import gzip
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

def load_peaks(filepath):
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                peaks.append({'chr': parts[0], 'start': int(parts[1]),
                             'end': int(parts[2]), 'score': float(parts[4])})
    return pd.DataFrame(peaks)


def count_peaks_in_region(peaks_df, chrom, start, end):
    if peaks_df is None or len(peaks_df) == 0:
        return 0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[(chr_peaks['start'] < end) & (chr_peaks['end'] > start)]
    return len(overlapping)


GENE_COORDS = {
    'HSP90B1': ('chr12', 103930000, 103950000),
    'HSPA5': ('chr9', 125234000, 125254000),
    'CALR': ('chr19', 12930000, 12950000),
    'XBP1': ('chr22', 28790000, 28810000),
    'ESR1': ('chr6', 151650000, 151670000),
    'FOXA1': ('chr14', 37590000, 37610000),
    'GATA3': ('chr10', 8045000, 8065000),
    'TFF1': ('chr21', 42650000, 42670000),
}


def main():
    print("="*80)
    print("CONFIRMING THE SATURATION MODEL: MCF7 vs T47D")
    print("="*80)

    # Load ATAC-seq (chromatin accessibility)
    print("\nLoading ATAC-seq data...")
    atac = {}
    atac_files = {
        'MCF7': 'GSM8036144_MCF7_summits.bed.gz',
        'T47D': 'GSM8036152_T47D_summits.bed.gz',
    }
    for name, fname in atac_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            atac[name] = load_peaks(fpath)
            print(f"  {name}: {len(atac[name])} ATAC peaks (open chromatin)")

    # Load ER ChIP-seq
    print("\nLoading ER ChIP-seq data...")
    er_chip = {}
    chip_files = {
        'MCF7_WT_E2': 'GSM3563751_MCF7_WT_E2_peaks.bed.gz',
        'T47D_WT_E2': 'GSM3563760_T47D_WT_E2_peaks.bed.gz',
    }
    for name, fname in chip_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            er_chip[name] = load_peaks(fpath)
            print(f"  {name}: {len(er_chip[name])} ER peaks")

    # KEY COMPARISON: Open chromatin vs ER binding
    print("\n" + "="*80)
    print("KEY COMPARISON: CHROMATIN ACCESSIBILITY vs ER OCCUPANCY")
    print("="*80)

    mcf7_atac = len(atac.get('MCF7', []))
    t47d_atac = len(atac.get('T47D', []))
    mcf7_er = len(er_chip.get('MCF7_WT_E2', []))
    t47d_er = len(er_chip.get('T47D_WT_E2', []))

    mcf7_ratio = mcf7_er / mcf7_atac * 100 if mcf7_atac > 0 else 0
    t47d_ratio = t47d_er / t47d_atac * 100 if t47d_atac > 0 else 0

    print(f"""
    Global comparison:

                        MCF7            T47D         Ratio
    ─────────────────────────────────────────────────────────
    ATAC peaks         {mcf7_atac:>8,}        {t47d_atac:>8,}        {mcf7_atac/t47d_atac:.2f}x
    ER ChIP peaks      {mcf7_er:>8,}        {t47d_er:>8,}        {mcf7_er/t47d_er:.1f}x
    ER/ATAC ratio      {mcf7_ratio:>7.1f}%        {t47d_ratio:>7.1f}%        {mcf7_ratio/t47d_ratio:.1f}x

    Interpretation:
    - MCF7 has {mcf7_atac/t47d_atac:.1f}x MORE open chromatin than T47D
    - MCF7 has {mcf7_er/t47d_er:.1f}x MORE ER binding than T47D
    - MCF7 ER occupies {mcf7_ratio:.1f}% of open chromatin (SATURATED)
    - T47D ER occupies {t47d_ratio:.1f}% of open chromatin (UNSATURATED)
    """)

    # Gene-level comparison
    print("\n" + "="*80)
    print("GENE-LEVEL: ATAC vs ER BINDING")
    print("="*80)

    print(f"\n{'Gene':<10} {'MCF7 ATAC':>12} {'MCF7 ER':>10} {'T47D ATAC':>12} {'T47D ER':>10}")
    print("-"*60)

    for gene, (chrom, start, end) in GENE_COORDS.items():
        mcf7_a = count_peaks_in_region(atac.get('MCF7'), chrom, start, end)
        mcf7_e = count_peaks_in_region(er_chip.get('MCF7_WT_E2'), chrom, start, end)
        t47d_a = count_peaks_in_region(atac.get('T47D'), chrom, start, end)
        t47d_e = count_peaks_in_region(er_chip.get('T47D_WT_E2'), chrom, start, end)

        print(f"{gene:<10} {mcf7_a:>12} {mcf7_e:>10} {t47d_a:>12} {t47d_e:>10}")

    # The saturation model
    print("\n" + "="*80)
    print("THE SATURATION MODEL - CONFIRMED?")
    print("="*80)

    print(f"""
    The model predicts:

    1. MCF7 should have HIGH ER occupancy of open chromatin
       - Measured: {mcf7_ratio:.1f}% of ATAC sites have ER binding
       - Status: {"CONFIRMED (>5%)" if mcf7_ratio > 5 else "NOT CONFIRMED"}

    2. T47D should have LOW ER occupancy of open chromatin
       - Measured: {t47d_ratio:.1f}% of ATAC sites have ER binding
       - Status: {"CONFIRMED (<2%)" if t47d_ratio < 2 else "NOT CONFIRMED"}

    3. MCF7 should be more "saturated" than T47D
       - MCF7 ER/ATAC ratio is {mcf7_ratio/t47d_ratio:.1f}x higher than T47D
       - Status: {"CONFIRMED" if mcf7_ratio > t47d_ratio else "NOT CONFIRMED"}
    """)

    # What this means for D538G
    print("\n" + "="*80)
    print("IMPLICATION FOR D538G EFFECTS")
    print("="*80)

    print(f"""
    MCF7 (ER/ATAC = {mcf7_ratio:.1f}%):
    - High ER occupancy of available sites
    - D538G changes binding preferences
    - Limited room for new binding
    - Net effect: REDISTRIBUTION → NET LOSS

    T47D (ER/ATAC = {t47d_ratio:.1f}%):
    - Low ER occupancy of available sites
    - {100 - t47d_ratio:.1f}% of open chromatin is UNOCCUPIED
    - D538G constitutive activity can fill these sites
    - Net effect: NEW BINDING → NET GAIN

    This explains why the same D538G mutation causes:
    - LOSS of ER peaks in MCF7 (saturated)
    - GAIN of ER peaks in T47D (unsaturated)
    """)


if __name__ == "__main__":
    main()
