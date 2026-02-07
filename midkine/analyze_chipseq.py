#!/usr/bin/env python
"""
Analyze ER ChIP-seq binding at secretory pathway genes.
"""

import os
import gzip
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

# Key genes to check for ER binding
# Format: gene -> (chr, tss_start, tss_end) - using hg19/GRCh37 coordinates
# Promoter region = TSS +/- 5kb
GENE_COORDS = {
    # Secretory pathway genes (from RNA-seq analysis)
    'XBP1': ('chr22', 29190548, 29200548),      # XBP1 promoter region
    'HSP90B1': ('chr12', 104323000, 104333000), # GRP94
    'HSPA5': ('chr9', 127997000, 128007000),    # BiP/GRP78
    'NCL': ('chr2', 232320000, 232330000),      # Nucleolin (MDK receptor)
    'SEC61B': ('chr9', 101985000, 101995000),
    'ATF6': ('chr1', 161732000, 161742000),
    'DDIT3': ('chr12', 57910000, 57920000),     # CHOP

    # MDK itself
    'MDK': ('chr11', 46400000, 46410000),

    # Classic ER targets (positive controls)
    'TFF1': ('chr21', 43782000, 43792000),      # pS2
    'GREB1': ('chr2', 11650000, 11660000),
    'PGR': ('chr11', 100900000, 100910000),

    # Contextual factors from spatial analysis
    'FOS': ('chr14', 75745000, 75755000),
    'MUC1': ('chr1', 155185000, 155195000),
    'FOXA1': ('chr14', 38059000, 38069000),
    'GATA3': ('chr10', 8096000, 8106000),
}


def load_peaks(filepath):
    """Load BED file with peaks."""
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                peaks.append({
                    'chr': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3],
                    'score': float(parts[4])
                })
    return pd.DataFrame(peaks)


def check_binding_at_gene(peaks_df, chrom, start, end):
    """Check if there's ER binding in the region. Return max score."""
    if peaks_df is None or len(peaks_df) == 0:
        return 0, 0

    # Filter to chromosome
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]

    # Find overlapping peaks
    overlapping = chr_peaks[
        (chr_peaks['start'] < end) &
        (chr_peaks['end'] > start)
    ]

    if len(overlapping) == 0:
        return 0, 0

    return len(overlapping), overlapping['score'].max()


def main():
    print("="*70)
    print("ER ChIP-seq BINDING ANALYSIS AT SECRETORY PATHWAY GENES")
    print("="*70)

    # Define samples
    samples = {
        'MCF7_WT_Veh': 'GSM3563750_MCF7_WT_Vehicle_peaks.bed.gz',
        'MCF7_WT_E2': 'GSM3563751_MCF7_WT_E2_peaks.bed.gz',
        'MCF7_D538G_Veh': 'GSM3563756_MCF7_D538G_Vehicle_peaks.bed.gz',
        'MCF7_D538G_E2': 'GSM3563757_MCF7_D538G_E2_peaks.bed.gz',
        'T47D_WT_Veh': 'GSM3563759_T47D_WT_Vehicle_peaks.bed.gz',
        'T47D_WT_E2': 'GSM3563760_T47D_WT_E2_peaks.bed.gz',
        'T47D_D538G_Veh': 'GSM3563765_T47D_D538G_Vehicle_peaks.bed.gz',
        'T47D_D538G_E2': 'GSM3563766_T47D_D538G_E2_peaks.bed.gz',
    }

    # Load all peak files
    peak_data = {}
    for name, filename in samples.items():
        filepath = os.path.join(DATA_DIR, filename)
        if os.path.exists(filepath):
            peak_data[name] = load_peaks(filepath)
            print(f"Loaded {name}: {len(peak_data[name])} peaks")
        else:
            print(f"Missing: {filename}")
            peak_data[name] = None

    # Analyze binding at each gene
    print("\n" + "="*70)
    print("ER BINDING AT KEY GENES (E2-stimulated conditions)")
    print("="*70)

    results = []
    for gene, (chrom, start, end) in GENE_COORDS.items():
        row = {'Gene': gene}

        for sample in ['MCF7_WT_E2', 'MCF7_D538G_E2', 'T47D_WT_E2', 'T47D_D538G_E2']:
            n_peaks, max_score = check_binding_at_gene(peak_data.get(sample), chrom, start, end)
            row[f'{sample}_peaks'] = n_peaks
            row[f'{sample}_score'] = max_score

        results.append(row)

    results_df = pd.DataFrame(results)

    # Calculate differential binding
    results_df['MCF7_D538G_vs_WT'] = results_df['MCF7_D538G_E2_score'] - results_df['MCF7_WT_E2_score']
    results_df['T47D_D538G_vs_WT'] = results_df['T47D_D538G_E2_score'] - results_df['T47D_WT_E2_score']
    results_df['Interaction'] = results_df['MCF7_D538G_vs_WT'] - results_df['T47D_D538G_vs_WT']

    # Display results
    print("\nER binding scores (higher = stronger binding):")
    print("-"*70)

    display_cols = ['Gene', 'MCF7_WT_E2_score', 'MCF7_D538G_E2_score', 'MCF7_D538G_vs_WT',
                    'T47D_WT_E2_score', 'T47D_D538G_E2_score', 'T47D_D538G_vs_WT', 'Interaction']

    print(results_df[display_cols].sort_values('Interaction', ascending=False).to_string(index=False))

    # Summary
    print("\n" + "="*70)
    print("SUMMARY: Does ChIP-seq support MCF7-specific secretory activation?")
    print("="*70)

    # Check secretory pathway genes specifically
    secretory_genes = ['XBP1', 'HSP90B1', 'HSPA5', 'NCL', 'SEC61B', 'ATF6']
    secretory_df = results_df[results_df['Gene'].isin(secretory_genes)]

    print("\nSecretory pathway genes:")
    for _, row in secretory_df.iterrows():
        gene = row['Gene']
        mcf7_diff = row['MCF7_D538G_vs_WT']
        t47d_diff = row['T47D_D538G_vs_WT']
        interaction = row['Interaction']

        mcf7_dir = "↑" if mcf7_diff > 0 else "↓" if mcf7_diff < 0 else "="
        t47d_dir = "↑" if t47d_diff > 0 else "↓" if t47d_diff < 0 else "="

        print(f"  {gene:10s}: MCF7 D538G{mcf7_dir}WT ({mcf7_diff:+.1f}), T47D D538G{t47d_dir}WT ({t47d_diff:+.1f}), Interaction={interaction:+.1f}")

    # Check if MCF7-specific
    mcf7_specific = secretory_df[secretory_df['Interaction'] > 0]
    print(f"\nGenes with MCF7-specific D538G ER binding: {mcf7_specific['Gene'].tolist()}")

    # Vehicle (ligand-independent) analysis
    print("\n" + "="*70)
    print("LIGAND-INDEPENDENT (Vehicle) ER BINDING")
    print("="*70)
    print("D538G is known for constitutive activity - check vehicle conditions")

    veh_results = []
    for gene, (chrom, start, end) in GENE_COORDS.items():
        row = {'Gene': gene}

        for sample in ['MCF7_WT_Veh', 'MCF7_D538G_Veh', 'T47D_WT_Veh', 'T47D_D538G_Veh']:
            n_peaks, max_score = check_binding_at_gene(peak_data.get(sample), chrom, start, end)
            row[f'{sample}_score'] = max_score

        row['MCF7_D538G_vs_WT_Veh'] = row['MCF7_D538G_Veh_score'] - row['MCF7_WT_Veh_score']
        row['T47D_D538G_vs_WT_Veh'] = row['T47D_D538G_Veh_score'] - row['T47D_WT_Veh_score']
        veh_results.append(row)

    veh_df = pd.DataFrame(veh_results)

    print("\nVehicle condition - D538G constitutive binding:")
    for gene in secretory_genes:
        row = veh_df[veh_df['Gene'] == gene].iloc[0]
        mcf7_d538g = row['MCF7_D538G_Veh_score']
        mcf7_wt = row['MCF7_WT_Veh_score']
        t47d_d538g = row['T47D_D538G_Veh_score']
        t47d_wt = row['T47D_WT_Veh_score']

        print(f"  {gene:10s}: MCF7-WT={mcf7_wt:.0f}, MCF7-D538G={mcf7_d538g:.0f} | T47D-WT={t47d_wt:.0f}, T47D-D538G={t47d_d538g:.0f}")

    # Save results
    output_dir = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk"
    results_df.to_csv(os.path.join(output_dir, "chipseq_er_binding.csv"), index=False)
    print(f"\nResults saved to {output_dir}/chipseq_er_binding.csv")


if __name__ == "__main__":
    main()
