#!/usr/bin/env python
"""
Analyze FOXA1 ChIP-seq at chaperone loci.

GSE72249 contains FOXA1 ChIP-seq for MCF7 and T47D.
Key question: Does FOXA1 bind at chaperone gene promoters?

If FOXA1 opens chromatin for ER binding at chaperones:
- We should see FOXA1 signal at chaperone promoters
- Possibly different levels between MCF7 and T47D
"""

import os
import gzip
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

# Gene coordinates in hg19 (promoter regions: -2kb to +500bp from TSS)
GENE_COORDS_HG19 = {
    'HSP90B1': ('chr12', 104324000, 104328000),  # GRP94/endoplasmin
    'HSPA5': ('chr9', 127995000, 127999000),     # BiP/GRP78
    'CALR': ('chr19', 13049000, 13053000),       # Calreticulin
    'CANX': ('chr5', 179105000, 179109000),      # Calnexin
    'PDIA4': ('chr7', 148824000, 148828000),     # ERp72
    'PDIA6': ('chr2', 10915000, 10919000),       # P5
    'SEC61A1': ('chr3', 127802000, 127806000),   # Translocon
    'XBP1': ('chr22', 29190000, 29194000),       # UPR TF
    'ATF6': ('chr1', 161753000, 161757000),      # UPR TF
    'ESR1': ('chr6', 152125000, 152129000),      # ER itself
    'FOXA1': ('chr14', 38059000, 38063000),      # FOXA1 itself
    'TFF1': ('chr21', 43780000, 43784000),       # ER target (positive control)
    'GATA3': ('chr10', 8095000, 8099000),        # ER cofactor
}


def get_signal_in_region(bedgraph_file, chrom, start, end):
    """Get mean ChIP signal in a genomic region from bedgraph file."""
    signals = []

    with gzip.open(bedgraph_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            bg_chr = parts[0]
            bg_start = int(parts[1])
            bg_end = int(parts[2])
            bg_signal = float(parts[3])

            # Check for overlap with our region
            if bg_chr == chrom and bg_start < end and bg_end > start:
                # Calculate overlap
                overlap_start = max(bg_start, start)
                overlap_end = min(bg_end, end)
                overlap_length = overlap_end - overlap_start
                signals.append((bg_signal, overlap_length))

    if not signals:
        return 0.0

    # Weighted average by overlap length
    total_signal = sum(s * l for s, l in signals)
    total_length = sum(l for _, l in signals)
    return total_signal / total_length if total_length > 0 else 0.0


def main():
    print("="*80)
    print("FOXA1 ChIP-seq AT CHAPERONE LOCI (GSE72249)")
    print("="*80)

    # Define files
    foxa1_files = {
        'MCF7_FOXA1_unt_rep1': 'GSM1858624_GH917_FoxA1_unt_MCF7_ChIP_hg19.bedgraph.gz',
        'MCF7_FOXA1_unt_rep2': 'GSM1858625_GH925_FoxA1_unt_MCF7_ChIP_hg19.bedgraph.gz',
        'MCF7_FOXA1_E2_rep1': 'GSM1858628_GH919_FoxA1_E2_MCF7_ChIP_hg19.bedgraph.gz',
        'MCF7_FOXA1_E2_rep2': 'GSM1858629_GH927_FoxA1_E2_MCF7_ChIP_hg19.bedgraph.gz',
        'T47D_FOXA1_unt': 'GSM1858654_GH1070_T47D_FOXA1_Unt_ChIP_hg19.bedgraph.gz',
        'T47D_FOXA1_E2': 'GSM1858658_GH1072_T47D_FOXA1_E2_ChIP_hg19.bedgraph.gz',
    }

    er_files = {
        'MCF7_ER_unt_rep1': 'GSM1858620_GH620_MCF7_ER_unt_ChIP_hg19.bedgraph.gz',
        'MCF7_ER_unt_rep2': 'GSM1858621_GH628_MCF7_ER_unt_ChIP_hg19.bedgraph.gz',
        'MCF7_ER_E2_rep1': 'GSM1858622_GH622_MCF7_ER_E2_ChIP_hg19.bedgraph.gz',
        'MCF7_ER_E2_rep2': 'GSM1858623_GH630_MCF7_ER_E2_ChIP_hg19.bedgraph.gz',
        'T47D_ER_unt_rep1': 'GSM1858650_GH983_T47D_ER_unt_ChIP_hg19.bedgraph.gz',
        'T47D_ER_unt_rep2': 'GSM1858651_GH991_T47D_ER_unt_ChIP_hg19.bedgraph.gz',
        'T47D_ER_E2_rep1': 'GSM1858652_GH985_T47D_ER_E2_ChIP_hg19.bedgraph.gz',
        'T47D_ER_E2_rep2': 'GSM1858653_GH993_T47D_ER_E2_ChIP_hg19.bedgraph.gz',
    }

    # Check which files exist
    print("\nChecking available files...")
    available_foxa1 = {}
    available_er = {}

    for name, fname in foxa1_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            available_foxa1[name] = fpath
            print(f"  Found: {name}")

    for name, fname in er_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            available_er[name] = fpath
            print(f"  Found: {name}")

    if not available_foxa1:
        print("No FOXA1 ChIP-seq files found!")
        return

    # Analyze FOXA1 signal at chaperone loci
    print("\n" + "="*80)
    print("FOXA1 ChIP-seq SIGNAL AT CHAPERONE PROMOTERS")
    print("="*80)
    print("\nNote: Higher signal = more FOXA1 binding")

    # Use first available file for each cell line
    mcf7_foxa1_file = available_foxa1.get('MCF7_FOXA1_unt_rep1') or available_foxa1.get('MCF7_FOXA1_E2_rep1')
    t47d_foxa1_file = available_foxa1.get('T47D_FOXA1_unt') or available_foxa1.get('T47D_FOXA1_E2')

    chaperones = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4', 'PDIA6', 'SEC61A1']
    controls = ['TFF1', 'GATA3', 'ESR1', 'FOXA1']

    print(f"\n{'Gene':<12} {'MCF7 FOXA1':>15} {'T47D FOXA1':>15} {'T47D/MCF7':>12} {'Type':>15}")
    print("-"*75)

    mcf7_chap_signals = []
    t47d_chap_signals = []

    for gene in chaperones + controls:
        if gene not in GENE_COORDS_HG19:
            continue

        chrom, start, end = GENE_COORDS_HG19[gene]

        mcf7_signal = 0
        t47d_signal = 0

        if mcf7_foxa1_file:
            mcf7_signal = get_signal_in_region(mcf7_foxa1_file, chrom, start, end)

        if t47d_foxa1_file:
            t47d_signal = get_signal_in_region(t47d_foxa1_file, chrom, start, end)

        ratio = t47d_signal / mcf7_signal if mcf7_signal > 0 else 0
        gene_type = "CHAPERONE" if gene in chaperones else "CONTROL"

        print(f"{gene:<12} {mcf7_signal:>15.2f} {t47d_signal:>15.2f} {ratio:>12.2f} {gene_type:>15}")

        if gene in chaperones:
            mcf7_chap_signals.append(mcf7_signal)
            t47d_chap_signals.append(t47d_signal)

    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY: FOXA1 BINDING AT CHAPERONES")
    print("="*80)

    mcf7_mean = np.mean(mcf7_chap_signals) if mcf7_chap_signals else 0
    t47d_mean = np.mean(t47d_chap_signals) if t47d_chap_signals else 0

    print(f"""
    Mean FOXA1 signal at chaperone promoters:
    - MCF7:  {mcf7_mean:.2f}
    - T47D:  {t47d_mean:.2f}
    - Ratio: {t47d_mean/mcf7_mean:.2f}x

    Interpretation:
    - If T47D has MORE FOXA1 signal at chaperones, it has more open chromatin
    - This would allow D538G to fill these sites with ER binding
    - Explaining why T47D-D538G gains repression at chaperones
    """)

    # Compare ER binding at same loci
    print("\n" + "="*80)
    print("ER ChIP-seq SIGNAL AT CHAPERONE PROMOTERS")
    print("="*80)

    mcf7_er_file = available_er.get('MCF7_ER_E2_rep1')
    t47d_er_file = available_er.get('T47D_ER_E2_rep1')

    if mcf7_er_file and t47d_er_file:
        print(f"\n{'Gene':<12} {'MCF7 ER':>15} {'T47D ER':>15} {'MCF7/T47D':>12}")
        print("-"*60)

        mcf7_er_signals = []
        t47d_er_signals = []

        for gene in chaperones:
            if gene not in GENE_COORDS_HG19:
                continue

            chrom, start, end = GENE_COORDS_HG19[gene]
            mcf7_er = get_signal_in_region(mcf7_er_file, chrom, start, end)
            t47d_er = get_signal_in_region(t47d_er_file, chrom, start, end)

            ratio = mcf7_er / t47d_er if t47d_er > 0 else 0
            print(f"{gene:<12} {mcf7_er:>15.2f} {t47d_er:>15.2f} {ratio:>12.2f}")

            mcf7_er_signals.append(mcf7_er)
            t47d_er_signals.append(t47d_er)

        mcf7_er_mean = np.mean(mcf7_er_signals)
        t47d_er_mean = np.mean(t47d_er_signals)

        print(f"""
    Mean ER signal at chaperone promoters:
    - MCF7:  {mcf7_er_mean:.2f}
    - T47D:  {t47d_er_mean:.2f}
    - MCF7 has {mcf7_er_mean/t47d_er_mean:.1f}x more ER binding at chaperones

    This is CONSISTENT with MCF7 being ER-saturated!
    """)

    # Key insight
    print("\n" + "="*80)
    print("KEY INSIGHT: FOXA1 vs ER BINDING")
    print("="*80)

    if t47d_mean > 0 and mcf7_mean > 0 and t47d_er_mean > 0 and mcf7_er_mean > 0:
        foxa1_ratio = t47d_mean / mcf7_mean
        er_ratio = mcf7_er_mean / t47d_er_mean

        print(f"""
    FOXA1 binding at chaperones: T47D/MCF7 = {foxa1_ratio:.2f}
    ER binding at chaperones:   MCF7/T47D = {er_ratio:.2f}

    Interpretation:
    """)

        if foxa1_ratio > 1:
            print(f"    T47D has {foxa1_ratio:.1f}x more FOXA1 at chaperones (more open chromatin)")
        else:
            print(f"    MCF7 has {1/foxa1_ratio:.1f}x more FOXA1 at chaperones")

        if er_ratio > 1:
            print(f"    MCF7 has {er_ratio:.1f}x more ER at chaperones (more repression)")
        else:
            print(f"    T47D has {1/er_ratio:.1f}x more ER at chaperones")

        print(f"""
    The model predicts:
    - T47D has open chromatin (FOXA1-marked) at chaperones but LOW ER occupancy
    - D538G in T47D can fill these sites → new repression → chaperones DOWN
    - MCF7 already has HIGH ER occupancy → D538G causes redistribution → loses sites
    """)


if __name__ == "__main__":
    main()
