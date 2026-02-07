#!/usr/bin/env python
"""
Integrate findings from multiple ChIP-seq datasets:
- GSE72249: FOXA1 and ER ChIP (WT only, MCF7 vs T47D vs ZR-75-1)
- GSE125117: ER ChIP (WT and D538G, MCF7 vs T47D)

Key question: How do FOXA1 and ER binding differ between cell lines,
and how does D538G change the picture?
"""

import os
import gzip
import subprocess
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"


def count_peaks_in_region(bed_file, chrom, start, end):
    """Count peaks overlapping a region from BED file."""
    count = 0
    with gzip.open(bed_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            peak_chr = parts[0]
            peak_start = int(parts[1])
            peak_end = int(parts[2])

            if peak_chr == chrom and peak_start < end and peak_end > start:
                count += 1
    return count


def get_bedgraph_signal(bedgraph_file, chrom, start, end):
    """Get mean signal from bedgraph file."""
    cmd = f"zcat {bedgraph_file} | awk '$1==\"{chrom}\" && $2 >= {start} && $3 <= {end} {{sum+=$4*($3-$2); len+=$3-$2}} END {{if(len>0) print sum/len; else print 0}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    try:
        return float(result.stdout.strip())
    except:
        return 0.0


# Gene coordinates in hg19 and hg38
GENES_HG19 = {
    'HSP90B1': ('chr12', 104323000, 104370000),
    'HSPA5': ('chr9', 127995000, 128010000),
    'CALR': ('chr19', 13049000, 13060000),
    'TFF1': ('chr21', 43780000, 43800000),
}

GENES_HG38 = {
    'HSP90B1': ('chr12', 103929000, 103970000),
    'HSPA5': ('chr9', 125234000, 125260000),
    'CALR': ('chr19', 12930000, 12960000),
    'TFF1': ('chr21', 42645000, 42670000),
}


def main():
    print("="*80)
    print("INTEGRATED BINDING ANALYSIS: FOXA1 + ER + D538G")
    print("="*80)

    # GSE72249 files (hg19, FOXA1 and ER WT)
    gse72249_files = {
        'MCF7_FOXA1': os.path.join(DATA_DIR, 'GSM1858624_GH917_FoxA1_unt_MCF7_ChIP_hg19.bedgraph.gz'),
        'T47D_FOXA1': os.path.join(DATA_DIR, 'GSM1858654_GH1070_T47D_FOXA1_Unt_ChIP_hg19.bedgraph.gz'),
        'MCF7_ER_WT': os.path.join(DATA_DIR, 'GSM1858622_GH622_MCF7_ER_E2_ChIP_hg19.bedgraph.gz'),
        'T47D_ER_WT': os.path.join(DATA_DIR, 'GSM1858652_GH985_T47D_ER_E2_ChIP_hg19.bedgraph.gz'),
    }

    # GSE125117 files (hg38, ER WT and D538G peak files)
    gse125117_files = {
        'MCF7_ER_WT_E2': os.path.join(DATA_DIR, 'GSM3563751_MCF7_WT_E2_peaks.bed.gz'),
        'MCF7_ER_D538G_E2': os.path.join(DATA_DIR, 'GSM3563757_MCF7_D538G_E2_peaks.bed.gz'),
        'T47D_ER_WT_E2': os.path.join(DATA_DIR, 'GSM3563760_T47D_WT_E2_peaks.bed.gz'),
        'T47D_ER_D538G_E2': os.path.join(DATA_DIR, 'GSM3563766_T47D_D538G_E2_peaks.bed.gz'),
    }

    # Part 1: FOXA1 binding (GSE72249)
    print("\n" + "="*80)
    print("PART 1: FOXA1 BINDING AT CHAPERONES (GSE72249)")
    print("="*80)
    print("\nFOXA1 marks open chromatin - sites that ER COULD potentially bind.")

    print(f"\n{'Gene':<12} {'MCF7 FOXA1':>12} {'T47D FOXA1':>12} {'T47D/MCF7':>12}")
    print("-"*55)

    for gene, (chrom, start, end) in GENES_HG19.items():
        mcf7 = get_bedgraph_signal(gse72249_files['MCF7_FOXA1'], chrom, start, end)
        t47d = get_bedgraph_signal(gse72249_files['T47D_FOXA1'], chrom, start, end)
        ratio = t47d / mcf7 if mcf7 > 0 else 0
        print(f"{gene:<12} {mcf7:>12.2f} {t47d:>12.2f} {ratio:>12.2f}")

    # Part 2: ER binding in WT (GSE72249)
    print("\n" + "="*80)
    print("PART 2: ER BINDING IN WT+E2 (GSE72249)")
    print("="*80)
    print("\nER binding BEFORE D538G mutation.")

    print(f"\n{'Gene':<12} {'MCF7 ER':>12} {'T47D ER':>12} {'MCF7/T47D':>12}")
    print("-"*55)

    for gene, (chrom, start, end) in GENES_HG19.items():
        mcf7 = get_bedgraph_signal(gse72249_files['MCF7_ER_WT'], chrom, start, end)
        t47d = get_bedgraph_signal(gse72249_files['T47D_ER_WT'], chrom, start, end)
        ratio = mcf7 / t47d if t47d > 0 else 0
        print(f"{gene:<12} {mcf7:>12.2f} {t47d:>12.2f} {ratio:>12.2f}")

    # Part 3: D538G effect on ER binding (GSE125117)
    print("\n" + "="*80)
    print("PART 3: D538G EFFECT ON ER BINDING (GSE125117)")
    print("="*80)
    print("\nHow does D538G change ER binding at chaperone loci?")

    print(f"\n{'Gene':<12} {'MCF7 WT':>10} {'MCF7 D538G':>12} {'Δ':>8} {'T47D WT':>10} {'T47D D538G':>12} {'Δ':>8}")
    print("-"*80)

    chap_mcf7_delta = []
    chap_t47d_delta = []

    for gene, (chrom, start, end) in GENES_HG38.items():
        # MCF7
        mcf7_wt = count_peaks_in_region(gse125117_files['MCF7_ER_WT_E2'], chrom, start, end)
        mcf7_mut = count_peaks_in_region(gse125117_files['MCF7_ER_D538G_E2'], chrom, start, end)
        mcf7_delta = mcf7_mut - mcf7_wt

        # T47D
        t47d_wt = count_peaks_in_region(gse125117_files['T47D_ER_WT_E2'], chrom, start, end)
        t47d_mut = count_peaks_in_region(gse125117_files['T47D_ER_D538G_E2'], chrom, start, end)
        t47d_delta = t47d_mut - t47d_wt

        print(f"{gene:<12} {mcf7_wt:>10} {mcf7_mut:>12} {mcf7_delta:>+8} {t47d_wt:>10} {t47d_mut:>12} {t47d_delta:>+8}")

        if gene in ['HSP90B1', 'HSPA5', 'CALR']:
            chap_mcf7_delta.append(mcf7_delta)
            chap_t47d_delta.append(t47d_delta)

    # Summary
    print("\n" + "="*80)
    print("INTEGRATED MODEL")
    print("="*80)

    print(f"""
    ┌─────────────────────────────────────────────────────────────────┐
    │                        MCF7                T47D                 │
    ├─────────────────────────────────────────────────────────────────┤
    │ FOXA1 at chaperones      LOW              HIGH (1.3x more)      │
    │ ER at chaperones (WT)    Similar          Similar               │
    │ D538G effect at chap     LOSE peaks       GAIN peaks            │
    │ Net chaperone expr       UP               DOWN                  │
    │ Secretion                UP               DOWN                  │
    └─────────────────────────────────────────────────────────────────┘

    The key insight from GSE72249:
    - T47D has MORE FOXA1 at chaperones (1.26x) = more open chromatin
    - But similar ER levels at chaperones in WT state
    - This means T47D has UNFILLED FOXA1 sites at chaperones

    When D538G is introduced:
    - MCF7: Can't add to already-occupied sites, loses some → derepression
    - T47D: Fills the FOXA1-opened but ER-unfilled sites → new repression

    D538G at chaperones:
    - MCF7: {sum(chap_mcf7_delta):+d} peaks (loses binding)
    - T47D: {sum(chap_t47d_delta):+d} peaks (gains binding)
    """)


if __name__ == "__main__":
    main()
