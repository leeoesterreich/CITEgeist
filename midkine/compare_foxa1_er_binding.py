#!/usr/bin/env python
"""
Compare FOXA1 and ER ChIP-seq between MCF7 and T47D at chaperone loci.
Uses correct hg19 coordinates.
"""

import os
import gzip
import subprocess

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

# Correct hg19 gene coordinates (promoter + gene body)
GENE_COORDS_HG19 = {
    'HSP90B1': ('chr12', 104323000, 104365000),
    'HSPA5': ('chr9', 127995000, 128005000),
    'CALR': ('chr19', 13049000, 13056000),
    'CANX': ('chr5', 179105000, 179150000),
    'XBP1': ('chr22', 29190000, 29200000),
    'TFF1': ('chr21', 43780000, 43790000),  # Known ER/FOXA1 target (positive control)
    'ESR1': ('chr6', 152125000, 152425000),
    'GREB1': ('chr2', 11650000, 11750000),  # Another ER target
}


def get_mean_signal(bedgraph_file, chrom, start, end):
    """Extract mean signal from bedgraph using awk."""
    cmd = f"zcat {bedgraph_file} | awk '$1==\"{chrom}\" && $2 >= {start} && $3 <= {end} {{sum+=$4*($3-$2); len+=$3-$2}} END {{if(len>0) print sum/len; else print 0}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    try:
        return float(result.stdout.strip())
    except:
        return 0.0


def main():
    print("="*80)
    print("FOXA1 AND ER ChIP-seq: MCF7 vs T47D AT KEY LOCI")
    print("="*80)

    # Files
    files = {
        'MCF7_FOXA1': os.path.join(DATA_DIR, 'GSM1858624_GH917_FoxA1_unt_MCF7_ChIP_hg19.bedgraph.gz'),
        'T47D_FOXA1': os.path.join(DATA_DIR, 'GSM1858654_GH1070_T47D_FOXA1_Unt_ChIP_hg19.bedgraph.gz'),
        'MCF7_ER_E2': os.path.join(DATA_DIR, 'GSM1858622_GH622_MCF7_ER_E2_ChIP_hg19.bedgraph.gz'),
        'T47D_ER_E2': os.path.join(DATA_DIR, 'GSM1858652_GH985_T47D_ER_E2_ChIP_hg19.bedgraph.gz'),
    }

    # Check files
    for name, fpath in files.items():
        if os.path.exists(fpath):
            print(f"  Found: {name}")
        else:
            print(f"  Missing: {name}")

    # FOXA1 comparison
    print("\n" + "="*80)
    print("FOXA1 ChIP-seq SIGNAL (untreated)")
    print("="*80)
    print(f"\n{'Gene':<12} {'MCF7':>12} {'T47D':>12} {'T47D/MCF7':>12}")
    print("-"*50)

    foxa1_mcf7 = []
    foxa1_t47d = []

    for gene, (chrom, start, end) in GENE_COORDS_HG19.items():
        mcf7 = get_mean_signal(files['MCF7_FOXA1'], chrom, start, end)
        t47d = get_mean_signal(files['T47D_FOXA1'], chrom, start, end)
        ratio = t47d / mcf7 if mcf7 > 0 else 0

        print(f"{gene:<12} {mcf7:>12.2f} {t47d:>12.2f} {ratio:>12.2f}")

        if gene in ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'XBP1']:
            foxa1_mcf7.append(mcf7)
            foxa1_t47d.append(t47d)

    # ER comparison
    print("\n" + "="*80)
    print("ER ChIP-seq SIGNAL (+E2)")
    print("="*80)
    print(f"\n{'Gene':<12} {'MCF7':>12} {'T47D':>12} {'MCF7/T47D':>12}")
    print("-"*50)

    er_mcf7 = []
    er_t47d = []

    for gene, (chrom, start, end) in GENE_COORDS_HG19.items():
        mcf7 = get_mean_signal(files['MCF7_ER_E2'], chrom, start, end)
        t47d = get_mean_signal(files['T47D_ER_E2'], chrom, start, end)
        ratio = mcf7 / t47d if t47d > 0 else 0

        print(f"{gene:<12} {mcf7:>12.2f} {t47d:>12.2f} {ratio:>12.2f}")

        if gene in ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'XBP1']:
            er_mcf7.append(mcf7)
            er_t47d.append(t47d)

    # Summary
    print("\n" + "="*80)
    print("SUMMARY: CHAPERONE LOCI")
    print("="*80)

    import numpy as np
    mean_foxa1_mcf7 = np.mean(foxa1_mcf7)
    mean_foxa1_t47d = np.mean(foxa1_t47d)
    mean_er_mcf7 = np.mean(er_mcf7)
    mean_er_t47d = np.mean(er_t47d)

    print(f"""
    FOXA1 binding at chaperones (mean signal):
    - MCF7:  {mean_foxa1_mcf7:.2f}
    - T47D:  {mean_foxa1_t47d:.2f}
    - T47D/MCF7 = {mean_foxa1_t47d/mean_foxa1_mcf7:.2f}x

    ER binding at chaperones (mean signal):
    - MCF7:  {mean_er_mcf7:.2f}
    - T47D:  {mean_er_t47d:.2f}
    - MCF7/T47D = {mean_er_mcf7/mean_er_t47d:.2f}x
    """)

    # Interpretation
    print("="*80)
    print("INTERPRETATION")
    print("="*80)

    foxa1_ratio = mean_foxa1_t47d / mean_foxa1_mcf7 if mean_foxa1_mcf7 > 0 else 0
    er_ratio = mean_er_mcf7 / mean_er_t47d if mean_er_t47d > 0 else 0

    print(f"""
    At chaperone gene loci:

    1. FOXA1 (chromatin accessibility marker):
       {"T47D has MORE FOXA1" if foxa1_ratio > 1 else "MCF7 has MORE FOXA1"} ({foxa1_ratio:.2f}x)
       → {"T47D has more open chromatin at chaperones" if foxa1_ratio > 1 else "MCF7 has more open chromatin"}

    2. ER (actual ER binding):
       {"MCF7 has MORE ER" if er_ratio > 1 else "T47D has MORE ER"} ({er_ratio:.2f}x)
       → {"MCF7 already occupies these sites" if er_ratio > 1 else "T47D has more ER at chaperones"}

    Model prediction:
    - If T47D has FOXA1-opened chromatin but LOW ER occupancy at chaperones:
      → D538G can fill these sites → NEW repression → chaperones DOWN
    - If MCF7 has HIGH ER occupancy at chaperones:
      → D538G redistributes → LOSES these sites → derepression → chaperones UP

    Evidence for saturation model:
    - T47D FOXA1/MCF7 FOXA1 = {foxa1_ratio:.2f} (open chromatin available in T47D)
    - MCF7 ER/T47D ER = {er_ratio:.2f} (MCF7 is more ER-occupied)
    """)


if __name__ == "__main__":
    main()
