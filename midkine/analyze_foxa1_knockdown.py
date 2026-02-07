#!/usr/bin/env python
"""
Analyze FOXA1 knockdown in T47D.

If FOXA1 opens chromatin for ER binding, then:
- FOXA1 KD → less open chromatin → less ER binding → less repression → chaperones UP?

Or if FOXA1 itself activates chaperones:
- FOXA1 KD → chaperones DOWN

Let's see what the data shows.
"""

import os
import gzip
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

def load_foxa1_rnaseq():
    """Load FOXA1 knockdown RNA-seq data."""
    df = pd.read_csv(os.path.join(DATA_DIR, "GSE254218_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
                     sep='\t', compression='gzip', index_col=0)

    # Map gene IDs to symbols
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(df.index.astype(str).tolist(), scopes='entrezgene',
                              fields='symbol', species='human', returnall=True)
        id_to_sym = {str(h['query']): h['symbol'] for h in results['out'] if 'symbol' in h}
        df.index = df.index.astype(str).map(lambda x: id_to_sym.get(x, x))
    except:
        pass

    # Rename columns based on GEO metadata
    df = df.rename(columns={
        'GSM8036191': 'T47D_siCtrl',
        'GSM8036192': 'T47D_siFOXA1',
        'GSM8036193': 'T47D_siGRHL2',
        'GSM8036194': 'HCC38_siCtrl',
        'GSM8036195': 'HCC38_siGRHL2',
    })

    return df


def load_atac_peaks(filepath):
    """Load ATAC-seq summit peaks."""
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


def count_peaks_at_gene(peaks_df, chrom, start, end):
    """Count ATAC peaks overlapping a gene region."""
    if peaks_df is None or len(peaks_df) == 0:
        return 0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[(chr_peaks['start'] < end) & (chr_peaks['end'] > start)]
    return len(overlapping)


# Gene coordinates (hg38)
GENE_COORDS = {
    'HSP90B1': ('chr12', 103930000, 103950000),
    'HSPA5': ('chr9', 125234000, 125254000),
    'CALR': ('chr19', 12930000, 12950000),
    'CANX': ('chr5', 179680000, 179700000),
    'PDIA4': ('chr7', 150030000, 150050000),
    'XBP1': ('chr22', 28790000, 28810000),
    'ATF6': ('chr1', 161120000, 161140000),
    'FOXA1': ('chr14', 37590000, 37610000),
    'ESR1': ('chr6', 151650000, 151670000),
    'GATA3': ('chr10', 8045000, 8065000),
}


def main():
    print("="*80)
    print("FOXA1 KNOCKDOWN ANALYSIS IN T47D")
    print("="*80)
    print("""
    Hypothesis: FOXA1 opens chromatin, allowing ER to bind and repress chaperones.

    Prediction: FOXA1 KD → closed chromatin → no ER binding → derepression → chaperones UP

    Alternative: FOXA1 directly activates chaperones → FOXA1 KD → chaperones DOWN
    """)

    # Load RNA-seq
    print("\nLoading RNA-seq data...")
    rna = load_foxa1_rnaseq()
    print(f"Columns: {rna.columns.tolist()}")

    # Verify FOXA1 knockdown worked
    print("\n" + "="*80)
    print("VERIFICATION: Did FOXA1 knockdown work?")
    print("="*80)

    if 'FOXA1' in rna.index:
        foxa1_ctrl = rna.loc['FOXA1', 'T47D_siCtrl']
        foxa1_kd = rna.loc['FOXA1', 'T47D_siFOXA1']
        fc = foxa1_kd / foxa1_ctrl if foxa1_ctrl > 0 else 0
        print(f"\nFOXA1 expression:")
        print(f"  siCtrl: {foxa1_ctrl:.1f} TPM")
        print(f"  siFOXA1: {foxa1_kd:.1f} TPM")
        print(f"  Fold change: {fc:.2f} ({(1-fc)*100:.0f}% knockdown)")

    # Analyze chaperone genes
    print("\n" + "="*80)
    print("CHAPERONE GENE EXPRESSION: siCtrl vs siFOXA1")
    print("="*80)

    chaperones = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4', 'PDIA6',
                  'SEC61A1', 'SEC61B', 'ERO1A', 'SSR1']

    print(f"\n{'Gene':<12} {'siCtrl':>12} {'siFOXA1':>12} {'FC':>10} {'Direction':>12}")
    print("-"*60)

    up_count = 0
    down_count = 0

    for gene in chaperones:
        if gene not in rna.index:
            continue

        ctrl = rna.loc[gene, 'T47D_siCtrl']
        kd = rna.loc[gene, 'T47D_siFOXA1']
        fc = kd / ctrl if ctrl > 0 else 0

        direction = "UP" if fc > 1.1 else "DOWN" if fc < 0.9 else "~"
        if fc > 1.1:
            up_count += 1
        elif fc < 0.9:
            down_count += 1

        print(f"{gene:<12} {ctrl:>12.1f} {kd:>12.1f} {fc:>10.2f} {direction:>12}")

    print(f"\nSummary: {up_count} UP, {down_count} DOWN with FOXA1 knockdown")

    # UPR transcription factors
    print("\n" + "="*80)
    print("UPR REGULATORS: siCtrl vs siFOXA1")
    print("="*80)

    upr_genes = ['XBP1', 'ATF6', 'ATF4', 'ATF3', 'DDIT3', 'ERN1', 'EIF2AK3']

    print(f"\n{'Gene':<12} {'siCtrl':>12} {'siFOXA1':>12} {'FC':>10} {'Direction':>12}")
    print("-"*60)

    for gene in upr_genes:
        if gene not in rna.index:
            continue

        ctrl = rna.loc[gene, 'T47D_siCtrl']
        kd = rna.loc[gene, 'T47D_siFOXA1']
        fc = kd / ctrl if ctrl > 0 else 0
        direction = "UP" if fc > 1.1 else "DOWN" if fc < 0.9 else "~"

        print(f"{gene:<12} {ctrl:>12.1f} {kd:>12.1f} {fc:>10.2f} {direction:>12}")

    # ER system genes
    print("\n" + "="*80)
    print("ER SYSTEM GENES: siCtrl vs siFOXA1")
    print("="*80)

    er_genes = ['ESR1', 'GATA3', 'PGR', 'TFF1', 'GREB1']

    print(f"\n{'Gene':<12} {'siCtrl':>12} {'siFOXA1':>12} {'FC':>10} {'Direction':>12}")
    print("-"*60)

    for gene in er_genes:
        if gene not in rna.index:
            continue

        ctrl = rna.loc[gene, 'T47D_siCtrl']
        kd = rna.loc[gene, 'T47D_siFOXA1']
        fc = kd / ctrl if ctrl > 0 else 0
        direction = "UP" if fc > 1.1 else "DOWN" if fc < 0.9 else "~"

        print(f"{gene:<12} {ctrl:>12.1f} {kd:>12.1f} {fc:>10.2f} {direction:>12}")

    # Load ATAC-seq if available
    print("\n" + "="*80)
    print("ATAC-seq CHROMATIN ACCESSIBILITY")
    print("="*80)

    atac_files = {
        'T47D_baseline': 'GSM8036152_T47D_summits.bed.gz',
        'T47D_siCtrl': 'GSM8036178_202340107A-T47D-siCtrl_S6_summits.bed.gz',
        'T47D_siFOXA1': 'GSM8036179_202340107A-T47D-siFOXA1_S7_summits.bed.gz',
    }

    atac = {}
    for name, fname in atac_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            atac[name] = load_atac_peaks(fpath)
            print(f"  {name}: {len(atac[name])} peaks")

    if atac:
        print(f"\n{'Gene':<12} {'siCtrl peaks':>15} {'siFOXA1 peaks':>15} {'Δ peaks':>12}")
        print("-"*60)

        for gene, (chrom, start, end) in GENE_COORDS.items():
            ctrl_peaks = count_peaks_at_gene(atac.get('T47D_siCtrl'), chrom, start, end)
            kd_peaks = count_peaks_at_gene(atac.get('T47D_siFOXA1'), chrom, start, end)
            delta = kd_peaks - ctrl_peaks

            print(f"{gene:<12} {ctrl_peaks:>15} {kd_peaks:>15} {delta:>+12}")

    # Interpretation
    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)

    print("""
    If chaperones go UP with FOXA1 knockdown:
    → FOXA1 is required for ER-mediated repression
    → Supports: FOXA1 opens chromatin → ER binds → repression
    → FOXA1 KD closes chromatin → ER can't bind → derepression

    If chaperones go DOWN with FOXA1 knockdown:
    → FOXA1 directly or indirectly activates chaperones
    → Does not support the ER repression model

    If chaperones don't change:
    → FOXA1 is not involved in chaperone regulation
    → The mechanism is independent of FOXA1
    """)


if __name__ == "__main__":
    main()
