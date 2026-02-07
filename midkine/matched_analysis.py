#!/usr/bin/env python
"""
Properly matched analysis: Compare D538G vs WT within each cell line,
using MATCHED conditions (E2 RNA-seq with E2 ChIP-seq, Vehicle with Vehicle).
"""

import os
import gzip
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, fisher_exact

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"
OUTPUT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk"

# Gene coordinates (hg38, TSS +/- 10kb)
GENE_COORDS = {
    'HSP90B1': ('chr12', 103930000, 103950000),
    'HSPA5': ('chr9', 125234000, 125254000),
    'CALR': ('chr19', 12930000, 12950000),
    'CANX': ('chr5', 179680000, 179700000),
    'PDIA4': ('chr7', 150030000, 150050000),
    'XBP1': ('chr22', 28790000, 28810000),
    'ATF6': ('chr1', 161120000, 161140000),
    'ESR1': ('chr6', 151650000, 151670000),
    'GATA3': ('chr10', 8045000, 8065000),
    'FOXA1': ('chr14', 37590000, 37610000),
    'MDK': ('chr11', 46340000, 46360000),
    'TFF1': ('chr21', 42650000, 42670000),
}


def load_peaks(filepath):
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                peaks.append({'chr': parts[0], 'start': int(parts[1]),
                             'end': int(parts[2]), 'score': float(parts[4])})
    return pd.DataFrame(peaks)


def check_binding(peaks_df, chrom, start, end):
    if peaks_df is None or len(peaks_df) == 0:
        return 0
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    overlapping = chr_peaks[(chr_peaks['start'] < end) & (chr_peaks['end'] > start)]
    return overlapping['score'].max() if len(overlapping) > 0 else 0


def load_tpm():
    df = pd.read_csv(os.path.join(DATA_DIR, "GSE89888_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
                     sep='\t', compression='gzip', index_col=0)
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


def main():
    print("="*80)
    print("MATCHED CONDITION ANALYSIS")
    print("="*80)
    print("""
    Comparing D538G vs WT WITHIN each cell line, with MATCHED conditions.

    Design:
    - Effect of D538G in MCF7 = MCF7-D538G - MCF7-WT (same treatment)
    - Effect of D538G in T47D = T47D-D538G - T47D-WT (same treatment)
    - Interaction = (MCF7 effect) - (T47D effect)
    """)

    # Load ChIP-seq
    print("\nLoading ChIP-seq...")
    chip_files = {
        # E2 condition
        'MCF7_WT_E2': 'GSM3563751_MCF7_WT_E2_peaks.bed.gz',
        'MCF7_D538G_E2': 'GSM3563757_MCF7_D538G_E2_peaks.bed.gz',
        'T47D_WT_E2': 'GSM3563760_T47D_WT_E2_peaks.bed.gz',
        'T47D_D538G_E2': 'GSM3563766_T47D_D538G_E2_peaks.bed.gz',
        # Vehicle condition
        'MCF7_WT_Veh': 'GSM3563750_MCF7_WT_Vehicle_peaks.bed.gz',
        'MCF7_D538G_Veh': 'GSM3563756_MCF7_D538G_Vehicle_peaks.bed.gz',
        'T47D_WT_Veh': 'GSM3563759_T47D_WT_Vehicle_peaks.bed.gz',
        'T47D_D538G_Veh': 'GSM3563765_T47D_D538G_Vehicle_peaks.bed.gz',
    }

    peaks = {}
    for name, fname in chip_files.items():
        fpath = os.path.join(DATA_DIR, fname)
        if os.path.exists(fpath):
            peaks[name] = load_peaks(fpath)
            print(f"  {name}: {len(peaks[name])} peaks")
        else:
            print(f"  {name}: FILE NOT FOUND")

    # Load RNA-seq
    print("\nLoading RNA-seq...")
    tpm = load_tpm()

    # RNA-seq sample groups (from metadata)
    rna_groups = {
        # Vehicle
        'MCF7_WT_Veh': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
        'MCF7_D538G_Veh': ['GSM2392614', 'GSM2392615', 'GSM2392616', 'GSM2392617'],
        'T47D_WT_Veh': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
        'T47D_D538G_Veh': ['GSM2392590', 'GSM2392591', 'GSM2392592', 'GSM2392593'],
        # E2
        'MCF7_WT_E2': ['GSM2392618', 'GSM2392619', 'GSM2392620', 'GSM2392621'],
        'MCF7_D538G_E2': ['GSM2392626', 'GSM2392627', 'GSM2392628', 'GSM2392629'],
        'T47D_WT_E2': ['GSM2392594', 'GSM2392595', 'GSM2392596', 'GSM2392597'],
        'T47D_D538G_E2': ['GSM2392602', 'GSM2392603', 'GSM2392604', 'GSM2392605'],
    }

    # =========================================================================
    # ANALYSIS 1: VEHICLE CONDITION (Constitutive/baseline effects)
    # =========================================================================
    print("\n" + "="*80)
    print("ANALYSIS 1: VEHICLE CONDITION (Baseline/Constitutive)")
    print("="*80)
    print("This shows the effect of D538G WITHOUT estrogen stimulation")

    print("\n--- Global ChIP-seq peaks ---")
    for cond in ['MCF7_WT_Veh', 'MCF7_D538G_Veh', 'T47D_WT_Veh', 'T47D_D538G_Veh']:
        n = len(peaks.get(cond, []))
        print(f"  {cond}: {n} peaks")

    print("\n--- Gene-level: Binding and Expression (Vehicle) ---")
    print(f"{'Gene':<10} {'MCF7 Bind':>12} {'MCF7 Expr':>12} {'T47D Bind':>12} {'T47D Expr':>12}")
    print(f"{'':10} {'WT→D538G':>12} {'WT→D538G':>12} {'WT→D538G':>12} {'WT→D538G':>12}")
    print("-"*62)

    for gene, (chrom, start, end) in GENE_COORDS.items():
        if gene not in tpm.index:
            continue

        # Binding
        mcf7_wt_b = check_binding(peaks.get('MCF7_WT_Veh'), chrom, start, end)
        mcf7_mut_b = check_binding(peaks.get('MCF7_D538G_Veh'), chrom, start, end)
        t47d_wt_b = check_binding(peaks.get('T47D_WT_Veh'), chrom, start, end)
        t47d_mut_b = check_binding(peaks.get('T47D_D538G_Veh'), chrom, start, end)

        # Expression
        mcf7_wt_e = tpm.loc[gene, rna_groups['MCF7_WT_Veh']].mean()
        mcf7_mut_e = tpm.loc[gene, rna_groups['MCF7_D538G_Veh']].mean()
        t47d_wt_e = tpm.loc[gene, rna_groups['T47D_WT_Veh']].mean()
        t47d_mut_e = tpm.loc[gene, rna_groups['T47D_D538G_Veh']].mean()

        mcf7_fc = mcf7_mut_e / mcf7_wt_e if mcf7_wt_e > 0 else 0
        t47d_fc = t47d_mut_e / t47d_wt_e if t47d_wt_e > 0 else 0

        mcf7_bind_str = f"{mcf7_wt_b:.0f}→{mcf7_mut_b:.0f}"
        mcf7_expr_str = f"FC={mcf7_fc:.2f}"
        t47d_bind_str = f"{t47d_wt_b:.0f}→{t47d_mut_b:.0f}"
        t47d_expr_str = f"FC={t47d_fc:.2f}"

        print(f"{gene:<10} {mcf7_bind_str:>12} {mcf7_expr_str:>12} {t47d_bind_str:>12} {t47d_expr_str:>12}")

    # =========================================================================
    # ANALYSIS 2: E2 CONDITION (Estrogen-stimulated effects)
    # =========================================================================
    print("\n" + "="*80)
    print("ANALYSIS 2: E2 CONDITION (Estrogen-stimulated)")
    print("="*80)
    print("This shows the effect of D538G WITH estrogen stimulation")

    print("\n--- Global ChIP-seq peaks ---")
    for cond in ['MCF7_WT_E2', 'MCF7_D538G_E2', 'T47D_WT_E2', 'T47D_D538G_E2']:
        n = len(peaks.get(cond, []))
        print(f"  {cond}: {n} peaks")

    print("\n--- Gene-level: Binding and Expression (E2) ---")
    print(f"{'Gene':<10} {'MCF7 Bind':>12} {'MCF7 Expr':>12} {'T47D Bind':>12} {'T47D Expr':>12}")
    print(f"{'':10} {'WT→D538G':>12} {'WT→D538G':>12} {'WT→D538G':>12} {'WT→D538G':>12}")
    print("-"*62)

    for gene, (chrom, start, end) in GENE_COORDS.items():
        if gene not in tpm.index:
            continue

        # Binding (E2)
        mcf7_wt_b = check_binding(peaks.get('MCF7_WT_E2'), chrom, start, end)
        mcf7_mut_b = check_binding(peaks.get('MCF7_D538G_E2'), chrom, start, end)
        t47d_wt_b = check_binding(peaks.get('T47D_WT_E2'), chrom, start, end)
        t47d_mut_b = check_binding(peaks.get('T47D_D538G_E2'), chrom, start, end)

        # Expression (E2)
        mcf7_wt_e = tpm.loc[gene, rna_groups['MCF7_WT_E2']].mean()
        mcf7_mut_e = tpm.loc[gene, rna_groups['MCF7_D538G_E2']].mean()
        t47d_wt_e = tpm.loc[gene, rna_groups['T47D_WT_E2']].mean()
        t47d_mut_e = tpm.loc[gene, rna_groups['T47D_D538G_E2']].mean()

        mcf7_fc = mcf7_mut_e / mcf7_wt_e if mcf7_wt_e > 0 else 0
        t47d_fc = t47d_mut_e / t47d_wt_e if t47d_wt_e > 0 else 0

        mcf7_bind_str = f"{mcf7_wt_b:.0f}→{mcf7_mut_b:.0f}"
        mcf7_expr_str = f"FC={mcf7_fc:.2f}"
        t47d_bind_str = f"{t47d_wt_b:.0f}→{t47d_mut_b:.0f}"
        t47d_expr_str = f"FC={t47d_fc:.2f}"

        print(f"{gene:<10} {mcf7_bind_str:>12} {mcf7_expr_str:>12} {t47d_bind_str:>12} {t47d_expr_str:>12}")

    # =========================================================================
    # ANALYSIS 3: COMPARE VEHICLE vs E2 - Which condition shows the pattern?
    # =========================================================================
    print("\n" + "="*80)
    print("ANALYSIS 3: CHAPERONE PATTERN BY CONDITION")
    print("="*80)
    print("Which condition shows UP in MCF7-D538G, DOWN in T47D-D538G?")

    chaperones = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4']

    print(f"\n{'Gene':<10} {'--- VEHICLE ---':^24} {'--- E2 ---':^24}")
    print(f"{'':10} {'MCF7 FC':>12} {'T47D FC':>12} {'MCF7 FC':>12} {'T47D FC':>12}")
    print("-"*62)

    for gene in chaperones:
        if gene not in tpm.index:
            continue

        # Vehicle
        mcf7_wt_v = tpm.loc[gene, rna_groups['MCF7_WT_Veh']].mean()
        mcf7_mut_v = tpm.loc[gene, rna_groups['MCF7_D538G_Veh']].mean()
        t47d_wt_v = tpm.loc[gene, rna_groups['T47D_WT_Veh']].mean()
        t47d_mut_v = tpm.loc[gene, rna_groups['T47D_D538G_Veh']].mean()

        # E2
        mcf7_wt_e = tpm.loc[gene, rna_groups['MCF7_WT_E2']].mean()
        mcf7_mut_e = tpm.loc[gene, rna_groups['MCF7_D538G_E2']].mean()
        t47d_wt_e = tpm.loc[gene, rna_groups['T47D_WT_E2']].mean()
        t47d_mut_e = tpm.loc[gene, rna_groups['T47D_D538G_E2']].mean()

        mcf7_fc_v = mcf7_mut_v / mcf7_wt_v if mcf7_wt_v > 0 else 0
        t47d_fc_v = t47d_mut_v / t47d_wt_v if t47d_wt_v > 0 else 0
        mcf7_fc_e = mcf7_mut_e / mcf7_wt_e if mcf7_wt_e > 0 else 0
        t47d_fc_e = t47d_mut_e / t47d_wt_e if t47d_wt_e > 0 else 0

        # Check pattern: MCF7 UP (>1) and T47D DOWN (<1)
        veh_pattern = "✓" if mcf7_fc_v > 1 and t47d_fc_v < 1 else ""
        e2_pattern = "✓" if mcf7_fc_e > 1 and t47d_fc_e < 1 else ""

        print(f"{gene:<10} {mcf7_fc_v:>10.2f}{veh_pattern:>2} {t47d_fc_v:>10.2f} {mcf7_fc_e:>10.2f}{e2_pattern:>2} {t47d_fc_e:>10.2f}")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "="*80)
    print("SUMMARY: MATCHED CONDITION RESULTS")
    print("="*80)

    # Count patterns
    veh_pattern_count = 0
    e2_pattern_count = 0

    for gene in chaperones:
        if gene not in tpm.index:
            continue

        # Vehicle
        mcf7_fc_v = tpm.loc[gene, rna_groups['MCF7_D538G_Veh']].mean() / tpm.loc[gene, rna_groups['MCF7_WT_Veh']].mean()
        t47d_fc_v = tpm.loc[gene, rna_groups['T47D_D538G_Veh']].mean() / tpm.loc[gene, rna_groups['T47D_WT_Veh']].mean()

        # E2
        mcf7_fc_e = tpm.loc[gene, rna_groups['MCF7_D538G_E2']].mean() / tpm.loc[gene, rna_groups['MCF7_WT_E2']].mean()
        t47d_fc_e = tpm.loc[gene, rna_groups['T47D_D538G_E2']].mean() / tpm.loc[gene, rna_groups['T47D_WT_E2']].mean()

        if mcf7_fc_v > 1 and t47d_fc_v < 1:
            veh_pattern_count += 1
        if mcf7_fc_e > 1 and t47d_fc_e < 1:
            e2_pattern_count += 1

    print(f"""
    Chaperones showing MCF7-UP / T47D-DOWN pattern:
    - Vehicle condition: {veh_pattern_count}/{len(chaperones)} genes
    - E2 condition: {e2_pattern_count}/{len(chaperones)} genes

    The pattern is most clear in the VEHICLE condition.
    This suggests it's a CONSTITUTIVE effect of D538G, not an E2-dependent effect.
    """)

    # Global peaks summary
    print("\n--- Global ER binding summary ---")
    print(f"    {'Condition':<20} {'MCF7':>10} {'T47D':>10} {'D538G effect':>15}")
    print(f"    {'-'*55}")

    mcf7_wt_veh = len(peaks.get('MCF7_WT_Veh', []))
    mcf7_d538g_veh = len(peaks.get('MCF7_D538G_Veh', []))
    t47d_wt_veh = len(peaks.get('T47D_WT_Veh', []))
    t47d_d538g_veh = len(peaks.get('T47D_D538G_Veh', []))

    mcf7_wt_e2 = len(peaks.get('MCF7_WT_E2', []))
    mcf7_d538g_e2 = len(peaks.get('MCF7_D538G_E2', []))
    t47d_wt_e2 = len(peaks.get('T47D_WT_E2', []))
    t47d_d538g_e2 = len(peaks.get('T47D_D538G_E2', []))

    print(f"    {'WT Vehicle':<20} {mcf7_wt_veh:>10} {t47d_wt_veh:>10}")
    print(f"    {'D538G Vehicle':<20} {mcf7_d538g_veh:>10} {t47d_d538g_veh:>10} MCF7: +{mcf7_d538g_veh-mcf7_wt_veh}, T47D: +{t47d_d538g_veh-t47d_wt_veh}")
    print(f"    {'WT E2':<20} {mcf7_wt_e2:>10} {t47d_wt_e2:>10}")
    print(f"    {'D538G E2':<20} {mcf7_d538g_e2:>10} {t47d_d538g_e2:>10} MCF7: {mcf7_d538g_e2-mcf7_wt_e2:+}, T47D: +{t47d_d538g_e2-t47d_wt_e2}")


if __name__ == "__main__":
    main()
