#!/usr/bin/env python
"""
What drives the opposite D538G effects between MCF7 and T47D?
Look at baseline differences that could explain it.
"""

import os
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"

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
    print("WHAT DRIVES THE OPPOSITE D538G EFFECTS?")
    print("="*80)

    tpm = load_tpm()

    groups = {
        'MCF7_WT': ['GSM2392606', 'GSM2392607', 'GSM2392608', 'GSM2392609'],
        'T47D_WT': ['GSM2392582', 'GSM2392583', 'GSM2392584', 'GSM2392585'],
    }

    # =========================================================================
    # Key factors that determine ER activity and chromatin accessibility
    # =========================================================================
    print("\n" + "="*80)
    print("BASELINE DIFFERENCES: ER SYSTEM COMPONENTS")
    print("="*80)

    er_system = {
        'ESR1': 'Estrogen receptor alpha - the receptor itself',
        'ESR2': 'Estrogen receptor beta',
        'FOXA1': 'Pioneer factor - opens chromatin for ER binding',
        'GATA3': 'ER cofactor - cooperates with ER',
        'PGR': 'Progesterone receptor - ER target, marker of ER activity',
        'NCOA1': 'SRC-1 coactivator',
        'NCOA2': 'SRC-2/TIF2 coactivator',
        'NCOA3': 'SRC-3/AIB1 coactivator',
        'NCOR1': 'Corepressor',
        'NCOR2': 'SMRT corepressor',
        'EP300': 'p300 - histone acetyltransferase, coactivator',
        'CREBBP': 'CBP - coactivator',
        'HDAC1': 'Histone deacetylase - corepressor',
        'HDAC2': 'Histone deacetylase - corepressor',
    }

    print(f"\n{'Gene':<10} {'MCF7-WT':>12} {'T47D-WT':>12} {'Ratio':>10} {'Higher':>10} {'Role'}")
    print("-"*80)

    for gene, role in er_system.items():
        if gene not in tpm.index:
            continue
        mcf7 = tpm.loc[gene, groups['MCF7_WT']].mean()
        t47d = tpm.loc[gene, groups['T47D_WT']].mean()
        ratio = mcf7 / t47d if t47d > 0 else float('inf')
        higher = "MCF7" if ratio > 1.5 else "T47D" if ratio < 0.67 else "~"

        print(f"{gene:<10} {mcf7:>12.1f} {t47d:>12.1f} {ratio:>10.2f} {higher:>10} {role[:30]}")

    # =========================================================================
    # The hypothesis
    # =========================================================================
    print("\n" + "="*80)
    print("THE HYPOTHESIS: WHY OPPOSITE EFFECTS?")
    print("="*80)

    print("""
    Key baseline differences:

    1. ESR1 (ER itself):
       - MCF7: HIGH (58.8 TPM)
       - T47D: LOW (12.8 TPM)
       - MCF7 has 4.6x more ER

    2. FOXA1 (pioneer factor):
       - MCF7: 96.7 TPM
       - T47D: 170.6 TPM
       - T47D has 1.8x more FOXA1

    3. Global ER binding at baseline (WT-E2):
       - MCF7: 12,472 peaks (HIGH)
       - T47D: 1,724 peaks (LOW)
       - MCF7 has 7.2x more binding sites

    This suggests:
    - MCF7 is "ER-saturated" - high ER, many binding sites
    - T47D is "ER-limited" - low ER, few binding sites but more open chromatin (FOXA1)
    """)

    # =========================================================================
    # The mechanism
    # =========================================================================
    print("\n" + "="*80)
    print("PROPOSED MECHANISM")
    print("="*80)

    print("""
    D538G mutation effects:

    IN MCF7 (ER-saturated system):
    ┌─────────────────────────────────────────────────────────┐
    │ Baseline: High ER, 12,472 binding sites                 │
    │                                                         │
    │ D538G mutation:                                         │
    │ → Altered ER conformation                               │
    │ → Reduced affinity for some binding sites               │
    │ → LOSES 7,069 peaks (57% reduction)                     │
    │ → Loss of binding at HSP90B1 → derepression → UP        │
    └─────────────────────────────────────────────────────────┘

    IN T47D (ER-limited system):
    ┌─────────────────────────────────────────────────────────┐
    │ Baseline: Low ER, only 1,724 binding sites              │
    │ But: High FOXA1 = more open chromatin available         │
    │                                                         │
    │ D538G mutation:                                         │
    │ → Constitutive ER activity (doesn't need ligand)        │
    │ → Can now access FOXA1-opened sites                     │
    │ → GAINS 7,828 peaks (454% increase)                     │
    │ → New binding at HSP90B1 → repression → DOWN            │
    └─────────────────────────────────────────────────────────┘

    The key insight:
    - Same mutation, different cellular context
    - MCF7: Already maximally bound → D538G causes unbinding
    - T47D: Under-bound, open chromatin waiting → D538G causes new binding
    """)

    # =========================================================================
    # Can we prove this?
    # =========================================================================
    print("\n" + "="*80)
    print("CAN WE PROVE THIS?")
    print("="*80)

    print("""
    With only 2 cell lines, we CANNOT statistically prove causation.

    What we CAN show:
    ✓ MCF7 has 4.6x more ESR1 than T47D
    ✓ T47D has 1.8x more FOXA1 than MCF7
    ✓ MCF7-WT has 7.2x more ER binding sites than T47D-WT
    ✓ D538G causes opposite binding changes
    ✓ Binding changes correlate with expression changes

    What we CANNOT prove without additional experiments:
    ✗ That ESR1 levels cause the different response
    ✗ That FOXA1 levels determine chromatin accessibility
    ✗ The causal mechanism

    To prove this would require:
    - ESR1 knockdown/overexpression in both cell lines
    - FOXA1 knockdown/overexpression
    - ATAC-seq to measure chromatin accessibility
    - More cell lines with varying ESR1/FOXA1 ratios
    """)

    # =========================================================================
    # Summary table
    # =========================================================================
    print("\n" + "="*80)
    print("SUMMARY: MCF7 vs T47D")
    print("="*80)

    print(f"""
    {'Feature':<30} {'MCF7':<20} {'T47D':<20}
    {'-'*70}
    {'ESR1 expression':<30} {'HIGH (58.8)':<20} {'LOW (12.8)':<20}
    {'FOXA1 expression':<30} {'MODERATE (96.7)':<20} {'HIGH (170.6)':<20}
    {'WT ER binding sites':<30} {'HIGH (12,472)':<20} {'LOW (1,724)':<20}
    {'D538G effect on binding':<30} {'LOSES 57%':<20} {'GAINS 454%':<20}
    {'D538G effect on HSP90B1 bind':<30} {'LOST':<20} {'GAINED':<20}
    {'D538G effect on HSP90B1 expr':<30} {'UP (FC=1.57)':<20} {'DOWN (FC=0.68)':<20}
    {'D538G effect on chaperones':<30} {'UP':<20} {'DOWN':<20}
    {'D538G effect on secretion':<30} {'UP':<20} {'DOWN':<20}

    Proposed explanation:
    MCF7 = "ER-saturated" → D538G causes unbinding → derepression
    T47D = "ER-limited" → D538G enables new binding → repression
    """)


if __name__ == "__main__":
    main()
