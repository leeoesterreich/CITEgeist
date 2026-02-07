#!/usr/bin/env python
"""
Analyze FOXA1 and ER knockdown/overexpression in MCF7.

Dataset GSE75329:
- MCF7L-P-siControl: Parental MCF7 control
- MCF7L-P-siFoxA1-Seq1/2: FOXA1 knockdown
- MCF7L-P-siER: ER knockdown
- MCF7L-FoxA1-Dox/+Dox: FOXA1 overexpression

Key questions:
1. Does FOXA1 KD in MCF7 affect chaperones same as T47D?
2. Does ER KD in MCF7 affect chaperones?
3. Does FOXA1 overexpression affect chaperones?
"""

import os
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"


def load_data():
    df = pd.read_csv(os.path.join(DATA_DIR, "GSE75329_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
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
        'GSM1953112': 'TamR_siCtrl',
        'GSM1953115': 'TamR_siFOXA1_1',
        'GSM1953117': 'TamR_siFOXA1_2',
        'GSM1953119': 'TamR_siER',
        'GSM1953120': 'TamR_siIL8',
        'GSM1953122': 'P_siCtrl',
        'GSM1953124': 'P_siFOXA1_1',
        'GSM1953126': 'P_siFOXA1_2',
        'GSM1953127': 'P_siER',
        'GSM1953129': 'FOXA1_noDox',
        'GSM1953131': 'FOXA1_Dox',
    })

    return df


def main():
    print("="*80)
    print("MCF7 FOXA1/ER PERTURBATION ANALYSIS (GSE75329)")
    print("="*80)

    df = load_data()
    print(f"\nSamples: {df.columns.tolist()}")

    # Verify knockdowns worked
    print("\n" + "="*80)
    print("VERIFICATION: Did knockdowns work?")
    print("="*80)

    if 'FOXA1' in df.index:
        print("\nFOXA1 expression:")
        print(f"  P_siCtrl:     {df.loc['FOXA1', 'P_siCtrl']:.1f} TPM")
        print(f"  P_siFOXA1_1:  {df.loc['FOXA1', 'P_siFOXA1_1']:.1f} TPM (FC={df.loc['FOXA1', 'P_siFOXA1_1']/df.loc['FOXA1', 'P_siCtrl']:.2f})")
        print(f"  P_siFOXA1_2:  {df.loc['FOXA1', 'P_siFOXA1_2']:.1f} TPM (FC={df.loc['FOXA1', 'P_siFOXA1_2']/df.loc['FOXA1', 'P_siCtrl']:.2f})")
        print(f"  FOXA1_noDox:  {df.loc['FOXA1', 'FOXA1_noDox']:.1f} TPM")
        print(f"  FOXA1_Dox:    {df.loc['FOXA1', 'FOXA1_Dox']:.1f} TPM (FC={df.loc['FOXA1', 'FOXA1_Dox']/df.loc['FOXA1', 'FOXA1_noDox']:.2f})")

    if 'ESR1' in df.index:
        print("\nESR1 expression:")
        print(f"  P_siCtrl:  {df.loc['ESR1', 'P_siCtrl']:.1f} TPM")
        print(f"  P_siER:    {df.loc['ESR1', 'P_siER']:.1f} TPM (FC={df.loc['ESR1', 'P_siER']/df.loc['ESR1', 'P_siCtrl']:.2f})")

    # Chaperones with FOXA1 knockdown
    print("\n" + "="*80)
    print("MCF7 PARENTAL: FOXA1 KNOCKDOWN → CHAPERONES")
    print("="*80)
    print("Prediction: If MCF7 is saturated, FOXA1 KD might have LESS effect than T47D")

    chaperones = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4', 'PDIA6', 'SEC61A1']

    print(f"\n{'Gene':<12} {'siCtrl':>10} {'siFOXA1_1':>12} {'FC1':>8} {'siFOXA1_2':>12} {'FC2':>8}")
    print("-"*70)

    foxa1_up = 0
    foxa1_down = 0

    for gene in chaperones:
        if gene not in df.index:
            continue
        ctrl = df.loc[gene, 'P_siCtrl']
        kd1 = df.loc[gene, 'P_siFOXA1_1']
        kd2 = df.loc[gene, 'P_siFOXA1_2']
        fc1 = kd1 / ctrl if ctrl > 0 else 0
        fc2 = kd2 / ctrl if ctrl > 0 else 0
        avg_fc = (fc1 + fc2) / 2

        if avg_fc > 1.1:
            foxa1_up += 1
        elif avg_fc < 0.9:
            foxa1_down += 1

        print(f"{gene:<12} {ctrl:>10.1f} {kd1:>12.1f} {fc1:>8.2f} {kd2:>12.1f} {fc2:>8.2f}")

    print(f"\nWith FOXA1 KD: {foxa1_up} UP, {foxa1_down} DOWN")

    # Chaperones with ER knockdown
    print("\n" + "="*80)
    print("MCF7 PARENTAL: ER KNOCKDOWN → CHAPERONES")
    print("="*80)
    print("Prediction: If ER represses chaperones, ER KD → chaperones UP")

    print(f"\n{'Gene':<12} {'siCtrl':>10} {'siER':>12} {'FC':>8} {'Direction':>12}")
    print("-"*60)

    er_up = 0
    er_down = 0

    for gene in chaperones:
        if gene not in df.index:
            continue
        ctrl = df.loc[gene, 'P_siCtrl']
        kd = df.loc[gene, 'P_siER']
        fc = kd / ctrl if ctrl > 0 else 0
        direction = "UP" if fc > 1.1 else "DOWN" if fc < 0.9 else "~"

        if fc > 1.1:
            er_up += 1
        elif fc < 0.9:
            er_down += 1

        print(f"{gene:<12} {ctrl:>10.1f} {kd:>12.1f} {fc:>8.2f} {direction:>12}")

    print(f"\nWith ER KD: {er_up} UP, {er_down} DOWN")

    # Chaperones with FOXA1 overexpression
    print("\n" + "="*80)
    print("MCF7: FOXA1 OVEREXPRESSION → CHAPERONES")
    print("="*80)
    print("Prediction: If FOXA1 enables ER repression, FOXA1 OE → more repression → chaperones DOWN")

    print(f"\n{'Gene':<12} {'noDox':>10} {'+Dox':>12} {'FC':>8} {'Direction':>12}")
    print("-"*60)

    oe_up = 0
    oe_down = 0

    for gene in chaperones:
        if gene not in df.index:
            continue
        nodox = df.loc[gene, 'FOXA1_noDox']
        dox = df.loc[gene, 'FOXA1_Dox']
        fc = dox / nodox if nodox > 0 else 0
        direction = "UP" if fc > 1.1 else "DOWN" if fc < 0.9 else "~"

        if fc > 1.1:
            oe_up += 1
        elif fc < 0.9:
            oe_down += 1

        print(f"{gene:<12} {nodox:>10.1f} {dox:>12.1f} {fc:>8.2f} {direction:>12}")

    print(f"\nWith FOXA1 OE: {oe_up} UP, {oe_down} DOWN")

    # Compare MCF7 vs T47D FOXA1 KD effects
    print("\n" + "="*80)
    print("COMPARISON: MCF7 vs T47D FOXA1 KNOCKDOWN")
    print("="*80)

    # Load T47D data for comparison
    t47d_df = pd.read_csv(os.path.join(DATA_DIR, "GSE254218_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
                          sep='\t', compression='gzip', index_col=0)
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        results = mg.querymany(t47d_df.index.astype(str).tolist(), scopes='entrezgene',
                              fields='symbol', species='human', returnall=True)
        id_to_sym = {str(h['query']): h['symbol'] for h in results['out'] if 'symbol' in h}
        t47d_df.index = t47d_df.index.astype(str).map(lambda x: id_to_sym.get(x, x))
    except:
        pass

    t47d_df = t47d_df.rename(columns={
        'GSM8036191': 'T47D_siCtrl',
        'GSM8036192': 'T47D_siFOXA1',
    })

    print(f"\n{'Gene':<12} {'MCF7 FC':>12} {'T47D FC':>12} {'Same direction?':>18}")
    print("-"*60)

    same_direction = 0
    for gene in chaperones:
        if gene not in df.index or gene not in t47d_df.index:
            continue

        mcf7_ctrl = df.loc[gene, 'P_siCtrl']
        mcf7_kd = (df.loc[gene, 'P_siFOXA1_1'] + df.loc[gene, 'P_siFOXA1_2']) / 2
        mcf7_fc = mcf7_kd / mcf7_ctrl if mcf7_ctrl > 0 else 1

        t47d_ctrl = t47d_df.loc[gene, 'T47D_siCtrl']
        t47d_kd = t47d_df.loc[gene, 'T47D_siFOXA1']
        t47d_fc = t47d_kd / t47d_ctrl if t47d_ctrl > 0 else 1

        same = "YES" if (mcf7_fc > 1) == (t47d_fc > 1) else "NO"
        if same == "YES":
            same_direction += 1

        print(f"{gene:<12} {mcf7_fc:>12.2f} {t47d_fc:>12.2f} {same:>18}")

    print(f"\nSame direction: {same_direction}/{len(chaperones)}")

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    print(f"""
    MCF7 FOXA1 knockdown:
    - Chaperones: {foxa1_up} UP, {foxa1_down} DOWN

    MCF7 ER knockdown:
    - Chaperones: {er_up} UP, {er_down} DOWN

    MCF7 FOXA1 overexpression:
    - Chaperones: {oe_up} UP, {oe_down} DOWN

    Interpretation:
    - If ER represses chaperones: ER KD should make them go UP
    - If FOXA1 enables ER repression: FOXA1 KD should make them go UP (derepression)
    - If FOXA1 enables ER repression: FOXA1 OE should make them go DOWN (more repression)
    """)


if __name__ == "__main__":
    main()
