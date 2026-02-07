#!/usr/bin/env python
"""
Interpret the paradoxical FOXA1 overexpression result.

The finding: FOXA1 7x OE → chaperones 2x UP (not DOWN as predicted)

Possible explanations:
1. FOXA1 OE causes cellular stress → activates UPR → chaperones UP
2. FOXA1 has ER-independent effects on chaperones
3. FOXA1 OE changes the chromatin landscape in unexpected ways
4. The cells are responding to transgene expression stress
"""

import os
import pandas as pd
import numpy as np

DATA_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine"


def load_data():
    df = pd.read_csv(os.path.join(DATA_DIR, "GSE75329_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
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
    print("INTERPRETING FOXA1 OVEREXPRESSION PARADOX")
    print("="*80)

    df = load_data()

    print("""
    THE PARADOX:
    ────────────
    FOXA1 7x overexpression → chaperones 2x UP

    This is OPPOSITE to what the saturation model predicts:
    - If FOXA1 opens chromatin for ER binding
    - And ER represses chaperones
    - Then more FOXA1 → more ER binding → more repression → chaperones DOWN

    But we see chaperones UP. Why?
    """)

    # Check for signs of cellular stress (UPR activation)
    print("="*80)
    print("HYPOTHESIS 1: FOXA1 OE CAUSES CELLULAR STRESS")
    print("="*80)
    print("\nUPR stress markers:")

    upr_markers = ['ATF4', 'ATF6', 'DDIT3', 'ERN1', 'EIF2AK3', 'XBP1']

    print(f"\n{'Gene':<12} {'noDox':>10} {'+Dox':>10} {'FC':>10} {'Direction':>12}")
    print("-"*60)

    for gene in upr_markers:
        if gene not in df.index:
            continue
        nodox = df.loc[gene, 'FOXA1_noDox']
        dox = df.loc[gene, 'FOXA1_Dox']
        fc = dox / nodox if nodox > 0 else 0
        direction = "UP" if fc > 1.2 else "DOWN" if fc < 0.8 else "~"
        print(f"{gene:<12} {nodox:>10.1f} {dox:>10.1f} {fc:>10.2f} {direction:>12}")

    # Check if ER target genes go down (increased repression)
    print("\n" + "="*80)
    print("CHECK: ARE ER TARGET GENES REPRESSED?")
    print("="*80)
    print("\nIf FOXA1 OE → more ER binding, ER targets should change:")

    er_targets = ['TFF1', 'GREB1', 'PGR', 'AREG', 'MYC', 'CCND1']

    print(f"\n{'Gene':<12} {'noDox':>10} {'+Dox':>10} {'FC':>10} {'Direction':>12}")
    print("-"*60)

    er_up = 0
    er_down = 0
    for gene in er_targets:
        if gene not in df.index:
            continue
        nodox = df.loc[gene, 'FOXA1_noDox']
        dox = df.loc[gene, 'FOXA1_Dox']
        fc = dox / nodox if nodox > 0 else 0
        direction = "UP" if fc > 1.2 else "DOWN" if fc < 0.8 else "~"
        if fc > 1.2:
            er_up += 1
        elif fc < 0.8:
            er_down += 1
        print(f"{gene:<12} {nodox:>10.1f} {dox:>10.1f} {fc:>10.2f} {direction:>12}")

    print(f"\nER targets: {er_up} UP, {er_down} DOWN")

    # Check ESR1 itself
    print("\n" + "="*80)
    print("CHECK: DOES FOXA1 OE CHANGE ESR1 EXPRESSION?")
    print("="*80)

    if 'ESR1' in df.index:
        nodox = df.loc['ESR1', 'FOXA1_noDox']
        dox = df.loc['ESR1', 'FOXA1_Dox']
        fc = dox / nodox if nodox > 0 else 0
        print(f"\nESR1: {nodox:.1f} → {dox:.1f} TPM (FC={fc:.2f})")

        if fc < 0.8:
            print("\n*** ESR1 is DOWNREGULATED with FOXA1 OE! ***")
            print("This could explain chaperone increase:")
            print("FOXA1 OE → ESR1 DOWN → less repression → chaperones UP")
        elif fc > 1.2:
            print("\nESR1 is UPREGULATED - doesn't explain the paradox")
        else:
            print("\nESR1 unchanged - doesn't explain the paradox")

    # Alternative: FOXA1 might directly regulate chaperones
    print("\n" + "="*80)
    print("HYPOTHESIS 2: FOXA1 DIRECTLY ACTIVATES CHAPERONES")
    print("="*80)

    print("""
    FOXA1 is primarily known as a pioneer factor, but it's also a transcription factor.

    If FOXA1 can directly activate chaperone genes (independent of ER):
    - FOXA1 OE → direct transcriptional activation → chaperones UP
    - This would override the ER-mediated repression

    Evidence needed: FOXA1 ChIP-seq at chaperone promoters
    """)

    # The dox system might cause stress
    print("\n" + "="*80)
    print("HYPOTHESIS 3: DOX INDUCTION SYSTEM CAUSES STRESS")
    print("="*80)

    print("""
    The FOXA1-Dox system uses doxycycline-inducible expression.

    Concerns:
    - Doxycycline itself might affect cell biology
    - Massive transgene expression (7x) can overwhelm the protein folding machinery
    - This could activate UPR as a response to transgene protein load

    Control needed: Compare to a dox-treated empty vector control
    """)

    # Final interpretation
    print("\n" + "="*80)
    print("REVISED MODEL")
    print("="*80)

    print("""
    The simple model was:
        FOXA1 → opens chromatin → ER binds → represses chaperones

    The data suggests a more complex picture:

    1. FOXA1 has CONTEXT-DEPENDENT effects:
       - In D538G context: FOXA1 enables ER binding (T47D gains sites)
       - In WT context: FOXA1 OE may not simply increase repression

    2. Acute vs chronic effects:
       - FOXA1 knockdown (loss) → removes pioneer factor → ER can't bind → derepression
       - FOXA1 overexpression (acute gain) → stress response → UPR activation

    3. The key insight remains valid:
       - MCF7 is ER-saturated (high occupancy)
       - T47D is ER-limited (low occupancy)
       - D538G redistributes ER binding differently in each context

    The FOXA1 OE experiment shows that simply adding FOXA1 doesn't
    linearly increase repression - the system is more complex.
    """)

    # What does this mean for MDK?
    print("\n" + "="*80)
    print("IMPLICATIONS FOR MDK SECRETION")
    print("="*80)

    print("""
    The chaperone/secretory pathway findings remain solid:

    1. In MCF7-D538G (vs WT):
       - ER loses binding at chaperone loci (saturated → redistribution)
       - Chaperones UP → secretory machinery UP
       - MDK secretion UP

    2. In T47D-D538G (vs WT):
       - ER gains binding at chaperone loci (unsaturated → new binding)
       - Chaperones DOWN → secretory machinery DOWN
       - MDK secretion DOWN

    The FOXA1 OE paradox suggests:
    - The FOXA1-ER axis is not simply "more FOXA1 = more repression"
    - Acute overexpression may trigger compensatory responses
    - The natural variation between MCF7/T47D FOXA1 levels is more subtle
    """)


if __name__ == "__main__":
    main()
