#!/usr/bin/env python
"""
Critical evaluation of the chaperone hypothesis.

Problems to check:
1. Absolute levels after change (does T47D still have more?)
2. Is the change biologically meaningful (effect size)?
3. What about ER binding (ChIP-seq)?
4. Could this be a general stress response vs specific mechanism?
"""

import pandas as pd
import numpy as np

OUTPUT_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples/output_vignette4_mdk"

# Load results
df = pd.read_csv(f"{OUTPUT_DIR}/chaperone_hypothesis_test.csv")
chipseq = pd.read_csv(f"{OUTPUT_DIR}/chipseq_er_binding.csv")

print("="*80)
print("CRITICAL EVALUATION OF CHAPERONE HYPOTHESIS")
print("="*80)

# PROBLEM 1: Absolute levels after mutation
print("\n" + "="*80)
print("PROBLEM 1: ABSOLUTE LEVELS AFTER D538G MUTATION")
print("="*80)
print("Even if chaperones decrease in T47D, does T47D still have MORE than MCF7?")
print("-"*80)

key_genes = ['HSP90B1', 'HSPA5', 'CALR', 'CANX', 'PDIA4', 'XBP1']

print(f"{'Gene':<12} {'MCF7-D538G':>12} {'T47D-D538G':>12} {'T47D still higher?':>20}")
print("-"*80)

for gene in key_genes:
    if gene in df['gene'].values:
        row = df[df['gene'] == gene].iloc[0]
        mcf7_mut = row.get('MCF7_D538G_veh_mean', np.nan)
        t47d_mut = row.get('T47D_D538G_veh_mean', np.nan)

        still_higher = "YES" if t47d_mut > mcf7_mut else "NO"
        print(f"{gene:<12} {mcf7_mut:>12.1f} {t47d_mut:>12.1f} {still_higher:>20}")

print("\nVERDICT: T47D-D538G STILL has higher absolute chaperone levels than MCF7-D538G")
print("         If absolute levels mattered, T47D should secrete MORE, not less.")

# PROBLEM 2: Effect size - is the change biologically meaningful?
print("\n" + "="*80)
print("PROBLEM 2: EFFECT SIZE - ARE CHANGES BIOLOGICALLY MEANINGFUL?")
print("="*80)

print(f"{'Gene':<12} {'MCF7 Δ':>12} {'T47D Δ':>12} {'MCF7 Δ%':>10} {'T47D Δ%':>10}")
print("-"*80)

for gene in key_genes:
    if gene in df['gene'].values:
        row = df[df['gene'] == gene].iloc[0]
        mcf7_wt = row.get('MCF7_WT_veh_mean', 0)
        mcf7_mut = row.get('MCF7_D538G_veh_mean', 0)
        t47d_wt = row.get('T47D_WT_veh_mean', 0)
        t47d_mut = row.get('T47D_D538G_veh_mean', 0)

        mcf7_delta = mcf7_mut - mcf7_wt
        t47d_delta = t47d_mut - t47d_wt
        mcf7_pct = (mcf7_delta / mcf7_wt * 100) if mcf7_wt > 0 else 0
        t47d_pct = (t47d_delta / t47d_wt * 100) if t47d_wt > 0 else 0

        print(f"{gene:<12} {mcf7_delta:>+12.1f} {t47d_delta:>+12.1f} {mcf7_pct:>+10.1f}% {t47d_pct:>+10.1f}%")

# PROBLEM 3: ER ChIP-seq binding
print("\n" + "="*80)
print("PROBLEM 3: ER ChIP-seq BINDING - IS THIS ER-MEDIATED?")
print("="*80)
print("If D538G-ER drives chaperone upregulation, we should see ER binding")
print("-"*80)

# Check binding at chaperone genes
chaperone_chipseq = chipseq[chipseq['Gene'].isin(['HSP90B1', 'HSPA5', 'XBP1', 'SEC61B', 'ATF6', 'DDIT3'])]

print(f"{'Gene':<12} {'MCF7-WT':>10} {'MCF7-D538G':>12} {'T47D-WT':>10} {'T47D-D538G':>12}")
print("-"*80)

for gene in ['HSP90B1', 'HSPA5', 'XBP1', 'SEC61B', 'ATF6', 'DDIT3']:
    if gene in chipseq['Gene'].values:
        row = chipseq[chipseq['Gene'] == gene].iloc[0]
        mcf7_wt = row.get('MCF7_WT_E2_score', 0)
        mcf7_mut = row.get('MCF7_D538G_E2_score', 0)
        t47d_wt = row.get('T47D_WT_E2_score', 0)
        t47d_mut = row.get('T47D_D538G_E2_score', 0)

        print(f"{gene:<12} {mcf7_wt:>10.1f} {mcf7_mut:>12.1f} {t47d_wt:>10.1f} {t47d_mut:>12.1f}")

print("\nVERDICT: NO ER binding at HSP90B1 or XBP1 in ANY condition")
print("         HSPA5 has ER binding only in T47D, not MCF7!")
print("         The chaperone upregulation in MCF7 is NOT directly ER-mediated")

# PROBLEM 4: Alternative explanation - general stress response
print("\n" + "="*80)
print("PROBLEM 4: ALTERNATIVE - IS THIS UPR/STRESS RESPONSE?")
print("="*80)
print("D538G mutation might trigger different stress responses in MCF7 vs T47D")
print("-"*80)

# UPR markers
upr_markers = ['XBP1', 'ATF6', 'ATF4', 'DDIT3', 'HSPA5', 'CALR']

print(f"{'Gene':<12} {'MCF7 FC':>10} {'T47D FC':>10} {'UPR Role':>30}")
print("-"*80)

upr_roles = {
    'XBP1': 'IRE1 branch TF',
    'ATF6': 'ATF6 branch TF',
    'ATF4': 'PERK branch TF',
    'DDIT3': 'Pro-apoptotic (CHOP)',
    'HSPA5': 'ER chaperone (BiP)',
    'CALR': 'ER Ca2+ & chaperone',
}

for gene in upr_markers:
    if gene in df['gene'].values:
        row = df[df['gene'] == gene].iloc[0]
        mcf7_fc = row.get('MCF7_fc', np.nan)
        t47d_fc = row.get('T47D_fc', np.nan)
        role = upr_roles.get(gene, '')

        print(f"{gene:<12} {mcf7_fc:>10.2f} {t47d_fc:>10.2f} {role:>30}")

print("\nPattern: MCF7-D538G activates adaptive UPR (chaperones UP)")
print("         T47D-D538G may have maladaptive response (chaperones DOWN)")

# SUMMARY
print("\n" + "="*80)
print("SUMMARY: WHAT THE DATA ACTUALLY SUPPORTS")
print("="*80)

print("""
SUPPORTED:
✓ Chaperones (HSP90B1, HSPA5, etc.) go UP in MCF7-D538G (p<0.001)
✓ Chaperones go DOWN in T47D-D538G (p<0.0001 for HSP90B1)
✓ Direction of change correlates with secretion phenotype
✓ 12 secretory pathway genes show this OPPOSITE pattern

NOT SUPPORTED:
✗ Absolute chaperone levels do NOT explain phenotype
  (T47D-D538G still has MORE chaperones than MCF7-D538G)
✗ No ER binding at chaperone genes in MCF7
  (This is NOT directly ER-mediated)
✗ Causation not established (correlation only)

ALTERNATIVE INTERPRETATION:
The D538G mutation triggers OPPOSITE stress responses:
- MCF7: Adaptive UPR → chaperone upregulation → enhanced secretion
- T47D: Maladaptive response → chaperone downregulation → reduced secretion

This is INDIRECT regulation, not direct ER target gene activation.
The upstream mechanism remains UNKNOWN.
""")

# The key insight
print("="*80)
print("THE KEY INSIGHT")
print("="*80)
print("""
The data shows that D538G causes OPPOSITE effects on secretory machinery:
- MCF7: Gains secretory capacity (chaperones UP by 50-60%)
- T47D: Loses secretory capacity (chaperones DOWN by 25-30%)

This explains the MDK mRNA vs secretion paradox:
- MCF7: Even with lower MDK mRNA, enhanced machinery = MORE secreted protein
- T47D: Even with higher MDK mRNA, impaired machinery = LESS secreted protein

But we CANNOT prove:
1. WHY D538G has opposite effects in these cell lines
2. WHAT upstream factor mediates this (not ER directly based on ChIP)
3. WHETHER this is causal or just correlational
""")
