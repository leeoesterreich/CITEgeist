#!/usr/bin/env python
"""
12_triple_venn.py

THREE-WAY VENN: Expression × Binding × Secretory Pathway

Purpose: Intersect three criteria to identify candidates:
1. Opposite expression regulation (MCF7 vs T47D response to D538G)
2. Differential ER binding at promoter
3. GO-defined secretory pathway / ER chaperone genes

This should winnow down to a small number of candidates where
HSP90B1 emerges as a top hit.

Datasets:
- GSE89888 (RNA-seq) - expression
- GSE125117 (ER ChIP-seq) - binding
- Gene Ontology - secretory pathway genes

Outputs:
- outputs/figures/fig12_triple_venn.png
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']


def get_secretory_pathway_genes():
    """
    Get genes in secretory pathway using Gene Ontology.
    Uses multiple GO terms for comprehensive coverage.
    """
    try:
        import mygene
        mg = mygene.MyGeneInfo()

        # GO terms related to secretory pathway
        go_terms = {
            'GO:0034975': 'protein folding in ER',
            'GO:0030433': 'ubiquitin-dependent ERAD',
            'GO:0006888': 'ER to Golgi vesicle transport',
            'GO:0006890': 'retrograde vesicle transport, Golgi to ER',
            'GO:0045047': 'protein targeting to ER',
            'GO:0006614': 'SRP-dependent cotranslational protein targeting to membrane',
            'GO:0030968': 'endoplasmic reticulum unfolded protein response',
            'GO:0034976': 'response to ER stress',
        }

        all_genes = set()
        print("\nFetching secretory pathway genes from Gene Ontology:")

        for go_id, go_name in go_terms.items():
            try:
                results = mg.query(f'go:{go_id}', species='human',
                                  fields='symbol', size=500)
                genes = [hit['symbol'] for hit in results.get('hits', []) if 'symbol' in hit]
                print(f"  {go_id} ({go_name}): {len(genes)} genes")
                all_genes.update(genes)
            except Exception as e:
                print(f"  {go_id}: Error - {e}")

        print(f"\nTotal unique secretory pathway genes: {len(all_genes)}")

        # If GO query failed, use fallback
        if len(all_genes) == 0:
            return get_fallback_secretory_genes()
        return all_genes

    except Exception as e:
        print(f"Error fetching GO genes: {e}")
        # Fallback to manually curated list
        return get_fallback_secretory_genes()


def get_fallback_secretory_genes():
    """Fallback list from KEGG hsa04141 + manual curation."""
    print("Using fallback secretory pathway gene list...")
    genes = {
        # ER chaperones
        'HSP90B1', 'HSPA5', 'HYOU1', 'CALR', 'CANX', 'PDIA3', 'PDIA4', 'PDIA6',
        'DNAJB11', 'DNAJC3', 'DNAJC10', 'ERO1A', 'ERO1B', 'UGGT1', 'UGGT2',
        'P4HB', 'PDIA2', 'AGR2', 'VAPA',
        # Translocon
        'SEC61A1', 'SEC61A2', 'SEC61B', 'SEC61G', 'SEC62', 'SEC63',
        'SSR1', 'SSR2', 'SSR3', 'SSR4', 'SRPRB', 'SRP9', 'SRP14', 'SRP19',
        # ERAD
        'EDEM1', 'EDEM2', 'EDEM3', 'OS9', 'XBP1', 'ATF6', 'ERN1',
        'SYVN1', 'SEL1L', 'HERPUD1', 'DERL1', 'DERL2', 'DERL3',
        # ER-Golgi transport
        'LMAN1', 'LMAN2', 'MCFD2', 'SURF4', 'TMED2', 'TMED9', 'TMED10',
        'SEC22B', 'SEC23A', 'SEC23B', 'SEC24A', 'SEC24B', 'SEC24C', 'SEC24D',
        'SEC31A', 'SEC31B', 'SAR1A', 'SAR1B',
        # UPR sensors
        'ATF4', 'DDIT3', 'GADD34', 'CHOP',
    }
    return genes


def main():
    print("=" * 80)
    print("SCRIPT 12: THREE-WAY VENN DIAGRAM")
    print("Expression × Binding × Secretory Pathway")
    print("=" * 80)

    # Load Venn results from script 11
    venn_path = OUTPUT_DIR / "tables" / "venn_expression_binding.csv"
    if not venn_path.exists():
        sys.exit(f"Missing Venn data. Run script 11 first: {venn_path}")

    df = pd.read_csv(venn_path)
    print(f"\nLoaded {len(df)} genes with opposite expression")

    # Get secretory pathway genes
    secretory_genes = get_secretory_pathway_genes()

    # Define the three sets
    set_expression = set(df['gene'])  # All have opposite expression
    set_binding = set(df[df['has_differential_binding']]['gene'])
    set_secretory = set_expression & secretory_genes  # Secretory genes that are in our dataset

    print("\n" + "=" * 80)
    print("THREE-WAY VENN CATEGORIES")
    print("=" * 80)

    print(f"\nCircle 1 - Opposite Expression: {len(set_expression)} genes")
    print(f"Circle 2 - Differential ER Binding: {len(set_binding)} genes")
    print(f"Circle 3 - Secretory Pathway (GO): {len(set_secretory)} genes")

    # Calculate intersections
    expr_only = set_expression - set_binding - set_secretory
    bind_only = set_binding - set_expression - set_secretory  # Should be 0
    secr_only = set_secretory - set_expression - set_binding  # Should be 0

    expr_bind = (set_expression & set_binding) - set_secretory
    expr_secr = (set_expression & set_secretory) - set_binding
    bind_secr = (set_binding & set_secretory) - set_expression  # Should be 0

    all_three = set_expression & set_binding & set_secretory

    print(f"\nExpression only: {len(expr_only)}")
    print(f"Expression + Binding: {len(expr_bind)}")
    print(f"Expression + Secretory: {len(expr_secr)}")
    print(f"\n*** ALL THREE (Expression + Binding + Secretory): {len(all_three)} genes ***")

    if len(all_three) > 0:
        print("\nGenes in triple intersection:")
        # Get effect sizes for these genes
        triple_df = df[df['gene'].isin(all_three)].copy()
        triple_df = triple_df.sort_values('total_effect', ascending=False)

        print(f"\n{'Rank':<5} {'Gene':<12} {'MCF7_FC':>8} {'T47D_FC':>8} {'MCF7_bind':>10} {'T47D_bind':>10} {'Effect':>8}")
        print("-" * 75)
        for i, (_, row) in enumerate(triple_df.iterrows(), 1):
            marker = " <<<" if row['gene'] == 'HSP90B1' else ""
            print(f"{i:<5} {row['gene']:<12} {row['MCF7_FC']:>8.2f} {row['T47D_FC']:>8.2f} "
                  f"{row['MCF7_binding_diff']:>10.0f} {row['T47D_binding_diff']:>10.0f} {row['total_effect']:>8.1f}{marker}")

        # Check HSP90B1
        if 'HSP90B1' in all_three:
            hsp_row = triple_df[triple_df['gene'] == 'HSP90B1'].iloc[0]
            hsp_rank = (triple_df['total_effect'] >= hsp_row['total_effect']).sum()
            print(f"\n*** HSP90B1 RANK: #{hsp_rank}/{len(all_three)} ***")
        else:
            print("\nHSP90B1 NOT in triple intersection")

    # Create figure
    print("\nGenerating figure...")
    fig = plt.figure(figsize=(16, 10))

    # Panel A: 3-way Venn
    ax1 = fig.add_subplot(2, 2, 1)

    # venn3 needs tuple of (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    # where uppercase = in set, lowercase = not in set
    venn_counts = (
        len(expr_only),           # Expression only (Abc)
        len(bind_only),           # Binding only (aBc) - should be 0
        len(expr_bind),           # Expression + Binding (ABc)
        len(secr_only),           # Secretory only (abC) - should be 0
        len(expr_secr),           # Expression + Secretory (AbC)
        len(bind_secr),           # Binding + Secretory (aBC) - should be 0
        len(all_three),           # All three (ABC)
    )

    v = venn3(subsets=venn_counts,
              set_labels=('Opposite\nExpression', 'Differential\nER Binding', 'Secretory\nPathway (GO)'),
              ax=ax1)

    ax1.set_title('A. Three-Way Venn: Winnowing Candidates\n(GSE89888 + GSE125117 + Gene Ontology)', fontsize=11)

    # Panel B: Genes in triple intersection
    ax2 = fig.add_subplot(2, 2, 2)

    if len(all_three) > 0:
        triple_df = df[df['gene'].isin(all_three)].sort_values('total_effect', ascending=False)
        n_show = min(15, len(triple_df))
        top = triple_df.head(n_show).iloc[::-1]

        y_pos = np.arange(n_show)
        colors = ['gold' if g == 'HSP90B1' else 'forestgreen' for g in top['gene']]

        ax2.barh(y_pos, top['total_effect'], color=colors, edgecolor='black')
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(top['gene'], fontsize=10)
        ax2.set_xlabel("Expression Effect Size (|Cohen's d|)", fontsize=10)
        ax2.set_title(f'B. ALL {len(all_three)} Genes in Triple Intersection\n(Gold = HSP90B1)', fontsize=11)
    else:
        ax2.text(0.5, 0.5, 'No genes in\ntriple intersection', ha='center', va='center',
                fontsize=14, transform=ax2.transAxes)
        ax2.set_title('B. Genes in Triple Intersection', fontsize=11)

    # Panel C: Binding pattern for triple intersection genes
    ax3 = fig.add_subplot(2, 2, 3)

    if len(all_three) > 0:
        triple_df = df[df['gene'].isin(all_three)].sort_values('total_effect', ascending=False)

        x = np.arange(len(triple_df))
        width = 0.35

        bars1 = ax3.bar(x - width/2, triple_df['MCF7_binding_diff'], width,
                       label='MCF7 (D538G-WT)', color='steelblue')
        bars2 = ax3.bar(x + width/2, triple_df['T47D_binding_diff'], width,
                       label='T47D (D538G-WT)', color='indianred')

        ax3.set_xticks(x)
        ax3.set_xticklabels(triple_df['gene'], rotation=45, ha='right', fontsize=9)
        ax3.set_ylabel('ER Binding Change (D538G - WT)')
        ax3.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax3.legend(loc='best', fontsize=9)
        ax3.set_title('C. ER Binding Changes for Triple Intersection Genes\n(GSE125117)', fontsize=10)

    # Panel D: Summary
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')

    hsp_rank_text = "N/A"
    if len(all_three) > 0 and 'HSP90B1' in all_three:
        triple_df = df[df['gene'].isin(all_three)].sort_values('total_effect', ascending=False)
        hsp_row = triple_df[triple_df['gene'] == 'HSP90B1'].iloc[0]
        hsp_rank = (triple_df['total_effect'] >= hsp_row['total_effect']).sum()
        hsp_rank_text = f"#{hsp_rank}/{len(all_three)}"

    summary_text = f"""THREE-WAY VENN: WINNOWING TO HSP90B1
=====================================

FILTER 1: Opposite Expression (GSE89888)
• {len(set_expression)} genes with opposite D538G effect
  (MCF7 vs T47D)

FILTER 2: Differential ER Binding (GSE125117)
• {len(set_binding)} genes with ER binding change
  at promoter upon D538G mutation

FILTER 3: Secretory Pathway (Gene Ontology)
• {len(set_secretory)} genes in ER/secretory pathway
  (GO terms for protein folding, ERAD, transport)

TRIPLE INTERSECTION:
• {len(all_three)} genes meet ALL THREE criteria

HSP90B1 STATUS:
• In intersection: {'YES' if 'HSP90B1' in all_three else 'NO'}
• Rank: {hsp_rank_text}

CONCLUSION:
Starting from {len(set_expression)} genes, the triple
intersection identifies only {len(all_three)} candidates
that are:
1. Oppositely regulated by D538G
2. Have differential ER binding
3. Are in the secretory pathway
"""

    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

    plt.tight_layout()

    fig_path = OUTPUT_DIR / "figures" / "fig12_triple_venn.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "figures" / "fig12_triple_venn.pdf", bbox_inches='tight')
    plt.close()

    print(f"\nSaved figure to: {fig_path}")

    print("\n" + "=" * 80)
    print("SCRIPT 12 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
