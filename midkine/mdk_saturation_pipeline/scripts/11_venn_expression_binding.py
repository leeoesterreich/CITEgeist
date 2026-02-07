#!/usr/bin/env python
"""
11_venn_expression_binding.py

VENN DIAGRAM APPROACH: Intersect expression and binding data genome-wide

Purpose:
1. Take the 666 genes with opposite expression regulation (from script 10)
2. Check genome-wide ER ChIP-seq binding at each gene's promoter
3. Identify genes with DIFFERENTIAL ER binding (D538G vs WT)
4. Create Venn diagram showing intersection

This approach winnows down candidates by requiring BOTH:
- Opposite expression regulation (MCF7 vs T47D response to D538G)
- Differential ER binding at the gene's locus

Datasets:
- GSE89888 (RNA-seq) - expression data
- GSE125117 (ER ChIP-seq) - binding data

Outputs:
- outputs/tables/venn_expression_binding.csv
- outputs/figures/fig11_venn_expression_binding.png
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import gzip
import yaml
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).parent
PIPELINE_DIR = SCRIPT_DIR.parent
CONFIG_PATH = PIPELINE_DIR / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

DATA_DIR = PIPELINE_DIR / config['paths']['data_dir']
OUTPUT_DIR = PIPELINE_DIR / config['paths']['output_dir']

# Promoter window around TSS
PROMOTER_WINDOW = 10000  # ±10kb from TSS


def load_peaks(filepath):
    """Load BED peak file with coordinates and signal."""
    peaks = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    peak = {
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'score': float(parts[4]) if len(parts) > 4 else 1.0
                    }
                    peaks.append(peak)
                except (ValueError, IndexError):
                    continue
    return pd.DataFrame(peaks)


def get_gene_coordinates(gene_list):
    """
    Get TSS coordinates for a list of genes using mygene.
    Returns dict: gene_symbol -> (chr, tss, strand)
    """
    try:
        import mygene
        mg = mygene.MyGeneInfo()

        # Query specific genes (batch query)
        print(f"Fetching gene coordinates for {len(gene_list)} genes from MyGene.info...")
        results = mg.querymany(gene_list, scopes='symbol',
                              fields='symbol,genomic_pos',
                              species='human', returnall=True)

        gene_coords = {}
        for hit in results['out']:
            if 'symbol' not in hit or 'genomic_pos' not in hit:
                continue

            symbol = hit['symbol']
            gpos = hit['genomic_pos']

            # Handle multiple genomic positions
            if isinstance(gpos, list):
                gpos = gpos[0]

            if isinstance(gpos, dict) and 'chr' in gpos:
                chrom = 'chr' + str(gpos['chr']) if not str(gpos['chr']).startswith('chr') else gpos['chr']
                strand = gpos.get('strand', 1)
                start = gpos.get('start', 0)
                end = gpos.get('end', 0)

                # TSS is start for + strand, end for - strand
                tss = start if strand == 1 else end
                gene_coords[symbol] = (chrom, tss, strand)

        print(f"Retrieved coordinates for {len(gene_coords)} genes")
        return gene_coords

    except Exception as e:
        print(f"Error fetching gene coordinates: {e}")
        return {}


def get_binding_at_promoter(peaks_df, chrom, tss, window=PROMOTER_WINDOW):
    """
    Get total ER ChIP-seq signal at a gene's promoter.
    Returns sum of peak scores within window of TSS.
    """
    if peaks_df is None or len(peaks_df) == 0:
        return 0.0

    # Filter to chromosome
    chr_peaks = peaks_df[peaks_df['chr'] == chrom]
    if len(chr_peaks) == 0:
        return 0.0

    # Find peaks overlapping promoter window
    prom_start = tss - window
    prom_end = tss + window

    overlapping = chr_peaks[
        (chr_peaks['start'] < prom_end) &
        (chr_peaks['end'] > prom_start)
    ]

    if len(overlapping) == 0:
        return 0.0

    return overlapping['score'].sum()


def main():
    print("=" * 80)
    print("SCRIPT 11: VENN DIAGRAM - EXPRESSION × BINDING INTERSECTION")
    print("=" * 80)

    # Load expression results from script 10
    expr_path = OUTPUT_DIR / "tables" / "unsupervised_hsp90b1_ranking.csv"
    if not expr_path.exists():
        sys.exit(f"Missing expression data. Run script 10 first: {expr_path}")

    expr_df = pd.read_csv(expr_path)
    print(f"\nLoaded {len(expr_df)} genes with significant opposite regulation")

    # Load ChIP-seq peak files (absolute paths)
    print("\nLoading ER ChIP-seq peak files (GSE125117)...")
    CHIPSEQ_DIR = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine')
    peak_files = {
        'MCF7_WT': CHIPSEQ_DIR / 'GSM3563751_MCF7_WT_E2_peaks.bed.gz',
        'MCF7_D538G': CHIPSEQ_DIR / 'GSM3563757_MCF7_D538G_E2_peaks.bed.gz',
        'T47D_WT': CHIPSEQ_DIR / 'GSM3563760_T47D_WT_E2_peaks.bed.gz',
        'T47D_D538G': CHIPSEQ_DIR / 'GSM3563766_T47D_D538G_E2_peaks.bed.gz',
    }

    peaks = {}
    for name, path in peak_files.items():
        if path.exists():
            peaks[name] = load_peaks(path)
            print(f"  {name}: {len(peaks[name])} peaks")
        else:
            print(f"  {name}: FILE NOT FOUND")
            peaks[name] = pd.DataFrame()

    # Get gene coordinates for genes in expression data
    gene_list = expr_df['gene'].tolist()
    gene_coords = get_gene_coordinates(gene_list)

    # Calculate binding at each gene's promoter
    print(f"\nCalculating ER binding at promoters (±{PROMOTER_WINDOW//1000}kb from TSS)...")

    binding_results = []
    genes_with_coords = 0

    for _, row in expr_df.iterrows():
        gene = row['gene']

        if gene not in gene_coords:
            continue

        genes_with_coords += 1
        chrom, tss, strand = gene_coords[gene]

        # Get binding in each condition
        mcf7_wt_bind = get_binding_at_promoter(peaks['MCF7_WT'], chrom, tss)
        mcf7_d538g_bind = get_binding_at_promoter(peaks['MCF7_D538G'], chrom, tss)
        t47d_wt_bind = get_binding_at_promoter(peaks['T47D_WT'], chrom, tss)
        t47d_d538g_bind = get_binding_at_promoter(peaks['T47D_D538G'], chrom, tss)

        # Calculate differential binding
        # MCF7: D538G vs WT
        mcf7_diff = mcf7_d538g_bind - mcf7_wt_bind
        # T47D: D538G vs WT
        t47d_diff = t47d_d538g_bind - t47d_wt_bind

        # Has any ER binding at promoter?
        has_binding = (mcf7_wt_bind > 0 or mcf7_d538g_bind > 0 or
                      t47d_wt_bind > 0 or t47d_d538g_bind > 0)

        # Has differential binding? (change > threshold)
        DIFF_THRESHOLD = 50  # Minimum signal change to count as differential
        has_diff_binding = (abs(mcf7_diff) > DIFF_THRESHOLD or
                           abs(t47d_diff) > DIFF_THRESHOLD)

        # Opposite binding change? (MCF7 goes one way, T47D goes other way)
        opposite_binding = ((mcf7_diff > DIFF_THRESHOLD and t47d_diff < -DIFF_THRESHOLD) or
                           (mcf7_diff < -DIFF_THRESHOLD and t47d_diff > DIFF_THRESHOLD))

        binding_results.append({
            'gene': gene,
            'chr': chrom,
            'tss': tss,
            'MCF7_WT_binding': mcf7_wt_bind,
            'MCF7_D538G_binding': mcf7_d538g_bind,
            'MCF7_binding_diff': mcf7_diff,
            'T47D_WT_binding': t47d_wt_bind,
            'T47D_D538G_binding': t47d_d538g_bind,
            'T47D_binding_diff': t47d_diff,
            'has_any_binding': has_binding,
            'has_differential_binding': has_diff_binding,
            'has_opposite_binding': opposite_binding,
            # Expression data
            'MCF7_FC': row['MCF7_FC'],
            'T47D_FC': row['T47D_FC'],
            'correlation_type': row['correlation_type'],
            'total_effect': row['total_effect'],
        })

    print(f"Found coordinates for {genes_with_coords}/{len(expr_df)} genes")

    bind_df = pd.DataFrame(binding_results)

    # Summary statistics
    print("\n" + "=" * 80)
    print("VENN DIAGRAM CATEGORIES")
    print("=" * 80)

    n_opposite_expr = len(bind_df)
    n_any_binding = bind_df['has_any_binding'].sum()
    n_diff_binding = bind_df['has_differential_binding'].sum()
    n_opposite_binding = bind_df['has_opposite_binding'].sum()

    print(f"\nCircle 1: Opposite Expression (D538G effect)")
    print(f"  {n_opposite_expr} genes")

    print(f"\nCircle 2: Any ER binding at promoter")
    print(f"  {n_any_binding} genes")

    print(f"\nCircle 3: Differential ER binding (D538G vs WT)")
    print(f"  {n_diff_binding} genes")

    # Intersection: opposite expression + differential binding
    intersection = bind_df[bind_df['has_differential_binding']]
    print(f"\n*** INTERSECTION (Opposite Expression + Differential Binding): {len(intersection)} genes ***")

    if len(intersection) > 0:
        print("\nTop candidates by effect size:")
        top_intersection = intersection.sort_values('total_effect', ascending=False).head(20)
        print(f"\n{'Rank':<5} {'Gene':<12} {'MCF7_FC':>8} {'T47D_FC':>8} {'MCF7_bind':>10} {'T47D_bind':>10} {'Effect':>8}")
        print("-" * 75)
        for i, (_, row) in enumerate(top_intersection.iterrows(), 1):
            print(f"{i:<5} {row['gene']:<12} {row['MCF7_FC']:>8.2f} {row['T47D_FC']:>8.2f} "
                  f"{row['MCF7_binding_diff']:>10.0f} {row['T47D_binding_diff']:>10.0f} {row['total_effect']:>8.1f}")

        # Check if HSP90B1 is in intersection
        if 'HSP90B1' in intersection['gene'].values:
            hsp_row = intersection[intersection['gene'] == 'HSP90B1'].iloc[0]
            hsp_rank = (intersection['total_effect'] >= hsp_row['total_effect']).sum()
            print(f"\n*** HSP90B1 is in the intersection! Rank: #{hsp_rank}/{len(intersection)} ***")
            print(f"    MCF7 binding change: {hsp_row['MCF7_binding_diff']:.0f}")
            print(f"    T47D binding change: {hsp_row['T47D_binding_diff']:.0f}")
        else:
            print("\nHSP90B1 is NOT in the intersection")
            # Check why
            if 'HSP90B1' in bind_df['gene'].values:
                hsp_row = bind_df[bind_df['gene'] == 'HSP90B1'].iloc[0]
                print(f"    HSP90B1 binding: MCF7_diff={hsp_row['MCF7_binding_diff']:.0f}, T47D_diff={hsp_row['T47D_binding_diff']:.0f}")
                print(f"    has_differential_binding: {hsp_row['has_differential_binding']}")

    # Save results
    bind_df.to_csv(OUTPUT_DIR / "tables" / "venn_expression_binding.csv", index=False)

    # Create figure
    print("\nGenerating Venn diagram figure...")
    fig = plt.figure(figsize=(16, 12))

    # Panel A: Venn diagram
    ax1 = fig.add_subplot(2, 2, 1)

    # Two-circle Venn: opposite expression vs differential binding
    set_expr = set(bind_df['gene'])
    set_binding = set(bind_df[bind_df['has_differential_binding']]['gene'])

    # Since all genes in bind_df have opposite expression, we show:
    # - How many have differential binding
    # - How many don't
    venn = venn2([set_expr - set_binding, set_binding],
                 set_labels=('Opposite Expression\nOnly', 'Differential Binding'),
                 ax=ax1)
    ax1.set_title('A. Venn: Opposite Expression ∩ Differential ER Binding\n(GSE89888 + GSE125117)', fontsize=11)

    # Panel B: Candidates in intersection
    ax2 = fig.add_subplot(2, 2, 2)

    if len(intersection) > 0:
        top_n = min(15, len(intersection))
        top = intersection.sort_values('total_effect', ascending=False).head(top_n)
        top = top.iloc[::-1]

        y_pos = np.arange(top_n)
        colors = ['gold' if g == 'HSP90B1' else 'steelblue' for g in top['gene']]

        ax2.barh(y_pos, top['total_effect'], color=colors, edgecolor='black')
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(top['gene'], fontsize=9)
        ax2.set_xlabel("Expression Effect Size (|Cohen's d|)")
        ax2.set_title(f'B. Top {top_n} Candidates in Intersection\n(Opposite Expression + Differential Binding)', fontsize=10)

    # Panel C: Binding changes for top candidates
    ax3 = fig.add_subplot(2, 2, 3)

    if len(intersection) > 0:
        top_bind = intersection.sort_values('total_effect', ascending=False).head(10)

        x = np.arange(len(top_bind))
        width = 0.35

        ax3.bar(x - width/2, top_bind['MCF7_binding_diff'], width, label='MCF7 (D538G-WT)', color='steelblue')
        ax3.bar(x + width/2, top_bind['T47D_binding_diff'], width, label='T47D (D538G-WT)', color='indianred')

        ax3.set_xticks(x)
        ax3.set_xticklabels(top_bind['gene'], rotation=45, ha='right', fontsize=9)
        ax3.set_ylabel('ER Binding Change (D538G - WT)')
        ax3.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax3.legend(loc='best', fontsize=9)
        ax3.set_title('C. ER Binding Changes at Candidate Promoters (GSE125117)', fontsize=10)

    # Panel D: Summary
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')

    # Find HSP90B1 rank in intersection
    hsp_in_intersection = 'HSP90B1' in intersection['gene'].values if len(intersection) > 0 else False
    hsp_rank_text = "N/A"
    if hsp_in_intersection:
        hsp_row = intersection[intersection['gene'] == 'HSP90B1'].iloc[0]
        hsp_rank = (intersection['total_effect'] >= hsp_row['total_effect']).sum()
        hsp_rank_text = f"#{hsp_rank}/{len(intersection)}"

    summary_text = f"""VENN DIAGRAM APPROACH: WINNOWING CANDIDATES
============================================

STARTING POINT:
• {n_opposite_expr} genes with opposite expression
  (MCF7 vs T47D response to D538G)

FILTER BY ER BINDING:
• {n_any_binding} have any ER binding at promoter
• {n_diff_binding} have DIFFERENTIAL ER binding
  (D538G changes binding vs WT)

INTERSECTION:
• {len(intersection)} genes have BOTH:
  1. Opposite expression regulation
  2. Differential ER binding at promoter

HSP90B1 STATUS:
• In intersection: {'YES' if hsp_in_intersection else 'NO'}
• Rank in intersection: {hsp_rank_text}

INTERPRETATION:
The intersection contains genes where:
• Expression is opposite in MCF7 vs T47D
• ER binding changes upon D538G mutation
These are the strongest candidates for
ER-mediated regulation of the phenotype."""

    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()

    fig_path = OUTPUT_DIR / "figures" / "fig11_venn_expression_binding.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / "figures" / "fig11_venn_expression_binding.pdf", bbox_inches='tight')
    plt.close()

    print(f"\nSaved figure to: {fig_path}")
    print(f"Saved table to: {OUTPUT_DIR / 'tables' / 'venn_expression_binding.csv'}")

    print("\n" + "=" * 80)
    print("SCRIPT 11 COMPLETE")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
