#!/usr/bin/env python3
"""
Quick visualization of per-cell-type correlation with αSMA-only fibroblast GT fix.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Apply clean style
plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.family': 'sans-serif',
})

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Old results (with Vimentin gate) - from design doc
OLD_RESULTS = {
    "B cells": 0.56,
    "CD4+ T cells": 0.42,
    "CD8+ T cells": 0.79,
    "Macrophages": 0.80,
    "Endothelial": 0.58,
    "Epithelial": 0.60,
    "Fibroblasts": 0.42,  # This was the problem!
}

# New results (αSMA only) - from evaluation
NEW_RESULTS = {
    "B cells": 0.556,
    "CD4+ T cells": 0.416,
    "CD8+ T cells": 0.793,
    "Macrophages": 0.802,
    "Endothelial": 0.592,
    "Epithelial": 0.598,
    "Fibroblasts": 0.535,  # Fixed!
}

NEW_STD = {
    "B cells": 0.418,
    "CD4+ T cells": 0.165,
    "CD8+ T cells": 0.062,
    "Macrophages": 0.019,
    "Endothelial": 0.161,
    "Epithelial": 0.090,
    "Fibroblasts": 0.133,
}

def main():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    cell_types = list(NEW_RESULTS.keys())
    x = np.arange(len(cell_types))
    width = 0.35

    # Panel A: Before vs After comparison
    ax = axes[0]
    old_vals = [OLD_RESULTS[ct] for ct in cell_types]
    new_vals = [NEW_RESULTS[ct] for ct in cell_types]
    new_errs = [NEW_STD[ct] for ct in cell_types]

    bars1 = ax.bar(x - width/2, old_vals, width, label='Old GT (αSMA + VIM)',
                   color='#d62728', alpha=0.7)
    bars2 = ax.bar(x + width/2, new_vals, width, yerr=new_errs,
                   label='New GT (αSMA only)', color='#2ca02c', alpha=0.8,
                   capsize=3, error_kw={'lw': 1})

    # Highlight Fibroblasts improvement
    fibroblast_idx = cell_types.index("Fibroblasts")
    bars1[fibroblast_idx].set_edgecolor('black')
    bars1[fibroblast_idx].set_linewidth(2)
    bars2[fibroblast_idx].set_edgecolor('black')
    bars2[fibroblast_idx].set_linewidth(2)

    # Add improvement annotation for Fibroblasts
    improvement = NEW_RESULTS["Fibroblasts"] - OLD_RESULTS["Fibroblasts"]
    ax.annotate(f'+{improvement:.2f}',
                xy=(fibroblast_idx + width/2, NEW_RESULTS["Fibroblasts"] + NEW_STD["Fibroblasts"] + 0.02),
                ha='center', fontsize=11, fontweight='bold', color='#2ca02c')

    ax.set_ylabel('Pearson Correlation (r)')
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45, ha='right')
    ax.set_ylim(0, 1.0)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, lw=0.8)
    ax.legend(loc='upper left', framealpha=0.9)
    ax.set_title('Per-Cell-Type Correlation: Before vs After αSMA Fix', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel B: Improvement delta
    ax2 = axes[1]
    deltas = [NEW_RESULTS[ct] - OLD_RESULTS[ct] for ct in cell_types]
    colors = ['#2ca02c' if d > 0 else '#d62728' for d in deltas]

    bars = ax2.bar(x, deltas, width=0.6, color=colors, alpha=0.8)

    # Highlight Fibroblasts
    bars[fibroblast_idx].set_edgecolor('black')
    bars[fibroblast_idx].set_linewidth(2)

    ax2.axhline(y=0, color='black', linestyle='-', lw=0.8)
    ax2.set_ylabel('Δ Pearson r (New - Old)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(cell_types, rotation=45, ha='right')
    ax2.set_ylim(-0.15, 0.15)
    ax2.set_title('Correlation Change After αSMA-Only GT Fix', fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Add value labels
    for i, (bar, delta) in enumerate(zip(bars, deltas)):
        if abs(delta) > 0.01:
            va = 'bottom' if delta > 0 else 'top'
            offset = 0.005 if delta > 0 else -0.005
            ax2.text(bar.get_x() + bar.get_width()/2, delta + offset,
                    f'{delta:+.2f}', ha='center', va=va, fontsize=9)

    plt.tight_layout()

    # Save
    out_path = OUTPUT_DIR / "asma_fix_comparison.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved to: {out_path}")

    # Also save PDF
    plt.savefig(OUTPUT_DIR / "asma_fix_comparison.pdf", bbox_inches='tight', facecolor='white')

    plt.close()

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY: αSMA-Only Fibroblast GT Fix")
    print("="*60)
    print(f"{'Cell Type':<20} {'Old r':<10} {'New r':<10} {'Δ':<10}")
    print("-"*50)
    for ct in cell_types:
        delta = NEW_RESULTS[ct] - OLD_RESULTS[ct]
        marker = "***" if ct == "Fibroblasts" else ""
        print(f"{ct:<20} {OLD_RESULTS[ct]:<10.3f} {NEW_RESULTS[ct]:<10.3f} {delta:+.3f} {marker}")
    print("="*60)
    print(f"Fibroblast improvement: +{NEW_RESULTS['Fibroblasts'] - OLD_RESULTS['Fibroblasts']:.3f}")
    print(f"  (from r={OLD_RESULTS['Fibroblasts']:.2f} to r={NEW_RESULTS['Fibroblasts']:.2f})")

if __name__ == '__main__':
    main()
