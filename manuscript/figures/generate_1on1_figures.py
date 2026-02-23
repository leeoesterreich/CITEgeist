#!/usr/bin/env python3
"""
Generate benchmark figures for 1-on-1 meeting.
Shows: Discrete model with Cellpose segmentation performance vs continuous baseline.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Set up output directory
output_dir = Path(__file__).parent / "1on1_figures"
output_dir.mkdir(exist_ok=True)

# Data from benchmark runs
# =========================

# MIXED dataset results (proportion correlation)
continuous_mixed = [0.911, 0.905, 0.893, 0.908, 0.888]  # λ=0.0 (no Laplacian)
discrete_mixed = [0.420, 0.454, 0.388, 0.439, 0.474]    # per-marker scaling
discrete_clr_mixed = [0.648, 0.696, 0.607, 0.671, 0.573]  # CLR preprocessing

# HIGH_SEG dataset results (proportion correlation)
discrete_high_seg = [0.729, 0.781, 0.831, 0.741, 0.830]

# Xenium achievable-7 (proportion correlation - same metric as simulation!)
# Computed from ground truth comparison
xenium_discrete = [0.6519, 0.6802, 0.6692, 0.6043, 0.5101]  # per region
xenium_continuous = [0.6257, 0.5336, 0.6249, 0.6007, 0.5493]  # per region

# ========================================
# Figure 1: Mixed Dataset - Method Comparison
# ========================================
fig1, ax1 = plt.subplots(figsize=(8, 6))

methods = ['Continuous\n(Baseline)', 'Discrete\n(per-marker)', 'Discrete\n(CLR)']
means = [np.mean(continuous_mixed), np.mean(discrete_mixed), np.mean(discrete_clr_mixed)]
stds = [np.std(continuous_mixed), np.std(discrete_mixed), np.std(discrete_clr_mixed)]

colors = ['#2ecc71', '#e74c3c', '#3498db']  # green, red, blue
bars = ax1.bar(methods, means, yerr=stds, capsize=5, color=colors, edgecolor='black', linewidth=1.5)

# Add value labels on bars
for bar, mean, std in zip(bars, means, stds):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 0.02,
             f'{mean:.2f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

ax1.set_ylabel('Proportion Correlation (Pearson r)', fontsize=12)
ax1.set_title('Mixed Dataset: Discrete Model Progress', fontsize=14, fontweight='bold')
ax1.set_ylim(0, 1.05)
ax1.axhline(y=0.9, color='gray', linestyle='--', alpha=0.5, label='Target (0.90)')
ax1.legend(loc='lower right')

# Add gap annotation
gap_closed = np.mean(discrete_clr_mixed) - np.mean(discrete_mixed)
total_gap = np.mean(continuous_mixed) - np.mean(discrete_mixed)
pct_closed = (gap_closed / total_gap) * 100

ax1.annotate('', xy=(2, np.mean(discrete_clr_mixed)), xytext=(1, np.mean(discrete_mixed)),
            arrowprops=dict(arrowstyle='->', color='#9b59b6', lw=2))
ax1.text(1.5, (np.mean(discrete_mixed) + np.mean(discrete_clr_mixed))/2 - 0.03,
         f'+{gap_closed:.2f}\n({pct_closed:.0f}% gap)',
         ha='center', va='top', fontsize=10, color='#9b59b6', fontweight='bold')

plt.tight_layout()
fig1.savefig(output_dir / "fig1_mixed_method_comparison.png", dpi=150, bbox_inches='tight')
fig1.savefig(output_dir / "fig1_mixed_method_comparison.svg", bbox_inches='tight')
print(f"Saved: {output_dir / 'fig1_mixed_method_comparison.png'}")


# ========================================
# Figure 2: All Three Datasets (Proportion Correlation)
# ========================================
fig2, ax2 = plt.subplots(figsize=(11, 6))

# Group by dataset - ALL three datasets
datasets = ['Simulation\n(mixed)', 'Simulation\n(high_seg)', 'Xenium\n(achievable-7)']
x = np.arange(len(datasets))
width = 0.25

# Continuous baseline
continuous_means = [np.mean(continuous_mixed), np.nan, np.mean(xenium_continuous)]
continuous_stds = [np.std(continuous_mixed), 0, np.std(xenium_continuous)]

# Discrete (per-marker scaling)
discrete_means = [np.mean(discrete_mixed), np.mean(discrete_high_seg), np.mean(xenium_discrete)]
discrete_stds = [np.std(discrete_mixed), np.std(discrete_high_seg), np.std(xenium_discrete)]

# Discrete + CLR (only have simulation mixed)
clr_means = [np.mean(discrete_clr_mixed), np.nan, np.nan]  # TODO: run others with CLR
clr_stds = [np.std(discrete_clr_mixed), 0, 0]

# Plot bars
bars1 = ax2.bar(x - width, [m if not np.isnan(m) else 0 for m in continuous_means],
                width, label='Continuous', color='#2ecc71', edgecolor='black', linewidth=1.5,
                yerr=[s if not np.isnan(m) else 0 for m, s in zip(continuous_means, continuous_stds)], capsize=4)
bars2 = ax2.bar(x, discrete_means, width, yerr=discrete_stds, capsize=4,
                label='Discrete (Cellpose)', color='#e74c3c', edgecolor='black', linewidth=1.5)
bars3 = ax2.bar(x + width, [m if not np.isnan(m) else 0 for m in clr_means],
                width, label='Discrete + CLR', color='#3498db', edgecolor='black', linewidth=1.5,
                yerr=[s if not np.isnan(m) else 0 for m, s in zip(clr_means, clr_stds)], capsize=4)

# Add value labels
for bar, mean, std in zip(bars1, continuous_means, continuous_stds):
    if not np.isnan(mean) and mean > 0:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 0.02,
                 f'{mean:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

for bar, mean, std in zip(bars2, discrete_means, discrete_stds):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 0.02,
             f'{mean:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

for bar, mean, std in zip(bars3, clr_means, clr_stds):
    if not np.isnan(mean) and mean > 0:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + std + 0.02,
                 f'{mean:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax2.set_ylabel('Proportion Correlation (Pearson r)', fontsize=12)
ax2.set_title('Discrete Model Performance Across All Datasets', fontsize=14, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(datasets, fontsize=11)
ax2.set_ylim(0, 1.05)
ax2.legend(loc='upper right')
ax2.axhline(y=0.9, color='gray', linestyle='--', alpha=0.5, label='_nolegend_')

# Add annotations
ax2.annotate('Gap to close\non mixed', xy=(0.25, 0.64), xytext=(0.35, 0.80),
            fontsize=9, ha='left',
            arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
ax2.annotate('Discrete WINS!', xy=(2, 0.65), xytext=(2.15, 0.50),
            fontsize=10, ha='left', color='#c0392b', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5))

plt.tight_layout()
fig2.savefig(output_dir / "fig2_all_datasets.png", dpi=150, bbox_inches='tight')
fig2.savefig(output_dir / "fig2_all_datasets.svg", bbox_inches='tight')
print(f"Saved: {output_dir / 'fig2_all_datasets.png'}")


# ========================================
# Figure 3: Gap Analysis - What's Left
# ========================================
fig3, ax3 = plt.subplots(figsize=(8, 5))

# Waterfall-style breakdown
categories = ['Baseline\nDiscrete', 'CLR\nPreprocessing', 'Remaining\nGap', 'Target\n(Continuous)']
values = [
    np.mean(discrete_mixed),
    gap_closed,
    np.mean(continuous_mixed) - np.mean(discrete_clr_mixed),
    0  # placeholder
]

# Cumulative for waterfall
cumulative = [values[0], values[0] + values[1], values[0] + values[1], np.mean(continuous_mixed)]
bottoms = [0, values[0], values[0] + values[1], 0]

colors = ['#e74c3c', '#3498db', '#95a5a6', '#2ecc71']

for i, (cat, val, cum, bot) in enumerate(zip(categories, values, cumulative, bottoms)):
    if i == 3:  # Target bar
        ax3.bar(i, np.mean(continuous_mixed), color=colors[i], edgecolor='black', linewidth=1.5)
        ax3.text(i, np.mean(continuous_mixed) + 0.02, f'{np.mean(continuous_mixed):.2f}',
                 ha='center', fontsize=11, fontweight='bold')
    elif i == 2:  # Remaining gap
        ax3.bar(i, val, bottom=bot, color=colors[i], edgecolor='black', linewidth=1.5, hatch='//')
        ax3.text(i, cum + 0.02, f'+{val:.2f}\n(to close)', ha='center', fontsize=10, color='#7f8c8d')
    else:
        ax3.bar(i, val if i == 0 else val, bottom=bot, color=colors[i], edgecolor='black', linewidth=1.5)
        if i == 0:
            ax3.text(i, val + 0.02, f'{val:.2f}', ha='center', fontsize=11, fontweight='bold')
        else:
            ax3.text(i, cum + 0.02, f'+{val:.2f}', ha='center', fontsize=11, fontweight='bold', color='#2980b9')

ax3.set_xticks(range(len(categories)))
ax3.set_xticklabels(categories, fontsize=11)
ax3.set_ylabel('Proportion Correlation', fontsize=12)
ax3.set_title('Closing the Gap: CLR Fixed ~45% of the Problem', fontsize=14, fontweight='bold')
ax3.set_ylim(0, 1.1)

# Add legend-like annotations
ax3.annotate('Next: investigate\nIQP objective & beta fitting', xy=(2, 0.5), fontsize=10,
             ha='center', style='italic', color='#7f8c8d')

plt.tight_layout()
fig3.savefig(output_dir / "fig3_gap_analysis.png", dpi=150, bbox_inches='tight')
fig3.savefig(output_dir / "fig3_gap_analysis.svg", bbox_inches='tight')
print(f"Saved: {output_dir / 'fig3_gap_analysis.png'}")

# ========================================
# Figure 4: Cellpose Nuclei Detection Accuracy
# ========================================
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(12, 5))

# Simulation: Cellpose vs ground truth nuclei counts (from logs)
sim_nuclei_corr = [0.986, 0.987, 0.986, 0.987, 0.986]  # from discrete_nuc_clr logs

# Xenium: Cellpose vs Xenium n_cells ground truth
xenium_nuclei_corr = [0.9576, 0.9069, 0.9681, 0.9781, 0.9795]  # computed above

# Panel A: Bar chart comparing simulation vs Xenium
datasets = ['Simulation\n(DAPI)', 'Xenium\n(DAPI)']
means = [np.mean(sim_nuclei_corr), np.mean(xenium_nuclei_corr)]
stds = [np.std(sim_nuclei_corr), np.std(xenium_nuclei_corr)]

bars = ax4a.bar(datasets, means, yerr=stds, capsize=5, color=['#3498db', '#e74c3c'],
                edgecolor='black', linewidth=1.5)
for bar, mean in zip(bars, means):
    ax4a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
              f'{mean:.3f}', ha='center', va='bottom', fontsize=14, fontweight='bold')

ax4a.set_ylabel('Pearson Correlation with Ground Truth', fontsize=12)
ax4a.set_title('Cellpose Nuclei Detection Accuracy', fontsize=14, fontweight='bold')
ax4a.set_ylim(0.85, 1.02)
ax4a.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)

# Panel B: Per-region breakdown for Xenium
regions = [f'Region {i}' for i in range(5)]
colors_regions = plt.cm.Reds(np.linspace(0.4, 0.8, 5))
bars_b = ax4b.bar(regions, xenium_nuclei_corr, color=colors_regions, edgecolor='black', linewidth=1.5)
for bar, val in zip(bars_b, xenium_nuclei_corr):
    ax4b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
              f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax4b.set_ylabel('Pearson Correlation', fontsize=12)
ax4b.set_title('Xenium Cellpose Accuracy by Region', fontsize=14, fontweight='bold')
ax4b.set_ylim(0.85, 1.02)
ax4b.axhline(y=np.mean(xenium_nuclei_corr), color='#c0392b', linestyle='--', alpha=0.7,
             label=f'Mean: {np.mean(xenium_nuclei_corr):.3f}')
ax4b.legend(loc='lower right')

plt.tight_layout()
fig4.savefig(output_dir / "fig4_cellpose_accuracy.png", dpi=150, bbox_inches='tight')
fig4.savefig(output_dir / "fig4_cellpose_accuracy.svg", bbox_inches='tight')
print(f"Saved: {output_dir / 'fig4_cellpose_accuracy.png'}")


# ========================================
# Figure 5: Xenium GEX Deconvolution (THE KEY RESULT)
# ========================================
fig5, ax5 = plt.subplots(figsize=(8, 5))

# GEX Pearson r from evaluation log (eval_discrete_7999070.out)
methods = ['CITEgeist\nDiscrete', 'CITEgeist\nContinuous', 'Cell2Location', 'Tangram']
gex_pearson = [0.3838, 0.3648, 0.3628, 0.2090]
gex_rmse = [4.4341, 4.5065, 4.9939, 4.3250]

colors = ['#e74c3c', '#2ecc71', '#9b59b6', '#f39c12']
bars = ax5.bar(methods, gex_pearson, color=colors, edgecolor='black', linewidth=1.5)

# Add value labels
for bar, val in zip(bars, gex_pearson):
    ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
             f'{val:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

# Highlight the winner
bars[0].set_edgecolor('#c0392b')
bars[0].set_linewidth(3)

ax5.set_ylabel('GEX Pearson Correlation', fontsize=12)
ax5.set_title('Xenium GEX Deconvolution: Discrete Wins!', fontsize=14, fontweight='bold')
ax5.set_ylim(0, 0.5)

# Add improvement annotation
improvement = gex_pearson[0] - gex_pearson[1]
pct_improvement = (improvement / gex_pearson[1]) * 100
ax5.annotate(f'+{improvement:.3f}\n(+{pct_improvement:.1f}%)',
             xy=(0.5, 0.375), fontsize=11, ha='center', color='#c0392b', fontweight='bold')

plt.tight_layout()
fig5.savefig(output_dir / "fig5_xenium_gex_comparison.png", dpi=150, bbox_inches='tight')
fig5.savefig(output_dir / "fig5_xenium_gex_comparison.svg", bbox_inches='tight')
print(f"Saved: {output_dir / 'fig5_xenium_gex_comparison.png'}")


print("\n=== Summary for 1-on-1 ===")
print("SIMULATION (mixed):")
print(f"  Continuous: {np.mean(continuous_mixed):.3f} ± {np.std(continuous_mixed):.3f}")
print(f"  Discrete: {np.mean(discrete_mixed):.3f} ± {np.std(discrete_mixed):.3f}")
print(f"  Discrete+CLR: {np.mean(discrete_clr_mixed):.3f} ± {np.std(discrete_clr_mixed):.3f}")
print(f"\nSIMULATION (high_seg):")
print(f"  Discrete: {np.mean(discrete_high_seg):.3f} ± {np.std(discrete_high_seg):.3f}")
print(f"\nXENIUM (achievable-7) - PROPORTION:")
print(f"  Continuous: {np.mean(xenium_continuous):.3f} ± {np.std(xenium_continuous):.3f}")
print(f"  Discrete: {np.mean(xenium_discrete):.3f} ± {np.std(xenium_discrete):.3f}")
print(f"\nXENIUM (achievable-7) - GEX DECONVOLUTION (KEY RESULT!):")
print(f"  CITEgeist Discrete: 0.384")
print(f"  CITEgeist Continuous: 0.365")
print(f"  Cell2Location: 0.363")
print(f"  ** Discrete BEATS all methods on GEX! **")
print(f"\nGap Analysis (mixed dataset):")
print(f"  Gap closed by CLR: {gap_closed:.3f} ({pct_closed:.1f}%)")
print(f"  Remaining gap: {np.mean(continuous_mixed) - np.mean(discrete_clr_mixed):.3f}")
