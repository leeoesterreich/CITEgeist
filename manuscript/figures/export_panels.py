#!/usr/bin/env python3
"""
Export individual figure panels as PNG (300 DPI) and SVG for Illustrator.

Creates per-figure folders under manuscript/figures/export/Figure{1-5}/
with individual panel files: panel_A.png, panel_A.svg, etc.

For schematic panels (SVG-based): copies existing rendered files.
For data panels: imports panel functions and renders each independently.
"""

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Make sure we can import figure modules from this directory
sys.path.insert(0, str(Path(__file__).parent))

from figure_style import apply_style, PALETTE

apply_style()

# Paths
SCRIPT_DIR = Path(__file__).parent
EXPORT_DIR = SCRIPT_DIR / "export"
SCHEMATIC_DIR = SCRIPT_DIR / "output" / "schematics"
SCHEMATIC_RENDERED = SCHEMATIC_DIR / "rendered"

DPI = 300


def export_panel(fig, ax, export_path_stem):
    """Save a single-panel figure as PNG and SVG."""
    for fmt in ("png", "svg"):
        path = f"{export_path_stem}.{fmt}"
        try:
            if fmt == "svg":
                fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
            else:
                fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
        except (ZeroDivisionError, ValueError):
            # Nested GridSpec can trigger ZeroDivisionError in constrained_layout
            # when using bbox_inches="tight"; fall back to saving without it
            if fmt == "svg":
                fig.savefig(path, format="svg", facecolor="white")
            else:
                fig.savefig(path, dpi=DPI, facecolor="white")
        print(f"    {Path(path).name}")
    plt.close(fig)


def copy_schematic(fig_num, panel_letter, svg_name):
    """Copy existing schematic SVG and rendered PNG."""
    dest_dir = EXPORT_DIR / f"Figure{fig_num}"
    dest_dir.mkdir(parents=True, exist_ok=True)
    stem = f"panel_{panel_letter}"

    # SVG
    svg_src = SCHEMATIC_DIR / svg_name
    if svg_src.exists():
        shutil.copy2(svg_src, dest_dir / f"{stem}.svg")
        print(f"    {stem}.svg (schematic)")

    # Rendered PNG
    png_name = svg_name.replace(".svg", ".png")
    png_src = SCHEMATIC_RENDERED / png_name
    if png_src.exists():
        shutil.copy2(png_src, dest_dir / f"{stem}.png")
        print(f"    {stem}.png (schematic)")


def render_data_panel(fig_num, panel_letter, panel_func, *args,
                      figsize=(5, 4), **kwargs):
    """Render a data panel function into its own figure and export."""
    dest_dir = EXPORT_DIR / f"Figure{fig_num}"
    dest_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    panel_func(ax, *args, **kwargs)
    try:
        fig.tight_layout()
    except RuntimeError:
        pass  # Colorbar layout conflict - skip tight_layout
    export_panel(fig, ax, str(dest_dir / f"panel_{panel_letter}"))


# =========================================================================
# Figure 1: Pipeline Overview (all schematics)
# =========================================================================
def export_figure1():
    print("\n  Figure 1: Pipeline Overview")
    copy_schematic(1, "A", "figure1_panel_a_pipeline.svg")
    copy_schematic(1, "B", "figure1_panel_b_spatial_stats.svg")
    copy_schematic(1, "C", "figure1_panel_c_resolution.svg")


# =========================================================================
# Figure 2: Profile Discovery (3 schematics + 1 data panel)
# =========================================================================
def export_figure2():
    print("\n  Figure 2: Profile Discovery")
    copy_schematic(2, "A", "figure2_panel_a_marker_interest.svg")
    copy_schematic(2, "B", "figure2_panel_b_profile_discovery.svg")
    copy_schematic(2, "C", "figure2_panel_c_xenium_demo.svg")

    # Panel D is a data panel from generate_figure2
    try:
        from generate_figure2 import panel_d_profile_table
        render_data_panel(2, "D", panel_d_profile_table, figsize=(5, 4))
    except Exception as e:
        print(f"    WARNING: Could not render 2D: {e}")


# =========================================================================
# Figure 3: Benchmarking (6 data panels)
# =========================================================================
def export_figure3():
    print("\n  Figure 3: Benchmarking")

    try:
        from generate_figure3 import (
            panel_a_profile_discovery,
            panel_b_xenium_proportions, load_xenium_summary,
            panel_c_simulated_proportions, load_simulation_props,
            panel_d_gex_benchmark, load_simulation_gex,
            panel_e_scatter, panel_f_spatial_comparison,
            load_scatter_data, load_spatial_coords,
        )

        # Load data
        xenium_summary = load_xenium_summary()
        sim_props = load_simulation_props()
        sim_gex = load_simulation_gex()
        gt, pred = load_scatter_data(region_id=0)
        coords = load_spatial_coords(region_id=0)

        # Panel A: Profile Discovery Accuracy (table, no data args)
        render_data_panel(3, "A", panel_a_profile_discovery, figsize=(6, 5))

        # Panel B: Xenium Proportion Benchmark
        render_data_panel(3, "B", panel_b_xenium_proportions, xenium_summary, figsize=(6, 4))

        # Panel C: Simulated Proportion Benchmark
        render_data_panel(3, "C", panel_c_simulated_proportions, sim_props, figsize=(6, 4))

        # Panel D: GEX Benchmark
        render_data_panel(3, "D", panel_d_gex_benchmark, sim_gex, figsize=(5, 4))

        # Panel E: Scatter (GT vs Predicted)
        render_data_panel(3, "E", panel_e_scatter, gt, pred, figsize=(5, 5))

        # Panel F: Spatial Comparison Maps
        render_data_panel(3, "F", panel_f_spatial_comparison, gt, pred, coords, figsize=(6, 5))

    except Exception as e:
        print(f"    WARNING: Could not render Fig3 data panels: {e}")
        import traceback
        traceback.print_exc()


# =========================================================================
# Figure 4: Midkine/ESR1 Case Study (8 data panels)
# =========================================================================
def render_panel_c_estrogene(fig_num, panel_letter, panel_func, gex_df,
                              d538g_labels, figsize=(6, 5)):
    """Special renderer for panel_c_estrogene which takes (fig, gs_cell, ...)
    instead of (ax, ...).  Creates a figure with a GridSpec and passes a
    GridSpec cell so the function can split it internally."""
    from matplotlib.gridspec import GridSpec

    dest_dir = EXPORT_DIR / f"Figure{fig_num}"
    dest_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=figsize)
    fig.set_layout_engine('none')  # Prevent ZeroDivisionError in nested GridSpec
    gs = GridSpec(1, 1, figure=fig,
                  left=0.10, right=0.92, top=0.92, bottom=0.08)
    panel_func(fig, gs[0, 0], gex_df, d538g_labels)
    export_panel(fig, None, str(dest_dir / f"panel_{panel_letter}"))


def export_figure4():
    print("\n  Figure 4: Midkine/ESR1 Case Study")

    try:
        from generate_figure4 import (
            panel_a_basal_ck, panel_b_mutation_spatial,
            panel_c_estrogene, panel_d_pathway_dotplot,
            panel_e_mdk_spatial, panel_f_commot,
            panel_g_elisa, panel_h_if,
            load_visium_adata, load_gene_expression,
            compute_basal_score, compute_d538g_score,
            classify_d538g_spots, compute_mdk_score,
        )

        # Load shared data
        adata = load_visium_adata()
        gex_df = load_gene_expression()
        basal_score = compute_basal_score(gex_df) if gex_df is not None else None
        d538g_score = compute_d538g_score(gex_df) if gex_df is not None else None
        d538g_labels = classify_d538g_spots(d538g_score) if d538g_score is not None else None
        mdk_score = compute_mdk_score(gex_df) if gex_df is not None else None

        # Panel A: Basal CK spatial (H&E)
        render_data_panel(4, "A", panel_a_basal_ck, adata, basal_score, figsize=(5, 5))

        # Panel B: D538G mutation signal (H&E)
        render_data_panel(4, "B", panel_b_mutation_spatial, adata, d538g_score, figsize=(5, 5))

        # Panel C: EstroGene split violin + heatmap (SPECIAL: fig + GridSpec)
        render_panel_c_estrogene(4, "C", panel_c_estrogene, gex_df, d538g_labels,
                                  figsize=(6, 5))

        # Panel D: Pathway dot plot (no data args)
        render_data_panel(4, "D", panel_d_pathway_dotplot, figsize=(5, 5))

        # Panel E: MDK spatial (H&E)
        render_data_panel(4, "E", panel_e_mdk_spatial, adata, mdk_score, figsize=(5, 5))

        # Panel F: COMMOT signaling (no data args)
        render_data_panel(4, "F", panel_f_commot, figsize=(5, 4))

        # Panel G: ELISA placeholder
        render_data_panel(4, "G", panel_g_elisa, figsize=(4, 3))

        # Panel H: IF placeholder
        render_data_panel(4, "H", panel_h_if, figsize=(4, 3))

    except Exception as e:
        print(f"    WARNING: Could not render Fig4 panels: {e}")
        import traceback
        traceback.print_exc()


# =========================================================================
# Figure 5: Cross-Sample Integration (8 data panels)
# =========================================================================
def export_figure5():
    print("\n  Figure 5: Cross-Sample Integration")

    try:
        from generate_figure5 import (
            panel_a_proportion_shift,
            panel_b_program_summary,
            panel_c_response_dotplot,
            panel_d_network,
            panel_e_colocalization,
            panel_f_exclusion,
            panel_g_morans_violin,
            panel_h_spatial_comparison,
        )

        # All Figure 5 panels load data internally -- just pass (ax)
        render_data_panel(5, "A", panel_a_proportion_shift, figsize=(6, 4))
        render_data_panel(5, "B", panel_b_program_summary, figsize=(6, 5))
        render_data_panel(5, "C", panel_c_response_dotplot, figsize=(5, 5))
        render_data_panel(5, "D", panel_d_network, figsize=(6, 6))
        render_data_panel(5, "E", panel_e_colocalization, figsize=(6, 4))
        render_data_panel(5, "F", panel_f_exclusion, figsize=(6, 4))
        render_data_panel(5, "G", panel_g_morans_violin, figsize=(6, 4))
        render_data_panel(5, "H", panel_h_spatial_comparison, figsize=(6, 5))

    except Exception as e:
        print(f"    WARNING: Could not render Fig5 panels: {e}")
        import traceback
        traceback.print_exc()


# =========================================================================
# Supplementary Figure S2: DE & Pathway Enrichment (2 data panels)
# =========================================================================
def export_supp_figure2():
    print("\n  Supp Figure S2: DE & Pathway Enrichment")

    try:
        from generate_supp_figure2 import panel_a_volcano, panel_b_pathway_dotplot

        # Panel A: Volcano plot
        render_data_panel("SuppS2", "A", panel_a_volcano, figsize=(6, 5))

        # Panel B: Pathway dot plot
        render_data_panel("SuppS2", "B", panel_b_pathway_dotplot, figsize=(6, 5))

    except Exception as e:
        print(f"    WARNING: Could not render Supp Fig S2 panels: {e}")
        import traceback
        traceback.print_exc()


# =========================================================================
# Main
# =========================================================================
def main():
    print("=" * 60)
    print("Exporting individual panels (PNG 300 DPI + SVG)")
    print(f"Output: {EXPORT_DIR}")
    print("=" * 60)

    # Create export directories for main figures
    for i in range(1, 6):
        (EXPORT_DIR / f"Figure{i}").mkdir(parents=True, exist_ok=True)

    # Create export directories for supplementary figures
    (EXPORT_DIR / "SuppS2").mkdir(parents=True, exist_ok=True)

    # Export main figures
    export_figure1()
    export_figure2()
    export_figure3()
    export_figure4()
    export_figure5()

    # Export supplementary figures
    export_supp_figure2()

    print("\n" + "=" * 60)
    print("Panel export complete")
    print("=" * 60)

    # Summary for main figures
    for i in range(1, 6):
        fig_dir = EXPORT_DIR / f"Figure{i}"
        files = sorted(fig_dir.glob("*"))
        print(f"\n  Figure{i}/  ({len(files)} files)")
        for f in files:
            print(f"    {f.name}")

    # Summary for supplementary figures
    for supp_dir in ["SuppS2"]:
        fig_dir = EXPORT_DIR / supp_dir
        if fig_dir.exists():
            files = sorted(fig_dir.glob("*"))
            print(f"\n  {supp_dir}/  ({len(files)} files)")
            for f in files:
                print(f"    {f.name}")


if __name__ == "__main__":
    main()
