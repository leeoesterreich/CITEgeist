
import os
import sys
import argparse
from datetime import datetime
from typing import Dict, Any, List, Optional

import numpy as np
import scanpy as sc
import pandas as pd
import squidpy as sq
import scipy.sparse

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Now import using the full package path
from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    # Module 2c: Profile selection
    select_profiles,
    # Module 4: Protein-anchored program discovery
    discover_anchored_programs,
    store_results_in_adata,
)

def main():

    cell_profiles = {
        "Cancer Cells": {
            "Major": ["EPCAM-1"],
            "Minor": ["SDC1-1",  "KRT5-1"]  # CD138 - possible cancer stem cell marker
        },
        "Macrophages": {
            "Major": ["CD68-1" ],  # General macrophage and M2-polarized macrophages
            "Minor": ["CD14-1"]  # Monocyte/macrophage lineage marker
        },
        "CD4 T Cells": {
            "Major": ["CD3E-1", "CD4-1"],  # General, Helper, and Cytotoxic T cells
        },
        "CD8 T Cells": {
            "Major": ["CD3E-1", "CD8A-1"],  # General, Helper, and Cytotoxic T cells
        },
        "B Cells": {
            "Major": ["MS4A1-1", "CD19-1"],  # General B cell markers and developmental marker
        },
        "Endothelial Cells": {
            "Major": ["PECAM1-1"],  # CD31 - endothelial cell marker
        },
        "Fibroblasts": {
            "Major": ["ACTA2-1"],  # α-SMA - myofibroblast marker, indicates activated stroma
        }
    }

    # Reload the model module if you're actively developing

    parser = argparse.ArgumentParser(description='Run CITEgeist on a single sample.')
    parser.add_argument('--path', type=str, required=True, help='Path to the sample')
    parser.add_argument('--radius', type=float, required=True, help='Radius for neighbor detection')
    parser.add_argument('--min_counts', type=float, default=100, required=True,
                       help='Minimum counts for filtering')
    parser.add_argument('--output_folder', type=str, default='output', help='Output folder')
    # Auto-profile discovery via spatial colocalization pipeline
    parser.add_argument('--auto-profiles', action='store_true', default=True,
                        help='Use auto-discovered profiles via spatial colocalization pipeline (default: True)')
    parser.add_argument('--manual-profiles', action='store_true', default=False,
                        help='Use manually defined profiles instead of auto-discovery')
    parser.add_argument('--top-k', type=int, default=3,
                        help='Mutual top-k for profile pair sparsification (default: 3, tested value)')
    parser.add_argument('--discovery-seed', type=int, default=1234,
                        help='Random seed for profile discovery reproducibility')
    parser.add_argument('--n-permutations', type=int, default=999,
                        help='Number of permutations for significance testing (default: 999)')
    parser.add_argument('--fdr-threshold', type=float, default=0.05,
                        help='FDR threshold for significant markers/pairs (default: 0.05)')
    # Module 2c: Spatial variance-based profile selection
    parser.add_argument('--variance-target', type=float, default=0.90,
                        help='Target fraction of spatial variance to explain (default: 0.90)')
    parser.add_argument('--min-marginal-gain', type=float, default=0.005,
                        help='Minimum marginal variance gain per profile (default: 0.5%%)')
    parser.add_argument('--validation-warn-only', action='store_true', default=False,
                        help='Only warn (do not error) on validation failures for auto-profiles')

    args = parser.parse_args()


    adata = sq.read.visium(args.path, counts_file='filtered_feature_bc_matrix.h5', load_images=True, gex_only=False)
    sample_name = args.path.split('/')[-3]

    

    print("Sample name: ", sample_name)

        # Initialize the model
    model = CitegeistModel(sample_name=sample_name, adata=adata, output_folder=args.output_folder)

    # Preprocess and run models
    # Split into gene expression and antibody capture datasets
    model.split_adata()

    model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=args.min_counts)

    model.copy_gex_to_protein_adata()

    # Save raw antibody data BEFORE preprocessing (Module 1-2c need raw data)
    # CLR transformation changes kurtosis/Moran's I characteristics
    X_antibody_raw = model.antibody_capture_adata.X.copy()
    if scipy.sparse.issparse(X_antibody_raw):
        X_antibody_raw = X_antibody_raw.toarray()

    # Preprocess datasets
    model.preprocess_gex()
    model.preprocess_antibody()  # Applies Winsorizing + CLR transformation
    print(model)

    # Load or discover cell profiles
    # Default to auto-profiles unless --manual-profiles is specified
    use_auto_profiles = args.auto_profiles and not args.manual_profiles

    # When using auto-profiles, default to warn-only validation since
    # rare cell types are valid and shouldn't cause failures
    validation_warn_only = args.validation_warn_only or use_auto_profiles

    if use_auto_profiles:
        # Use spatial colocalization pipeline for auto profile discovery
        print("\n=== Running Spatial Colocalization Pipeline ===")

        # Get spatial coordinates and marker names
        coords = model.antibody_capture_adata.obsm.get('spatial', None)
        if coords is None:
            raise ValueError("Spatial coordinates required for auto-profile discovery")

        # Check for NaN/inf in coordinates and filter if needed
        finite_mask = np.all(np.isfinite(coords), axis=1)
        if not np.all(finite_mask):
            n_invalid = (~finite_mask).sum()
            print(f"  Warning: Filtering {n_invalid} spots with non-finite coordinates")
            coords = coords[finite_mask]
            X_antibody_raw = X_antibody_raw[finite_mask]
            # Also filter the model's adata objects to keep them in sync
            valid_obs = model.antibody_capture_adata.obs_names[finite_mask]
            model.antibody_capture_adata = model.antibody_capture_adata[valid_obs, :].copy()
            model.gene_expression_adata = model.gene_expression_adata[valid_obs, :].copy()

        marker_names = list(model.antibody_capture_adata.var_names)

        # NOTE: Use RAW antibody data (X_antibody_raw) for Module 1-2c
        # The CLR transformation changes kurtosis/Moran's I characteristics
        # significantly and causes adaptive thresholds to fail.

        # Module 1: Identify spatially interesting markers (RAW DATA)
        print("\nModule 1: Identifying interesting markers (using raw data)...")
        marker_result = identify_interesting_markers(
            X=X_antibody_raw,
            coords=coords,
            marker_names=marker_names,
            morans_k=8,
            morans_n_perm=args.n_permutations,
            seed=args.discovery_seed,
            verbose=True,
        )
        interesting_markers = marker_result.interesting_markers
        print(f"  Found {len(interesting_markers)} spatially interesting markers")

        if len(interesting_markers) < 2:
            print("  Warning: Not enough interesting markers. Using manual profiles.")
            model.load_cell_profile_dict(cell_profiles)
        else:
            # Module 2a: Analyze marker colocalization (RAW DATA)
            print("\nModule 2a: Analyzing marker colocalization (using raw data)...")
            coloc_result = analyze_marker_colocalization(
                X=X_antibody_raw,
                coords=coords,
                marker_names=marker_names,
                markers_to_analyze=interesting_markers,
                neighbor_k=8,
                n_permutations=args.n_permutations,
                seed=args.discovery_seed,
                verbose=True,
            )
            print(f"  Found {len(coloc_result.pairs)} significant colocalization pairs")

            # Module 2b: Discover profiles
            print("\nModule 2b: Discovering profiles...")
            discovery_result = discover_profiles(
                colocalization_result=coloc_result,
                fdr_alpha=args.fdr_threshold,
                top_k=args.top_k,
                seed=args.discovery_seed,
                verbose=True,
            )
            print(f"  Discovered {len(discovery_result.profiles)} candidate profiles")

            for i, profile in enumerate(discovery_result.profiles):
                print(f"  {i+1}. {profile}")

            # Module 2c: Select profiles by spatial variance
            print("\nModule 2c: Selecting profiles by spatial variance...")
            selection_result = select_profiles(
                X=X_antibody_raw,
                coords=coords,
                marker_names=marker_names,
                profiles=discovery_result.profiles,
                interesting_markers=interesting_markers,
                colocalization_result=coloc_result,
                min_spatial_explained=args.variance_target,
                min_marginal_gain=args.min_marginal_gain,
                verbose=True,
            )
            selected_profiles = selection_result.selected_profiles
            n_selected = selection_result.optimal_n
            total_ve = float(selection_result.variance_explained[-1]) if n_selected > 0 else 0.0
            print(f"\nModule 2c: Selected {n_selected} profiles (explains {total_ve:.1%} spatial variance)")
            print(f"  Stopping reason: {selection_result.stopping_reason}")

            for i, profile in enumerate(selected_profiles):
                print(f"  {i+1}. {profile}")

            # Convert to cell_profile_dict format
            auto_cell_profiles = {}
            for i, profile in enumerate(selected_profiles):
                profile_name = f"Profile_{i+1}"
                markers_list = list(profile) if not isinstance(profile, list) else profile
                auto_cell_profiles[profile_name] = {"Major": markers_list}

            model.load_cell_profile_dict(auto_cell_profiles)
    else:
        # Use manual cell profile dictionary
        model.load_cell_profile_dict(cell_profiles)

    # Skip explicit Gurobi registration - module load sets GRB_LICENSE_FILE env var
    # model.register_gurobi("/ihome/crc/install/gurobi/gurobi1102/linux64/lic/gurobi.lic")

    global_cell_type_proportions_df, finetuned_cell_type_proportions_df = model.run_cell_proportion_model(
        radius=args.radius,
        validation_warn_only=validation_warn_only,
    )

    # Plot cell proportions (Append cell proportions) 
    model.append_proportions_to_adata(key='global')

    # Run gene expression model (if needed)
    # model.run_gex_model()
    prop_adata = model.get_adata()

    # Plot cell proportions (Append cell proportions) 
    model.append_proportions_to_adata(key='finetuned')

    prop_adata = model.get_adata()

    # Plot cell proportions for the first few profiles
    # Use profile names from the model's loaded cell profile dict
    profile_names = list(model.cell_profile_dict.keys()) if model.cell_profile_dict else []
    plot_profiles = profile_names[:min(3, len(profile_names))]  # Plot first 3 profiles
    if plot_profiles:
        try:
            sc.pl.spatial(prop_adata, color=plot_profiles, ncols=3,
                         save=f'_{sample_name}_proportions.png')
        except Exception as e:
            print(f"Could not plot proportions: {e}")


    pass1_results = model.run_cell_expression_pass1(
                radius=args.radius, 
                max_workers=None, 
                checkpoint_interval=100, 
                output_dir="checkpoints", 
                rerun=True
            )

    # Plot cell proportions (Append cell proportions) 
    model.append_gex_to_adata(pass_number=1)


    prop_gex_adata = model.get_adata()

    print(prop_gex_adata)

    # =========================================================================
    # Module 4: Protein-Anchored Spatial Program Discovery
    # =========================================================================
    print("\n=== Module 4: Protein-Anchored Program Discovery ===")

    # Use finetuned proportions from Module 3 as anchor weights
    cell_type_proportions = finetuned_cell_type_proportions_df

    # Get the profile discovery result if we used auto-profiles
    profile_result = discovery_result if use_auto_profiles and 'discovery_result' in dir() else None

    # Run Module 4
    program_result = discover_anchored_programs(
        adata=model.gene_expression_adata,
        cell_type_proportions=cell_type_proportions,
        profile_discovery_result=profile_result,
        protein_adata=model.antibody_capture_adata,
        K_programs=3,  # Discover 3 programs per cell type
        lambda_spatial=0.1,
        lambda_sparsity=0.01,
        min_proportion_threshold=0.1,
        validate_with_proteins=True,
        top_n_genes=50,
        random_state=42,
    )

    # Print summary
    print(program_result.summary())

    # Store results in AnnData
    store_results_in_adata(prop_gex_adata, program_result)

    # Print detailed results for each anchor
    print("\n=== Discovered Programs ===")
    for anchor_name, anchor_result in program_result.results_by_anchor.items():
        print(f"\n{anchor_name} ({anchor_result.n_spots_used} spots):")

        # Print protein validation correlations if available
        if not anchor_result.protein_correlations.empty:
            print(f"  Anchor proteins: {anchor_result.anchor_proteins}")
            if 'anchor_correlation' in anchor_result.protein_correlations.columns:
                print(f"  Anchor validation scores:")
                for k, row in anchor_result.protein_correlations.iterrows():
                    validated = "✓" if row.get('validated', False) else "✗"
                    print(f"    Program {k}: r={row['anchor_correlation']:.3f} {validated}")

        # Print top genes per program
        for prog in anchor_result.programs:
            moran_sig = "*" if prog.spatial_moran_pvalue < 0.05 else ""
            print(f"  Program {prog.program_id}: Moran's I={prog.spatial_moran_i:.3f}{moran_sig}")
            print(f"    Top genes: {', '.join(prog.top_genes[:10])}")

    # Save program results to CSV
    output_path = os.path.join(args.output_folder, f"{sample_name}_module4_programs.csv")
    program_result.to_dataframe().to_csv(output_path, index=False)
    print(f"\nSaved program summary to: {output_path}")

    # Visualize program activities spatially (optional)
    for anchor_name in list(program_result.results_by_anchor.keys())[:3]:  # First 3 anchors
        key = f'X_anchored_programs_{anchor_name}'
        if key in prop_gex_adata.obsm:
            # Add program activities to obs for plotting
            H = prop_gex_adata.obsm[key]
            for k in range(min(3, H.shape[1])):  # First 3 programs
                col_name = f'{anchor_name}_program_{k}'
                prop_gex_adata.obs[col_name] = H[:, k]

            # Plot
            try:
                prog_cols = [f'{anchor_name}_program_{k}' for k in range(min(3, H.shape[1]))]
                sc.pl.spatial(prop_gex_adata, color=prog_cols, ncols=3,
                             title=[f'{anchor_name} P{k}' for k in range(len(prog_cols))],
                             save=f'_{sample_name}_{anchor_name}_programs.png')
            except Exception as e:
                print(f"Could not plot {anchor_name}: {e}")

    print("\n=== Module 4 Complete ===")


if __name__ == "__main__":
    main()
