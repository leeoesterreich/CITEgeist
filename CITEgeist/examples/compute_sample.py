"""
This script processes spatial transcriptomics data using the CITEgeist model.

It performs the following main steps:
1. Loads and preprocesses Visium spatial transcriptomics data
2. Applies cell type deconvolution using predefined cell profiles
3. Runs the CITEgeist model to analyze spatial gene expression patterns
4. Generates visualizations and outputs results

The script requires a Gurobi license for optimization computations.

Example usage:
    python compute_sample.py --path /path/to/visium/data --radius 400 --min_counts 100 
        --output_folder output --gurobi_license /path/to/gurobi.lic

Requirements:
    - Gurobi optimizer with valid license
    - Python packages: scanpy, squidpy, numpy, pandas
    - CITEgeist package installed in environment
"""

import os
import sys
import argparse
from datetime import datetime
from typing import Dict, Any, List, Optional

import numpy as np
import scanpy as sc
import pandas as pd
import squidpy as sq

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

# Now import using the full package path
from model.citegeist_model import CitegeistModel

def verify_gurobi_license(license_path: str) -> None:
    """
    Verify that the Gurobi license file exists at the specified path.
    
    Args:
        license_path (str): Path to the Gurobi license file
        
    Raises:
        FileNotFoundError: If the license file does not exist at the specified path
    """
    if not os.path.isfile(license_path):
        raise FileNotFoundError(f"Gurobi license file not found at: {license_path}")

def main():
    """
    Main function to run the CITEgeist analysis pipeline.
    
    This function:
    1. Parses command line arguments
    2. Loads and preprocesses the spatial transcriptomics data
    3. Initializes and runs the CITEgeist model
    4. Generates visualizations and saves results
    """

    # Define cell type profiles with their characteristic markers
    # Each cell type has Major (strong) and optional Minor (supporting) markers
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

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run CITEgeist on a single sample.')
    parser.add_argument('--path', type=str, required=True, help='Path to the Visium data folder')
    parser.add_argument('--radius', type=float, required=True, help='Spatial radius (in μm) for neighbor detection')
    parser.add_argument('--min_counts', type=float, default=100, required=True, 
                       help='Minimum count threshold for filtering low-quality spots')
    parser.add_argument('--output_folder', type=str, default='output', help='Directory to save output files')
    parser.add_argument('--gurobi_license', type=str, required=True, help='Path to Gurobi license file (gurobi.lic)')

    args = parser.parse_args()

    # Verify that the Gurobi license file exists
    verify_gurobi_license(args.gurobi_license)

    # Load the Visium dataset
    adata = sq.read.visium(args.path, counts_file='filtered_feature_bc_matrix.h5', load_images=True, gex_only=False)
    sample_name = args.path.split('/')[-3]
    print("Sample name: ", sample_name)

    # Initialize the CITEgeist model
    model = CitegeistModel(sample_name=sample_name, adata=adata, output_folder='output')

    # Load predefined cell type profiles
    model.load_cell_profile_dict(cell_profiles)

    # Preprocess the data
    model.split_adata()  # Split into gene expression and antibody capture datasets
    model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=args.min_counts)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody()
    print(model)

    # Initialize Gurobi optimizer
    model.register_gurobi(args.gurobi_license)

    # Run cell type proportion estimation
    global_cell_type_proportions_df, finetuned_cell_type_proportions_df = model.run_cell_proportion_model(radius=args.radius)

    # Add results to AnnData object and visualize
    model.append_proportions_to_adata(key='global')
    model.append_proportions_to_adata(key='finetuned')
    prop_adata = model.get_adata()
    
    # Generate spatial plot of Cancer Cell proportions
    sc.pl.spatial(prop_adata, color="Cancer Cells")

    # Run detailed cell expression analysis
    pass1_results = model.run_cell_expression_pass1(
                radius=args.radius, 
                max_workers=None, 
                checkpoint_interval=100, 
                output_dir="checkpoints", 
                rerun=True
            )

    # Add gene expression results to AnnData object
    model.append_gex_to_adata(pass_number=1)
    prop_gex_adata = model.get_adata()
    print(prop_gex_adata)


if __name__ == "__main__":
    main()
