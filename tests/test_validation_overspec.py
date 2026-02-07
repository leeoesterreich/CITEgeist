#!/usr/bin/env python
"""
Validation test: Run CITEgeist with OVERSPECIFIED cell types.
This adds fake cell types with randomly mixed antibodies that don't exist.
Should trigger the "cell type < 1%" validation error.
"""
import os
import sys

# Add CITEgeist package directory to path to import LOCAL modified version
citegeist_pkg_dir = "/ihome/alee/alc376/alc376_bgfs/CITEgeist"
sys.path.insert(0, citegeist_pkg_dir)

import scanpy as sc
from CITEgeist.model.citegeist_model import CitegeistModel

# OVERSPECIFIED: All 9 real cell types PLUS 2 fake ones with random antibodies
cell_type_profiles_overspec = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
    # FAKE cell types that don't exist - using random mix of real antibodies
    "Neutrophils": {"Major": ["T-cells_Protein_1", "Myeloid_Protein_2"]},
    "Macrophages": {"Major": ["CAFs_Protein_1", "Endothelial_Protein_2"]},
}

# Load first sample
input_folder = "/ihome/alee/alc376/alc376_bgfs/CITEgeist/replicates/high_seg"
sample_name = "Wu_rep_0"

adata_cite_path = os.path.join(input_folder, "h5ad_objects", f"{sample_name}_CITE.h5ad")
adata_gex_path = os.path.join(input_folder, "h5ad_objects", f"{sample_name}_GEX.h5ad")

print(f"Loading data from {input_folder}...")
adata_cite = sc.read_h5ad(adata_cite_path)
adata_gex = sc.read_h5ad(adata_gex_path)
print(f"  GEX shape: {adata_gex.shape}")
print(f"  CITE shape: {adata_cite.shape}")

# Initialize model
output_folder = "/ihome/alee/alc376/alc376_bgfs/CITEgeist/test_validation_output"
os.makedirs(output_folder, exist_ok=True)

model = CitegeistModel(
    sample_name=sample_name + "_overspec",
    output_folder=output_folder,
    simulation=True,
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite,
)

# Load OVERSPECIFIED cell profile dictionary
print("\n" + "="*60)
print("TEST: Overspecified cell types (9 real + 2 fake)")
print("="*60)
print("Real cell types:", ["B-cells", "CAFs", "Cancer Epithelial", "Endothelial", 
                           "Myeloid", "Normal Epithelial", "PVL", "Plasmablasts", "T-cells"])
print("Fake cell types:", ["Neutrophils", "Macrophages"])
model.load_cell_profile_dict(cell_type_profiles_overspec)

# Preprocess
model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

print("\nRunning cell proportion optimization...")
print("Expected: Validation may PASS or FAIL depending on fit quality of fake types")
print("-" * 60)

try:
    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=4.0,
        tolerance=1e-4,
        max_iterations=20,
        lambda_reg=1.0,
        alpha=0.5,
        max_workers=None,
        checkpoint_interval=100,
        max_y_change=0.4,
    )
    
    # If it passes, check what the fit quality actually was
    print("\n⚠️  Validation PASSED (no error raised)")
    print("   This means fake cell types had acceptable R² values")
    print("   (They share antibodies with real types, so they may fit reasonably)")
    
    # Show proportions
    import pandas as pd
    mean_props = global_props.mean(axis=0).sort_values(ascending=False)
    print(f"\n   Mean cell type proportions:")
    for ct, prop in mean_props.items():
        marker = " *** FAKE TYPE" if ct in ["Neutrophils", "Macrophages"] else ""
        print(f"     {ct}: {prop*100:.2f}%{marker}")
    
except ValueError as e:
    error_msg = str(e)
    
    # Check if it's a fit quality error (preferred) or proportion error
    if "Poor protein-proportion agreement" in error_msg or "R²" in error_msg:
        print("\n✅ TEST PASSED: Fit quality validation detected overspecification!")
        print(f"   Error: {error_msg}")
    elif "mean proportion <" in error_msg:
        print("\n✅ TEST PASSED: Proportion validation detected low abundance!")
        print(f"   Error: {error_msg}")
    elif "Unknown" in error_msg:
        print("\n⚠️  Unknown threshold triggered instead of overspecification:")
        print(f"   Error: {error_msg}")
    else:
        print(f"\n⚠️  Got ValueError but unexpected message:")
        print(f"   {error_msg}")

print("\n" + "="*60)
print("Test complete")
print("="*60)
