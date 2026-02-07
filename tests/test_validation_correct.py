#!/usr/bin/env python
"""
Validation test: Run CITEgeist with CORRECT cell types (baseline).
All 9 real cell types should be present. Validation should PASS.
"""
import os
import sys

# Add CITEgeist package directory to path to import LOCAL modified version
citegeist_pkg_dir = "/ihome/alee/alc376/alc376_bgfs/CITEgeist"
sys.path.insert(0, citegeist_pkg_dir)

import scanpy as sc
from CITEgeist.model.citegeist_model import CitegeistModel

# CORRECT: All 9 real cell types
cell_type_profiles_correct = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
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
    sample_name=sample_name + "_correct",
    output_folder=output_folder,
    simulation=True,
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite,
)

# Load CORRECT cell profile dictionary
print("\n" + "="*60)
print("TEST: Correct cell types (all 9 real types)")
print("="*60)
print("Cell types defined:", list(cell_type_profiles_correct.keys()))
model.load_cell_profile_dict(cell_type_profiles_correct)

# Preprocess
model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

print("\nRunning cell proportion optimization...")
print("Expected: Validation should PASS (Unknown < 5%, all cell types > 1%)")
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
    print("\n✅ TEST PASSED: Validation succeeded with correct cell types!")
    print(f"   Shape of global proportions: {global_props.shape}")
    print(f"   Cell types: {list(global_props.columns)}")
    
    # Show mean proportions
    import pandas as pd
    mean_props = global_props.mean(axis=0).sort_values(ascending=False)
    print(f"\n   Mean cell type proportions:")
    for ct, prop in mean_props.items():
        print(f"     {ct}: {prop*100:.2f}%")
    
except ValueError as e:
    error_msg = str(e)
    print(f"\n❌ TEST FAILED: Validation raised unexpected error:")
    print(f"   {error_msg}")

print("\n" + "="*60)
print("Test complete")
print("="*60)
