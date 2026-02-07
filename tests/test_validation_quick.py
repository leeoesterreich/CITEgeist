#!/usr/bin/env python
"""
Quick validation test: Run CITEgeist with underspecified cell types.
This should trigger the Unknown > 5% validation error.
"""
import os
import sys

# Add CITEgeist package directory to path to import LOCAL modified version
citegeist_pkg_dir = "/ihome/alee/alc376/alc376_bgfs/CITEgeist"
sys.path.insert(0, citegeist_pkg_dir)

import scanpy as sc
from CITEgeist.model.citegeist_model import CitegeistModel

# UNDERSPECIFIED: Only 3 cell types (should trigger Unknown > 5% error)
cell_type_profiles_underspec = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
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
    sample_name=sample_name,
    output_folder=output_folder,
    simulation=True,
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite,
)

# Load UNDERSPECIFIED cell profile dictionary (Unknown should auto-add)
print("\n" + "="*60)
print("TEST: Underspecified cell types (only 3 types)")
print("="*60)
print("Cell types defined:", list(cell_type_profiles_underspec.keys()))
model.load_cell_profile_dict(cell_type_profiles_underspec)

# Preprocess
model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

# Skip Gurobi registration - module load sets GRB_LICENSE_FILE env var
# model.register_gurobi("/ihome/crc/install/gurobi/gurobi1203/linux64/lic/gurobi.lic")

print("\nRunning cell proportion optimization...")
print("Expected: Validation should FAIL with 'Unknown > 5%' error")
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
    print("\n❌ TEST FAILED: Validation did not trigger!")
    print("   Expected Unknown > 5% error, but optimization completed successfully")
    
except ValueError as e:
    error_msg = str(e)
    if "Unknown cell type has mean proportion" in error_msg and "too few cell types" in error_msg:
        print("\n✅ TEST PASSED: Validation correctly detected underspecified cell types!")
        print(f"   Error message: {error_msg}")
    else:
        print(f"\n⚠️  Got ValueError but unexpected message:")
        print(f"   {error_msg}")

print("\n" + "="*60)
print("Test complete")
print("="*60)
