#!/usr/bin/env python
"""
Test: Cell types with markers that DON'T EXIST in the data should fail validation.
This tests the marker-specific fit quality approach.
"""
import os
import sys

# Add CITEgeist package directory to path
citegeist_pkg_dir = "/ihome/alee/alc376/alc376_bgfs/CITEgeist"
sys.path.insert(0, citegeist_pkg_dir)

import scanpy as sc
from CITEgeist.model.citegeist_model import CitegeistModel

# OVERSPECIFIED with NONEXISTENT markers
cell_type_profiles_truly_fake = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
    # FAKE with markers that DON'T EXIST
    "FakeCellType1": {"Major": ["DOES_NOT_EXIST_1", "DOES_NOT_EXIST_2"]},
    "FakeCellType2": {"Major": ["ALSO_NOT_REAL_1", "ALSO_NOT_REAL_2"]},
}

# Load data
input_folder = "/ihome/alee/alc376/alc376_bgfs/CITEgeist/replicates/high_seg"
sample_name = "Wu_rep_0"

adata_cite_path = os.path.join(input_folder, "h5ad_objects", f"{sample_name}_CITE.h5ad")
adata_gex_path = os.path.join(input_folder, "h5ad_objects", f"{sample_name}_GEX.h5ad")

print(f"Loading data...")
adata_cite = sc.read_h5ad(adata_cite_path)
adata_gex = sc.read_h5ad(adata_gex_path)

# Initialize model
output_folder = "/ihome/alee/alc376/alc376_bgfs/CITEgeist/test_validation_output"
os.makedirs(output_folder, exist_ok=True)

model = CitegeistModel(
    sample_name=sample_name + "_truly_fake",
    output_folder=output_folder,
    simulation=True,
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite,
)

print("\n" + "="*60)
print("TEST: Overspecified with NONEXISTENT markers")
print("="*60)
print("Real cell types have real markers")
print("Fake types: FakeCellType1, FakeCellType2 with NONEXISTENT markers")
model.load_cell_profile_dict(cell_type_profiles_truly_fake)

# Preprocess
model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

print("\nRunning cell proportion optimization...")
print("Expected: Fit quality validation should FAIL (fake types have R²=0.0)")
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
    print("   Fake cell types with nonexistent markers should have R²=0.0 < 0.25")
    
except ValueError as e:
    error_msg = str(e)
    
    if "Poor protein-proportion agreement" in error_msg or "R²" in error_msg:
        if "FakeCellType1" in error_msg or "FakeCellType2" in error_msg:
            print("\n✅ TEST PASSED: Marker-specific validation caught fake cell types!")
            print(f"   Error: {error_msg}")
        else:
            print("\n⚠️  R² error but didn't mention fake types:")
            print(f"   {error_msg}")
    else:
        print(f"\n⚠️  Got ValueError but unexpected message:")
        print(f"   {error_msg}")

print("\n" + "="*60)
print("Test complete")
print("="*60)
