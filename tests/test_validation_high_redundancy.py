#!/usr/bin/env python
"""
Test: Overspecified with redundancy > 20% threshold.
More aggressive overspecification: 7 real + 3 fake types for 30% redundancy.
"""
import os
import sys

# Add CITEgeist package directory to path
citegeist_pkg_dir = "/ihome/alee/alc376/alc376_bgfs/CITEgeist"
sys.path.insert(0, citegeist_pkg_dir)

import scanpy as sc
from CITEgeist.model.citegeist_model import CitegeistModel

# OVERSPECIFIED: 7 real types + 3 fake with shared markers
# 3/10 = 30% redundancy > 20% threshold
cell_type_profiles_high_redundancy = {
    # Real types (7 of 9)
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
    # Fake types (3) - cannibalize real types
    "FakeType1": {"Major": ["T-cells_Protein_1", "Myeloid_Protein_2"]},
    "FakeType2": {"Major": ["B-cells_Protein_1", "Endothelial_Protein_2"]},
    "FakeType3": {"Major": ["Cancer Epithelial_Protein_1", "PVL_Protein_2"]},
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
    sample_name=sample_name + "_high_redundancy",
    output_folder=output_folder,
    simulation=True,
    gene_expression_adata=adata_gex,
    antibody_capture_adata=adata_cite,
)

print("\n" + "="*60)
print("TEST: High redundancy (7 real + 3 fake = 30% redundancy)")
print("="*60)
print("Expected redundancy: 3/(7+3) = 30% > 20% threshold")
print("Should FAIL with redundancy error")
model.load_cell_profile_dict(cell_type_profiles_high_redundancy)

# Preprocess
model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
model.preprocess_gex(target_sum=10000)
model.preprocess_antibody()

print("\nRunning cell proportion optimization...")
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
    print("   Expected redundancy > 20% error")
    
except ValueError as e:
    error_msg = str(e)
    
    if "redundant" in error_msg.lower():
        print("\n✅ TEST PASSED: Redundancy validation detected overspecification!")
        print(f"   Error: {error_msg}")
    else:
        print(f"\n⚠️  Got ValueError but unexpected message:")
        print(f"   {error_msg}")

print("\n" + "="*60)
print("Test complete")
print("="*60)
