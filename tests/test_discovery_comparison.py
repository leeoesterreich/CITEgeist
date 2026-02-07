#!/usr/bin/env python
"""Compare discover_profiles vs discover_profiles_continuous on simulated data."""
import sys
sys.path.insert(0, '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist')

import numpy as np
import scanpy as sc
from CITEgeist.model import (
    identify_interesting_markers,
    analyze_marker_colocalization,
)
from CITEgeist.model.spatial_colocalization import (
    discover_profiles,
    discover_profiles_continuous,
)

# Load simulated data
adata = sc.read_h5ad('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_0_CITE.h5ad')

X = adata.X
if hasattr(X, 'toarray'):
    X = X.toarray()
coords = adata.obsm['spatial']
marker_names = list(adata.var_names)

# Ground truth profiles (2 markers per cell type)
GROUND_TRUTH = {
    "B-cells": ["B-cells_Protein_1", "B-cells_Protein_2"],
    "CAFs": ["CAFs_Protein_1", "CAFs_Protein_2"],
    "Cancer Epithelial": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"],
    "Endothelial": ["Endothelial_Protein_1", "Endothelial_Protein_2"],
    "Myeloid": ["Myeloid_Protein_1", "Myeloid_Protein_2"],
    "Normal Epithelial": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"],
    "PVL": ["PVL_Protein_1", "PVL_Protein_2"],
    "Plasmablasts": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"],
    "T-cells": ["T-cells_Protein_1", "T-cells_Protein_2"],
}

print("="*70)
print("TEST: discover_profiles vs discover_profiles_continuous")
print("="*70)
print(f"Data: {adata.shape[0]} spots, {len(marker_names)} markers")
print(f"Ground truth: {len(GROUND_TRUTH)} cell types, 2 markers each")
print()

# Module 1: Identify interesting markers
print("Module 1: Identifying interesting markers...")
m1_result = identify_interesting_markers(
    X=X, coords=coords, marker_names=marker_names,
    morans_k=8, morans_n_perm=99, seed=42, verbose=False
)
interesting = m1_result.interesting_markers
print(f"  Found {len(interesting)} interesting markers")

# Module 2a: Colocalization
print("Module 2a: Analyzing colocalization...")
coloc_result = analyze_marker_colocalization(
    X=X, coords=coords, marker_names=marker_names,
    markers_to_analyze=interesting,
    neighbor_k=8, n_permutations=199, seed=42, verbose=False
)
print(f"  Analyzed {len(coloc_result.pairs)} pairs")

# Test 1: Original discover_profiles (FDR-based)
print()
print("-"*70)
print("TEST 1: discover_profiles (FDR-based)")
print("-"*70)
result_fdr = discover_profiles(
    coloc_result, fdr_alpha=0.05, top_k=3, seed=42, verbose=False
)

# Count correctly grouped pairs
def count_correct_pairs(profiles, gt):
    correct = 0
    for ct, gt_markers in gt.items():
        gt_set = set(gt_markers)
        for prof in profiles:
            prof_set = set(prof)
            if gt_set.issubset(prof_set) or prof_set.issubset(gt_set):
                if len(gt_set & prof_set) == 2:
                    correct += 1
                    break
    return correct

fdr_correct = count_correct_pairs(result_fdr.profiles, GROUND_TRUTH)
print(f"Profiles discovered: {len(result_fdr.profiles)}")
print(f"Singletons: {len(result_fdr.singletons)}")
print(f"Correctly grouped cell types: {fdr_correct}/9")

# Show profiles
print("Profiles (cell-type markers only):")
for i, prof in enumerate(result_fdr.profiles):
    ct_markers = [m for m in prof if not m.startswith('Nonspecific')]
    if ct_markers:
        print(f"  {i}: {ct_markers}")

# Test 2: discover_profiles_continuous
print()
print("-"*70)
print("TEST 2: discover_profiles_continuous (no FDR gate)")
print("-"*70)
result_cont = discover_profiles_continuous(
    coloc_result, top_k=3, distance_metric="colocalization_score", seed=42, verbose=False
)

cont_correct = count_correct_pairs(result_cont.profiles, GROUND_TRUTH)
print(f"Profiles discovered: {len(result_cont.profiles)}")
print(f"Singletons: {len(result_cont.singletons)}")
print(f"Correctly grouped cell types: {cont_correct}/9")

# Show profiles
print("Profiles (cell-type markers only):")
for i, prof in enumerate(result_cont.profiles):
    ct_markers = [m for m in prof if not m.startswith('Nonspecific')]
    if ct_markers:
        print(f"  {i}: {ct_markers}")

print()
print("="*70)
print("SUMMARY")
print("="*70)
print(f"discover_profiles (FDR):        {fdr_correct}/9 cell types correctly grouped")
print(f"discover_profiles_continuous:   {cont_correct}/9 cell types correctly grouped")
