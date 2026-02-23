#!/usr/bin/env python
"""Analyze antibody data distribution for rare cell types."""
import sys
sys.path.insert(0, '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist')

import numpy as np
import pandas as pd
import scanpy as sc

def main():
    # Load CITE-seq data
    cite_path = '/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/scCube_12k/replicates/mixed/h5ad_objects/Wu_rep_0_CITE.h5ad'
    adata = sc.read_h5ad(cite_path)

    print('=== ANTIBODY DATA ANALYSIS (raw) ===')
    print(f'Shape: {adata.shape}')

    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Load ground truth proportions
    gt_path = '/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/scCube_12k/replicates/mixed/ST_sim/Wu_ST_0_prop.csv'
    gt_df = pd.read_csv(gt_path, index_col=0).drop(columns=['spot_x', 'spot_y'], errors='ignore')

    # For each cell type specific marker, check distribution
    markers_of_interest = {
        'Endothelial': ['Endothelial_Protein_1', 'Endothelial_Protein_2'],
        'Normal Epithelial': ['Normal Epithelial_Protein_1', 'Normal Epithelial_Protein_2'],
        'PVL': ['PVL_Protein_1', 'PVL_Protein_2'],
        'Plasmablasts': ['Plasmablasts_Protein_1', 'Plasmablasts_Protein_2'],
        'CAFs': ['CAFs_Protein_1', 'CAFs_Protein_2'],
        'Cancer Epithelial': ['Cancer Epithelial_Protein_1', 'Cancer Epithelial_Protein_2'],
    }

    for ct, markers in markers_of_interest.items():
        gt_prop = gt_df[ct].mean() * 100 if ct in gt_df.columns else 0
        print(f'\n{ct} (GT: {gt_prop:.1f}%):')
        for m in markers:
            if m in adata.var_names:
                idx = list(adata.var_names).index(m)
                vals = X[:, idx]
                print(f'  {m}:')
                print(f'    min={vals.min():.1f}, max={vals.max():.1f}, mean={vals.mean():.2f}')
                print(f'    p5={np.percentile(vals, 5):.1f}, p50={np.percentile(vals, 50):.1f}, p95={np.percentile(vals, 95):.1f}')
                print(f'    nonzero={100*np.sum(vals > 0)/len(vals):.1f}%')

    # Compare rare vs dominant markers
    print('\n\n=== RARE vs DOMINANT MARKER COMPARISON ===')
    rare_markers = ['Endothelial_Protein_1', 'PVL_Protein_1', 'Normal Epithelial_Protein_1', 'Plasmablasts_Protein_1']
    dom_markers = ['CAFs_Protein_1', 'Cancer Epithelial_Protein_1']

    print('\nRare cell type markers (GT < 2%):')
    for m in rare_markers:
        if m in adata.var_names:
            idx = list(adata.var_names).index(m)
            vals = X[:, idx]
            print(f'  {m}: mean={vals.mean():.2f}, p50={np.percentile(vals, 50):.1f}, p95={np.percentile(vals, 95):.1f}')

    print('\nDominant cell type markers (GT > 25%):')
    for m in dom_markers:
        if m in adata.var_names:
            idx = list(adata.var_names).index(m)
            vals = X[:, idx]
            print(f'  {m}: mean={vals.mean():.2f}, p50={np.percentile(vals, 50):.1f}, p95={np.percentile(vals, 95):.1f}')

    # Now show what happens after preprocessing
    print('\n\n=== AFTER DISCRETE PREPROCESSING (per_marker scaling) ===')
    from CITEgeist.model.citegeist_model import CitegeistModel

    gex_path = '/ix1/alee/LO_LAB/Personal/Brent_Schlegel/bts76/Projects/CITEgeist/Wu_Visium/Simulations/scCube_12k/replicates/mixed/h5ad_objects/Wu_rep_0_GEX.h5ad'
    adata_gex = sc.read_h5ad(gex_path)

    model = CitegeistModel(
        sample_name='debug',
        output_folder='/tmp/citegeist_debug',
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata,
    )
    model.preprocess_antibody_discrete(scale_mode='per_marker')

    X_proc = model.antibody_matrix_discrete
    marker_names_proc = model.discrete_marker_names

    print('Rare markers after preprocessing:')
    for m in rare_markers:
        if m in marker_names_proc:
            idx = marker_names_proc.index(m)
            vals = X_proc[:, idx]
            print(f'  {m}: mean={vals.mean():.3f}, min={vals.min():.3f}, max={vals.max():.3f}')

    print('\nDominant markers after preprocessing:')
    for m in dom_markers:
        if m in marker_names_proc:
            idx = marker_names_proc.index(m)
            vals = X_proc[:, idx]
            print(f'  {m}: mean={vals.mean():.3f}, min={vals.min():.3f}, max={vals.max():.3f}')

if __name__ == '__main__':
    main()
