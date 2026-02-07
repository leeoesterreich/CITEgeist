#!/usr/bin/env python
"""
Analyze protein marker profiles for each RNA cluster to create accurate cell type annotations.
"""

import sys
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / "src"))
from load_xenium import load_xenium_data, split_gex_protein
from rna_cell_types import load_rna_clusters

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
OUTPUT_FILE = Path(__file__).parent / "cluster_profiles.json"

def main():
    logger.info("Loading Xenium data...")
    adata = load_xenium_data(DATA_DIR)
    adata_gex, adata_protein = split_gex_protein(adata)

    logger.info("Loading RNA clusters...")
    clusters_df = load_rna_clusters(DATA_DIR, n_clusters=10)

    proteins = list(adata_protein.var_names)

    # Efficient alignment using pandas index
    logger.info("Aligning clusters with protein data...")
    common_cells = clusters_df.index.intersection(adata_protein.obs_names)
    logger.info(f"Found {len(common_cells)} common cells")

    # Subset adata to common cells only (much faster than np.isin)
    adata_subset = adata_protein[list(common_cells)]
    clusters_aligned = clusters_df.loc[common_cells]

    # Convert to dense matrix (only for the subset)
    logger.info("Converting to dense matrix...")
    X_aligned = adata_subset.X.toarray() if hasattr(adata_subset.X, 'toarray') else adata_subset.X

    logger.info(f"Analyzing {len(common_cells)} cells across 10 clusters")

    # Key markers for cell typing
    key_markers = [
        'CD3E', 'CD4', 'CD8A', 'CD20', 'CD138', 'CD68', 'CD163', 'HLA-DR',
        'CD11c', 'CD16', 'GranzymeB', 'PanCK', 'CD31', 'alphaSMA', 'Vimentin',
        'CD45', 'CD45RA', 'CD45RO', 'E-Cadherin', 'PD-1', 'PD-L1', 'Ki-67', 'PCNA',
        'Beta-catenin', 'PTEN', 'LAG-3', 'VISTA'
    ]

    results = {}

    # Get protein names from subset
    subset_proteins = list(adata_subset.var_names)

    for cluster in range(1, 11):
        cluster_mask = clusters_aligned['Cluster'] == cluster
        n_cells = cluster_mask.sum()

        cluster_data = {'n_cells': int(n_cells), 'markers': {}}

        for marker in key_markers:
            if marker in subset_proteins:
                idx = subset_proteins.index(marker)
                values = X_aligned[cluster_mask.values, idx]
                cluster_data['markers'][marker] = {
                    'mean': float(np.mean(values)),
                    'median': float(np.median(values)),
                    'pct_positive': float(100 * np.sum(values > 0) / len(values)) if len(values) > 0 else 0,
                    'max': float(np.max(values)) if len(values) > 0 else 0,
                }

        results[f'cluster_{cluster}'] = cluster_data

    # Save results
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(results, f, indent=2)

    logger.info(f"Saved cluster profiles to {OUTPUT_FILE}")

    # Print summary table
    print("\n" + "=" * 120)
    print("Mean protein expression per RNA cluster:")
    print("=" * 120)

    # Build summary DataFrame
    summary_data = {}
    for cluster in range(1, 11):
        cluster_key = f'cluster_{cluster}'
        summary_data[cluster] = {}
        summary_data[cluster]['n_cells'] = results[cluster_key]['n_cells']
        for marker in key_markers:
            if marker in results[cluster_key]['markers']:
                summary_data[cluster][marker] = results[cluster_key]['markers'][marker]['mean']

    df = pd.DataFrame(summary_data).T
    df.index.name = 'Cluster'
    print(df.round(1).to_string())

    print("\n\nCell counts per cluster:")
    for c in range(1, 11):
        print(f"  Cluster {c}: {results[f'cluster_{c}']['n_cells']:,} cells")

    print("\n\nTop markers per cluster (mean > 50):")
    for cluster in range(1, 11):
        markers = results[f'cluster_{cluster}']['markers']
        high_markers = [(m, v['mean']) for m, v in markers.items() if v['mean'] > 50]
        high_markers.sort(key=lambda x: -x[1])
        print(f"  Cluster {cluster}: {high_markers[:6]}")

    # Suggest cell type annotations
    print("\n\n" + "=" * 80)
    print("SUGGESTED CELL TYPE ANNOTATIONS:")
    print("=" * 80)

    for cluster in range(1, 11):
        markers = results[f'cluster_{cluster}']['markers']
        n_cells = results[f'cluster_{cluster}']['n_cells']

        # Determine cell type based on marker profile
        cd3e = markers.get('CD3E', {}).get('mean', 0)
        cd4 = markers.get('CD4', {}).get('mean', 0)
        cd8a = markers.get('CD8A', {}).get('mean', 0)
        cd20 = markers.get('CD20', {}).get('mean', 0)
        cd138 = markers.get('CD138', {}).get('mean', 0)
        cd68 = markers.get('CD68', {}).get('mean', 0)
        cd163 = markers.get('CD163', {}).get('mean', 0)
        hladr = markers.get('HLA-DR', {}).get('mean', 0)
        cd11c = markers.get('CD11c', {}).get('mean', 0)
        cd16 = markers.get('CD16', {}).get('mean', 0)
        granzymeb = markers.get('GranzymeB', {}).get('mean', 0)
        panck = markers.get('PanCK', {}).get('mean', 0)
        cd31 = markers.get('CD31', {}).get('mean', 0)
        alphasma = markers.get('alphaSMA', {}).get('mean', 0)
        vimentin = markers.get('Vimentin', {}).get('mean', 0)

        # Decision logic
        if cd3e > 200 and cd8a > 100:
            cell_type = "CD8+ T cells"
            rationale = f"CD3E={cd3e:.0f}, CD8A={cd8a:.0f}"
        elif cd3e > 100 and cd4 > 50 and cd8a < 50:
            cell_type = "CD4+ T cells"
            rationale = f"CD3E={cd3e:.0f}, CD4={cd4:.0f}, CD8A={cd8a:.0f}"
        elif cd3e > 100:
            cell_type = "T cells (mixed)"
            rationale = f"CD3E={cd3e:.0f}, CD4={cd4:.0f}, CD8A={cd8a:.0f}"
        elif cd20 > 200:
            cell_type = "B cells"
            rationale = f"CD20={cd20:.0f}"
        elif cd138 > 100 and cd20 < 50:
            cell_type = "Plasma cells"
            rationale = f"CD138={cd138:.0f}, CD20={cd20:.0f}"
        elif cd68 > 300 and cd163 > 100:
            cell_type = "M2-like Macrophages"
            rationale = f"CD68={cd68:.0f}, CD163={cd163:.0f}"
        elif cd68 > 200:
            cell_type = "Macrophages"
            rationale = f"CD68={cd68:.0f}, HLA-DR={hladr:.0f}"
        elif cd11c > 100 and hladr > 100:
            cell_type = "Dendritic cells"
            rationale = f"CD11c={cd11c:.0f}, HLA-DR={hladr:.0f}"
        elif cd16 > 100 and granzymeb > 50 and cd3e < 50:
            cell_type = "NK cells"
            rationale = f"CD16={cd16:.0f}, GranzymeB={granzymeb:.0f}, CD3E={cd3e:.0f}"
        elif panck > 30:
            cell_type = "Epithelial"
            rationale = f"PanCK={panck:.0f}"
        elif cd31 > 100:
            cell_type = "Endothelial"
            rationale = f"CD31={cd31:.0f}"
        elif alphasma > 50:
            cell_type = "Fibroblasts/CAFs"
            rationale = f"alphaSMA={alphasma:.0f}, Vimentin={vimentin:.0f}"
        elif vimentin > 200 and alphasma < 30:
            cell_type = "Mesenchymal/Stromal"
            rationale = f"Vimentin={vimentin:.0f}, alphaSMA={alphasma:.0f}"
        else:
            cell_type = "Unknown"
            rationale = "No clear marker profile"

        print(f"  Cluster {cluster:2d} ({n_cells:,} cells): {cell_type:25s} | {rationale}")

if __name__ == "__main__":
    main()
