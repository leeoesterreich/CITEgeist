#!/usr/bin/env python
"""
Analyze concordance between protein-gated and RNA-defined cell types.

Compares cell type assignments from:
1. Protein hierarchical gating (create_protein_gt.py logic)
2. RNA k-means clustering (analysis.tar.gz)

Output: Confusion matrix, agreement rates, spot-level correlations.
"""

import argparse
import logging
import tarfile
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.metrics import confusion_matrix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
XENIUM_DIR = Path('/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma')
PSEUDOVISIUM_DIR = Path('/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium')
OUTPUT_DIR = PSEUDOVISIUM_DIR / 'analysis'

# Cell type order - using combined "T cells" since RNA can't distinguish CD4/CD8
CELL_TYPE_ORDER = [
    "B cells",
    "T cells",  # Combined CD4+ and CD8+ (RNA can't separate them reliably)
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def get_threshold(expr_df: pd.DataFrame, marker: str, percentile: float = 50) -> float:
    """Get threshold as percentile of non-zero values."""
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0


def classify_cells_by_protein(expr_df: pd.DataFrame) -> pd.Series:
    """
    Classify cells using hierarchical protein gating.
    Matches logic in create_protein_gt.py.
    """
    cell_ids = expr_df.index
    cell_types = pd.Series('Unknown', index=cell_ids)

    # Calculate thresholds
    CD3E_thresh = get_threshold(expr_df, 'CD3E', 50)
    CD4_thresh = get_threshold(expr_df, 'CD4', 50)
    CD8A_thresh = get_threshold(expr_df, 'CD8A', 50)
    CD20_thresh = get_threshold(expr_df, 'CD20', 25)
    CD68_thresh = get_threshold(expr_df, 'CD68', 50)
    CD31_thresh = get_threshold(expr_df, 'CD31', 50)
    PanCK_thresh = get_threshold(expr_df, 'PanCK', 25)
    ECad_thresh = get_threshold(expr_df, 'E-Cadherin', 90)
    alphaSMA_thresh = get_threshold(expr_df, 'alphaSMA', 75)

    # Hierarchical gating
    # 1. B cells (CD20+)
    b_cells = expr_df['CD20'] > CD20_thresh
    cell_types[b_cells] = 'B cells'

    # 2. T cells (CD3E+) - combined, since RNA can't distinguish CD4/CD8
    t_cell_base = expr_df['CD3E'] > CD3E_thresh
    t_cells = t_cell_base & (cell_types == 'Unknown')
    cell_types[t_cells] = 'T cells'

    # 3. Macrophages (CD68+ CD3E-)
    cd68_pos = expr_df['CD68'] > CD68_thresh
    cd3e_neg = expr_df['CD3E'] < CD3E_thresh
    macrophages = cd68_pos & cd3e_neg & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'

    # 4. Endothelial (CD31+ CD68- CD3E-)
    cd31_pos = expr_df['CD31'] > CD31_thresh
    cd68_neg = expr_df['CD68'] < CD68_thresh
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'

    # 5. Epithelial (PanCK+ or E-Cadherin high)
    panck_pos = expr_df['PanCK'] > PanCK_thresh
    ecad_high = expr_df['E-Cadherin'] > ECad_thresh
    epithelial = (panck_pos | ecad_high) & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'

    # 6. Fibroblasts (alphaSMA high)
    asma_high = expr_df['alphaSMA'] > alphaSMA_thresh
    myofib = asma_high & ~cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[myofib] = 'Fibroblasts'

    return cell_types


def load_protein_data() -> Tuple[pd.DataFrame, pd.Series]:
    """Load Xenium protein data and apply gating."""
    logger.info("Loading Xenium protein data...")
    adata = sc.read_10x_h5(XENIUM_DIR / 'cell_feature_matrix.h5', gex_only=False)

    protein_mask = adata.var['feature_types'] == 'Protein Expression'
    adata_protein = adata[:, protein_mask].copy()

    X = adata_protein.X.toarray() if hasattr(adata_protein.X, 'toarray') else adata_protein.X
    proteins = list(adata_protein.var_names)
    cell_ids = [str(x) for x in adata_protein.obs_names]

    expr_df = pd.DataFrame(X, index=cell_ids, columns=proteins)
    logger.info(f"Loaded {len(expr_df)} cells x {len(proteins)} proteins")

    logger.info("Applying protein gating...")
    protein_cell_types = classify_cells_by_protein(expr_df)

    type_counts = protein_cell_types.value_counts()
    logger.info(f"Protein-gated cell types:\n{type_counts}")

    return expr_df, protein_cell_types


# RNA cluster to cell type mapping (from rna_cell_types.py)
RNA_CLUSTER_MAP = {
    1: "CD8+ T cells",
    2: "Macrophages",
    3: "Mixed Immune",
    4: "Epithelial",
    5: "Fibroblasts",
    6: "Stromal",
    7: "Endothelial",
    8: "B cells",
    9: "Proliferating T",
    10: "Vascular Stromal",
}

# Simplified mapping - RNA can't distinguish CD4/CD8 T cells reliably
# (CD4 gene is lowly expressed, only ~34% detection in protein-gated CD4+ T cells)
RNA_SIMPLIFIED_MAP = {
    "CD8+ T cells": "T cells",      # Combine all T cells
    "Macrophages": "Macrophages",
    "Mixed Immune": "Macrophages",  # Closest match
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
    "Stromal": "Fibroblasts",
    "Endothelial": "Endothelial",
    "B cells": "B cells",
    "Proliferating T": "T cells",   # Combine all T cells
    "Vascular Stromal": "Endothelial",
}


def load_rna_clusters() -> pd.Series:
    """Load RNA k-means clusters from Xenium analysis output."""
    logger.info("Loading RNA k-means clusters...")

    analysis_tar = XENIUM_DIR / 'analysis.tar.gz'
    cluster_path = 'analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv'

    with tarfile.open(analysis_tar, 'r:gz') as tar:
        f = tar.extractfile(cluster_path)
        clusters_df = pd.read_csv(f)

    clusters_df['cell_id'] = clusters_df['Barcode'].astype(str)
    clusters_df = clusters_df.set_index('cell_id')

    # Map cluster number to cell type
    rna_cell_types = clusters_df['Cluster'].map(RNA_CLUSTER_MAP)

    # Simplify to match protein-gated categories
    rna_cell_types_simplified = rna_cell_types.map(RNA_SIMPLIFIED_MAP)

    type_counts = rna_cell_types_simplified.value_counts()
    logger.info(f"RNA-defined cell types:\n{type_counts}")

    return rna_cell_types_simplified


def build_concordance_matrix(
    protein_types: pd.Series,
    rna_types: pd.Series,
    output_dir: Path,
) -> Tuple[pd.DataFrame, Dict[str, float], float]:
    """Build and visualize confusion matrix between protein and RNA cell types."""
    logger.info("Building concordance matrix...")

    # Align indices
    common_cells = protein_types.index.intersection(rna_types.index)
    logger.info(f"Common cells: {len(common_cells)}")

    protein_aligned = protein_types.loc[common_cells]
    rna_aligned = rna_types.loc[common_cells]

    # Get all cell types (excluding Unknown)
    all_types = CELL_TYPE_ORDER + ['Unknown']

    # Build confusion matrix
    cm = confusion_matrix(
        protein_aligned,
        rna_aligned,
        labels=all_types,
    )

    cm_df = pd.DataFrame(cm, index=all_types, columns=all_types)
    cm_df.index.name = 'Protein-gated'
    cm_df.columns.name = 'RNA-defined'

    # Calculate agreement metrics
    total = len(common_cells)
    diagonal = np.diag(cm).sum()
    overall_agreement = diagonal / total * 100

    logger.info(f"Overall concordance: {overall_agreement:.1f}%")

    # Per-cell-type concordance
    per_type_concordance = {}
    for i, ct in enumerate(all_types):
        if cm[i, :].sum() > 0:
            per_type_concordance[ct] = cm[i, i] / cm[i, :].sum() * 100

    logger.info("Per-cell-type concordance (protein → same RNA type):")
    for ct, pct in per_type_concordance.items():
        logger.info(f"  {ct}: {pct:.1f}%")

    # Save raw matrix
    cm_df.to_csv(output_dir / 'concordance_matrix.csv')

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(10, 8))

    # Normalize by row for visualization
    cm_norm = cm.astype(float)
    row_sums = cm_norm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    cm_norm = cm_norm / row_sums * 100

    sns.heatmap(
        cm_norm,
        annot=True,
        fmt='.0f',
        cmap='Blues',
        xticklabels=all_types,
        yticklabels=all_types,
        ax=ax,
        cbar_kws={'label': '% of protein-gated type'}
    )

    ax.set_xlabel('RNA-defined cell type')
    ax.set_ylabel('Protein-gated cell type')
    ax.set_title(f'Cell Type Concordance (Overall: {overall_agreement:.1f}%)')

    plt.tight_layout()
    plt.savefig(output_dir / 'concordance_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved heatmap to {output_dir / 'concordance_heatmap.png'}")

    return cm_df, per_type_concordance, overall_agreement


def calculate_spot_correlations(
    protein_types: pd.Series,
    rna_types: pd.Series,
    output_dir: Path,
) -> pd.DataFrame:
    """Calculate spot-level proportion correlations between protein and RNA GT."""
    from scipy import stats as scipy_stats

    logger.info("Calculating spot-level correlations...")

    # Load cell-to-spot mapping
    cell_to_spot_path = PSEUDOVISIUM_DIR / 'data_protein_gt' / 'cell_to_spot_mapping.csv'
    cell_to_spot = pd.read_csv(cell_to_spot_path, index_col=0)
    cell_to_spot.index = cell_to_spot.index.astype(str)

    # Align with cell types
    common_cells = protein_types.index.intersection(rna_types.index).intersection(cell_to_spot.index)
    logger.info(f"Cells with spot mapping: {len(common_cells)}")

    protein_aligned = protein_types.loc[common_cells]
    rna_aligned = rna_types.loc[common_cells]
    spots = cell_to_spot.loc[common_cells, 'spot_id']

    # Calculate proportions per spot for each GT
    spot_ids = spots.unique()

    protein_props = pd.DataFrame(0.0, index=spot_ids, columns=CELL_TYPE_ORDER)
    rna_props = pd.DataFrame(0.0, index=spot_ids, columns=CELL_TYPE_ORDER)

    for spot_id in spot_ids:
        spot_mask = spots == spot_id
        spot_protein = protein_aligned[spot_mask]
        spot_rna = rna_aligned[spot_mask]

        n_cells = spot_mask.sum()
        if n_cells == 0:
            continue

        # Protein proportions
        for ct in CELL_TYPE_ORDER:
            protein_props.loc[spot_id, ct] = (spot_protein == ct).sum() / n_cells
            rna_props.loc[spot_id, ct] = (spot_rna == ct).sum() / n_cells

    # Calculate correlation per cell type
    correlations = {}
    for ct in CELL_TYPE_ORDER:
        r, p = scipy_stats.pearsonr(protein_props[ct], rna_props[ct])
        correlations[ct] = {'pearson_r': r, 'p_value': p}

    corr_df = pd.DataFrame(correlations).T
    logger.info(f"Spot-level correlations:\n{corr_df}")

    # Save
    corr_df.to_csv(output_dir / 'spot_level_correlations.csv')

    # Plot scatter for each cell type
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, ct in enumerate(CELL_TYPE_ORDER):
        ax = axes[i]
        ax.scatter(protein_props[ct], rna_props[ct], alpha=0.3, s=10)
        ax.plot([0, 1], [0, 1], 'r--', label='y=x')
        ax.set_xlabel('Protein-GT proportion')
        ax.set_ylabel('RNA-GT proportion')
        ax.set_title(f'{ct}\nr={correlations[ct]["pearson_r"]:.2f}')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    # Hide extra subplot
    axes[7].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / 'spot_correlations.png', dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved scatter plots to {output_dir / 'spot_correlations.png'}")

    return corr_df


def generate_report(
    overall_agreement: float,
    per_type_concordance: Dict[str, float],
    spot_corr: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Generate markdown summary report."""
    logger.info("Generating summary report...")

    report = f"""# Protein vs RNA Cell Type Concordance Analysis

**Date**: {pd.Timestamp.now().strftime('%Y-%m-%d')}

## Summary

Overall cell-level concordance: **{overall_agreement:.1f}%**

## Per-Cell-Type Concordance

| Cell Type | Concordance (%) |
|-----------|-----------------|
"""

    for ct, pct in per_type_concordance.items():
        report += f"| {ct} | {pct:.1f} |\n"

    report += """
## Spot-Level Correlations

| Cell Type | Pearson r | p-value |
|-----------|-----------|---------|
"""

    for ct in spot_corr.index:
        r = spot_corr.loc[ct, 'pearson_r']
        p = spot_corr.loc[ct, 'p_value']
        report += f"| {ct} | {r:.3f} | {p:.2e} |\n"

    # Interpretation
    if overall_agreement >= 80:
        interpretation = "Current protein-GT benchmark is **defensible** (>80% concordance)."
    elif overall_agreement >= 60:
        interpretation = "Consider reporting **both benchmarks** in supplement (60-80% concordance)."
    else:
        interpretation = "**Dual benchmarks required** for fair comparison (<60% concordance)."

    report += f"""
## Interpretation

{interpretation}

## Figures

- `concordance_heatmap.png` - Cell-level confusion matrix
- `spot_correlations.png` - Spot-level proportion scatter plots
- `concordance_matrix.csv` - Raw confusion matrix
- `spot_level_correlations.csv` - Correlation values
"""

    report_path = output_dir / 'concordance_report.md'
    with open(report_path, 'w') as f:
        f.write(report)

    logger.info(f"Saved report to {report_path}")


def main():
    """Run concordance analysis."""
    parser = argparse.ArgumentParser(description='Analyze protein vs RNA cell type concordance')
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Protein vs RNA Cell Type Concordance Analysis")
    logger.info("=" * 60)

    # Load protein data and apply gating
    expr_df, protein_cell_types = load_protein_data()

    # Load RNA clusters
    rna_cell_types = load_rna_clusters()

    # Build concordance matrix
    cm_df, per_type, overall = build_concordance_matrix(
        protein_cell_types, rna_cell_types, output_dir
    )

    # Calculate spot-level correlations
    spot_corr = calculate_spot_correlations(
        protein_cell_types, rna_cell_types, output_dir
    )

    # Generate report
    generate_report(overall, per_type, spot_corr, output_dir)

    logger.info("=" * 60)
    logger.info("Analysis complete!")
    logger.info(f"Results in: {output_dir}")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
