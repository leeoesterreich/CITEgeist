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

# Cell type order (matches create_protein_gt.py)
CELL_TYPE_ORDER = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
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

    # 2. CD4+ T cells (CD3E+ CD4+ CD8A-)
    t_cell_base = expr_df['CD3E'] > CD3E_thresh
    cd4_pos = expr_df['CD4'] > CD4_thresh
    cd8_neg = expr_df['CD8A'] < CD8A_thresh
    cd4_tcells = t_cell_base & cd4_pos & cd8_neg & (cell_types == 'Unknown')
    cell_types[cd4_tcells] = 'CD4+ T cells'

    # 3. CD8+ T cells (CD3E+ CD8A+)
    cd8_pos = expr_df['CD8A'] > CD8A_thresh
    cd8_tcells = t_cell_base & cd8_pos & (cell_types == 'Unknown')
    cell_types[cd8_tcells] = 'CD8+ T cells'

    # 4. Macrophages (CD68+ CD3E-)
    cd68_pos = expr_df['CD68'] > CD68_thresh
    cd3e_neg = expr_df['CD3E'] < CD3E_thresh
    macrophages = cd68_pos & cd3e_neg & (cell_types == 'Unknown')
    cell_types[macrophages] = 'Macrophages'

    # 5. Endothelial (CD31+ CD68- CD3E-)
    cd31_pos = expr_df['CD31'] > CD31_thresh
    cd68_neg = expr_df['CD68'] < CD68_thresh
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == 'Unknown')
    cell_types[endothelial] = 'Endothelial'

    # 6. Epithelial (PanCK+ or E-Cadherin high)
    panck_pos = expr_df['PanCK'] > PanCK_thresh
    ecad_high = expr_df['E-Cadherin'] > ECad_thresh
    epithelial = (panck_pos | ecad_high) & (cell_types == 'Unknown')
    cell_types[epithelial] = 'Epithelial'

    # 7. Fibroblasts (alphaSMA high)
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

    # TODO: Implement remaining steps
    logger.info("Analysis complete!")


if __name__ == '__main__':
    main()
