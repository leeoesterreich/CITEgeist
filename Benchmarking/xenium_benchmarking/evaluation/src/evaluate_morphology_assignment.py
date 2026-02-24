# evaluate_morphology_assignment.py
"""
Evaluate morphology-guided cell type assignment accuracy.

Compares Module 3b's morphology-based nucleus assignment against:
1. Ground truth (Protein-gated and RNA clustering)
2. Baseline methods (random, uniform, spot-proportion-only)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Protein-gated ground truth cell types (7 types + Unknown)
PROTEIN_GT_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

# RNA clustering ground truth cell types (6 types, merged T cells)
RNA_GT_CELL_TYPES = [
    "B cells",
    "T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def load_ground_truth(gt_type: str, gt_dir: Path) -> pd.DataFrame:
    """
    Load ground truth cell type assignments.

    Args:
        gt_type: "protein" or "rna"
        gt_dir: Path to xenium_pseudovisium directory

    Returns:
        DataFrame with cell_id index and cell_type column
    """
    gt_dir = Path(gt_dir)

    if gt_type == "protein":
        gt_path = gt_dir / "data_protein_gt" / "cell_type_assignments.csv"
    elif gt_type == "rna":
        gt_path = gt_dir / "data_rna_gt" / "cell_type_assignments.csv"
    else:
        raise ValueError(f"Unknown gt_type: {gt_type}. Use 'protein' or 'rna'.")

    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth file not found: {gt_path}")

    gt_df = pd.read_csv(gt_path)

    # Handle different index column formats
    # Protein GT has unnamed first column (index), RNA GT has explicit cell_id
    if "cell_id" in gt_df.columns:
        gt_df = gt_df.set_index("cell_id")
    elif gt_df.columns[0] == "Unnamed: 0" or gt_df.columns[0] == "":
        # First column is unnamed index
        gt_df = gt_df.set_index(gt_df.columns[0])
        gt_df.index.name = "cell_id"

    logger.info(f"Loaded {gt_type} GT: {len(gt_df)} cells, {gt_df['cell_type'].nunique()} types")

    return gt_df
