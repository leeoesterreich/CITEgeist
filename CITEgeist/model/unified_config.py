"""Shared configuration for the unified PC-MIL pipeline."""

import os
from pathlib import Path

# === Cell Type Profile (9 types, K=9) ===
# Nested format for Module 3 (validated by validate_cell_profile_dict)
CELL_PROFILES_NESTED = {
    "Epithelial": {"Major": ["EPCAM-1"]},
    "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    "CD8_T_Cells": {"Major": ["CD3E-1", "CD8A-1"]},
    "CD4_T_Cells": {"Major": ["CD4-1"]},
    "B_Cells": {"Major": ["CD19-1"]},
    "Endothelial": {"Major": ["PECAM1-1"]},
    "Fibroblasts": {"Major": ["ACTA2-1"]},
    "Monocytes": {"Major": ["CD14-1"]},
    "Dendritic_Cells": {"Major": ["ITGAX-1", "HLA-DRA-1"]},
}

CELL_TYPE_NAMES = list(CELL_PROFILES_NESTED.keys())
K = len(CELL_TYPE_NAMES)  # 9

# === RNA Marker Validation Dictionary ===
RNA_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18"],
    "Macrophages": ["CD68", "CD163", "CSF1R"],
    "CD8_T_Cells": ["CD8A", "CD8B", "GZMB"],
    "CD4_T_Cells": ["CD4", "IL7R", "CCR7"],
    "B_Cells": ["CD19", "MS4A1", "CD79A"],
    "Fibroblasts": ["COL1A1", "DCN", "VIM"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Monocytes": ["CD14", "FCGR3A", "S100A8"],
    "Dendritic_Cells": ["ITGAX", "HLA-DRA", "CLEC10A"],
}

# === Patient Sample List ===
# Placeholder sample IDs — replace with your own (e.g. from the GEO accession in the manuscript).
PATIENT_SAMPLES = [
    "sample_P1_S1",
    "sample_P1_S2",
    "sample_P2_S1",
    "sample_P2_S2",
    "sample_P3_S1",
    "sample_P3_S2",
    "sample_P4_S1",
    "sample_P4_S2",
    "sample_P5_S1",
    "sample_P5_S2",
    "sample_P6_S1",
    "sample_P6_S2",
]

# === Paths ===
DATA_DIR = Path(os.environ.get("CITEGEIST_DATA_DIR", "data/processed_files"))
OUTPUT_BASE = Path("output/unified_pipeline")

# === Training Defaults ===
MAX_EPOCHS = 200
PATIENCE = 30
LAMBDA_RECON = 1.0
LAMBDA_ENTROPY = 0.1
LAMBDA_DIVERSITY = 0.5
LAMBDA_HUNGARIAN = 0.0  # Disabled for unified pipeline
RECON_WARMUP_EPOCHS = 20
PROTEIN_DROPOUT = 0.3

# === ViT Config ===
VIT_MODEL = "vit_small_patch16_224"
VIT_DIM = 384
PATCH_SIZE = 224
