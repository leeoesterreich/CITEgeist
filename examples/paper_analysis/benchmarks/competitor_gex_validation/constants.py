"""Shared constants for competitor GEX biology validation."""
from pathlib import Path

# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT_ROOT = Path("/path/to/CITEgeist_analysis")
OUTPUT_ROOT = PROJECT_ROOT / "output" / "competitor_gex"

WU_DATA_DIR = Path("/path/to/Wu_2021_BRCA_scRNASeq")
SPACERANGER_ROOT = Path("/path/to/CITEgeist_public_data/processed_files")
CITEGEIST_OUTPUT = PROJECT_ROOT / "output"

ERPOS_CIDS = [
    "CID3941",
    "CID3948",
    "CID4040",
    "CID4067",
    "CID4290A",
    "CID4398",
    "CID4461",
    "CID4463",
    "CID4471",
    "CID4530N",
    "CID4535",
]

SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

WU_TO_BREAST_MAPPING = {
    "Cancer Epithelial": "Cancer Epithelial",
    "Myeloid": "Macrophages",
    "CAFs": "Fibroblasts",
    "PVL": "Perivascular",
    "Endothelial": "Endothelial",
    "T-cells": "T cells",
    "Normal Epithelial": "Normal Epithelial",
    "Plasmablasts": "Plasma cells",
    "B-cells": "B-cells",
}

CELL_TYPES = sorted(set(WU_TO_BREAST_MAPPING.values()))

CELLS_PER_SPOT = 5
SPOT_RADIUS_UM = 27.5

CANDIDATE_GENES = ["MDK", "PTN", "SDC4", "NCL", "CD274"]

# Secretory-program signature genes (must match mdk_analysis/scripts/spatial_morans.py)
SECRETORY_GENES = [
    "HSP90B1",
    "HSPA5",
    "CALR",
    "CANX",
    "PDIA4",
    "PDIA6",
    "SEC23A",
    "SEC61B",
    "ATF6",
    "MAN1A1",
    "XBP1",
]

# Minimum nonzero cells for MDK to be evaluable (matches CITEgeist's threshold)
MDK_NONZERO_THRESHOLD = 10

COMMOT_PATHWAYS = ["s-MDK-SDC4", "s-MDK-NCL", "s-PTN-SDC4", "s-PTN-NCL", "s-MIF-CD74_CD44"]
COMMOT_DISTANCE_THR = 500
MORANS_N_PERMS = 999
MORANS_K_NEIGH = 6

NMF_K_PROGRAMS = 5
NMF_LAMBDA_SPARSITY = 0.01
