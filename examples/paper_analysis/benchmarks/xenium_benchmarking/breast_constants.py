"""
Breast cancer tissue configuration for Xenium benchmarking.

Adds breast cancer support alongside existing RCC benchmark.
Cell type taxonomy: 8 types collapsed from 16 Xenium vendor annotations.
Reference: Wu et al. 2021 (Nature Genetics) breast cancer atlas.
"""

from pathlib import Path
from typing import Dict, List, Optional

# Base directory for xenium benchmarking
_BENCH_DIR = Path(__file__).parent

# =============================================================================
# BREAST CANCER CELL TYPE DEFINITIONS
# =============================================================================

BREAST_GT_MAPPING: Dict[str, str] = {
    "0_Proliferative_Tumor_Cells": "Cancer Epithelial",
    "4_Tumor_Cells": "Cancer Epithelial",
    "7_Tumor_Cells": "Cancer Epithelial",
    "8_Tumor_Cells": "Cancer Epithelial",
    "9_Tumor_Cells": "Cancer Epithelial",
    "12_Tumor_Cells": "Cancer Epithelial",
    "15_Tumor_Cells": "Cancer Epithelial",
    "1_Macrophages": "Macrophages",
    "6_T_Lymphocytes": "T cells",
    "13_Plasma_Cells": "Plasma cells",
    "2_Stromal_Cells": "Fibroblasts",
    "3_Perivascular_Cells": "Perivascular",
    "5_Endothelial_Cells": "Endothelial",
    "11_Normal_Epithelial_Cells": "Normal Epithelial",
    "14_Normal_Myoepithelial_Cells": "Normal Epithelial",
    "10_Tumor_Adjacent_Myoepithelial_Cells": "Normal Epithelial",
}

BREAST_8_CELL_TYPES: List[str] = [
    "Cancer Epithelial",
    "Macrophages",
    "Fibroblasts",
    "Perivascular",
    "Endothelial",
    "T cells",
    "Normal Epithelial",
    "Plasma cells",
]

WU_TO_BREAST_MAPPING: Dict[str, Optional[str]] = {
    "Cancer Epithelial": "Cancer Epithelial",
    "Myeloid": "Macrophages",
    "CAFs": "Fibroblasts",
    "PVL": "Perivascular",
    "Endothelial": "Endothelial",
    "T-cells": "T cells",
    "Normal Epithelial": "Normal Epithelial",
    "Plasmablasts": "Plasma cells",
    "B-cells": "B-cells",  # kept in reference, not in evaluated_types
}

# =============================================================================
# DATA PATHS
# =============================================================================

XENIUM_BREAST_DIR = Path("/path/to/Xenium_BreastCancer")
XENIUM_BREAST_ZIP = XENIUM_BREAST_DIR / "Human_Breast_Biomarkers_S1_Top_outs.zip"
XENIUM_BREAST_CELL_GROUPS = XENIUM_BREAST_DIR / "Human_Breast_Biomarkers_S1_Top_cell_groups.csv"

WU_2021_DIR = Path("/path/to/Wu_2021_BRCA_scRNASeq")

# =============================================================================
# TISSUE CONFIGS
# =============================================================================

TISSUE_CONFIGS = {
    "breast": {
        "xenium_data": str(XENIUM_BREAST_DIR),
        "scrna_reference": str(_BENCH_DIR / "reference_data" / "Wu_2021"),
        "gt_source": "xenium_clusters",
        "gt_mapping": BREAST_GT_MAPPING,
        "cell_types": BREAST_8_CELL_TYPES,
        "n_regions": 5,
        "reference_label_col": "celltype_major",
        "reference_type_mapping": WU_TO_BREAST_MAPPING,
        "evaluated_types": BREAST_8_CELL_TYPES,
        "data_dir": str(_BENCH_DIR / "data_breast"),
        "prefix": "Xenium",
    },
}
