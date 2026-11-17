"""Single source of truth for the Xenium GT rebuild: labels, recipe values and input/output paths."""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import os
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]  # Benchmarking/xenium_benchmarking

CT6 = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]
CT7 = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]

CD8_GENES = ["CD8A", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "CCL5", "KLRC1", "LAG3", "PDCD1", "TIGIT"]
CD4_GENES = ["CD4", "IL7R", "CCR7", "SELL", "FOXP3", "IL2RA", "CTLA4"]
AMBIGUOUS_ABS_SCORE = 0.7

# Approved fixed marker rules (all genes are on the panel)
MARKER_SETS = {
    "T cells": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2"],
    "B cells": ["MS4A1", "CD79A", "CD19", "BANK1", "MZB1", "TNFRSF17"],
    "Macrophages": ["CD68", "CD163", "CD14", "AIF1", "MRC1", "MS4A6A", "FCGR1A"],
    "Fibroblasts": ["PDGFRA", "DPT", "FBLN1", "MFAP5", "OGN", "SFRP4", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "CLEC14A", "EGFL7", "ERG", "MMRN2", "ADGRL4"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH16", "MUC1", "FXYD2", "GATM", "DPEP1"],
}
MARKER_AMBIGUOUS_MARGIN = 0.25  # top minus second type score below this -> flagged ambiguous (label still kept)

RECIPE = {
    "feature_type": "Gene Expression",
    "n_genes_expected": 405,
    "remove_control_probes": True,
    "drop_zero_count_cells": True,
    "normalize_total_target_sum": 1e4,
    "log1p": True,
    "level1": {
        "n_pcs": 30,
        "n_neighbors": 15,
        "leiden_resolution": 2.0,
        "leiden_flavor": "igraph",
        "leiden_n_iterations": 2,
        "leiden_directed": False,
    },
    "tcell": {"n_pcs": 20, "n_neighbors": 15, "leiden_resolution": 1.0},
    "cd8_genes": CD8_GENES,
    "cd4_genes": CD4_GENES,
    "score": "mean z-scored CD8 genes minus mean z-scored CD4 genes",
    "cd8_if_score_gt": 0.0,
    "ambiguous_abs_score": AMBIGUOUS_ABS_SCORE,
    "seed": 0,
    "seed_stability": 1,
    "pca_solver": "arpack",
    "neighbors_metric": "euclidean",
    "singler_version": None,  # legacy slot (cluster naming no longer uses SingleR)
    "marker_sets": MARKER_SETS,
    "marker_ambiguous_margin": MARKER_AMBIGUOUS_MARGIN,
    "exclude_if_top_score_le": 0.0,  # cluster top type score <= this -> cluster excluded from the GT (cells dropped)
}

# Stability reruns: defined once; the stage runner and the shell loop both read this.
RECIPE["stability_variants"] = {
    "leiden_1.5": {"resolution": 1.5, "seed": RECIPE["seed"]},
    "leiden_2.5": {"resolution": 2.5, "seed": RECIPE["seed"]},
    "seed_1": {"resolution": RECIPE["level1"]["leiden_resolution"], "seed": RECIPE["seed_stability"]},
}

OUT_ROOT = Path(os.environ.get("GT_REBUILD_OUT_ROOT") or BASE / "ground_truth_xenium_clusters")  # env: tests only
XENIUM_DIR = str(_cg.root("DATA") / "10x_PublicData" / "Xenium_RNA_Proteomic_RenalCellCarcinoma")
MAPPING_CSV = BASE.parent / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv"
REF_LABELS_6 = BASE / "reference_data" / "GSE156632" / "processed" / "rctd" / "cell_types.csv"
REF_MTX_DIR = BASE / "reference_data" / "GSE156632" / "processed" / "seurat" / "mtx"
OLD_LABELS_DIR = BASE / "ground_truth_singler" / "singler_labels"
