#!/usr/bin/env python
"""
Generate granular 10-cell-type pseudo-Visium dataset.

This script creates a benchmark dataset using the UNSIMPLIFIED RNA clustering
(10 cell types instead of 6), providing maximum granularity to highlight
CITEgeist's proteomic integration advantage.

Cell types (10 total):
1. CD8+ T cells     - CD3E=378, CD8A=210
2. Macrophages      - CD68=430, CD163=88
3. Mixed Immune     - CD3E=142, CD8A=118, HLA-DR=142
4. Epithelial       - PanCK=39, Vimentin=311
5. Myofibroblasts   - alphaSMA=108, Vimentin=374
6. Stromal          - Mixed low markers
7. Endothelial      - CD31=168
8. B cells          - CD20=293, CD45RA=398
9. Proliferating T  - CD3E=679, PCNA=83
10. Vascular Stromal - CD31=53, Vimentin=209

Reference:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import sys
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))
from split_regions import create_all_region_datasets

DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
REPO_ROOT = Path(__file__).parent.parent.parent.parent
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_granular_gt"


def main():
    logger.info("=" * 70)
    logger.info("GRANULAR 10-CELL-TYPE GROUND TRUTH GENERATION")
    logger.info("=" * 70)
    logger.info("")
    logger.info("This script creates the UNSIMPLIFIED (10 cell type) ground truth.")
    logger.info("Cell types: CD8+ T cells, Macrophages, Mixed Immune, Epithelial,")
    logger.info("            Myofibroblasts, Stromal, Endothelial, B cells,")
    logger.info("            Proliferating T, Vascular Stromal")
    logger.info("")
    logger.info(f"Data directory: {DATA_DIR}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info("")

    # Create granular dataset with simplify_cell_types=False
    summary = create_all_region_datasets(
        data_dir=DATA_DIR,
        output_dir=str(OUTPUT_DIR),
        n_regions=5,
        spot_diameter=55.0,
        center_spacing=100.0,
        min_cells=3,
        use_rna_clusters=True,  # Use RNA-based clustering (recommended)
        n_clusters=10,
        simplify_cell_types=False,  # KEY: Keep 10 cell types unsimplified
    )

    logger.info("")
    logger.info("=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Ground truth method: {summary['ground_truth_method']}")
    logger.info(f"Total spots: {summary['total_spots']}")
    logger.info(f"Total cells: {summary['total_cells']}")
    logger.info(f"Genes: {summary['n_genes']}")
    logger.info(f"Proteins: {summary['n_proteins']}")
    logger.info(f"Cell types: {summary['cell_types']}")
    logger.info(f"Number of regions: {summary['n_regions']}")

    logger.info("")
    logger.info("Per-region breakdown:")
    for region in summary["regions"]:
        logger.info(f"  Region {region['region_id']}: {region['n_spots']} spots")

    logger.info("")
    logger.info("=" * 70)
    logger.info("NEXT STEPS")
    logger.info("=" * 70)
    logger.info("1. Run CITEgeist benchmark with --use-granular-profiles flag")
    logger.info("2. Compare results with simplified 6-cell-type benchmark")
    logger.info("3. Evaluate metrics with evaluation scripts")
    logger.info("")


if __name__ == "__main__":
    main()
