#!/bin/bash
#SBATCH --job-name=sace_gex_12patient
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=output/morphology_assignment_v2/logs/sace_gex_%j.out
#SBATCH --error=output/morphology_assignment_v2/logs/sace_gex_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

set -u

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd "$REPO_ROOT"

# Run SACE GEX for each sample using existing cell_assignments
python -u -c "
import sys, logging
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, '.')

PATIENT_DATA_ROOT = Path('/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files')
QP_DIR = Path('output/module3_cuopt_qp')
V2_DIR = Path('output/morphology_assignment_v2')

SAMPLES = [
    'HCC22-088-P1-S1', 'HCC22-088-P1-S2',
    'HCC22-088-P2-S1', 'HCC22-088-P2-S2',
    'HCC22-088-P3-S1_A', 'HCC22-088-P3-S2',
    'HCC22-088-P4-S1', 'HCC22-088-P4-S2_1i_rep',
    'HCC22-088-P5-S1', 'HCC22-088-P5-S2_F_rep',
    'HCC22-088-P6-S1', 'HCC22-088-P6-S2_D',
]
SAMPLE_METADATA = {s: {'min_counts': 100 if 'S1' in s else 25} for s in SAMPLES}

MODEL_PROFILE_DICT = {
    'Endothelial': {'Major': ['PECAM1-1']},
    'Fibroblasts': {'Major': ['ACTA2-1']},
    'B_Cells': {'Major': ['CD19-1']},
    'Macrophages': {'Major': ['CD68-1', 'CD163-1']},
    'Monocytes': {'Major': ['CD14-1']},
    'CD8_T_Cells': {'Major': ['CD8A-1', 'CD3E-1']},
    'CD4_T_Cells': {'Major': ['CD4-1', 'CD3E-1']},
    'Cancer_Luminal': {'Major': ['EPCAM-1']},
    'Cancer_Basal': {'Major': ['KRT5-1', 'SDC1-1', 'EPCAM-1']},
    'Dendritic_Cells': {'Major': ['ITGAX-1', 'HLA-DRA-1']},
}

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

import squidpy as sq
from CITEgeist.model import CitegeistModel

for sample in SAMPLES:
    sample_dir = V2_DIR / sample
    assignments_csv = sample_dir / 'cell_assignments.csv'
    h5ad_path = sample_dir / f'{sample}_single_cell.h5ad'

    if h5ad_path.exists():
        logger.info('SKIP %s — h5ad already exists', sample)
        continue

    if not assignments_csv.exists():
        logger.warning('SKIP %s — no cell_assignments.csv', sample)
        continue

    logger.info('=== SACE GEX for %s ===', sample)

    # Load cell assignments
    assign_df = pd.read_csv(assignments_csv)
    cell_assignments = dict(zip(assign_df['nucleus_id'].astype(str), assign_df['assigned_type']))

    # Load ensemble proportions
    ensemble_df = pd.read_csv(sample_dir / 'ensemble_proportions.csv', index_col=0)

    # Load nucleus-spot mapping
    nuc_map = pd.read_csv(sample_dir / 'nucleus_spot_mapping.csv')
    cell_spot_map = nuc_map.rename(columns={'nucleus_id': 'cell_id', 'centroid_x': 'x', 'centroid_y': 'y'})
    assigned_ids = set(assign_df['nucleus_id'].values)
    cell_spot_map = cell_spot_map[cell_spot_map['cell_id'].isin(assigned_ids)].copy()
    cell_spot_map['cell_id'] = cell_spot_map['cell_id'].astype(str)

    # Load patient data
    patient_dir = PATIENT_DATA_ROOT / sample / 'outs'
    adata = sq.read.visium(str(patient_dir), counts_file='filtered_feature_bc_matrix.h5',
                           load_images=True, gex_only=False)

    meta = SAMPLE_METADATA[sample]
    model = CitegeistModel(sample_name=sample, adata=adata, output_folder=str(sample_dir))
    model.split_adata()

    # Filter NaN spatial coords
    for attr in ('gene_expression_adata', 'antibody_capture_adata'):
        ad = getattr(model, attr)
        if ad is not None and 'spatial' in ad.obsm:
            coords = ad.obsm['spatial']
            valid = np.all(np.isfinite(coords), axis=1)
            n_invalid = (~valid).sum()
            if n_invalid > 0:
                logger.warning('Filtering %d NaN spatial coords from %s', n_invalid, attr)
                setattr(model, attr, ad[valid].copy())

    model.filter_gex(min_counts=meta['min_counts'])
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(MODEL_PROFILE_DICT)

    # Align ensemble proportions to model's filtered spots
    model_spots = list(model.gene_expression_adata.obs_names)
    common = [s for s in model_spots if s in ensemble_df.index]
    logger.info('Aligned proportions: %d model spots, %d ensemble spots, %d common',
                len(model_spots), len(ensemble_df), len(common))
    model.results['cell_prop'] = ensemble_df.loc[common]
    # Also filter GEX adata to common spots
    model.gene_expression_adata = model.gene_expression_adata[common].copy()
    if model.antibody_capture_adata is not None:
        ab_common = [s for s in common if s in model.antibody_capture_adata.obs_names]
        model.antibody_capture_adata = model.antibody_capture_adata[ab_common].copy()

    # Filter cell data to common spots only
    common_set = set(common)
    cell_spot_map = cell_spot_map[cell_spot_map['spot_barcode'].isin(common_set)].copy()
    valid_cells = set(cell_spot_map['cell_id'].values)
    cell_assignments = {k: v for k, v in cell_assignments.items() if k in valid_cells}
    logger.info('Filtered to %d cells in %d common spots', len(cell_assignments), len(common))

    # Run SACE
    spot_type_gex, cell_adata, diagnostics = model.run_sace_allocation(
        cell_assignments=cell_assignments,
        cell_spot_map=cell_spot_map,
        max_iter=1,
    )

    cell_adata.write_h5ad(str(h5ad_path))
    logger.info('Saved %s: %d cells x %d genes', h5ad_path.name, *cell_adata.shape)

logger.info('=== SACE GEX complete for all samples ===')
"
