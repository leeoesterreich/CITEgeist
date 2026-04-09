#!/bin/bash
#SBATCH --job-name=test_module6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=output/module6_test_%j.out
#SBATCH --error=output/module6_test_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

set -u

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python -u -c "
import sys, logging, json
import numpy as np
import pandas as pd
sys.path.insert(0, '.')

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

import squidpy as sq
from CITEgeist.model import CitegeistModel

SAMPLE = 'HCC22-088-P1-S1'
PATIENT_DIR = '/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files'
QP_DIR = 'output/module3_cuopt_qp'
OUT_DIR = 'output/module6_test'

import os
os.makedirs(OUT_DIR, exist_ok=True)

# Load patient data
logger.info('Loading %s...', SAMPLE)
adata = sq.read.visium(
    f'{PATIENT_DIR}/{SAMPLE}/outs',
    counts_file='filtered_feature_bc_matrix.h5',
    load_images=True,
    gex_only=False,
)

model = CitegeistModel(sample_name=SAMPLE, adata=adata, output_folder=OUT_DIR)
model.split_adata()

# Filter NaN spatial coords
valid = np.all(np.isfinite(model.gene_expression_adata.obsm['spatial']), axis=1)
n_invalid = (~valid).sum()
if n_invalid > 0:
    logger.warning('Filtering %d NaN spatial coords', n_invalid)
    model.gene_expression_adata = model.gene_expression_adata[valid].copy()
    model.antibody_capture_adata = model.antibody_capture_adata[valid].copy()

model.filter_gex(min_counts=100)
model.preprocess_gex()
model.preprocess_antibody()

# Load existing QP proportions instead of re-running Module 3
qp_path = f'{QP_DIR}/{SAMPLE}/{SAMPLE}_cell_prop_global_results.csv'
logger.info('Loading QP proportions from %s', qp_path)
qp_props = pd.read_csv(qp_path, index_col=0)

# Align to model's filtered spots
common = [s for s in model.gene_expression_adata.obs_names if s in qp_props.index]
logger.info('Aligned %d/%d spots', len(common), len(model.gene_expression_adata))

# Drop recon_error column if present
prop_cols = [c for c in qp_props.columns if c != 'recon_error']
model.results['cell_prop'] = qp_props.loc[common, prop_cols]
model.gene_expression_adata = model.gene_expression_adata[common].copy()
model.antibody_capture_adata = model.antibody_capture_adata[common].copy()

# Load profile dict (needed for model state but not for Module 6)
PROFILE_DICT = {
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
model.load_cell_profile_dict(PROFILE_DICT)

# === Run Module 6 ===
logger.info('=== Running Module 6: Functional Protein Annotation ===')
results = model.run_functional_annotation(
    max_iter=200,
    lr=0.01,
    device='cpu',
)

# === Validate ===
logger.info('')
logger.info('=== VALIDATION ===')

lam_df = results['functional_lambda']
summary = results['functional_summary']

# 1. Check lambda learned
logger.info('Lambda matrix shape: %s', lam_df.shape)
logger.info('Top 10 emission rates:')
lam_flat = lam_df.stack().sort_values(ascending=False)
for (ct, marker), val in lam_flat.head(10).items():
    logger.info('  %s × %s: λ=%.1f', ct, marker, val)

# 2. Check gates
n_gated = sum(1 for v in summary.values() if isinstance(v, dict) and v.get('gmm_converged', False))
n_total = len(summary)
logger.info('GMM gating: %d/%d pairs converged', n_gated, n_total)

# 3. Check spatial significance
n_sig = sum(1 for v in summary.values() if isinstance(v, dict) and v.get('morans_p', 1.0) < 0.05)
logger.info('Spatially significant pairs (Morans I p<0.05): %d', n_sig)

# 4. Print biological highlights
logger.info('')
logger.info('=== BIOLOGICAL HIGHLIGHTS ===')
for (ct, marker), stats in sorted(summary.items()):
    if not isinstance(stats, dict):
        continue
    if stats.get('gmm_converged', False):
        logger.info('  %s × %s: %.1f%% positive, λ=%.1f, Morans I=%.3f (p=%.2e)',
                    ct, marker,
                    stats.get('pct_positive', 0),
                    stats.get('lambda', 0),
                    stats.get('morans_i', 0),
                    stats.get('morans_p', 1))

# 5. Check output files
for suffix in ['_functional_lambda.csv', '_functional_intensity.csv',
               '_functional_gates.csv', '_functional_summary.json']:
    fpath = f'{OUT_DIR}/{SAMPLE}{suffix}'
    exists = os.path.exists(fpath)
    logger.info('Output %s: %s', suffix, 'OK' if exists else 'MISSING')

logger.info('')
logger.info('=== Module 6 integration test COMPLETE ===')
"
