#!/bin/bash
#SBATCH --job-name=vae_pipeline_v2
#SBATCH --output=slurm/logs/vae_pipeline_v2_%j.out
#SBATCH --error=slurm/logs/vae_pipeline_v2_%j.err
#SBATCH --time=24:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# VAE + Prototype Training Pipeline with v2 Preprocessing
#
# Uses the new preprocessing with:
# 1. Global percentile normalization (preserves intensity across patches)
# 2. Size features (log1p(w, h, area) for each nucleus)
#
# Requires: sbatch_preprocess_v2.sh to have run first
#
# Stages:
# 1. Combine patches from all regions
# 2. Train VAE (2-channel: DAPI + boundary)
# 3. Prepare proportions file
# 4. Reorganize spot patches with size features
# 5. Train prototypes with DirectSoftmax + size features

set -e  # Exit on error

# Environment setup
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load cuda/12.1.1

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Paths
BASE_DIR="Benchmarking/xenium_benchmarking"
PATCHES_INPUT="${BASE_DIR}/CITEgeist/output/patches_v2"
PROP_BASE="${BASE_DIR}/CITEgeist/output/hybrid_cellpose"
OUTPUT_DIR="${BASE_DIR}/CITEgeist/output/vae_sinkhorn_v2"

mkdir -p ${OUTPUT_DIR}/{patches_combined,vae,prototypes,spot_patches,evaluation}

echo "============================================"
echo "VAE Pipeline v2 (with global normalization)"
echo "============================================"
echo "Input patches: ${PATCHES_INPUT}"
echo "Output: ${OUTPUT_DIR}"
echo "============================================"

# Step 1: Verify preprocessing was done
echo ""
echo "Step 1: Verifying v2 preprocessing..."
echo "============================================"

for region in 0 1 2 3 4; do
    GLOBAL_STATS="${PATCHES_INPUT}/region_${region}/global_stats.npz"
    if [ ! -f "${GLOBAL_STATS}" ]; then
        echo "ERROR: Missing global_stats.npz for region ${region}"
        echo "Run sbatch_preprocess_v2.sh first!"
        exit 1
    fi
done
echo "All regions have global_stats.npz - v2 preprocessing verified"

# Step 2: Combine patches for VAE training
echo ""
echo "Step 2: Combining patches for VAE training..."
echo "============================================"

python -c "
import numpy as np
from pathlib import Path
from tqdm import tqdm

patches_dir = Path('${PATCHES_INPUT}')
output_dir = Path('${OUTPUT_DIR}/patches_combined')
all_patches = []

for region_dir in sorted(patches_dir.glob('region_*')):
    for patch_file in tqdm(list(region_dir.glob('spot_*_patches.npy')), desc=f'Loading {region_dir.name}'):
        patches = np.load(patch_file)
        all_patches.append(patches)

if all_patches:
    combined = np.concatenate(all_patches, axis=0)
    print(f'Combined {len(combined)} patches from all regions')
    print(f'Shape: {combined.shape}')

    # Verify normalization (should be in [0, 1] range)
    print(f'Value range: [{combined.min():.4f}, {combined.max():.4f}]')

    # Save in chunks for VAE training
    chunk_size = 10000
    for i in range(0, len(combined), chunk_size):
        chunk = combined[i:i+chunk_size]
        np.save(output_dir / f'patches_chunk_{i//chunk_size:04d}.npy', chunk)

    print(f'Saved {len(combined)//chunk_size + 1} chunk files')
else:
    print('ERROR: No patches found!')
    exit(1)
"

# Step 3: Train VAE with 2 channels
echo ""
echo "Step 3: Training VAE (2-channel)..."
echo "============================================"

python -m CITEgeist.model.train_vae \
    --patches-dir "${OUTPUT_DIR}/patches_combined" \
    --output-dir "${OUTPUT_DIR}/vae" \
    --in-channels 2 \
    --latent-dim 128 \
    --batch-size 256 \
    --epochs 50 \
    --lr 1e-4 \
    --beta 0.5 \
    --device cuda \
    --num-workers 8

# Verify VAE checkpoint has preprocessing_version
echo "Verifying VAE checkpoint..."
python -c "
import torch
ckpt = torch.load('${OUTPUT_DIR}/vae/vae_final.pt', map_location='cpu')
if 'preprocessing_version' not in ckpt:
    print('ERROR: VAE checkpoint missing preprocessing_version!')
    exit(1)
print(f'VAE preprocessing_version: {ckpt[\"preprocessing_version\"]}')
"

# Step 4: Prepare proportions file
echo ""
echo "Step 4: Preparing proportions file..."
echo "============================================"

python -c "
import pandas as pd
from pathlib import Path

prop_base = Path('${PROP_BASE}')
output = []

for region in range(5):
    prop_file = prop_base / f'Xenium_region_{region}' / f'Xenium_region_{region}_cell_prop_finetuned_results.csv'
    if not prop_file.exists():
        print(f'Skipping region {region} - no proportions')
        continue

    df = pd.read_csv(prop_file, index_col=0)
    df['spot_id'] = df.index
    df['region'] = region
    df = df.reset_index(drop=True)
    output.append(df)

combined = pd.concat(output, ignore_index=True)
combined.to_csv('${OUTPUT_DIR}/proportions_combined.csv', index=False)
print(f'Combined proportions for {len(combined)} spots from {len(output)} regions')
"

# Step 5: Reorganize spot patches with size features
echo ""
echo "Step 5: Reorganizing patches with size features..."
echo "============================================"

python -c "
import shutil
from pathlib import Path

src_dir = Path('${PATCHES_INPUT}')
dst_dir = Path('${OUTPUT_DIR}/spot_patches')
dst_dir.mkdir(exist_ok=True)

patch_count = 0
size_count = 0

for region_dir in sorted(src_dir.glob('region_*')):
    region = region_dir.name.replace('region_', '')

    # Copy patches
    for patch_file in region_dir.glob('spot_*_patches.npy'):
        spot_id = patch_file.stem.replace('_patches', '').replace('spot_', '', 1)
        new_name = f'r{region}_{spot_id}_patches.npy'
        shutil.copy(patch_file, dst_dir / new_name)
        patch_count += 1

    # Copy size features
    for size_file in region_dir.glob('spot_*_sizes.npy'):
        spot_id = size_file.stem.replace('_sizes', '').replace('spot_', '', 1)
        new_name = f'r{region}_{spot_id}_sizes.npy'
        shutil.copy(size_file, dst_dir / new_name)
        size_count += 1

    # Also copy global_stats.npz for validation
    global_stats = region_dir / 'global_stats.npz'
    if global_stats.exists():
        shutil.copy(global_stats, dst_dir / f'r{region}_global_stats.npz')

print(f'Copied {patch_count} patch files')
print(f'Copied {size_count} size feature files')

if patch_count != size_count:
    print(f'WARNING: Mismatch between patches ({patch_count}) and sizes ({size_count})')

# Copy one global_stats.npz as the canonical one
# (they should all be the same since preprocessing uses the same image)
shutil.copy(src_dir / 'region_0' / 'global_stats.npz', dst_dir / 'global_stats.npz')
"

# Update proportions to match new spot IDs
python -c "
import pandas as pd

df = pd.read_csv('${OUTPUT_DIR}/proportions_combined.csv')
df['spot_id'] = 'r' + df['region'].astype(str) + '_' + df['spot_id']
df = df.drop('region', axis=1)
df.to_csv('${OUTPUT_DIR}/proportions_for_training.csv', index=False)
print(f'Updated {len(df)} spot IDs')
"

# Step 6: Train Prototypes with DirectSoftmax + size features
echo ""
echo "Step 6: Training Prototypes (DirectSoftmax + size features)..."
echo "============================================"

python -m CITEgeist.model.train_prototypes \
    --vae-checkpoint "${OUTPUT_DIR}/vae/vae_final.pt" \
    --patches-dir "${OUTPUT_DIR}/spot_patches" \
    --proportions-file "${OUTPUT_DIR}/proportions_for_training.csv" \
    --output-dir "${OUTPUT_DIR}/prototypes" \
    --n-types 7 \
    --latent-dim 128 \
    --projection-dim 32 \
    --epochs 50 \
    --lr 1e-3 \
    --use-direct-softmax \
    --softmax-temperature 0.1 \
    --repulsion-weight 1.0 \
    --repulsion-margin 0.5 \
    --entropy-weight 0.01 \
    --device cuda

echo ""
echo "============================================"
echo "VAE Pipeline v2 Complete!"
echo "============================================"
echo "VAE checkpoint: ${OUTPUT_DIR}/vae/vae_final.pt"
echo "Prototype checkpoint: ${OUTPUT_DIR}/prototypes/prototypes_final.pt"
echo ""
echo "Key improvements:"
echo "  - Global percentile normalization (preserves intensity)"
echo "  - Size features (log1p(w, h, area))"
echo "  - Preprocessing version tracked in checkpoints"
echo "============================================"
