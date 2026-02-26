#!/bin/bash
#SBATCH --job-name=vae_pipeline
#SBATCH --output=logs/vae_pipeline_%j.out
#SBATCH --error=logs/vae_pipeline_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# VAE + Sinkhorn Single-Cell Assignment Pipeline
# Runs: patch preparation -> VAE training -> prototype training

set -e  # Exit on error

# Environment setup
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Paths
BASE_DIR="Benchmarking/xenium_benchmarking"
IMAGE_BASE="${BASE_DIR}/scResolve/images/morphology_hires"
SC_OUTPUT="${BASE_DIR}/CITEgeist/output/single_cell_resolution"
PROP_BASE="${BASE_DIR}/CITEgeist/output/hybrid_cellpose"
OUTPUT_DIR="${BASE_DIR}/CITEgeist/output/vae_sinkhorn"
PATCHES_DIR="${OUTPUT_DIR}/patches"

mkdir -p ${OUTPUT_DIR}/{patches,vae,prototypes,evaluation}
mkdir -p ${PATCHES_DIR}

echo "============================================"
echo "VAE + Sinkhorn Pipeline Starting"
echo "============================================"

# Step 1: Prepare patches for all regions
echo ""
echo "Step 1: Preparing patches for all regions..."
echo "============================================"

for region in 0 1 2 3 4; do
    echo "Processing Xenium_region_${region}..."

    IMAGE_PATH="${IMAGE_BASE}/Xenium_region_${region}/morphology.png"
    MASK_PATH="${SC_OUTPUT}/Xenium_region_${region}/Xenium_region_${region}_cellpose_mask.npy"
    SPOT_MAP="${SC_OUTPUT}/Xenium_region_${region}/Xenium_region_${region}_nuclei_spot_map.csv"
    REGION_OUTPUT="${PATCHES_DIR}/region_${region}"

    if [ ! -f "${MASK_PATH}" ]; then
        echo "  Skipping region ${region} - no mask found"
        continue
    fi

    python -m Benchmarking.xenium_benchmarking.CITEgeist.src.prepare_patches \
        --image "${IMAGE_PATH}" \
        --mask "${MASK_PATH}" \
        --nuclei_spot_map "${SPOT_MAP}" \
        --output_dir "${REGION_OUTPUT}" \
        --expansion 0.75 \
        --patch_size 96

    echo "  Done with region ${region}"
done

# Step 2: Combine patches from all regions for VAE training
echo ""
echo "Step 2: Combining patches for VAE training..."
echo "============================================"

python -c "
import numpy as np
from pathlib import Path
from tqdm import tqdm

patches_dir = Path('${PATCHES_DIR}')
all_patches = []

for region_dir in sorted(patches_dir.glob('region_*')):
    for patch_file in tqdm(list(region_dir.glob('spot_*_patches.npy')), desc=f'Loading {region_dir.name}'):
        patches = np.load(patch_file)
        all_patches.append(patches)

if all_patches:
    combined = np.concatenate(all_patches, axis=0)
    print(f'Combined {len(combined)} patches from all regions')

    # Save in chunks for VAE training
    chunk_size = 10000
    combined_dir = patches_dir / 'combined'
    combined_dir.mkdir(exist_ok=True)

    for i in range(0, len(combined), chunk_size):
        chunk = combined[i:i+chunk_size]
        np.save(combined_dir / f'patches_chunk_{i//chunk_size:04d}.npy', chunk)

    print(f'Saved {len(combined)//chunk_size + 1} chunk files')
"

# Step 3: Train VAE (Stage 1)
echo ""
echo "Step 3: Training VAE (Stage 1)..."
echo "============================================"

python -m CITEgeist.model.train_vae \
    --patches-dir "${PATCHES_DIR}/combined" \
    --output-dir "${OUTPUT_DIR}/vae" \
    --latent-dim 128 \
    --batch-size 64 \
    --epochs 50 \
    --lr 1e-4 \
    --beta 0.5 \
    --device cpu

# Step 4: Prepare proportions file for prototype training
echo ""
echo "Step 4: Preparing proportions file..."
echo "============================================"

python -c "
import pandas as pd
from pathlib import Path

patches_dir = Path('${PATCHES_DIR}')
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

# Step 5: Prepare spot-level patches for prototype training
echo ""
echo "Step 5: Reorganizing patches for prototype training..."
echo "============================================"

python -c "
import shutil
from pathlib import Path

src_dir = Path('${PATCHES_DIR}')
dst_dir = Path('${OUTPUT_DIR}/spot_patches')
dst_dir.mkdir(exist_ok=True)

count = 0
for region_dir in sorted(src_dir.glob('region_*')):
    region = region_dir.name.replace('region_', '')
    for patch_file in region_dir.glob('spot_*_patches.npy'):
        spot_id = patch_file.stem.replace('_patches', '')
        # Add region prefix to avoid collisions
        new_name = f'r{region}_{spot_id}_patches.npy'
        shutil.copy(patch_file, dst_dir / new_name)
        count += 1

print(f'Copied {count} spot patch files')
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

# Step 6: Train Prototypes (Stage 2)
echo ""
echo "Step 6: Training Prototypes (Stage 2)..."
echo "============================================"

python -m CITEgeist.model.train_prototypes \
    --vae-checkpoint "${OUTPUT_DIR}/vae/vae_final.pt" \
    --patches-dir "${OUTPUT_DIR}/spot_patches" \
    --proportions-file "${OUTPUT_DIR}/proportions_for_training.csv" \
    --output-dir "${OUTPUT_DIR}/prototypes" \
    --n-types 7 \
    --latent-dim 128 \
    --projection-dim 32 \
    --epochs 30 \
    --lr 1e-3 \
    --sinkhorn-temp 0.1 \
    --device cpu

echo ""
echo "============================================"
echo "VAE + Sinkhorn Pipeline Complete!"
echo "============================================"
echo "VAE checkpoint: ${OUTPUT_DIR}/vae/vae_final.pt"
echo "Prototype checkpoint: ${OUTPUT_DIR}/prototypes/prototypes_final.pt"
