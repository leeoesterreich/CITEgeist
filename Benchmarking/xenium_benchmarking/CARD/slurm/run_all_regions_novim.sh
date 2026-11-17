#!/bin/bash
#SBATCH --job-name=card_novim
#SBATCH --array=0-4
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Run CARD deconvolution with VIM excluded (fairer benchmark)
# VIM is a mesenchymal marker expressed in RCC cells, not fibroblast-specific

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
set -e

REGION_ID=${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "CARD Benchmark (VIM excluded) - Region ${REGION_ID}"
echo "=============================================="
echo "Start time: $(date)"

# Paths
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
REF_DIR="${BASE_DIR}/reference_data/GSE156632/processed_protein7/rctd"
SPATIAL_DIR="${BASE_DIR}/CARD/output_protein_gt/temp_csvs"
OUTPUT_DIR="${BASE_DIR}/CARD/output_protein_gt_novim"

mkdir -p "${OUTPUT_DIR}/Xenium_region_${REGION_ID}"

# Activate CARD environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CARD_env

echo "R: $(which R)"

# Run CARD with VIM excluded
Rscript -e "
library(CARD)
library(Matrix)

cat('Loading data...\n')
spatial_counts <- as.matrix(read.csv('${SPATIAL_DIR}/Xenium_region_${REGION_ID}_counts.csv', row.names=1, check.names=FALSE))
spatial_coords <- read.csv('${SPATIAL_DIR}/Xenium_region_${REGION_ID}_coords.csv', row.names=1, check.names=FALSE)
ref_counts <- as.matrix(read.csv('${REF_DIR}/counts.csv', row.names=1, check.names=FALSE))
ref_cell_types <- read.csv('${REF_DIR}/cell_types.csv', row.names=1, check.names=FALSE)

cat(sprintf('  Spatial: %d genes x %d spots\n', nrow(spatial_counts), ncol(spatial_counts)))
cat(sprintf('  Reference: %d genes x %d cells\n', nrow(ref_counts), ncol(ref_counts)))

# EXCLUDE VIM
cat('Excluding VIM from both datasets...\n')
spatial_counts <- spatial_counts[rownames(spatial_counts) != 'VIM', ]
ref_counts <- ref_counts[rownames(ref_counts) != 'VIM', ]
cat(sprintf('  After exclusion - Spatial: %d genes, Reference: %d genes\n', nrow(spatial_counts), nrow(ref_counts)))

# Prepare metadata
spatial_location <- as.data.frame(spatial_coords)
colnames(spatial_location) <- c('x', 'y')

sc_meta <- data.frame(
  cellID = rownames(ref_cell_types),
  cellType = ref_cell_types\$cell_type,
  sampleInfo = 'sample1',
  row.names = rownames(ref_cell_types)
)

# Run CARD
cat('Creating CARD object...\n')
CARD_obj <- createCARDObject(
  sc_count = ref_counts,
  sc_meta = sc_meta,
  spatial_count = spatial_counts,
  spatial_location = spatial_location,
  ct.varname = 'cellType',
  ct.select = unique(sc_meta\$cellType),
  sample.varname = 'sampleInfo',
  minCountGene = 100,
  minCountSpot = 5
)

cat('Running CARD deconvolution...\n')
CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)

# Save results
proportions <- as.data.frame(CARD_obj@Proportion_CARD)
output_file <- '${OUTPUT_DIR}/Xenium_region_${REGION_ID}/Xenium_region_${REGION_ID}_cell_prop_finetuned_results.csv'
write.csv(proportions, output_file, row.names=TRUE)

cat(sprintf('\nResults: %d spots x %d cell types\n', nrow(proportions), ncol(proportions)))
cat('Mean proportions:\n')
means <- colMeans(proportions)
for (ct in names(sort(means, decreasing=TRUE))) {
  cat(sprintf('  %s: %.1f%%\n', ct, means[ct] * 100))
}
cat(sprintf('\nSaved to: %s\n', output_file))
"

echo ""
echo "=============================================="
echo "Region ${REGION_ID} Complete!"
echo "=============================================="
echo "End time: $(date)"
