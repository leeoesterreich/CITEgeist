#!/bin/bash
#SBATCH --job-name=singler_postprocess
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

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
set -euo pipefail

PY_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"
export PATH="${PY_ENV}/bin:$PATH"
eval "$(conda shell.bash hook)"
set +u  # conda activate.d binutils script references unbound ADDR2LINE under `set -u`
conda activate "${PY_ENV}"
set -u

BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
GT_DIR="${BASE_DIR}/ground_truth_singler"
LABELS_DIR="${GT_DIR}/singler_labels"
MAPPING_CSV="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt/cell_to_spot_mapping.csv"

echo "=== Step 1: CD4/CD8 protein split ==="
python "${GT_DIR}/src/split_cd4_cd8.py" \
    --singler-labels "${LABELS_DIR}/singler_all_cells.csv" \
    --protein-csv "${LABELS_DIR}/xenium_cell_proteins.csv" \
    --output-dir "${LABELS_DIR}"

echo "=== Step 2: Aggregate to spots (6-type) ==="
python "${GT_DIR}/src/aggregate_to_spots.py" \
    --labels-csv "${LABELS_DIR}/singler_all_cells.csv" \
    --mapping-csv "${MAPPING_CSV}" \
    --output-dir "${GT_DIR}/ground_truth_6type" \
    --n-types 6

echo "=== Step 3: Aggregate to spots (7-type) ==="
python "${GT_DIR}/src/aggregate_to_spots.py" \
    --labels-csv "${LABELS_DIR}/singler_labels_7type.csv" \
    --mapping-csv "${MAPPING_CSV}" \
    --output-dir "${GT_DIR}/ground_truth_7type" \
    --n-types 7

echo "=== Step 4: Validate annotations + concordance ==="
python "${GT_DIR}/src/validate_annotations.py" \
    --singler-labels "${LABELS_DIR}/singler_all_cells.csv" \
    --xenium-counts "${LABELS_DIR}/xenium_cell_gex.csv" \
    --output-dir "${LABELS_DIR}"

echo "=== Done. Check validation_report.json and cd4_cd8_split_report.json ==="
