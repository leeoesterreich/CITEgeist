#!/bin/bash
#SBATCH --job-name=regen_protein_gt
#SBATCH --output=regen_protein_gt_%j.out
#SBATCH --error=regen_protein_gt_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Regenerate protein-based ground truth with αSMA-only fibroblast definition
# (Vimentin gate removed due to lack of specificity in RCC)

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../_paths.sh"
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
echo "=== Regenerating Protein GT (αSMA-only fibroblasts) ==="
echo "Start time: $(date)"

# Activate CITEgeist conda environment
export PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env/bin:$PATH"
source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

cd ${CITEGEIST_REPO_ROOT}

# Backup old ground truth
BACKUP_DIR="Benchmarking/xenium_pseudovisium/data_protein_gt_backup_$(date +%Y%m%d_%H%M%S)"
echo "Backing up old ground truth to: ${BACKUP_DIR}"
cp -r Benchmarking/xenium_pseudovisium/data_protein_gt "${BACKUP_DIR}"

# Regenerate ground truth
python Benchmarking/xenium_pseudovisium/src/create_protein_gt.py

echo "=== Done ==="
echo "End time: $(date)"
