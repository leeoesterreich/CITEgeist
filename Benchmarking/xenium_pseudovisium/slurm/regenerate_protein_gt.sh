#!/bin/bash
#SBATCH --job-name=regen_protein_gt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

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
set -e

export PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env/bin:$PATH"
source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

REPO="${CITEGEIST_REPO_ROOT}"

echo "Re-generating protein GT with achievable-7 types..."
echo "Start time: $(date)"

python "${REPO}/Benchmarking/xenium_pseudovisium/src/create_protein_gt.py"

echo ""
echo "Verifying output..."
python -c "
import pandas as pd
gt = pd.read_csv('${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth/Xenium_region_0_prop.csv', index_col=0)
print('Columns:', list(gt.columns))
print('Shape:', gt.shape)
ct_cols = [c for c in gt.columns if c not in ['n_cells','spot_x','spot_y']]
print(f'Cell type columns ({len(ct_cols)}):', ct_cols)
for c in ct_cols:
    print(f'  {c}: mean={gt[c].mean():.4f}')
"

echo ""
echo "Done! End time: $(date)"
