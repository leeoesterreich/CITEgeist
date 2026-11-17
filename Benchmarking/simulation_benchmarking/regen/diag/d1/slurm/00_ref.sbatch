#!/bin/bash
#SBATCH --job-name=d1_ref
#SBATCH --cluster=htc
#SBATCH --partition=htc

#SBATCH --time=06:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
# A10 diagnostic D1: full-gene ER+ reference; fixed script copies, only the reference changes.
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../../../_paths.sh"
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
set -eo pipefail
B=${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking
D=${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/regen/diag/d1
ENVS=${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate $ENVS/CITEgeist_env
cd ${CITEGEIST_REPO_ROOT}
T0=$(date +%s)
python -u $D/src/build_full_reference.py $D/reference/Wu_scRNA_ref_ERpos_fullgene.h5ad
md5sum $D/reference/Wu_scRNA_ref_ERpos_fullgene.h5ad | tee $D/reference/fullgene_h5ad.md5
conda deactivate; conda activate $ENVS/r_deconv_env
Rscript $B/models/convert_reference_formats.R $D/reference/Wu_scRNA_ref_ERpos_fullgene.h5ad
md5sum $D/reference/Wu_scRNA_ref_ERpos_fullgene.Rds | tee $D/reference/fullgene_Rds.md5
echo "RDS_DONE $(( $(date +%s)-T0 ))"
conda deactivate; conda activate $ENVS/CITEgeist_env
LOC=$(mktemp -d); cp $D/src/convert_reference_to_csv_fullgene.py $LOC/; cd $LOC
python -u convert_reference_to_csv_fullgene.py
mkdir -p $D/reference/reference_csv
for f in reference_counts.csv reference_cell_types.csv; do cp $LOC/reference_csv/$f $D/reference/reference_csv/$f.tmp && mv $D/reference/reference_csv/$f.tmp $D/reference/reference_csv/$f; cmp $LOC/reference_csv/$f $D/reference/reference_csv/$f; done
md5sum $D/reference/reference_csv/*.csv | tee $D/reference/fullgene_card_csv.md5
rm -rf $LOC; echo "ALL_DONE $(( $(date +%s)-T0 ))"
