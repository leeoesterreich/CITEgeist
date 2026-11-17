#!/bin/bash
#SBATCH --job-name=regen_seurat
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-9
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=2
# A10 Step 3: rerun seurat with the ORIGINAL per-condition script and settings; outputs only under regen/results.
# Array task t: t<5 -> high_seg rep t ; t>=5 -> mixed rep t-5.
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
set -eo pipefail
B=${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking
if [ "$SLURM_ARRAY_TASK_ID" -lt 5 ]; then COND=high_seg; REP=$SLURM_ARRAY_TASK_ID; else COND=mixed; REP=$((SLURM_ARRAY_TASK_ID-5)); fi
OUT=$B/regen/results/seurat/$COND/Wu_rep_$REP
mkdir -p $OUT
eval "$(conda shell.bash hook 2>/dev/null)"
T0=$(date +%s)
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env
test -s $B/models/Wu_scRNA_ref_ERpos.Rds
# Seurat_regen.R converts the h5ad to h5seurat NEXT TO the input; work on a copy so replicates/ is untouched.
SCR=$OUT/scratch; mkdir -p $SCR
cp $B/replicates/$COND/h5ad_objects/Wu_rep_${REP}_GEX.h5ad $SCR/
md5sum $B/replicates/$COND/h5ad_objects/Wu_rep_${REP}_GEX.h5ad
cd $B/regen/src
# predictions CSV is written before the (non-essential) h5Seurat/PDF steps; require the CSV, tolerate a late plot error
Rscript --vanilla Seurat_regen.R $SCR/Wu_rep_${REP}_GEX.h5ad $OUT Wu_rep_$REP || echo "WARNING: Seurat script exited non-zero after predictions"
test -s $OUT/Wu_rep_${REP}_Seurat_deconv_predictions.csv
rm -rf $SCR
echo "RUNTIME_SEC $(( $(date +%s) - T0 )) method=Seurat cond=$COND rep=$REP"
