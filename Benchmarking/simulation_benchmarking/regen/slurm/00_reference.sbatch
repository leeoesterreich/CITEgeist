#!/bin/bash
#SBATCH --job-name=regen_ref
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --time=03:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
# A10 Step 3a: rebuild RCTD/Seurat Rds (r_deconv_env, no module load r) and CARD reference_csv (CITEgeist_env).
# Rds is built in regen/reference from a copy of the h5ad (models/*.h5seurat is NOT touched), then copied to models/.
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
ENVS=${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}
eval "$(conda shell.bash hook 2>/dev/null)"
T0=$(date +%s)
md5sum $B/models/Wu_scRNA_ref_ERpos.h5ad | tee $B/regen/reference/input_h5ad.md5
mkdir -p $B/regen/reference
cp $B/models/Wu_scRNA_ref_ERpos.h5ad $B/regen/reference/Wu_scRNA_ref_ERpos.h5ad
cmp $B/models/Wu_scRNA_ref_ERpos.h5ad $B/regen/reference/Wu_scRNA_ref_ERpos.h5ad
conda activate $ENVS/r_deconv_env
echo "R: $(which Rscript)"; Rscript --version
cd ${CITEGEIST_REPO_ROOT}
Rscript $B/models/convert_reference_formats.R $B/regen/reference/Wu_scRNA_ref_ERpos.h5ad
md5sum $B/regen/reference/Wu_scRNA_ref_ERpos.Rds | tee $B/regen/reference/Wu_scRNA_ref_ERpos.Rds.md5
cp -n $B/regen/reference/Wu_scRNA_ref_ERpos.Rds $B/models/Wu_scRNA_ref_ERpos.Rds
cmp $B/regen/reference/Wu_scRNA_ref_ERpos.Rds $B/models/Wu_scRNA_ref_ERpos.Rds
T1=$(date +%s); echo "RDS_RUNTIME_SEC $((T1-T0))"
# CARD reference_csv: run the ORIGINAL script from a node-local copy (NFS wide-CSV corruption), then verified copy.
conda deactivate; conda activate $ENVS/CITEgeist_env
LOC=$(mktemp -d)
cp $B/CARD/convert_reference_to_csv.py $LOC/
cd $LOC
python -u convert_reference_to_csv.py
ls -la $LOC/reference_csv
mkdir -p $B/CARD/reference_csv
for f in reference_counts.csv reference_cell_types.csv; do
  cp $LOC/reference_csv/$f $B/CARD/reference_csv/$f.tmp && mv $B/CARD/reference_csv/$f.tmp $B/CARD/reference_csv/$f
  cmp $LOC/reference_csv/$f $B/CARD/reference_csv/$f
done
md5sum $B/CARD/reference_csv/*.csv | tee $B/regen/reference/card_reference_csv.md5
rm -rf $LOC
T2=$(date +%s); echo "CARD_CSV_RUNTIME_SEC $((T2-T1))"
echo DONE
