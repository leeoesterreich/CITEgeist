#!/bin/bash
#SBATCH --job-name=v2off
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --exclude=htc-1024-n2,htc-n37
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
# One official v2 Xenium stage (v42 plan Task 8). Env: STAGE (x1 x2 x3 x4 x5 x5_ddof1 x6 b15 x9) and SHA (the commit that holds this code).
# Submitted by .commands/submit-v2-official.sh. Runs from the main checkout (the stage scripts resolve untracked data and imports from their
# own checkout path, so a detached worktree cannot run them) after proving that every tracked file of the code tree (v2_official.CODE_TREE:
# CITEgeist/, Benchmarking/xenium_benchmarking/ and the imported scripts) equals commit SHA; the stamp step proves it again.
# Every stage writes NEW dirs, runs its acceptance check against the experiment record where one exists, stamps its outputs
# (scripts/v42_stamp.py) and last writes the stage-complete marker.
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
set -euo pipefail
: "${STAGE:?set STAGE}"; : "${SHA:?set SHA to the commit that holds this code}"
PY="$(cg_env_python CITEgeist_env)"
REPO=${CITEGEIST_REPO_ROOT}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID PYTHONPATH=$REPO
cd "$REPO"
XB=Benchmarking/xenium_benchmarking
RF=$XB/rna_refinement
"$PY" -u $RF/v2_official.py check-code "$SHA"
STARTED=$("$PY" -c "import sys; sys.path.insert(0,'scripts'); import v42_stamp; print(v42_stamp.utc_now())")
eval "$("$PY" -u $RF/v2_official.py paths)"
GT_GEX=$XB/ground_truth_xenium_rnaprot/ground_truth_7type/ground_truth_gex
stamp() { "$PY" -u $RF/v2_official.py stamp --started "$STARTED" "$@"; }
echo "stage=$STAGE sha=$SHA host=$(hostname) start=$STARTED"
case "$STAGE" in
  x1)  "$PY" -u $RF/v2_official.py x1 ;;
  x2)  "$PY" -u $XB/spots_baseline/score_p138_strict.py --labelset rnaprot --arm v2
       "$PY" -u $RF/v2_official.py x2-check
       for N in 6 7; do
         P=PRED$N
         stamp --stage X2 --producer $XB/spots_baseline/score_p138_strict.py --role x2_strict_scoring \
           --artifact "$RES_RNAPROT/method_comparison_singler_${N}type/CITEgeist_${N}type_results.json" --input-glob "${!P}" baseline.csv.gz
       done ;;
  x3)  for R in 0 1 2 3 4; do
         "$PY" -u $XB/spots_baseline/sace_props_ablation.py marker_share "$R" "$SACE" "$PRED7" ms_v41 no_round
       done
       "$PY" -u $RF/v2_official.py x3-manifest ;;
  x4)  export XEN_LABELSET=rnaprot_v41 XEN_REGISTER_V2_SACE=1
       for G in matched own; do
         "$PY" -u $XB/evaluation/src/compare_xenium_gex_methods.py --spots all --genes "$G" --methods CITEgeist_v2_SACE Cell2Location Tangram \
           --method_dir "CITEgeist_v2_SACE=$SACE" --output_dir "$FAIR/$G"
       done
       "$PY" -u $RF/v2_official.py x4-check "$FAIR"
       "$PY" -u $XB/evaluation/src/matched_gene_r.py --all-spots --v2-sace-dir "$SACE" --summary "$FAIR/own/xenium_gex_comparison_7type_summary.csv" --out "$FAIR/matched_gene_r.json"
       for G in matched own; do
         for F in xenium_gex_comparison_7type_summary.csv xenium_gex_comparison_7type.csv; do
           stamp --stage X4 --producer $XB/evaluation/src/compare_xenium_gex_methods.py --role x4_fair_gex_$G --artifact "$FAIR/$G/$F" \
             --input-glob "$SACE" '*_layer.csv' --input-glob $XB/Cell2Location/output_singler_7type '*_layer.csv' \
             --input-glob $XB/Tangram/output_singler_7type_native '*_layer.csv' --input-glob $GT_GEX '*_GT.csv'
         done
       done
       stamp --stage X4 --producer $XB/evaluation/src/matched_gene_r.py --role x4_matched_gene_r --artifact "$FAIR/matched_gene_r.json" \
         --inputs-from-json "$FAIR/matched_gene_r.json" input_sha256 ;;
  x5)  # fair standard for CITEgeist, C2L and Tangram alike: all GT spots, an absent spot is zero (RMSE_SIGMA_ALL_SPOTS)
       export XEN_LABELSET=rnaprot_v41 XEN_V2_SACE_DIR="$SACE" RMSE_SIGMA_OUT_DIR="$FAIR" RMSE_SIGMA_ALL_SPOTS=1
       "$PY" -u scripts/compute_gex_rmse_sd_ratio.py
       "$PY" -u scripts/compute_gex_rmse_sd_ratio_matched_c2l.py
       stamp --stage X5 --producer scripts/compute_gex_rmse_sd_ratio.py --role x5_rmse_over_gene_sd --artifact "$FAIR/rmse_sd_ratio_provenance.json" \
         --inputs-from-json "$FAIR/rmse_sd_ratio_provenance.json" input_sha256
       stamp --stage X5 --producer scripts/compute_gex_rmse_sd_ratio_matched_c2l.py --role x5_rmse_over_gene_sd_matched_c2l \
         --artifact "$FAIR/rmse_sd_ratio_matched_c2l.json" --inputs-from-json "$FAIR/rmse_sd_ratio_matched_c2l.json" input_sha256 ;;
  x5_ddof1)  # corrected X5 rerun (ddof=1): same inputs, NEW dir + run id + marker; $FAIR (the original X5) is only read
       eval "$("$PY" -u $RF/v2_official.py paths --x5-ddof1 | grep -E '^(X5_OUT|RID_X5_DDOF1)=')"
       export XEN_LABELSET=rnaprot_v41 XEN_V2_SACE_DIR="$SACE" RMSE_SIGMA_OUT_DIR="$X5_OUT" RMSE_SIGMA_ALL_SPOTS=1
       "$PY" -u scripts/compute_gex_rmse_sd_ratio.py
       "$PY" -u scripts/compute_gex_rmse_sd_ratio_matched_c2l.py
       stamp --stage X5 --run-id "$RID_X5_DDOF1" --producer scripts/compute_gex_rmse_sd_ratio.py --role x5_rmse_over_gene_sd --artifact "$X5_OUT/rmse_sd_ratio_provenance.json" \
         --inputs-from-json "$X5_OUT/rmse_sd_ratio_provenance.json" input_sha256
       stamp --stage X5 --run-id "$RID_X5_DDOF1" --producer scripts/compute_gex_rmse_sd_ratio_matched_c2l.py --role x5_rmse_over_gene_sd_matched_c2l \
         --artifact "$X5_OUT/rmse_sd_ratio_matched_c2l.json" --inputs-from-json "$X5_OUT/rmse_sd_ratio_matched_c2l.json" input_sha256 ;;
  x6)  mkdir "$SC"
       "$PY" -u $XB/CITEgeist/src/eval_single_cell_assignment.py --stamp --out "$SC/eval_results_v2_rnaprot.json" \
         --pred-template "$PRED7/Xenium_region_{region}/predictions/baseline.csv.gz" \
         --gt-labels $XB/ground_truth_xenium_rnaprot/labels/cell_labels_7type.csv
       "$PY" -u $RF/v2_official.py x6-check "$SC/eval_results_v2_rnaprot.json"
       stamp --stage X6 --producer $XB/CITEgeist/src/eval_single_cell_assignment.py --role x6_single_cell_official \
         --artifact "$SC/eval_results_v2_rnaprot.json" --input-glob "$PRED7" baseline.csv.gz --input $XB/ground_truth_xenium_rnaprot/labels/cell_labels_7type.csv ;;
  b15) "$PY" -u $XB/spots_baseline/score_p138_strict.py --labelset rnagt --arm v2
       (cd $XB && XEN_LABELSET=rnaprot "$PY" -u -m gt_rebuild.sensitivity_table --with-v2 --rescore-v2)
       stamp --stage X2 --producer $XB/gt_rebuild/sensitivity_table.py --role b15_annotation_sensitivity \
         --artifact "$RES_RNAPROT/sensitivity_supplement_6type.json" --inputs-from-json "$RES_RNAPROT/sensitivity_supplement_6type.json" inputs_sha256 \
         --input-glob "$PRED6" baseline.csv.gz ;;
  x9)  "$PY" -u $RF/v2_official.py x9 "$X9_OUT" ;;
  *)   echo "unknown STAGE $STAGE" >&2; exit 2 ;;
esac
"$PY" -u $RF/v2_official.py complete "$STAGE"
echo "stage=$STAGE done $(date -Is)"
