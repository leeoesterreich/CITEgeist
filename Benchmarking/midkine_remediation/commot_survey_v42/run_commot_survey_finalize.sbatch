#!/bin/bash
#SBATCH --job-name=cmsurvey_b26_fin
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
# Dependent job (afterok on the array): asserts all 24 task outputs + stamps, then the check 19 ranking, final stamps, _COMPLETE.json.
# Submitted by .commands/submit-commot-survey-stage.sh with env SHA, CODE_DIR (clean detached worktree), MAIN. Needs pandas only (no commot).
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
MAIN=${CITEGEIST_REPO_ROOT}
cd "$MAIN"
export MPLCONFIGDIR=/tmp/mpl_${SLURM_JOB_ID}
PY="$(cg_env_python CITEgeist_env)"
: "${CODE_DIR:?}" "${SHA:?}"
export V2_NAMESPACE_ROOT=$MAIN PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$CODE_DIR:$CODE_DIR/scripts
cd "$CODE_DIR"
[ "$(git rev-parse HEAD)" = "$SHA" ] && [ -z "$(git status --porcelain -uno)" ] || { echo "worktree is not clean at $SHA" >&2; exit 1; }
echo "stage=B26 finalize sha=$SHA host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/commot_survey_v42/commot_survey.py --finalize
echo "done end=$(date -Is)"
