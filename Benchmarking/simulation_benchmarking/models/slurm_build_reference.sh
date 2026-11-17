#!/bin/bash
#SBATCH --job-name=build_sim_reference
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --output=Benchmarking/simulation_benchmarking/models/build_reference_%j.out
#SBATCH --error=Benchmarking/simulation_benchmarking/models/build_reference_%j.err

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
set -eo pipefail

eval "$(conda shell.bash hook 2>/dev/null)"

# Stage 1: Build h5ad (Python)
echo "=== Stage 1: build_reference.py ==="
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env
set -u

cd ${CITEGEIST_REPO_ROOT}

python -u Benchmarking/simulation_benchmarking/models/build_reference.py

# Stage 2: Convert formats (R)
echo "=== Stage 2: convert_reference_formats.R ==="
set +u
conda deactivate
module load gcc/12.2.0 r/4.4.0
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env
set -u

Rscript Benchmarking/simulation_benchmarking/models/convert_reference_formats.R \
    Benchmarking/simulation_benchmarking/models/Wu_scRNA_ref_ERpos.h5ad

echo "=== All stages complete ==="
ls -lh Benchmarking/simulation_benchmarking/models/Wu_scRNA_ref_ERpos.*
