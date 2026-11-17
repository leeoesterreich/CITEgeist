#!/bin/bash
#SBATCH --job-name=mil_3ch_s4
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/mil_3ch_s4_%j.out
#SBATCH --error=slurm/logs/mil_3ch_s4_%j.err

# Stage 4 ONLY: MIL benchmark using pre-trained SimCLR from regen run
# Uses GT proportions as prior (deconv results deleted in disk cleanup)

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

export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v "ondemand-jupyter" | tr '\n' ':' | sed 's/:$//')

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

BENCH="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
OUT_ROOT="${BENCH}/output/regen_3ch_fair"

cd "${BENCH}"

echo "============================================"
echo "Stage 4: MIL Benchmark (GT prior)"
echo "Started: $(date)"
echo "============================================"

for REGION in 1 2 4; do
    echo ""
    echo "Region ${REGION} (2ch)..."
    python src/benchmark_mil_assignment.py \
        --region ${REGION} \
        --simclr-checkpoint "${OUT_ROOT}/simclr_2ch/simclr_best.pt" \
        --patches-dir "${OUT_ROOT}/patches_2ch" \
        --output-dir "${OUT_ROOT}/mil_2ch" \
        --in-channels 2 \
        --use-gt-prior \
        --lambda-prior 1.0 --n-epochs 100 --device cuda --batch-size 64

    echo "Region ${REGION} (3ch)..."
    python src/benchmark_mil_assignment.py \
        --region ${REGION} \
        --simclr-checkpoint "${OUT_ROOT}/simclr_3ch/simclr_best.pt" \
        --patches-dir "${OUT_ROOT}/patches_3ch" \
        --output-dir "${OUT_ROOT}/mil_3ch" \
        --in-channels 3 \
        --use-gt-prior \
        --lambda-prior 1.0 --n-epochs 100 --device cuda --batch-size 64

    echo "Region ${REGION} complete: $(date)"
done

echo ""
echo "============================================"
echo "=== FAIR 3ch vs 2ch RESULTS (same cells, GT prior) ==="
python3 -c "
import json, os, statistics
regions = [1, 2, 4]
r2s, r3s = [], []
for r in regions:
    f2 = '${OUT_ROOT}' + f'/mil_2ch/region_{r}/metrics.json'
    f3 = '${OUT_ROOT}' + f'/mil_3ch/region_{r}/metrics.json'
    d2 = json.load(open(f2))
    d3 = json.load(open(f3))
    v2 = d2.get('proportion_pearson_r', 0)
    v3 = d3.get('proportion_pearson_r', 0)
    r2s.append(v2); r3s.append(v3)
    sc2 = d2.get('sc_accuracy', 0)
    sc3 = d3.get('sc_accuracy', 0)
    print(f'Region {r}: 2ch r={v2:.4f} sc={sc2:.1%}  |  3ch r={v3:.4f} sc={sc3:.1%}  |  delta={v3-v2:+.4f}')
m2 = statistics.mean(r2s)
m3 = statistics.mean(r3s)
print(f'Mean:     2ch r={m2:.4f}              |  3ch r={m3:.4f}              |  delta={m3-m2:+.4f}')
print()
if m3 > m2 + 0.02:
    print('VERDICT: 3ch IMPROVES over 2ch (adopt 3ch)')
elif m3 < m2 - 0.02:
    print('VERDICT: 3ch HURTS vs 2ch (stay with 2ch)')
else:
    print('VERDICT: No significant difference (stay with 2ch for simplicity)')
"
echo "============================================"
echo "All complete: $(date)"
