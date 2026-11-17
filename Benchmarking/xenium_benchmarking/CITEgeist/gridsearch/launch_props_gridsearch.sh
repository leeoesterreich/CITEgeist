#!/bin/bash
# Launch all proportions grid search jobs: 5 regions x 5 chunks = 25 submissions
# Throttled: waits for queue to drain below MAX_QUEUED before submitting next chunk
# Resumable: skips (combo, region) pairs that already have a result JSON
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

PARTITION="${1:-l40s}"
CHUNK_SIZE=500
N_COMBOS=2500
N_CHUNKS=$((N_COMBOS / CHUNK_SIZE))
MAX_QUEUED=400  # stay under GPU cluster MaxSubmitJobsPerAccount (~500)
RESULTS_DIR="${SCRIPT_DIR}/results/props"

mkdir -p results/props/logs

echo "Submitting props grid search: ${N_COMBOS} combos x 5 regions"
echo "Partition: ${PARTITION}, ${N_CHUNKS} chunks of ${CHUNK_SIZE}"
echo "Throttle: wait until queue < ${MAX_QUEUED} (gs_props jobs only)"
echo ""

wait_for_slots() {
    while true; do
        QUEUED=$(squeue -h -r -u "$USER" --cluster=gpu -n gs_props 2>/dev/null | wc -l)
        if [ "$QUEUED" -lt "$MAX_QUEUED" ]; then
            return
        fi
        echo "    [throttle] ${QUEUED} gs_props jobs queued, waiting for < ${MAX_QUEUED}..."
        sleep 30
    done
}

build_skip_list() {
    local region_id=$1
    local offset=$2
    local chunk_size=$3
    local skip_ids=()
    for i in $(seq 0 $((chunk_size - 1))); do
        local combo_id=$((offset + i))
        if [ -f "${RESULTS_DIR}/${combo_id}/region_${region_id}.json" ]; then
            skip_ids+=("$i")
        fi
    done
    echo "${skip_ids[*]:-}"
}

TOTAL_SUBMITTED=0
TOTAL_SKIPPED=0
for REGION_ID in 0 1 2 3 4; do
    for CHUNK in $(seq 0 $((N_CHUNKS - 1))); do
        OFFSET=$((CHUNK * CHUNK_SIZE))
        MAX_IDX=$((CHUNK_SIZE - 1))

        SKIP=$(build_skip_list "$REGION_ID" "$OFFSET" "$CHUNK_SIZE")
        N_SKIP=$(echo "$SKIP" | wc -w)
        N_REMAINING=$((CHUNK_SIZE - N_SKIP))

        if [ "$N_REMAINING" -eq 0 ]; then
            echo "  region=${REGION_ID} chunk=${CHUNK} offset=${OFFSET} -> ALL DONE (skipped ${N_SKIP})"
            TOTAL_SKIPPED=$((TOTAL_SKIPPED + CHUNK_SIZE))
            continue
        fi

        # Build array spec excluding completed tasks
        if [ "$N_SKIP" -gt 0 ]; then
            # Build comma-separated list of remaining task IDs
            ARRAY_SPEC=""
            for i in $(seq 0 "$MAX_IDX"); do
                if ! echo " $SKIP " | grep -q " $i "; then
                    if [ -z "$ARRAY_SPEC" ]; then
                        ARRAY_SPEC="$i"
                    else
                        ARRAY_SPEC="${ARRAY_SPEC},$i"
                    fi
                fi
            done
        else
            ARRAY_SPEC="0-${MAX_IDX}"
        fi

        wait_for_slots

        JOB_ID=$(sbatch \
            --cluster=gpu \
            --partition="${PARTITION}" \
            --array="${ARRAY_SPEC}" \
            --export="REGION_ID=${REGION_ID},OFFSET=${OFFSET}" \
            sbatch_props_gridsearch.sbatch \
            | awk '{print $4}')

        echo "  region=${REGION_ID} chunk=${CHUNK} offset=${OFFSET} -> job ${JOB_ID} (${N_REMAINING} tasks, ${N_SKIP} skipped)"
        TOTAL_SUBMITTED=$((TOTAL_SUBMITTED + N_REMAINING))
        TOTAL_SKIPPED=$((TOTAL_SKIPPED + N_SKIP))
    done
done

echo ""
echo "Submitted ${TOTAL_SUBMITTED} tasks, skipped ${TOTAL_SKIPPED} already-done"
