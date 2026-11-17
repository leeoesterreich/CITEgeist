#!/bin/bash
# Launch all SACE grid search jobs: 1125 combos x 5 regions
# Resumable: skips (combo, region) pairs that already have a result JSON
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

PROPS_DIR="${1:?Usage: $0 <props-dir-for-best-config>}"
CHUNK_SIZE=375
N_COMBOS=1125
N_CHUNKS=$((N_COMBOS / CHUNK_SIZE))
# HTC QOS MaxSubmitJobsPerUser=1250; keep well below with 375-task chunks
MAX_QUEUED=800

mkdir -p results/sace/logs

echo "Submitting SACE grid search: ${N_COMBOS} combos x 5 regions"
echo "Props dir: ${PROPS_DIR}"
echo "${N_CHUNKS} chunks of ${CHUNK_SIZE}"
echo "Throttle: wait until queue < ${MAX_QUEUED} (gs_sace jobs only)"
echo ""

build_skip_list() {
    local REGION_ID=$1
    local OFFSET=$2
    local CHUNK_SIZE=$3
    local skip_ids=()
    for i in $(seq 0 $((CHUNK_SIZE - 1))); do
        local COMBO_ID=$((OFFSET + i))
        if [ -f "results/sace/${COMBO_ID}/region_${REGION_ID}.json" ]; then
            skip_ids+=("$i")
        fi
    done
    echo "${skip_ids[*]:-}"
}

submit_chunk() {
    local ARRAY_SPEC=$1
    local REGION_ID=$2
    local OFFSET=$3
    local MAX_RETRIES=30
    local attempt=0

    while [ "$attempt" -lt "$MAX_RETRIES" ]; do
        SBATCH_OUT=$(sbatch \
            --cluster=htc \
            --array="${ARRAY_SPEC}" \
            --export="REGION_ID=${REGION_ID},OFFSET=${OFFSET},PROPS_DIR=${PROPS_DIR}" \
            sbatch_sace_gridsearch.sbatch 2>&1) || true
        JOB_ID=$(echo "$SBATCH_OUT" | awk '/Submitted/{print $4}')
        if [ -n "$JOB_ID" ]; then
            echo "$JOB_ID"
            return 0
        fi
        attempt=$((attempt + 1))
        echo "    [submit attempt ${attempt}/${MAX_RETRIES}] $(echo "$SBATCH_OUT" | grep -i error | tail -1)" >&2
        # Wait for queue to drain — use -r to expand array elements (MaxSubmit counts elements, not job lines)
        while true; do
            QUEUED=$(squeue -u "$USER" --cluster=htc -h -n gs_sace -r 2>/dev/null | wc -l)
            if [ "$QUEUED" -lt "$MAX_QUEUED" ]; then
                break
            fi
            echo "    [throttle] ${QUEUED} gs_sace elements queued (expanded), waiting for < ${MAX_QUEUED}..." >&2
            sleep 120
        done
        sleep 10
    done
    return 1
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
            TOTAL_SKIPPED=$((TOTAL_SKIPPED + N_SKIP))
            continue
        fi

        # Throttle: wait until queue is below threshold BEFORE submitting (-r expands array elements)
        while true; do
            QUEUED=$(squeue -u "$USER" --cluster=htc -h -n gs_sace -r 2>/dev/null | wc -l)
            if [ "$QUEUED" -lt "$MAX_QUEUED" ]; then
                break
            fi
            echo "    [throttle] ${QUEUED} gs_sace elements queued (expanded), waiting for < ${MAX_QUEUED}..."
            sleep 120
        done

        # Build array spec excluding completed tasks
        if [ "$N_SKIP" -gt 0 ]; then
            ARRAY_SPEC=""
            for i in $(seq 0 "$MAX_IDX"); do
                if ! echo " $SKIP " | grep -q " $i "; then
                    if [ -n "$ARRAY_SPEC" ]; then
                        ARRAY_SPEC="${ARRAY_SPEC},$i"
                    else
                        ARRAY_SPEC="$i"
                    fi
                fi
            done
        else
            ARRAY_SPEC="0-${MAX_IDX}"
        fi

        JOB_ID=$(submit_chunk "$ARRAY_SPEC" "$REGION_ID" "$OFFSET")
        if [ -z "$JOB_ID" ]; then
            echo "  region=${REGION_ID} chunk=${CHUNK} offset=${OFFSET} -> ALL RETRIES FAILED, skipping"
            continue
        fi

        echo "  region=${REGION_ID} chunk=${CHUNK} offset=${OFFSET} -> job ${JOB_ID} (${N_REMAINING} tasks, ${N_SKIP} skipped)"
        TOTAL_SUBMITTED=$((TOTAL_SUBMITTED + N_REMAINING))
        TOTAL_SKIPPED=$((TOTAL_SKIPPED + N_SKIP))
    done
done

echo ""
echo "Submitted ${TOTAL_SUBMITTED} tasks, skipped ${TOTAL_SKIPPED} already-done"
