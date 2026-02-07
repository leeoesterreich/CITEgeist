#!/bin/bash
# Monitor Module 1-2 discovery jobs and run comparison when complete

JOB_ID=7873527
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/output/module12_discovery"
EXAMPLES_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/examples"
EXPECTED_FILES=14

echo "Monitoring job ${JOB_ID} for Module 1-2 discovery..."
echo "Expected output files: ${EXPECTED_FILES}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Started: $(date)"
echo "==========================================="

while true; do
    # Check if any jobs are still running (exclude header line)
    RUNNING=$(squeue -u alc376 --cluster=htc 2>/dev/null | grep -c "${JOB_ID}")

    # Count completed JSON files
    COMPLETED=$(ls -1 ${OUTPUT_DIR}/*_module12_discovery.json 2>/dev/null | wc -l)

    echo "[$(date '+%H:%M:%S')] Jobs running: ${RUNNING}, Files completed: ${COMPLETED}/${EXPECTED_FILES}"

    if [ ${RUNNING} -eq 0 ]; then
        echo ""
        echo "==========================================="
        echo "All jobs finished at $(date)"
        echo "Completed files: ${COMPLETED}"
        echo ""

        # List any failed jobs by checking for missing files
        echo "Checking for missing outputs..."
        MISSING=0
        for sample in HCC22-088-P1-S1 HCC22-088-P1-S2 HCC22-088-P2-S1 HCC22-088-P2-S2 \
                      HCC22-088-P3-S1_A HCC22-088-P3-S2 HCC22-088-P4-S1 HCC22-088-P4-S2 \
                      HCC22-088-P4-S2_1i_rep HCC22-088-P5-S1 HCC22-088-P5-S2 \
                      HCC22-088-P5-S2_F_rep HCC22-088-P6-S1 HCC22-088-P6-S2_D; do
            if [ ! -f "${OUTPUT_DIR}/${sample}_module12_discovery.json" ]; then
                echo "  MISSING: ${sample}"
                MISSING=$((MISSING+1))
            fi
        done

        if [ ${MISSING} -eq 0 ]; then
            echo "All ${EXPECTED_FILES} files present!"
        fi

        if [ ${COMPLETED} -ge 1 ]; then
            echo ""
            echo "Running comparison analysis on ${COMPLETED} files..."
            echo "==========================================="

            # Activate conda and run comparison
            source /ihome/alee/alc376/miniconda3/etc/profile.d/conda.sh
            conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
            cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

            python examples/compare_profiles.py \
                --input-dir output/module12_discovery \
                --output-dir output/profile_comparison

            echo ""
            echo "==========================================="
            echo "Comparison complete!"
            echo "Results in: output/profile_comparison/"
            ls -la output/profile_comparison/
        else
            echo ""
            echo "WARNING: No files completed."
            echo "Check SLURM logs in examples/slurm_log/ for errors."
        fi

        break
    fi

    sleep 60
done

echo ""
echo "Monitoring finished at $(date)"
