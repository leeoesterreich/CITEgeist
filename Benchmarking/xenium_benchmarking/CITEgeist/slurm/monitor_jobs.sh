#!/bin/bash
# Monitor benchmark jobs and report completion status

ACHIEVABLE_JOB=7487089
AUTODISCOVERY_JOB=7487094

echo "Monitoring jobs: Achievable_7 ($ACHIEVABLE_JOB), Autodiscovery ($AUTODISCOVERY_JOB)"
echo "Started at: $(date)"
echo ""

while true; do
    # Check job status
    ACHIEVABLE_STATUS=$(squeue -j $ACHIEVABLE_JOB --cluster=htc -h 2>/dev/null | wc -l)
    AUTODISCOVERY_STATUS=$(squeue -j $AUTODISCOVERY_JOB --cluster=htc -h 2>/dev/null | wc -l)
    
    echo "[$(date '+%H:%M:%S')] Achievable_7: $ACHIEVABLE_STATUS tasks, Autodiscovery: $AUTODISCOVERY_STATUS tasks"
    
    # Check if both complete
    if [ "$ACHIEVABLE_STATUS" -eq 0 ] && [ "$AUTODISCOVERY_STATUS" -eq 0 ]; then
        echo ""
        echo "=========================================="
        echo "BOTH JOBS COMPLETE at $(date)"
        echo "=========================================="
        
        # Check for errors
        echo ""
        echo "Checking for errors..."
        ACHIEVABLE_ERRORS=$(grep -l "Error\|FAILED\|Traceback" /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/xenium_achievable_7_*.err 2>/dev/null | wc -l)
        AUTODISCOVERY_ERRORS=$(grep -l "Error\|FAILED\|Traceback" /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/autodiscovery_bench_*.err 2>/dev/null | wc -l)
        
        echo "  Achievable_7 errors: $ACHIEVABLE_ERRORS"
        echo "  Autodiscovery errors: $AUTODISCOVERY_ERRORS"
        
        # Check output files
        echo ""
        echo "Checking output files..."
        ACHIEVABLE_OUTPUTS=$(ls /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7/Xenium_region_*_gene_expression_pass1.parquet 2>/dev/null | wc -l)
        AUTODISCOVERY_OUTPUTS=$(ls /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_autodiscovery/Xenium_region_*_gene_expression_pass1.parquet 2>/dev/null | wc -l)
        
        echo "  Achievable_7 GEX outputs: $ACHIEVABLE_OUTPUTS/5"
        echo "  Autodiscovery GEX outputs: $AUTODISCOVERY_OUTPUTS/5"
        
        if [ "$ACHIEVABLE_OUTPUTS" -eq 5 ] && [ "$AUTODISCOVERY_OUTPUTS" -eq 5 ]; then
            echo ""
            echo "SUCCESS: All outputs present!"
            echo ""
            echo "Ready to run Module 4/5 validation:"
            echo "  cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm"
            echo "  sbatch run_module4_validation.sh"
        else
            echo ""
            echo "WARNING: Some outputs may be missing"
        fi
        
        exit 0
    fi
    
    sleep 300  # Check every 5 minutes
done
