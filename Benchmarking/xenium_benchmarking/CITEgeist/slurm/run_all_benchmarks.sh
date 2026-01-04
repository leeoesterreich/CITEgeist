#!/bin/bash
# ============================================================================
# Master Script: Run All CITEgeist Benchmarking Modes
# ============================================================================
# This script submits all benchmarking jobs in the correct order:
#
# Phase 1: Profile Discovery & Ground Truth Generation (in parallel)
#   - Mode A: Profile Discovery Evaluation
#   - Mode C: Manual Profiles Benchmark (granular 10-cell-type)
#
# Phase 2: Auto-Discovery Benchmark (depends on Phase 1)
#   - Mode B: Full Auto-Discovery Pipeline
#
# Phase 3: Comparison (depends on Phase 2)
#   - Generate comparison figures
#
# Usage:
#   ./run_all_benchmarks.sh
# ============================================================================

GREEN="\033[1;32m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RESET="\033[0m"

SLURM_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm"

echo -e "${BLUE}======================================${RESET}"
echo -e "${BLUE}CITEgeist Benchmarking Suite${RESET}"
echo -e "${BLUE}======================================${RESET}"

# Create slurm_log directory
mkdir -p "${SLURM_DIR}/slurm_log"

# ============================================================================
# Phase 1: Profile Discovery & Manual Profiles (parallel)
# ============================================================================
echo -e "\n${YELLOW}Phase 1: Submitting Profile Discovery and Manual Profiles benchmarks...${RESET}"

# Mode A: Profile Discovery Evaluation
JOB_A=$(sbatch "${SLURM_DIR}/evaluate_profile_discovery.sh" | awk '{print $4}')
echo -e "${GREEN}Mode A (Profile Discovery): Job ${JOB_A}${RESET}"

# Mode C: Manual Profiles (granular)
JOB_C=$(sbatch "${SLURM_DIR}/xenium_benchmark_granular.sh" | awk '{print $4}')
echo -e "${GREEN}Mode C (Manual Profiles): Job ${JOB_C}${RESET}"

# ============================================================================
# Phase 2: Auto-Discovery (depends on Phase 1)
# ============================================================================
echo -e "\n${YELLOW}Phase 2: Submitting Auto-Discovery benchmark (depends on Phase 1)...${RESET}"

JOB_B=$(sbatch --dependency=afterok:${JOB_A}:${JOB_C} "${SLURM_DIR}/run_autodiscovery_benchmark.sh" | awk '{print $4}')
echo -e "${GREEN}Mode B (Auto-Discovery): Job ${JOB_B} (depends on ${JOB_A}, ${JOB_C})${RESET}"

# ============================================================================
# Phase 3: Comparison (depends on Phase 2)
# ============================================================================
echo -e "\n${YELLOW}Phase 3: Submitting comparison figures job (depends on Phase 2)...${RESET}"

JOB_COMPARE=$(sbatch --dependency=afterok:${JOB_B}:${JOB_C} "${SLURM_DIR}/compare_benchmarks.sh" | awk '{print $4}')
echo -e "${GREEN}Comparison: Job ${JOB_COMPARE} (depends on ${JOB_B}, ${JOB_C})${RESET}"

# ============================================================================
# Summary
# ============================================================================
echo -e "\n${BLUE}======================================${RESET}"
echo -e "${BLUE}All Jobs Submitted${RESET}"
echo -e "${BLUE}======================================${RESET}"
echo -e "Mode A (Profile Discovery):  ${JOB_A}"
echo -e "Mode C (Manual Profiles):    ${JOB_C}"
echo -e "Mode B (Auto-Discovery):     ${JOB_B}"
echo -e "Comparison:                  ${JOB_COMPARE}"
echo -e ""
echo -e "Monitor with: squeue -u \$USER --cluster=htc"
echo -e "View logs in: ${SLURM_DIR}/slurm_log/"
echo -e "${BLUE}======================================${RESET}"
