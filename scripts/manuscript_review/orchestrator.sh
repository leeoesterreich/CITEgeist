#!/bin/bash
# Manuscript Review Agent Orchestrator
#
# This script prepares the review workspace and documents launch commands.
# Due to hpc-isolate architecture, agents are launched manually in phases.
#
# Usage:
#   1. Run this script to set up workspace
#   2. Launch Phase 1 agents (6 section reviewers)
#   3. Wait for completion, then launch Phase 2 (4 cross-cutting)
#   4. Wait for completion, then launch Phase 3 (synthesis)

set -e

# Configuration
WORKSPACE="${HOME}/alc376_bgfs/manuscript_review_$(date +%Y-%m-%d)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CITEGEIST_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Data paths for --data mounts
MANUSCRIPT_PATH="${CITEGEIST_ROOT}/manuscript"
BENCHMARKING_PATH="${CITEGEIST_ROOT}/Benchmarking"
CODEBASE_PATH="${CITEGEIST_ROOT}/CITEgeist"

echo "=== Manuscript Review Agent Orchestrator ==="
echo ""
echo "Workspace: ${WORKSPACE}"
echo "Manuscript: ${MANUSCRIPT_PATH}"
echo "Benchmarking: ${BENCHMARKING_PATH}"
echo "Codebase: ${CODEBASE_PATH}"
echo ""

# Create workspace
mkdir -p "${WORKSPACE}/reviews"
mkdir -p "${WORKSPACE}/prompts"

# Copy prompts to workspace
cp -r "${SCRIPT_DIR}/prompts/"* "${WORKSPACE}/prompts/"
cp -r "${SCRIPT_DIR}/sections/"* "${WORKSPACE}/prompts/"

echo "Workspace created at: ${WORKSPACE}"
echo ""

# Generate filled section prompts
echo "=== Generating Section Agent Prompts ==="

for section in 2.1 2.2 2.3 2.4 2.5 2.6; do
    section_file="${WORKSPACE}/prompts/section_${section}.md"
    template="${WORKSPACE}/prompts/section_template.md"
    section_text="${WORKSPACE}/prompts/section_${section}.md"

    # Extract section title
    case $section in
        2.1) title="Framework Overview" ;;
        2.2) title="Modules 1-2: Automated Profile Discovery" ;;
        2.3) title="Module 3: Protein-Anchored Deconvolution" ;;
        2.4) title="Module 4: Spatially Coherent Program Discovery" ;;
        2.5) title="Module 5: Cross-Sample Integration" ;;
        2.6) title="Interpretable Outputs Enable Downstream Analysis" ;;
    esac

    output_file="section-${section}-review.md"

    # Create filled prompt
    filled_prompt="${WORKSPACE}/prompts/section_${section}_filled.md"

    # Read template and section text
    section_content=$(cat "${WORKSPACE}/prompts/section_${section}.md")

    cat > "${filled_prompt}" << PROMPT_EOF
# Section Review Agent Prompt

You are reviewing **Section ${section}: ${title}** of a manuscript being submitted to Cell Reports Methods / Patterns (Cell Press).

## Critical Context

CITEgeist has been **desk-rejected before** for "insufficient impact" when only showing Module 3. The full 5-module pipeline is now complete. **Your job is to ensure this section makes maximum impact.**

This is not a gentle review. Be ruthlessly critical. Find every weakness. Identify every missed opportunity. The goal is to make this manuscript rejection-proof.

## CITEgeist's Key Advantages (these MUST be hammered home)

1. **Reference-free**: Uses same-slide protein measurements, no external scRNA-seq atlas required
2. **Spatially-native**: Spatial statistics integrated at every module (Moran's I, Laplacian regularization, bivariate spatial analysis), not bolted on as an afterthought
3. **End-to-end pipeline**: Raw spatial multi-omics data → cross-patient conserved programs in one unified framework
4. **Competitive accuracy**: Matches Cell2Location/RCTD benchmarks WITHOUT requiring their reference data dependencies

## Your Section

${section_content}

## Your Task

Evaluate this section against **6 criteria**. For each, be specific and cite file paths where relevant.

### 1. Referenced Figures Assessment
- Which figures does this section reference?
- Do the figures effectively support the claims made in the text?
- Are the figures publication-quality?
- Could better visualizations be generated from available outputs in /data/outputs/?
- Are there panels that should be added or removed?
- Does every figure panel get referenced in the text?

### 2. Most Impressive Result Available?
- Given all available outputs (not just what's in the paper), is this section showing the strongest evidence?
- Search /data/Benchmarking/ and /data/manuscript/figures/ for potentially stronger results
- If better evidence exists, specify the file path and what it shows

### 3. Scientific Edge Communicated?
- Does this section explicitly state what CITEgeist does that competitors cannot?
- Is the reference-free advantage mentioned?
- Is the spatially-native design highlighted?
- Could the advantage be stated more directly or forcefully?

### 4. Unused Capabilities or Results?
- Review /data/codebase/ for capabilities not demonstrated
- Check /data/Benchmarking/ for analyses that were run but not included
- What work has been done that isn't being shown?

### 5. Quantitative Backing?
- Are all claims supported by specific numbers?
- Flag any vague statements like "competitive", "comparable", "similar"
- Every performance claim should have metrics attached

### 6. Anticipated Reviewer Objections?
- What would a skeptical reviewer attack?
- What are the weak points in the argument?
- Is the section proactively addressing likely criticisms?

## Output Format

Write your review to: /workspace/reviews/${output_file}

Use the structure specified in the template, with:
- Executive Summary
- Findings (6 sections)
- Recommendations table (CRITICAL/MAJOR/MINOR)

## Data Access

You have read-only access to:
- /data/manuscript/ - Manuscript text and figures
- /data/Benchmarking/ - All benchmarking outputs
- /data/codebase/ - Full CITEgeist source code

Be thorough. Check everything. This manuscript cannot afford another rejection.
PROMPT_EOF

    echo "Created: ${filled_prompt}"
done

echo ""
echo "=== Launch Commands ==="
echo ""
echo "Phase 1: Section Agents (launch all 6 in parallel)"
echo "---------------------------------------------------"

for section in 2.1 2.2 2.3 2.4 2.5 2.6; do
    echo ""
    echo "# Section ${section}"
    echo "/hpc-isolate:launch ${WORKSPACE} \\"
    echo "    --data ${MANUSCRIPT_PATH} \\"
    echo "    --data ${BENCHMARKING_PATH} \\"
    echo "    --data ${CODEBASE_PATH} \\"
    echo "    --cpus 4 --mem 8G --time 04:00:00 \\"
    echo "    --prompt \"Read prompts/section_${section}_filled.md and execute the review. Write output to reviews/section-${section}-review.md\""
done

echo ""
echo ""
echo "Phase 2: Cross-Cutting Agents (launch after Phase 1 completes)"
echo "--------------------------------------------------------------"

for agent in benchmarking_rigor figure_text_alignment narrative_flow competitor_positioning; do
    output_file="${agent//_/-}-review.md"
    echo ""
    echo "# ${agent}"
    echo "/hpc-isolate:launch ${WORKSPACE} \\"
    echo "    --data ${MANUSCRIPT_PATH} \\"
    echo "    --data ${BENCHMARKING_PATH} \\"
    echo "    --data ${CODEBASE_PATH} \\"
    echo "    --cpus 4 --mem 8G --time 04:00:00 \\"
    echo "    --prompt \"Read prompts/crosscutting/${agent}.md and execute the review. Also read reviews/ for Phase 1 results. Write output to reviews/${output_file}\""
done

echo ""
echo ""
echo "Phase 3: Synthesis Agent (launch after Phase 2 completes)"
echo "----------------------------------------------------------"
echo ""
echo "/hpc-isolate:launch ${WORKSPACE} \\"
echo "    --data ${MANUSCRIPT_PATH} \\"
echo "    --cpus 4 --mem 8G --time 02:00:00 \\"
echo "    --prompt \"Read prompts/synthesis.md and all files in reviews/. Produce CONSOLIDATED_ISSUES.md\""
echo ""
echo "=== Setup Complete ==="
echo ""
echo "Workspace ready at: ${WORKSPACE}"
echo "Next: Launch Phase 1 agents using commands above"
