# Cross-Cutting Review: Benchmarking Rigor

You are reviewing the **benchmarking rigor** of a manuscript being submitted to Cell Reports Methods / Patterns.

## Critical Context

CITEgeist has been desk-rejected before for "insufficient impact." Methods journals live and die by benchmarking quality. If reviewers perceive unfair comparisons, cherry-picked metrics, or poorly tuned baselines, the paper is dead.

Your job: **Find every potential benchmarking weakness before reviewers do.**

## Your Scope

Review the ENTIRE manuscript for benchmarking issues, with special focus on:
- Section 2.3 (Module 3 deconvolution benchmarking)
- Figure 3 (main benchmarking results)
- Any comparative claims throughout

## Evaluation Criteria

### 1. Fair Comparisons
- Were all methods given equal opportunity to succeed?
- Same preprocessing applied to all methods?
- Same input data (gene sets, cell types, spatial resolution)?
- Were method-specific optimizations applied fairly?

### 2. Baseline Tuning
- Were competitor methods (Cell2Location, RCTD, Tangram, Seurat) run with optimal parameters?
- Check `/data/outputs/` for parameter sweep results
- Were default parameters used when optimized parameters would be fairer?
- Is there documentation of parameter choices?

### 3. Appropriate Metrics
- Are the metrics appropriate for the claims being made?
- Pearson r, JSD, RMSE, MAE - are these the right choices?
- Are there other metrics that would tell a different story?
- Is statistical significance properly assessed?

### 4. Ground Truth Validity
- Is the Xenium pseudo-Visium approach a valid ground truth?
- Are there limitations to aggregating single-cell data into spots?
- Is this acknowledged in the manuscript?

### 5. Cherry-Picking Concerns
- Are results shown for all regions, or just favorable ones?
- Check `/data/outputs/` for per-region results that weren't shown
- Is there variability across regions that's being hidden?
- Would a reviewer suspect selective reporting?

### 6. Competitor Positioning Fairness
- Are Tangram and Seurat being treated fairly?
- The manuscript shows r=0.14 for Tangram, r=0.17 for Seurat - is this legitimate or are they being set up to fail?
- Were these methods used appropriately for their intended purpose?

## Section Agent Reviews

Read the Phase 1 section reviews in `/workspace/reviews/` for additional context, especially `section-2.3-review.md`.

## Output Format

Write your review to: `/workspace/reviews/benchmarking-rigor-review.md`

```markdown
# Cross-Cutting Review: Benchmarking Rigor

## Executive Summary
[2-3 sentences: Is the benchmarking defensible? What are the main risks?]

## Findings

### Fair Comparisons
[Your analysis]

### Baseline Tuning
[Your analysis - check for parameter documentation]

### Appropriate Metrics
[Your analysis]

### Ground Truth Validity
[Your analysis]

### Cherry-Picking Concerns
[Your analysis - check for hidden variability]

### Competitor Positioning Fairness
[Your analysis - especially Tangram/Seurat treatment]

## Recommendations

| Priority | Issue | Recommendation | Evidence Location |
|----------|-------|----------------|-------------------|
| CRITICAL | [issue] | [fix] | [file path] |
| MAJOR | [issue] | [fix] | [file path] |
| MINOR | [issue] | [fix] | [file path] |
```

## Data Access

- `/data/manuscript/` - Manuscript and figures
- `/data/outputs/` - Benchmarking results, per-region breakdowns
- `/data/codebase/Benchmarking/` - Benchmark implementation details
- `/workspace/reviews/` - Phase 1 section reviews

A hostile reviewer will attack the benchmarking first. Make it bulletproof.
