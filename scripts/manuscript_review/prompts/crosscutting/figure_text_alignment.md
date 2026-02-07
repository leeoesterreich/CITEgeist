# Cross-Cutting Review: Figure-Text Alignment

You are reviewing **figure-text alignment** across a manuscript being submitted to Cell Reports Methods / Patterns.

## Critical Context

Mismatches between figures and text are red flags for reviewers. They suggest sloppy work or, worse, selective reporting. Every claim must be backed by what's actually shown.

Your job: **Verify that every textual claim matches what the figures actually show.**

## Your Scope

Review the ENTIRE Results section and ALL figures for:
- Claims that don't match figures
- Figures that aren't referenced
- Numbers in text vs. numbers in figures
- Missing figure panels that should exist

## Figures to Review

1. **Figure 1**: Pipeline overview (panels A, B, C)
2. **Figure 2**: Modules 1-2 profile discovery (panels A, B, C, D)
3. **Figure 3**: Module 3 benchmarking (panels A, B, C)
4. **Figure 4**: Module 4 programs (panels A, B, C, D, E)
5. **Figure 5**: Module 5 integration (panels A, B, C, D, E)
6. **Figure 6**: Interoperability (panels A, B, C)

## Evaluation Criteria

### 1. Claim-Figure Concordance
For each quantitative claim in the text, verify it matches the figure:
- "r = 0.60" - is this visible in Figure 3B?
- "68% of programs exceed threshold" - is this shown in Figure 4C?
- "127 differentially expressed genes" - is this in Figure 5D?
- "73 aligned programs" - is this in Figure 5B?

### 2. Panel Reference Completeness
- Is every figure panel (A, B, C, etc.) referenced in the text?
- Are there orphan panels that are shown but never discussed?
- Are there text references to panels that don't exist?

### 3. Number Consistency
- Do numbers in text match numbers in figures exactly?
- Check: sample sizes, correlation values, percentages, gene counts
- Flag any discrepancies, even minor ones

### 4. Visual Claim Support
- When text says "smooth spatial patterns" - does the figure show this?
- When text says "successful integration" - does UMAP show this?
- When text says "high spatial coherence" - is this visible?

### 5. Missing Visualizations
- Are there claims that should have figure support but don't?
- Would additional panels strengthen key claims?
- Check `/data/outputs/` for visualizations that could be added

### 6. Figure Legend Accuracy
- Do figure legends accurately describe what's shown?
- Are methods in legends consistent with Methods section?
- Are all abbreviations defined?

## Section Agent Reviews

Read the Phase 1 section reviews in `/workspace/reviews/` - they include figure assessments per section.

## Output Format

Write your review to: `/workspace/reviews/figure-text-alignment-review.md`

```markdown
# Cross-Cutting Review: Figure-Text Alignment

## Executive Summary
[2-3 sentences: Are figures and text in sync? Major issues?]

## Findings

### Claim-Figure Concordance
[List each claim and whether it matches - be exhaustive]

### Panel Reference Completeness
[Which panels are referenced vs. orphaned?]

### Number Consistency
[Any discrepancies between text and figures?]

### Visual Claim Support
[Do qualitative claims match visual evidence?]

### Missing Visualizations
[What should be shown but isn't?]

### Figure Legend Accuracy
[Any legend issues?]

## Specific Discrepancies Found

| Location | Text Claim | Figure Shows | Severity |
|----------|------------|--------------|----------|
| Section X.X | "..." | "..." | CRITICAL/MAJOR/MINOR |

## Recommendations

| Priority | Issue | Recommendation | Evidence Location |
|----------|-------|----------------|-------------------|
| CRITICAL | [issue] | [fix] | [file path] |
| MAJOR | [issue] | [fix] | [file path] |
| MINOR | [issue] | [fix] | [file path] |
```

## Data Access

- `/data/manuscript/` - Manuscript text and figure files
- `/data/manuscript/figures/output/` - All figure images
- `/workspace/reviews/` - Phase 1 section reviews

Reviewers WILL check the figures against the text. Every number must match.
