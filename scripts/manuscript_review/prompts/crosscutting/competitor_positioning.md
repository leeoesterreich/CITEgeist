# Cross-Cutting Review: Competitor Positioning

You are reviewing **competitor positioning** in a manuscript being submitted to Cell Reports Methods / Patterns.

## Critical Context

Methods papers must position against the competition. Too harsh = reviewers defend the attacked methods. Too soft = reviewers ask "why do we need this?" The balance is critical.

Your job: **Ensure CITEgeist's advantages are clearly stated while treating competitors fairly.**

## Competitors Mentioned

1. **Cell2Location** - Reference-based deconvolution (Nat Biotechnol 2022)
2. **RCTD** - Reference-based deconvolution
3. **Tangram** - Deep learning spatial mapping (Nat Methods 2021)
4. **Seurat** - Label transfer approach (Nat Biotechnol 2018)

## Your Scope

Review the ENTIRE manuscript for competitor positioning:
- Introduction (why existing methods are insufficient)
- Results (benchmarking comparisons)
- Discussion (positioning and limitations)

## Evaluation Criteria

### 1. Fair Characterization
- Are competitor methods described accurately?
- Are their intended use cases acknowledged?
- Are their strengths mentioned, not just weaknesses?
- Would the developers of these methods feel fairly represented?

### 2. CITEgeist Edge Clarity
For each advantage, is it:
- Clearly stated?
- Supported by evidence?
- Differentiated from what competitors offer?

Advantages to check:
- Reference-free operation
- Spatially-native design
- End-to-end pipeline
- Same-slide proteomics anchoring

### 3. Apples-to-Apples Comparisons
- Is the benchmarking comparison appropriate?
- Are we comparing CITEgeist's strengths to competitors' weaknesses?
- Would a fairer comparison tell the same story?

### 4. Tangram/Seurat Treatment
The manuscript shows Tangram r=0.14, Seurat r=0.17 vs. CITEgeist r=0.60.
- Is this a fair comparison or are these methods being misused?
- Are these methods designed for this task?
- Will reviewers familiar with these methods object?

### 5. Limitation Acknowledgment
- Does CITEgeist acknowledge where competitors are better?
- Are there use cases where Cell2Location would be preferred?
- Is the "requires protein panel" limitation discussed?
- Is this honest self-assessment present?

### 6. Claim Hedging
- Are comparative claims appropriately hedged?
- "Competitive" vs "superior" vs "comparable"
- Is the language defensive enough without underselling?

## Section Agent Reviews

Read Phase 1 section reviews, especially:
- `section-2.3-review.md` (benchmarking section)
- Phase 1 reviews may flag competitor issues

## Output Format

Write your review to: `/workspace/reviews/competitor-positioning-review.md`

```markdown
# Cross-Cutting Review: Competitor Positioning

## Executive Summary
[2-3 sentences: Is CITEgeist's edge clear? Are competitors treated fairly?]

## Findings

### Fair Characterization
| Method | Treatment | Fair? | Concerns |
|--------|-----------|-------|----------|
| Cell2Location | [how described] | Yes/No | [issues] |
| RCTD | ... | ... | ... |
| Tangram | ... | ... | ... |
| Seurat | ... | ... | ... |

### CITEgeist Edge Clarity
| Advantage | Clearly Stated? | Evidence? | Location |
|-----------|-----------------|-----------|----------|
| Reference-free | Yes/No | Yes/No | Section X.X |
| Spatially-native | ... | ... | ... |
| End-to-end | ... | ... | ... |
| Same-slide anchoring | ... | ... | ... |

### Apples-to-Apples Comparisons
[Your analysis]

### Tangram/Seurat Treatment
[Specific assessment - will their users object?]

### Limitation Acknowledgment
[Does CITEgeist admit its weaknesses?]

### Claim Hedging
[Is language appropriately calibrated?]

## Recommendations

| Priority | Issue | Recommendation | Location |
|----------|-------|----------------|----------|
| CRITICAL | [issue] | [fix] | Section X.X |
| MAJOR | [issue] | [fix] | Section X.X |
| MINOR | [issue] | [fix] | Section X.X |
```

## Data Access

- `/data/manuscript/` - Full manuscript
- `/data/codebase/Benchmarking/` - How comparisons were run
- `/workspace/reviews/` - Phase 1 section reviews

Reviewers often ARE the developers of competing methods. Treat them fairly.
