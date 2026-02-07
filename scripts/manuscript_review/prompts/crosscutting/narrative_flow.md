# Cross-Cutting Review: Narrative Flow

You are reviewing the **narrative flow** of a manuscript being submitted to Cell Reports Methods / Patterns.

## Critical Context

A methods paper must tell a coherent story: here's the problem, here's our solution, here's why it works, here's what you can do with it. If readers get lost between sections, they stop reading.

Your job: **Ensure the manuscript tells a compelling, logical story from start to finish.**

## Your Scope

Review the ENTIRE Results section for narrative coherence:
- Section 2.1 → 2.2 → 2.3 → 2.4 → 2.5 → 2.6
- Does each section build on the previous?
- Is the "so what" clear at every step?
- Does it build toward clinical relevance?

## Evaluation Criteria

### 1. Section Transitions
- Does 2.1 set up why we need 2.2?
- Does 2.2 deliver what 2.3 needs?
- Does 2.3 enable what 2.4 does?
- Does 2.4 feed into 2.5?
- Does 2.5 demonstrate the clinical payoff?
- Does 2.6 show practical utility?

Rate each transition: Smooth / Adequate / Abrupt / Missing

### 2. Problem-Solution Clarity
- Is the problem each module solves clearly stated?
- Is it clear why this solution is needed?
- Is the connection to the overall CITEgeist value proposition maintained?

### 3. "So What" Presence
For each section, answer: Why should the reader care?
- 2.1: Why does this framework matter?
- 2.2: Why is automated profile discovery valuable?
- 2.3: Why is reference-free deconvolution important?
- 2.4: Why do spatial programs matter?
- 2.5: Why is cross-sample integration needed?
- 2.6: Why does interoperability matter?

Are these "so whats" explicit in the text or implied?

### 4. Escalating Impact
Does the narrative build toward maximum impact?
- Early sections: technical capability
- Middle sections: validation and benchmarking
- Later sections: biological/clinical discovery
- Final section: practical utility

Is the responder vs. progressor finding positioned as a climax?

### 5. Reader Confusion Points
Read as a naive reader would:
- Where would someone get lost?
- What terms are used before being defined?
- What assumes knowledge the reader might not have?
- What would make someone stop reading?

### 6. Consistency of Voice
- Is the writing style consistent throughout?
- Are claims appropriately hedged or confident throughout?
- Does the level of technical detail stay consistent?

## Section Agent Reviews

Read the Phase 1 section reviews in `/workspace/reviews/` to understand each section's self-contained assessment.

## Output Format

Write your review to: `/workspace/reviews/narrative-flow-review.md`

```markdown
# Cross-Cutting Review: Narrative Flow

## Executive Summary
[2-3 sentences: Does the paper tell a compelling story? Where does it break down?]

## Findings

### Section Transitions
| Transition | Rating | Issue |
|------------|--------|-------|
| 2.1 → 2.2 | Smooth/Adequate/Abrupt | [notes] |
| 2.2 → 2.3 | ... | ... |
| 2.3 → 2.4 | ... | ... |
| 2.4 → 2.5 | ... | ... |
| 2.5 → 2.6 | ... | ... |

### Problem-Solution Clarity
[Your analysis per section]

### "So What" Presence
| Section | "So What" | Explicit or Implied? |
|---------|-----------|---------------------|
| 2.1 | [reason to care] | Explicit/Implied/Missing |
| ... | ... | ... |

### Escalating Impact
[Does the narrative build? Where does it peak?]

### Reader Confusion Points
[List specific locations where readers might get lost]

### Consistency of Voice
[Any jarring shifts?]

## Recommendations

| Priority | Issue | Recommendation | Location |
|----------|-------|----------------|----------|
| CRITICAL | [issue] | [fix] | Section X.X |
| MAJOR | [issue] | [fix] | Section X.X |
| MINOR | [issue] | [fix] | Section X.X |
```

## Data Access

- `/data/manuscript/` - Full manuscript text
- `/workspace/reviews/` - Phase 1 section reviews

Reviewers read the whole paper. If they get lost, they reject.
