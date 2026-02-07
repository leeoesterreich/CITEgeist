# Manuscript Review Agent Swarm Design

**Date**: 2026-02-07
**Target**: Cell Reports Methods / Patterns (Cell Press multi-journal submission)
**Goal**: Systematic agent-driven review of Results section to maximize impact before submission

---

## Background

CITEgeist has been desk-rejected before for "insufficient impact" when only Module 3 was shown. The full 5-module pipeline is now complete. This review ensures every section:
1. Shows the most impressive available result
2. Clearly communicates scientific edge and contribution
3. Takes full advantage of all work done
4. Anticipates and preempts reviewer objections

---

## Agent Architecture

### Tier 1: Section Agents (6 agents, run in parallel)

| Agent | Section | Focus |
|-------|---------|-------|
| `section-2.1-reviewer` | Framework Overview | Is the pipeline introduction compelling? Does it set up the value prop? |
| `section-2.2-reviewer` | Modules 1-2 (Profile Discovery) | Is autodiscovery being sold hard enough? Best validation shown? |
| `section-2.3-reviewer` | Module 3 (Deconvolution) | Is the benchmarking the best we have? Reference-free advantage clear? |
| `section-2.4-reviewer` | Module 4 (Program Discovery) | Are the program examples compelling? Moran's I validation sufficient? |
| `section-2.5-reviewer` | Module 5 (Cross-Sample Integration) | Is the responder/progressor story maximally impactful? |
| `section-2.6-reviewer` | Downstream Tools | Is interoperability demonstrated convincingly? |

### Tier 2: Cross-Cutting Agents (4 agents, run after Tier 1)

| Agent | Scope | Focus |
|-------|-------|-------|
| `benchmarking-rigor` | Whole manuscript | Fair comparisons? Proper baselines? Appropriate metrics? |
| `figure-text-alignment` | Whole manuscript | Every claim backed by figure? Orphan results? |
| `narrative-flow` | Whole manuscript | Does each section build? Is "so what" clear? |
| `competitor-positioning` | Whole manuscript | Fair to other methods? CITEgeist edge articulated? |

### Tier 3: Synthesis Agent (1 agent, runs last)

Reads all 10 reviews and produces `CONSOLIDATED_ISSUES.md` with prioritized action items.

---

## Evaluation Criteria

### Section Agent Questions (6 per section)

1. **Referenced Figures Assessment**
   - Which figures does this section reference?
   - Do they effectively support the claims?
   - Are they publication-quality?
   - Could better visualizations be generated from available outputs?
   - Are there panels that should be added or removed?

2. **Most Impressive Result Available?**
   - Given all available outputs (not just what's in the paper), is this section showing the strongest evidence?
   - Are there better examples, stronger numbers, or more compelling visualizations in the outputs folder?

3. **Scientific Edge Communicated?**
   - Does this section explicitly state what CITEgeist does that others can't?
   - Is the "reference-free", "spatially-native", or "end-to-end" advantage hammered home?

4. **Unused Capabilities or Results?**
   - Are there capabilities in the codebase that aren't being demonstrated?
   - Results that were generated but not included?
   - Analyses that would strengthen the section?

5. **Quantitative Backing?**
   - Are claims supported by numbers?
   - Vague statements like "competitive accuracy" should have specific metrics attached.

6. **Anticipated Reviewer Objections?**
   - What would a skeptical reviewer attack?
   - Is the section proactively addressing likely criticisms?

### Cross-Cutting Agent Criteria

- **Benchmarking Rigor**: Methods tuned fairly? Same preprocessing? Appropriate metrics for the claim? Statistical tests where needed?
- **Figure-Text Alignment**: Every figure panel referenced? Text claims match visual? Results mentioned but not shown?
- **Narrative Flow**: Section transitions logical? Reader knows why each section matters? Builds toward clinical relevance?
- **Competitor Positioning**: Other methods' limitations stated fairly? Direct comparisons where possible? Acknowledging when others are better?

---

## Data Access

Each agent gets read-only access to:

```
/data/manuscript/          # Manuscript text + figures
  ├── CITEgeist_Patterns_v4.md
  ├── CITEgeist_Patterns_v4_Figures.md
  └── figures/             # All figure files

/data/outputs/             # Analysis outputs
  ├── xenium_benchmarking/ # Benchmark results
  ├── patient_analyses/    # Per-patient Module 3-5 outputs
  └── figures_generated/   # Scripts + generated figures

/data/codebase/            # Full CITEgeist repo (read-only)
  ├── CITEgeist/model/     # Core implementation
  ├── Benchmarking/        # Benchmark framework
  └── docs/plans/          # Design docs showing what was built
```

Each agent writes to:

```
/workspace/reviews/
  ├── section-2.1-review.md
  ├── section-2.2-review.md
  ├── section-2.3-review.md
  ├── section-2.4-review.md
  ├── section-2.5-review.md
  ├── section-2.6-review.md
  ├── benchmarking-rigor-review.md
  ├── figure-text-alignment-review.md
  ├── narrative-flow-review.md
  ├── competitor-positioning-review.md
  └── CONSOLIDATED_ISSUES.md
```

---

## Execution Workflow

### Phase 1: Section Agents (parallel)

Six isolated Claude instances launched simultaneously via `/hpc-isolate:launch`:
- Same read-only data mounts (manuscript, outputs, codebase)
- Unique prompt specifying their section and evaluation criteria
- 4 CPUs, 8GB RAM, 4-hour time limit

```
┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│ Section 2.1 │ │ Section 2.2 │ │ Section 2.3 │
└─────────────┘ └─────────────┘ └─────────────┘
┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│ Section 2.4 │ │ Section 2.5 │ │ Section 2.6 │
└─────────────┘ └─────────────┘ └─────────────┘
        ↓ all complete ↓
```

### Phase 2: Cross-Cutting Agents (parallel)

Four agents launched after Phase 1 completes. They read both the manuscript AND the section reviews from Phase 1.

```
┌──────────────┐ ┌─────────────┐ ┌──────────────┐ ┌────────────┐
│ Benchmarking │ │ Figure-Text │ │ Narrative    │ │ Competitor │
└──────────────┘ └─────────────┘ └──────────────┘ └────────────┘
        ↓ all complete ↓
```

### Phase 3: Synthesis Agent

Single agent reads all 10 reviews, produces `CONSOLIDATED_ISSUES.md`.

### Timeline

- Phase 1: ~2-3 hours each, parallel
- Phase 2: ~2 hours each, parallel
- Phase 3: ~30 minutes
- **Total wall time: ~5-6 hours**

---

## Output Format

### Individual Review Structure

```markdown
# Section X.X Review: [Title]

## Executive Summary
[2-3 sentence verdict: is this section landing its punch?]

## Findings

### Referenced Figures Assessment
[Which figures does this section reference? Do they effectively support the claims?
Are the figures publication-quality? Could better visualizations be generated from
available outputs? Are there panels that should be added or removed?]

### Most Impressive Result Available?
[What's shown vs. what exists in outputs. Specific file paths if better options found.]

### Scientific Edge Communicated?
[Is reference-free / spatially-native / end-to-end advantage explicit?]

### Unused Capabilities or Results?
[Code features not demonstrated. Analyses not included.]

### Quantitative Backing?
[Claims without numbers. Vague statements that need metrics.]

### Anticipated Reviewer Objections?
[What a skeptic would attack. Suggested preemptive responses.]

## Recommendations

| Priority | Issue | Recommendation | Evidence Location |
|----------|-------|----------------|-------------------|
| CRITICAL | ... | ... | /data/outputs/... |
| MAJOR | ... | ... | ... |
| MINOR | ... | ... | ... |
```

### Consolidated Issues Format

```markdown
# Manuscript Review: Consolidated Issues

## Critical (Must Fix Before Submission)
- [ ] [Section 2.3] Benchmarking shows r=0.60 but outputs have r=0.65 with updated radius
- [ ] [Figure-Text] Figure 4D claims "68% spatial coherence" but figure shows 62%

## Major (Strongly Recommended)
- [ ] [Section 2.5] Response analysis undersells - PyDESeq2 found 127 DE genes, only 7 mentioned
- [ ] [Competitor] Tangram comparison may be unfair - did we use optimal parameters?

## Minor (Nice to Have)
- [ ] [Narrative] Transition from 2.3→2.4 is abrupt, needs bridging sentence
```

---

## Agent Prompt Template

### Section Agent Prompt

```markdown
You are reviewing Section [X.X] of a manuscript being submitted to Cell Reports Methods / Patterns.

## Context
CITEgeist is a spatial transcriptomics framework. It has been desk-rejected before for
"insufficient impact" when only showing Module 3. The full 5-module pipeline is now complete.
Your job is to ensure this section makes maximum impact.

## Your Section
[Section title and full text pasted here]

## CITEgeist's Key Advantages (hammer these)
1. Reference-free: Uses same-slide protein, no external scRNA-seq atlas
2. Spatially-native: Spatial statistics at every module, not bolted on
3. End-to-end pipeline: Raw data → cross-patient programs in one framework
4. Competitive accuracy: Matches Cell2Location/RCTD without their reference requirements

## Your Task
Evaluate this section against 6 criteria. For each, be specific and cite file paths.

1. **Referenced Figures**: Which figures support this section? Are they compelling?
2. **Most Impressive Result**: Is there stronger evidence in /data/outputs/?
3. **Scientific Edge**: Is the advantage over other methods explicit?
4. **Unused Work**: Are there capabilities/results not being shown?
5. **Quantitative Backing**: Are claims supported by numbers?
6. **Reviewer Objections**: What would a skeptic attack?

Write your review to: /workspace/reviews/section-[X.X]-review.md
```

### Cross-Cutting Agent Prompts

Similar structure but scoped across all sections, focusing on their specific dimension (benchmarking rigor, figure-text alignment, narrative flow, or competitor positioning).

### Synthesis Agent Prompt

```markdown
You are synthesizing 10 manuscript reviews into a prioritized action list.

Read all reviews in /workspace/reviews/ and produce CONSOLIDATED_ISSUES.md with:
1. Critical issues (must fix before submission)
2. Major issues (strongly recommended)
3. Minor issues (nice to have)

For each issue:
- Tag the source section or cross-cutting dimension
- Include specific recommendation
- Reference evidence location if applicable

Format as a checkbox list for easy tracking.
```

---

## Implementation Files

```
scripts/
  └── manuscript_review_orchestrator.sh    # Launches all agents in phases

prompts/
  ├── section_agent_template.md            # Template filled per section
  ├── crosscutting/
  │   ├── benchmarking_rigor.md
  │   ├── figure_text_alignment.md
  │   ├── narrative_flow.md
  │   └── competitor_positioning.md
  └── synthesis_agent.md
```

---

## Execution Command

```bash
/hpc-isolate:launch /path/to/review_workspace \
    --data /ix1/.../CITEgeist/manuscript \
    --data /ix1/.../CITEgeist/Benchmarking \
    --data /ix1/.../CITEgeist/CITEgeist \
    --prompt "Execute manuscript review orchestration per scripts/manuscript_review_orchestrator.sh"
```

---

## Success Criteria

1. All 10 agents complete successfully
2. `CONSOLIDATED_ISSUES.md` produced with prioritized items
3. At least one "critical" or "major" issue identified per section (if none found, section is likely underselling)
4. Actionable recommendations with specific file paths for evidence
