# Manuscript Review Agent Swarm - Implementation Plan

**Date**: 2026-02-07
**Design Doc**: `docs/plans/2026-02-07-manuscript-review-agents-design.md`
**Execution Mode**: Interactive (building together, then launch)

---

## Phase 1: Create Directory Structure

```bash
mkdir -p scripts/manuscript_review
mkdir -p scripts/manuscript_review/prompts/crosscutting
```

**Files to create:**
- `scripts/manuscript_review/orchestrator.sh`
- `scripts/manuscript_review/prompts/section_template.md`
- `scripts/manuscript_review/prompts/crosscutting/benchmarking_rigor.md`
- `scripts/manuscript_review/prompts/crosscutting/figure_text_alignment.md`
- `scripts/manuscript_review/prompts/crosscutting/narrative_flow.md`
- `scripts/manuscript_review/prompts/crosscutting/competitor_positioning.md`
- `scripts/manuscript_review/prompts/synthesis.md`

---

## Phase 2: Extract Section Text

Parse `manuscript/CITEgeist_Patterns_v4.md` to extract each Results subsection:
- Section 2.1: Lines 43-57
- Section 2.2: Lines 59-69
- Section 2.3: Lines 71-84
- Section 2.4: Lines 86-98
- Section 2.5: Lines 100-111
- Section 2.6: Lines 113-126

Store as individual files for agent prompts:
- `scripts/manuscript_review/sections/section_2.1.md`
- `scripts/manuscript_review/sections/section_2.2.md`
- etc.

---

## Phase 3: Create Prompt Templates

### 3.1 Section Agent Template

Create `prompts/section_template.md` with:
- Context about CITEgeist and prior rejections
- The 4 key advantages to hammer
- 6 evaluation criteria
- Output format specification
- Placeholder for section text: `{{SECTION_TEXT}}`

### 3.2 Cross-Cutting Agent Prompts

Create 4 focused prompts:
- `benchmarking_rigor.md` - Fair comparisons, proper tuning, appropriate metrics
- `figure_text_alignment.md` - Claims match figures, no orphan results
- `narrative_flow.md` - Section transitions, "so what" clarity
- `competitor_positioning.md` - Fair to other methods, clear CITEgeist edge

### 3.3 Synthesis Agent Prompt

Create `synthesis.md` with:
- Instructions to read all 10 reviews
- Priority categorization (Critical/Major/Minor)
- Output format for CONSOLIDATED_ISSUES.md

---

## Phase 4: Create Orchestrator Script

`orchestrator.sh` will:
1. Create `/workspace/reviews/` directory
2. Generate filled prompts for each section agent
3. Document launch commands for all 11 agents (6 + 4 + 1)

Note: Due to hpc-isolate architecture, we'll launch agents manually in phases rather than having a script spawn sub-jobs. The orchestrator prepares everything.

---

## Phase 5: Identify Data Paths

Determine exact paths for `--data` mounts:
- Manuscript: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript`
- Outputs: Need to identify where analysis outputs live
- Codebase: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist`

Also need:
- Figures directory location
- Benchmarking results location
- Patient analysis outputs location

---

## Phase 6: Create Review Workspace

Create a dedicated workspace for the review:
```bash
mkdir -p ~/alc376_bgfs/manuscript_review_2026-02-07
```

Copy prompts and orchestrator to workspace.

---

## Phase 7: Launch Phase 1 Agents

Launch 6 section agents in parallel via `/hpc-isolate:launch`:
```bash
# Each agent gets same data mounts, unique prompt
/hpc-isolate:launch ~/alc376_bgfs/manuscript_review_2026-02-07 \
    --data /ix1/.../manuscript \
    --data /ix1/.../CITEgeist \
    --prompt "Review section 2.1 per prompts/section_2.1_filled.md, write to reviews/"
```

Monitor with `/hpc-isolate:list` until all complete.

---

## Phase 8: Launch Phase 2 Agents

After Phase 1 completes, launch 4 cross-cutting agents.
They read both manuscript AND Phase 1 reviews (now in workspace).

---

## Phase 9: Launch Synthesis Agent

After Phase 2 completes, launch synthesis agent to produce `CONSOLIDATED_ISSUES.md`.

---

## Phase 10: Review Results

Read `CONSOLIDATED_ISSUES.md` and triage:
- Critical: Address immediately
- Major: Address before submission
- Minor: Address if time permits

---

## Execution Checklist

- [ ] Create directory structure
- [ ] Extract section text from manuscript
- [ ] Create section agent template
- [ ] Create 4 cross-cutting prompts
- [ ] Create synthesis prompt
- [ ] Identify all data paths
- [ ] Create review workspace
- [ ] Launch 6 section agents (Phase 1)
- [ ] Wait for Phase 1 completion
- [ ] Launch 4 cross-cutting agents (Phase 2)
- [ ] Wait for Phase 2 completion
- [ ] Launch synthesis agent (Phase 3)
- [ ] Review CONSOLIDATED_ISSUES.md
