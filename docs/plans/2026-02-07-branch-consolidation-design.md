# Branch Consolidation and Main Merge Design

**Date:** 2026-02-07
**Status:** Approved
**Purpose:** Consolidate development branches and merge to main for publication

## Background

The CITEgeist repository has accumulated multiple branches during active development:
- `hierarchical_approach`: 383 commits ahead of main (primary development branch)
- `docs`: 172 commits ahead of main (fully contained in hierarchical_approach)
- Multiple feature branches with varying states of completion
- Stale experimental branches

This design consolidates the work into a clean state for publication.

## Branch Analysis

### Primary Branches
| Branch | Commits from main | Relationship |
|--------|-------------------|--------------|
| `hierarchical_approach` | 383 | Primary development branch |
| `docs` | 172 | Ancestor of hierarchical_approach |

### Feature Branches with Unique Commits
| Branch | Unique Commits | Decision |
|--------|---------------|----------|
| `feature/singlecell-demonstration` | 14 | Keep (WIP) |
| `feature/autodiscovery-enhancement` | 6 | Archive (reverted changes) |
| `feature/em-joint-optimization` | 17 | Archive (failed experiment) |
| `dev` | 2 | Archive (old, superseded) |
| `feature/hierarchical-profiles` | 0 | Delete (fully merged) |

### Stale Branches (No Unique Work)
- `phase1`, `phase2`, `phase3`, `phase3_iterative`, `phase1_locoregional`
- `cleanup`, `cuOPT`, `nmf_approach`
- `claude/*` (automated session branches)
- `overnight/*` (automated overnight branches)

## Execution Plan

### Phase 1: Commit Untracked Files

Before merging, commit essential outputs on `hierarchical_approach`:

**Files to commit (~691 files):**
- Figures: `figures/`, `manuscript/figures/`, benchmark figures
- Scripts: Benchmark src/, evaluation src/, example scripts
- Evaluation results: `evaluation/results/`, metrics CSVs
- Design docs: `docs/plans/*.md`
- Midkine analysis: `mdk_saturation_pipeline/outputs/`
- Ground truth data: `xenium_pseudovisium/data/ground_truth_gex/`

**Files excluded via .gitignore:**
- Large data files: `*.h5ad`, `*.h5seurat`, `*.h5`
- Reference data: `GSE*` directories
- Large outputs: `*_layer_pass1.csv`, `*.parquet`
- Checkpoints: `*.npy`, `*.npz`
- Transient: SLURM logs, `.worktrees/`, tool configs

### Phase 2: Archive Stale Branches

Create archive tags before deletion:
```bash
git tag archive/dev dev
git tag archive/autodiscovery-enhancement feature/autodiscovery-enhancement
git tag archive/em-joint-optimization feature/em-joint-optimization
git push origin --tags
```

Delete branches:
```bash
# Local
git branch -D dev feature/autodiscovery-enhancement feature/em-joint-optimization

# Remote
git push origin --delete dev feature/autodiscovery-enhancement feature/em-joint-optimization
```

### Phase 3: Merge to Main

```bash
git checkout main
git merge docs                    # Fast-forward
git merge hierarchical_approach   # Merge commit
git push origin main
```

### Phase 4: Post-Merge Cleanup

Rename hierarchical_approach to dev:
```bash
git branch -m hierarchical_approach dev
git push origin :hierarchical_approach  # Delete old remote
git push origin dev -u                   # Push new dev
```

Delete merged/stale branches:
```bash
# Delete docs (merged)
git branch -D docs
git push origin --delete docs

# Delete stale branches
git branch -D phase1 phase2 phase3 phase3_iterative phase1_locoregional
git branch -D cleanup cuOPT nmf_approach feature/hierarchical-profiles

# Delete from remote
git push origin --delete phase1 phase2 phase3 phase3_iterative phase1_locoregional
git push origin --delete cleanup cuOPT nmf_approach feature/hierarchical-profiles

# Delete claude/* and overnight/* branches
git branch -D $(git branch --list 'claude/*')
git branch -D $(git branch --list 'overnight/*')
git push origin --delete $(git branch -r --list 'origin/claude/*' | sed 's|origin/||')
git push origin --delete $(git branch -r --list 'origin/overnight/*' | sed 's|origin/||')
```

## Final Repository State

### Branches
| Branch | Purpose |
|--------|---------|
| `main` | Publication-ready, contains all consolidated work |
| `dev` | Active development (renamed from hierarchical_approach) |
| `feature/singlecell-demonstration` | WIP for future single-cell work |

### Archive Tags
| Tag | Original Branch |
|-----|-----------------|
| `archive/dev` | Old dev branch |
| `archive/autodiscovery-enhancement` | Failed autodiscovery experiments |
| `archive/em-joint-optimization` | Failed EM optimization experiments |

## Rollback Plan

If issues arise after merge:
```bash
# Reset main to pre-merge state
git checkout main
git reset --hard <commit-before-merge>
git push origin main --force

# Restore archived branches from tags
git checkout -b dev archive/dev
```

## Verification Checklist

- [ ] All 691 untracked files committed
- [ ] No files > 100MB in commit
- [ ] Archive tags created and pushed
- [ ] main contains all work from hierarchical_approach
- [ ] dev branch exists and tracks main
- [ ] Stale branches deleted locally and remotely
- [ ] GitHub repository shows correct branch structure
