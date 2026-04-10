#!/usr/bin/env bash
# sync_to_main.sh — Copy allowlisted paths from dev to main, remove excluded paths.
# Usage: ./scripts/sync_to_main.sh [--dry-run]
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
DATE=$(date +%Y-%m-%d)
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# ── Preflight checks ──────────────────────────────────────────────────────
CURRENT_BRANCH=$(git branch --show-current)
if [[ "$CURRENT_BRANCH" != "dev" ]]; then
  echo "ERROR: must be on dev branch (currently: $CURRENT_BRANCH)" >&2
  exit 1
fi
if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "ERROR: working tree is not clean — commit or stash first" >&2
  exit 1
fi

# ── Fetch both remotes ────────────────────────────────────────────────────
echo "Fetching remotes..."
git fetch origin
git fetch dev-repo

# ── Switch to main ────────────────────────────────────────────────────────
echo "Switching to main..."
git checkout main
git pull origin main

# ── Copy allowlisted paths from dev ───────────────────────────────────────
ALLOWLIST=(
  "CITEgeist/"
  "tests/"
  "examples/"
  "docs/CITEgeist/"
  "docs/Benchmarking/"
  "docs/codebase_index.md"
  "docs/codebase_index_summary.md"
  "docs/installation.md"
  "docs/analysis/"
  "docs/detection_estimation_algorithm.md"
  "docs/detection_post_filtering.md"
  "README.md"
  "CONTRIBUTING.md"
  "LICENSE"
  "pyproject.toml"
  "pytest.ini"
  "requirements.txt"
  "requirements-dev.txt"
  ".pre-commit-config.yaml"
  ".gitignore"
  "scripts/sync_to_main.sh"
  ".github/"
)

echo "Copying allowlisted paths from dev..."
for path in "${ALLOWLIST[@]}"; do
  # Strip trailing slash for git checkout
  clean_path="${path%/}"
  if git show "dev:${clean_path}" &>/dev/null 2>&1; then
    git checkout dev -- "${clean_path}" 2>/dev/null && echo "  ✓ $path" || echo "  ✗ FAILED: $path"
  else
    echo "  - SKIP (not in dev): $path"
  fi
done

# ── Remove excluded paths tracked on main ─────────────────────────────────
EXCLUDELIST=(
  "docs/superpowers/"
  "docs/plans/"
  "docs/handoffs/"
  "docs/review_2026-02-07/"
  "docs/review_2026-04-01/"
  "docs/PROGRESS_REPORT.md"
  "docs/Design_AutoProfileDiscovery.md"
  "docs/diagram.md"
  "docs/mdk_evidence_ledger.md"
  "docs/module-1-2-leiden-comparison.md"
  "Benchmarking/"
  "midkine/"
  "output/"
  "manuscript/"
  "powerpoints/"
  "slurm/"
  ".babysit/"
  ".superpowers/"
  "AGENTS.md"
  "CLAUDE.md"
)

echo "Removing excluded paths from main..."
for path in "${EXCLUDELIST[@]}"; do
  if git ls-files --error-unmatch "${path%/}" &>/dev/null 2>&1 || \
     git ls-files "${path%/}" | grep -q .; then
    git rm -r --cached --ignore-unmatch "${path%/}"
    rm -rf "${REPO_ROOT}/${path%/}"
    echo "  ✓ removed: $path"
  fi
done

# ── Show diff ─────────────────────────────────────────────────────────────
echo ""
echo "=== Changes staged for commit ==="
git diff --staged --stat
echo ""

if $DRY_RUN; then
  echo "DRY RUN — reverting all changes"
  git checkout -- .
  git clean -fd
  git checkout dev
  exit 0
fi

# ── Confirm + commit + push ───────────────────────────────────────────────
read -p "Proceed with commit and push to origin main? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  echo "Aborted — reverting changes"
  git checkout -- .
  git clean -fd
  git checkout dev
  exit 1
fi

git commit -m "chore: sync dev→main $DATE"
git push origin main
git checkout dev
echo ""
echo "Done! Synced dev→main on $DATE"
