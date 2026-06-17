#!/usr/bin/env bash
# sync_to_main.sh — Copy allowlisted paths from dev to main, remove excluded paths.
# Usage: ./scripts/sync_to_main.sh [--simulate|--dry-run]
#   --simulate  safe read-only audit of the would-ship public tree (no branch mutation)
#   --dry-run   deprecated alias for --simulate
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
DATE=$(date +%Y-%m-%d)
DRY_RUN=false
SIMULATE=false
for _arg in "$@"; do
  [[ "$_arg" == "--dry-run"  ]] && DRY_RUN=true
  [[ "$_arg" == "--simulate" ]] && SIMULATE=true
done

# ── Safe read-only simulation (MUST come before any working-tree mutation) ──
# --dry-run is a deprecated alias for --simulate: routes to the safe read-only
# simulation instead of the hazardous checkout+clean dry-run (retired 2026-06-17).
if $SIMULATE || $DRY_RUN; then
  if $DRY_RUN && ! $SIMULATE; then
    echo "[deprecated] --dry-run now runs the safe read-only simulation (--simulate); the old checkout+clean dry-run is retired."
  fi
  echo "Fetching origin/main (read-only, safe)..."
  git fetch origin main
  exec python "${REPO_ROOT}/scripts/audit_public_tree.py" \
    --simulate-sync \
    --ref origin/main \
    --dev-ref dev \
    --allow-from "${REPO_ROOT}/scripts/sync_to_main.sh"
fi

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
# Reject untracked (non-ignored) files. The revert path runs `git clean -fd`,
# which deletes ALL untracked files — including pre-existing ones unrelated to
# the sync. Catching them here (still on dev) prevents data loss.
UNTRACKED=$(git ls-files --others --exclude-standard)
if [[ -n "$UNTRACKED" ]]; then
  echo "ERROR: untracked files present — commit, stash (git stash -u), or remove them first." >&2
  echo "       The revert path runs 'git clean -fd' and would DELETE these:" >&2
  echo "$UNTRACKED" | sed 's/^/         /' >&2
  exit 1
fi

# Safe revert back to a clean dev checkout. Untracked files are rejected at
# preflight, so anything untracked here was created by this script (dev copies)
# and is safe to clean.
revert_and_return_to_dev() {
  git reset --hard HEAD
  git clean -fd
  git checkout dev
}

# ── Fetch both remotes ────────────────────────────────────────────────────
echo "Fetching remotes..."
git fetch origin
git fetch dev-repo

# ── Switch to main ────────────────────────────────────────────────────────
echo "Switching to main..."
git checkout main
# ERR trap: from here on we are on `main`. Any unexpected failure must return to
# dev WITHOUT git clean. Installed BEFORE `git pull` so a pull failure is caught,
# and kept active through staging/audit/commit/push/final-checkout — cleared only
# once we are safely back on dev (normal completion).
trap 'echo "ERROR: sync failed unexpectedly — returning to dev (no clean)." >&2; git checkout -f dev 2>/dev/null || true' ERR
git pull origin main

# ── Copy allowlisted paths from dev ───────────────────────────────────────
ALLOWLIST=(
  "CITEgeist/"
  "tests/"
  "examples/"
  "docs/CITEgeist/"
  "docs/Benchmarking/"
  "docs/installation.md"
  "docs/analysis/"
  "docs/detection_estimation_algorithm.md"
  "docs/detection_post_filtering.md"
  "docs/quickstart_real_visium.md"
  "docs/source/"
  "CITEgeist_env.yml"
  "postBuild"
  ".pylintrc"
  ".readthedocs.yaml"
  "README.md"
  "CONTRIBUTING.md"
  "LICENSE"
  "pyproject.toml"
  "pytest.ini"
  "requirements.txt"
  "requirements-dev.txt"
  ".gitignore"
  "scripts/sync_to_main.sh"
  ".github/"
)

# ── Wipe .github/workflows/ so stale-on-main-only workflow files don't survive ──
# git checkout only copies/updates; it never deletes files present on main but absent from dev.
if git ls-files .github/workflows/ | grep -q .; then
  git rm -r --ignore-unmatch .github/workflows/ >/dev/null  # tracked files only; leaves gitignored local data
  echo "Cleared .github/workflows/ (will be restored from dev)"
fi

# ── Wipe stale directories that should not be on main ──
# These are untracked on dev but survive sync because git checkout copies but never deletes.
STALE_WIPE_DIRS=(
  "CITEgeist/model/_archive/"
  "CITEgeist/junk/"
  "CITEgeist/examples/_archive/"
  "CITEgeist/tests/"
  "examples/_archive/"
)
for dir in "${STALE_WIPE_DIRS[@]}"; do
  if git ls-files "${dir}" | grep -q .; then
    git ls-files "${dir}" | xargs git rm -f
    echo "Cleared stale: ${dir}"
  fi
done

# ── Wipe root-level model .py files that were git-mv'd into subpackages on dev ──
# Uses git ls-files (tracked-only) to avoid failing on untracked working-tree files.
KEEP_MODEL_FILES=("__init__.py" "citegeist_model.py" "checkpoints.py" "utils.py" "unified_config.py" "module2_proposal_builder.py")
for f in $(git ls-files 'CITEgeist/model/*.py'); do
  base=$(basename "$f")
  dir=$(dirname "$f")
  [[ "$dir" != "CITEgeist/model" ]] && continue
  keep=false
  for k in "${KEEP_MODEL_FILES[@]}"; do [[ "$base" == "$k" ]] && keep=true && break; done
  if ! $keep; then
    git rm -f "$f"
    echo "Cleared stale root model file: $f"
  fi
done

echo "Copying allowlisted paths from dev..."
for path in "${ALLOWLIST[@]}"; do
  # Strip trailing slash for git checkout
  clean_path="${path%/}"
  # Directory entries: drop main's version of the dir from the index FIRST, so
  # stale main-orphans (files present on main but removed from dev — e.g. old
  # examples/, retired model/test files from a prior sync) do NOT survive the
  # `git checkout dev --` overlay (checkout only adds/updates, never deletes).
  # This makes each allowlisted directory contain EXACTLY dev's content.
  if [[ "$path" == */ ]]; then
    if git ls-files "${clean_path}/" | grep -q .; then
      git ls-files "${clean_path}/" | xargs git rm -f >/dev/null
    fi
  fi
  if git show "dev:${clean_path}" &>/dev/null 2>&1; then
    if git checkout dev -- "${clean_path}" 2>/dev/null; then
      echo "  ✓ $path"
    else
      echo "  ✗ FAILED to refresh from dev: $path" >&2
      exit 1
    fi
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
  "repro/"
  "figures/"
  "docs/codebase_index.md"
  "docs/codebase_index_summary.md"
  "docs/core_scripts/"
  "CITEgeist/examples/"
  "CITEgeist/slurm/"
  "examples/slurm/"
  "examples/_archive/"
  "tests/test_figures/"
  "tests/test_figure_tools/"
  "tests/repo_audit/"
  "tests/benchmarking/"
  "tests/logs/"
  "tests/slurm_log/"
  "tests/slurm/"
  "tests/test_v33_manuscript.py"
  "tests/test_v34_manuscript.py"
  "tests/test_v34_comment_restoration.py"
  "tests/test_manuscript_data.py"
  "tests/test_validate_docx_package.py"
  "tests/test_swap_docx_media.py"
  "tests/test_smoke_benchmark_outputs.py"
  "tests/test_integrate_summary.py"
  "tests/test_if_analysis.py"
  "tests/test_elisa_analysis.py"
  "tests/test_compute_gex_rmse_sd_ratio.py"
  "tests/test_panel_label_role.py"
  "tests/sbatch_module2c_slow.sh"
  "tests/sbatch_test_sweep.sh"
  "docs/source/notebooks/"         # patient vignette notebooks with private /bgfs paths + executed outputs
  "tests/test_v22_regressions.py"  # manuscript-regression test referencing private patient data
  "tests/test_smoke_patient_pipeline.py"  # private-patient-data integration test
  "docs/CITEgeist/REPRODUCIBILITY.md"    # references removed repro/ hub
  "examples/sample_paths.txt"      # bare patient specimen-ID list
  "output/"
  "manuscript/"
  "powerpoints/"
  "slurm/"
  ".babysit/"
  ".superpowers/"
  "AGENTS.md"
  "CLAUDE.md"
  ".pre-commit-config.yaml"
)

echo "Removing excluded paths from main..."
for path in "${EXCLUDELIST[@]}"; do
  if git ls-files --error-unmatch "${path%/}" &>/dev/null 2>&1 || \
     git ls-files "${path%/}" | grep -q .; then
    git rm -rf --ignore-unmatch "${path%/}" >/dev/null  # -f: excluded subpaths of allowlisted dirs are already staged from dev (differ from HEAD); force removal. tracked files only; leaves gitignored local data
    echo "  ✓ removed: $path"
  fi
done

# ── Show diff ─────────────────────────────────────────────────────────────
echo ""
echo "=== Changes staged for commit ==="
git diff --staged --stat
echo ""

# ── Audit the staged public tree (path-management + PII + excluded-ref/log) ──
echo ""
echo "=== Auditing staged public tree ==="
# BLOCKING #1: the auditor + the allow/deny source live on `dev`, NOT on `main`.
# We are on `main` now (post-checkout), so scripts/audit_public_tree.py and the
# up-to-date scripts/sync_to_main.sh are NOT in the working tree. Materialize both
# from `dev` into a temp dir and run from there. (audit_ref only shells out to
# `git ls-tree`/`git show` on the tree SHA, so its CWD/branch is irrelevant.)
AUDIT_DIR=$(mktemp -d)
git show dev:scripts/audit_public_tree.py > "${AUDIT_DIR}/audit_public_tree.py"
git show dev:scripts/sync_to_main.sh       > "${AUDIT_DIR}/sync_to_main.sh"
# HEAD is still the last main commit; audit the WORKING TREE + index instead.
# Stage everything so git ls-tree of a temp tree reflects what will be committed:
AUDIT_TREE=$(git write-tree)
if ! python "${AUDIT_DIR}/audit_public_tree.py" \
      --ref "$AUDIT_TREE" --allow-from "${AUDIT_DIR}/sync_to_main.sh"; then
  echo "AUDIT FAILED — staged tree has leakage/management violations (see above)." >&2
  echo "Reverting to dev." >&2
  rm -rf "${AUDIT_DIR}"
  revert_and_return_to_dev
  exit 2
fi
rm -rf "${AUDIT_DIR}"
echo "Audit clean."

# ── Confirm + commit + push ───────────────────────────────────────────────
read -p "Proceed with commit and push to origin main? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  echo "Aborted — reverting changes"
  revert_and_return_to_dev
  exit 1
fi

git commit -m "chore: sync dev→main $DATE"
git push origin main
git checkout dev
trap - ERR  # back on dev — clear the ERR trap (normal completion)
echo ""
echo "Done! Synced dev→main on $DATE"
