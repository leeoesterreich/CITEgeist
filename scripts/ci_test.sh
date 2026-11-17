#!/usr/bin/env bash
# The CI `test` job's pytest command. Run from the repository root; used by CI
# and by the acceptance gate so both execute the same thing.
set -euo pipefail

# pytest imports every collected module before marker filtering, so a
# module-scope import of a package outside the dependency closure
# aborts the run with RC 2. tests/ci_collect_skip.txt lists any such
# paths; it is optional, and paths in it that do not exist are
# skipped, so this works unchanged on a tree that has none of them.
IGNORES=""
if [ -f tests/ci_collect_skip.txt ]; then
  while read -r line; do
    case "$line" in ''|'#'*) continue ;; esac
    [ -e "$line" ] && IGNORES="$IGNORES --ignore=$line"
  done < tests/ci_collect_skip.txt
fi
echo "collection ignores:${IGNORES:- (none)}"
# shellcheck disable=SC2086  # IGNORES is intentionally word-split
pytest tests/ \
  -m "unit and not slow and not requires_data and not requires_cuopt and not requires_hpc" \
  $IGNORES \
  --cov=CITEgeist/model \
  --cov-report=xml \
  --cov-report=term-missing \
  --cov-fail-under=15 \
  --tb=short \
  -v
