#!/bin/bash
# Consolidate duplicate folders in CITEgeist project
# This script follows git-tracked canonical folders and archives duplicates

set -e

PROJ_ROOT="/ihome/alee/alc376/alc376_bgfs/CITEgeist"
ARCHIVE_DIR="${PROJ_ROOT}/_archive_duplicates/$(date +%Y%m%d_%H%M%S)"

cd "$PROJ_ROOT"

echo "=== CITEgeist Folder Consolidation ===" echo ""
echo "Canonical (git-tracked) locations:"
echo "  - tests/ (project root)"
echo "  - CITEgeist/model/"
echo ""
echo "Duplicates to consolidate:"
echo "  - CITEgeist/tests/"
echo ""

# Create archive directory
mkdir -p "$ARCHIVE_DIR"

# ============================================================================
# 1. Handle CITEgeist/tests/ duplicate
# ============================================================================
echo "Step 1: Analyzing CITEgeist/tests/ vs tests/"

# Copy unique files from CITEgeist/tests/ to tests/ if they don't exist
if [ -d "CITEgeist/tests" ]; then
    for file in CITEgeist/tests/*; do
        basename_file=$(basename "$file")
        target="tests/$basename_file"
        
        if [ ! -e "$target" ]; then
            echo "  [COPY] Unique file: $basename_file → tests/"
            cp -r "$file" "$target"
        elif ! diff -q "$file" "$target" > /dev/null 2>&1; then
            echo "  [DIFF] File differs: $basename_file (keeping tests/ version, archiving CITEgeist/tests/ version)"
        fi
    done
    
    # Archive the duplicate directory
    echo "  [ARCHIVE] Moving CITEgeist/tests/ → $ARCHIVE_DIR/"
    mv CITEgeist/tests "$ARCHIVE_DIR/"
else
    echo "  CITEgeist/tests/ does not exist, skipping"
fi

# ============================================================================
# 2. Check for any other duplicate 'model' directories
# ============================================================================
echo ""
echo "Step 2: Checking for duplicate 'model' directories"

if [ -d "model" ] && [ -d "CITEgeist/model" ]; then
    echo "  WARNING: Both /model and /CITEgeist/model exist"
    echo "  Git-tracked canonical: CITEgeist/model/"
    
    # Check if /model has unique content
    has_unique=false
    for file in model/*; do
        basename_file=$(basename "$file")
        target="CITEgeist/model/$basename_file"
        
        if [ ! -e "$target" ]; then
            echo "  [UNIQUE] $basename_file only in /model"
            has_unique=true
        fi
    done
    
    if [ "$has_unique" = true ]; then
        echo "  [ACTION NEEDED] /model has unique files. Manual review recommended."
    else
        echo "  [ARCHIVE] /model appears to be duplicate of CITEgeist/model/"
        mv model "$ARCHIVE_DIR/"
    fi
else
    echo "  Only one model directory exists (expected)"
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=== Consolidation Complete ===" 
echo" "
echo "Archived duplicates in: $ARCHIVE_DIR"
echo ""
echo "Canonical structure:"
echo "  /ihome/alee/alc376/alc376_bgfs/CITEgeist/"
echo "  ├── CITEgeist/model/    (source code)"
echo "  └── tests/              (test files)"
echo ""
echo "To prevent future confusion, consider adding to .gitignore:"
echo "  CITEgeist/tests/"
echo ""
