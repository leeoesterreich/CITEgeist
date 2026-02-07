# Archived Test Code

This directory contains archived test files for deprecated or archived CITEgeist functionality.

## Archived Files

### RNA-Only Pipeline Tests (Archived: February 2026)

**Files:**
- `test_rna_marker_selection.py` - Unit tests for RNA marker selection
- `test_rna_vs_protein_xenium.py` - Comparison tests: RNA vs protein-based approaches
- `sbatch_rna_marker_tests.sh` - SLURM batch script for RNA marker tests
- `sbatch_rna_vs_protein_validation.sh` - SLURM batch script for RNA vs protein validation

**Reason for archiving:**
These tests are for the RNA-only pipeline which was archived due to poor performance caused by RNA noise and dropout issues. See `CITEgeist/model/_archive/README.md` for more details.

## Note

These test files are NOT maintained and reference archived modules. They will not run without the corresponding archived model files.
