# Archived Model Code

This directory contains archived CITEgeist model code that is no longer used in production.

## Why Archive Instead of Delete?

These files represent completed development efforts that may be useful for:
- Reference for future development
- Understanding design decisions
- Potential reuse of specific algorithms
- Historical documentation

## Archived Files

### RNA-Only Pipeline (Archived: February 2026)

**Files:**
- `rna_pipeline.py` - Main RNA-only pipeline orchestration
- `rna_marker_selection.py` - RNA marker detection and selection algorithms

**Reason for archiving:**
The RNA-only pipeline was developed as an alternative to the protein-based (CITE-seq) approach. However, testing revealed that RNA-only marker detection is significantly less reliable due to:

1. **RNA dropout** - Single-cell RNA-seq data suffers from technical dropout, where expressed genes fail to be detected, leading to false negatives
2. **Noise** - RNA expression is inherently noisier than protein expression, making spatial colocalization patterns harder to detect
3. **Lower signal-to-noise ratio** - Protein markers (from CITE-seq) provide cleaner signals for cell type identification

The protein-based pipeline (`marker_interest.py`, `spatial_colocalization.py`) remains the recommended approach for CITEgeist analysis.

### Previously Archived Files

- `auto_profile_discovery.py` - Early automated profile discovery approach
- `background_correction.py` - Background correction routines
- `cell_assignment.py` - Cell assignment algorithms
- `cell_classification.py` - Cell classification code
- `joint_optimization.py` - Joint optimization approach
- `marker_quality.py` - Marker quality assessment
- `profile_matching.py` - Profile matching utilities
- `spot_level_enhancement.py` - Spot-level enhancement methods

## Note

These files are NOT maintained and may have outdated dependencies or API incompatibilities with the current codebase. Use at your own risk for reference purposes only.
