# Folder Consolidation Log

## Benchmarking Consolidation - December 27, 2025

### Issue
Duplicate folder structures existed:
- `benchmarking/` (root) and `Benchmarking/benchmarking/` (nested)
- `simulation/` (root) and `Benchmarking/simulation/` (nested)

Files had diverged with `Benchmarking/` versions being ~0.5-3 seconds newer, while root `benchmarking/` had unique Runtime log files.

### Resolution

#### New Structure
```
Benchmarking/
├── README.md
├── simulation_generation/       (renamed from simulation/)
├── simulation_benchmarking/     (renamed from benchmarking/)
│   ├── Runtime/                 (merged logs from root benchmarking/)
│   └── ...
└── xenium_benchmarking/         (renamed from xenium/)
```

#### Changes Made
1. Archived root-level folders to `_archive_duplicates/20251227_benchmarking_consolidation/`
2. Copied unique Runtime logs from root `benchmarking/` to `Benchmarking/benchmarking/`
3. Renamed folders within `Benchmarking/`:
   - `simulation` → `simulation_generation`
   - `benchmarking` → `simulation_benchmarking`
   - `xenium` → `xenium_benchmarking`
4. Removed root-level `benchmarking/` and `simulation/`
5. Updated `.gitignore` to prevent future duplication

#### Files Merged
- `Runtime/high_seg/*.log` (5 files)
- `Runtime/mixed/*.log` (5 files)

---

# Folder Consolidation - December 3, 2024

## Issue
Git pull created duplicate folder structure causing import confusion:
- `/tests/` and `/CITEgeist/tests/`  
- Both contained similar but not identical files

## Resolution

### Consolidated Structure
```
/ihome/alee/alc376/alc376_bgfs/CITEgeist/
├── CITEgeist/model/     # Source code (git-tracked canonical)
└── tests/                # Test files (git-tracked canonical)
```

### Archived Duplicates
- `CITEgeist/tests/` → `_archive_duplicates/20251203_141743/tests/`

### Changes Made
1. Identified git-tracked canonical locations (`tests/` and `CITEgeist/model/`)
2. Compared duplicate test directories
3. All unique files from tests/ were already newer/canonical
4. Archived `CITEgeist/tests/` to prevent confusion

### Prevention
Consider adding to `.gitignore`:
```
CITEgeist/tests/
_archive_duplicates/
```

## Validation
After consolidation, integration test with real simulated data **PASSED**:
- Tested underspecified cell types (3 instead of 9)
- Correctly triggered validation error: "Unknown cell type has 25.07% > 5.00%"
- Validation implementation working correctly with actual data
