# Auto Profile Discovery Implementation Plan & Progress

**Created:** 2025-12-05
**Status:** In Progress - Environment Setup Phase
**Primary Goal:** Validate mathematical formulation on Wu_rep_0 simulated h5ad data

---

## Quick Resume Instructions

To resume this work:

```bash
# 1. Load environment
module load gurobi/12.0.3
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# 2. Check Python environment (should be Python 3.10.x from CITEgeist_env)
which python
python --version

# 3. Install pytest if needed
pip install pytest

# 4. Test imports
python -c "from CITEgeist.model import discover_profiles; print('Import OK')"

# 5. Continue with Phase 2 (fix paths in test file)
```

---

## Current Status Summary

### ✅ Completed
- [x] Comprehensive exploration of codebase (auto_profile_discovery.py, tests, design doc)
- [x] Implementation plan created
- [x] Gurobi module loaded (11.0.2 available)
- [x] `discover_profiles` imports successfully
- [x] CITEgeist package structure verified

### ⚠️ In Progress
- [ ] **Environment Setup** - Python 3.12 detected (should be 3.10 from CITEgeist_env)
- [ ] pytest installation needed

### 🔜 Next Steps
1. Activate proper conda environment (CITEgeist_env with Python 3.10)
2. Install pytest in environment
3. Fix path in test_profile_discovery.py:385
4. Run unit tests
5. Manual exploration of Wu_rep_0 data
6. Run simulated data test (MAIN GOAL)

---

## Implementation Phases

### Phase 1: Environment Setup ⏳ IN PROGRESS

**Goal:** Set up Python 3.10 environment with all dependencies

**Critical Commands:**
```bash
# Load Gurobi
module load gurobi/12.0.3

# Activate correct conda environment
# NOTE: Having issues with conda activate, may need to use alternative method
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/

# OR try direct path
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"

# Verify Python 3.10
python --version  # Should show 3.10.x

# Install pytest if needed
pip install pytest

# Verify imports
python -c "import gurobipy; print(f'Gurobi {gurobipy.gurobi.version()}')"
python -c "from CITEgeist.model import discover_profiles; print('Import OK')"
```

**Current Issue:**
- Using Python 3.12 from base miniconda, need Python 3.10 from CITEgeist_env
- conda activate giving segmentation faults
- pytest not available in current PATH

**Resolution:**
- Try direct PATH export to CITEgeist_env/bin
- Or manually activate environment
- Install pytest once in correct environment

---

### Phase 2: Fix Critical Path Blockers 📋 PENDING

**File to Edit:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests/test_profile_discovery.py`

**Change at Line 385:**
```python
# FROM:
cite_path = Path("replicates/high_seg/h5ad_objects/Wu_rep_0_CITE.h5ad")

# TO:
cite_path = Path(__file__).parent.parent / "replicates" / "high_seg" / "h5ad_objects" / "Wu_rep_0_CITE.h5ad"
```

**Verify data exists:**
```bash
ls -lh /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_0_CITE.h5ad
```

---

### Phase 3: Validate Core Components 🧪 PENDING

**Unit Tests to Run:**
```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Test 1: Standardization
pytest tests/test_profile_discovery.py::TestStandardizeMarkers -v

# Test 2: Scoring
pytest tests/test_profile_discovery.py::TestScoreCandidate -v

# Test 3: Profile matrices
pytest tests/test_profile_discovery.py::TestBuildProfileMatrix -v

# Test 4: EM refinement
pytest tests/test_profile_discovery.py::TestEMRefine -v

# Test 5: Basic discovery
pytest tests/test_profile_discovery.py::TestDiscoverProfiles::test_single_dominant_population -v
pytest tests/test_profile_discovery.py::TestDiscoverProfiles::test_reproducibility -v
```

**Expected:** All pass without errors

**If failures:**
- Import errors → Check scipy, numpy versions
- Assertion failures → May need to adjust fixture signal strength
- Type errors → Check Python 3.10 compatibility

---

### Phase 4: Simulated Data Test (MAIN GOAL) 🎯 PENDING

#### Step 4.1: Manual Exploration

**Run this script first:**
```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
python3 << 'EOF'
import scanpy as sc
import numpy as np
from CITEgeist.model import discover_profiles

# Load data
print("Loading Wu_rep_0_CITE.h5ad...")
adata = sc.read_h5ad("replicates/high_seg/h5ad_objects/Wu_rep_0_CITE.h5ad")

print(f"\n=== Data Structure ===")
print(f"Shape: {adata.shape} (spots × markers)")
print(f"Data type: {type(adata.X)}")

print(f"\n=== Markers ({len(adata.var_names)}) ===")
for i, m in enumerate(adata.var_names):
    print(f"  {i:2d}: {m}")

# Check for spatial coordinates
coords = None
if 'spatial' in adata.obsm:
    coords = adata.obsm['spatial']
    print(f"\n=== Spatial Coordinates ===")
    print(f"From: adata.obsm['spatial'], Shape: {coords.shape}")
elif 'spot_x' in adata.obs and 'spot_y' in adata.obs:
    coords = adata.obs[['spot_x', 'spot_y']].values
    print(f"\n=== Spatial Coordinates ===")
    print(f"From: adata.obs[['spot_x', 'spot_y']], Shape: {coords.shape}")
else:
    print("\n=== No Spatial Coordinates Found ===")

# Basic statistics
print(f"\n=== Marker Statistics ===")
X_dense = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
for i, name in enumerate(adata.var_names):
    mean = X_dense[:, i].mean()
    std = X_dense[:, i].std()
    max_val = X_dense[:, i].max()
    print(f"  {name:20s}: mean={mean:6.2f}, std={std:6.2f}, max={max_val:6.2f}")

print("\n=== Running Discovery (verbose mode) ===")
result = discover_profiles(
    adata.X,
    list(adata.var_names),
    max_k=2,          # Start with pairs
    seed=1234,
    n_perm=100,       # Faster for testing
    alpha=0.1,        # More permissive
    coords=coords,
    verbose=True
)

print(f"\n=== Discovery Results ===")
print(f"Found {len(result.profiles)} profiles:")
for name, profile in result.profiles.items():
    markers = profile['Major']
    beta_str = ', '.join([f"{result.beta[m]:.2f}" for m in markers])
    print(f"  {name}: {markers}")
    print(f"    β: [{beta_str}]")

print(f"\nBIC trace: {result.bic_trace}")
print(f"Iterations: {result.n_iterations}")

# Check for background markers
background = [name for name in result.profiles.keys() if "Nonspecific" in name]
if background:
    print(f"\n⚠️  WARNING: Found {len(background)} background markers:")
    for bg in background:
        print(f"    - {bg}")
else:
    print(f"\n✅ No 'Nonspecific' markers in discovered profiles")
EOF
```

**Expected Output:**
- Data loads successfully
- Discovers 2-8 profiles
- No "Nonspecific_*" markers in profiles
- Clear biological interpretation
- BIC trace shows convergence

#### Step 4.2: Run Pytest

Once manual exploration succeeds:
```bash
pytest tests/test_profile_discovery.py::TestDiscoverProfiles::test_simulated_h5ad_background_filtered -v --tb=long
```

**Success Criteria:**
- ✅ Test passes
- ✅ Discovers ≥1 profile
- ✅ No background markers
- ✅ Assertion passes

---

### Phase 5: Comprehensive Testing (OPTIONAL)

After main goal achieved:
```bash
# All discovery tests
pytest tests/test_profile_discovery.py::TestDiscoverProfiles -v

# Edge cases
pytest tests/test_profile_discovery.py -v -k "edge or sparse or reproducibility"
```

---

## File Locations

### Implementation Files (Already Complete)
- `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/model/auto_profile_discovery.py` (652 lines)
- `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/model/__init__.py` (exports configured)

### Test Files
- `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests/test_profile_discovery.py` (456 lines)
- Needs path fix at line 385

### Data Files
- `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/replicates/high_seg/h5ad_objects/Wu_rep_0_CITE.h5ad`
- Multiple replicates available (Wu_rep_0 through Wu_rep_4)

### Documentation
- `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/docs/Design_AutoProfileDiscovery.md` (545 lines)
- `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CLAUDE.md` (needs environment setup section added)

---

## Known Issues & Solutions

### Issue: Python 3.12 instead of 3.10
**Solution:** Ensure CITEgeist_env conda environment is properly activated
```bash
# Try direct PATH modification
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
```

### Issue: pytest not found
**Solution:** Install in environment
```bash
pip install pytest
# Or if in conda env
conda install pytest
```

### Issue: Found 0 profiles
**Solution:** Increase signal strength or lower alpha
```python
# In fixtures (line ~50): Change += 3.0 to += 5.0
# In discover_profiles call: alpha=0.1 instead of 0.05
```

### Issue: Remote I/O errors on /ix1
**Solution:** This is known filesystem issue, retry commands or use /ihome/ mount

### Issue: Conda activation segfault
**Solution:** Use direct PATH export or check conda installation

---

## Success Metrics

### Primary Goal (Must Achieve)
- [x] Environment loads successfully
- [x] Imports work: `from CITEgeist.model import discover_profiles`
- [ ] Unit tests pass
- [ ] Manual exploration discovers profiles on Wu_rep_0_CITE.h5ad
- [ ] `test_simulated_h5ad_background_filtered` PASSES
- [ ] No background markers in discovered profiles
- [ ] **Math formulation validated on simulated data**

### Secondary Goals (Nice to Have)
- [ ] All integration tests pass
- [ ] Edge case tests pass
- [ ] Reproducibility verified
- [ ] Sparse matrix handling confirmed

---

## Troubleshooting Quick Reference

| Symptom | Likely Cause | Solution |
|---------|-------------|----------|
| pytest not found | Wrong environment | `pip install pytest` |
| Import error | Package not installed | `pip install -e .` |
| Wrong Python version | Conda env not active | Export PATH to CITEgeist_env/bin |
| Gurobi license error | Module not loaded | `module load gurobi/12.0.3` |
| File not found | Wrong working directory | cd to /ix1/.../CITEgeist |
| Found 0 profiles | Weak signal | Increase signal or lower alpha |
| Too many profiles | Overfitting | Check BIC trace |
| Remote I/O error | Filesystem issue | Retry or use /ihome mount |

---

## Next Session Checklist

When resuming:

1. ✅ Navigate to project directory
   ```bash
   cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
   ```

2. ✅ Load Gurobi
   ```bash
   module load gurobi/12.0.3
   ```

3. ✅ Verify Python environment
   ```bash
   python --version  # Should be 3.10.x
   which python      # Should be in CITEgeist_env
   ```

4. ✅ Test imports
   ```bash
   python -c "from CITEgeist.model import discover_profiles; print('OK')"
   ```

5. 🔜 Continue with Phase 2 (fix path in test file)

6. 🔜 Run unit tests (Phase 3)

7. 🔜 Manual exploration (Phase 4.1)

8. 🔜 Run simulated data test (Phase 4.2) - **MAIN GOAL**

---

## Contact & Resources

- **Implementation Plan:** /ihome/alee/alc376/.claude/plans/recursive-skipping-pancake.md
- **This File:** /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/IMPLEMENTATION_PLAN_AutoProfileDiscovery.md
- **Design Doc:** docs/Design_AutoProfileDiscovery.md
- **Main Code:** CITEgeist/model/auto_profile_discovery.py
- **Tests:** tests/test_profile_discovery.py

---

**Last Updated:** 2025-12-05
**Session Status:** Paused after Phase 1 environment verification
**Ready to Resume:** Phase 2 path fixes
