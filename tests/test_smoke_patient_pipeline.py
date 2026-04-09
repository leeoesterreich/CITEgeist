"""End-to-end smoke test for the CITEgeist patient pipeline.

Requires a real patient sample. Set CITEGEIST_SMOKE_SAMPLE to the path of
a combined .h5ad file (output of sq.read.visium + split_adata), e.g.:

    export CITEGEIST_SMOKE_SAMPLE=/ix1/alee/.../output/HCC22-088-P4-S2_combined.h5ad

Or create tests/fixtures/smoke_sample_path.txt with the path on one line.

Marked slow + requires_data; excluded from default CI (-m "not slow").
"""
import os
import pytest
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH_FILE = REPO_ROOT / "tests" / "fixtures" / "smoke_sample_path.txt"
N_SPOTS_SUBSET = 200


def _discover_sample_path():
    """Return path to smoke test sample, or None."""
    env_path = os.environ.get("CITEGEIST_SMOKE_SAMPLE")
    if env_path and Path(env_path).exists():
        return Path(env_path)
    if FIXTURE_PATH_FILE.exists():
        candidate = Path(FIXTURE_PATH_FILE.read_text().strip())
        if candidate.exists():
            return candidate
    return None


@pytest.fixture(scope="module")
def sample_path():
    path = _discover_sample_path()
    if path is None:
        pytest.skip(
            "Smoke test sample not found. Set CITEGEIST_SMOKE_SAMPLE env var "
            "or create tests/fixtures/smoke_sample_path.txt with the path."
        )
    return path


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.requires_data
def test_pipeline_m1_to_m35_runs_without_error(sample_path, tmp_path):
    """Full pipeline M1→M3→M3-post(Bayesian)→M3-gex(SACE)→M3.5 on 200-spot subset."""
    import scanpy as sc
    from CITEgeist.model.citegeist_model import CitegeistModel

    # Load + subset
    adata = sc.read_h5ad(sample_path)
    if adata.n_obs > N_SPOTS_SUBSET:
        adata = adata[:N_SPOTS_SUBSET].copy()

    # Instantiate model
    model = CitegeistModel(adata, simulation=False)

    # M1: marker interest detection
    m1_result = model.run_marker_interest()
    assert m1_result is not None, "M1 returned None"
    assert len(m1_result.interesting_markers) > 0, "M1 found no interesting markers"

    # M2: profile discovery (use top 2 markers only for speed)
    top_markers = m1_result.interesting_markers[:2]
    m2_result = model.run_profile_discovery(interesting_markers=top_markers)
    assert m2_result is not None, "M2 returned None"

    # M3: QP proportions (CPU fallback expected in test env; skip if cuOPT unavailable)
    try:
        prop_result = model.run_cell_proportion_model(profiles=m2_result.profiles)
    except RuntimeError as e:
        if "cuOPT" in str(e) or "GPU" in str(e):
            pytest.skip(f"cuOPT not available in test env: {e}")
        raise

    assert prop_result is not None, "M3 returned None"
    prop_df = prop_result.proportions
    assert prop_df.shape[0] == N_SPOTS_SUBSET
    # Proportions must sum to <= 1 per spot (with small tolerance)
    row_sums = prop_df.values.sum(axis=1)
    assert np.all(row_sums <= 1.05), f"Proportions sum > 1.05 in {(row_sums > 1.05).sum()} spots"
    assert np.all(prop_df.values >= -1e-6), "Negative proportions found"

    # M3-post: Bayesian cell assignment
    assign_result = model.run_cell_assignment(
        prop_result=prop_result, assignment_method="bayesian"
    )
    # Bayesian may fall back to hungarian if embeddings unavailable — both are OK
    assert assign_result is not None, "M3-post returned None"
    assert len(assign_result.cell_assignments) > 0, "No cells assigned"

    # M3-gex: SACE per-cell GEX
    sace_result = model.run_sace_gex(
        prop_result=prop_result,
        assign_result=assign_result,
        max_iter=1,
    )
    assert sace_result is not None, "M3-gex returned None"
    assert sace_result.adata is not None, "SACE adata is None"
    assert "cell_type" in sace_result.adata.obs.columns, "cell_type missing from SACE adata.obs"

    # M3.5: SACE protein + functional gating
    m35_result = model.run_sace_protein(prop_result=prop_result, assign_result=assign_result)
    assert m35_result is not None, "M3.5 returned None"
