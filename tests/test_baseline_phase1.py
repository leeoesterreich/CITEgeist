"""PHASE1 switch in protein_sensitivity.baseline: unset/qp keeps the pinned QP call byte-for-byte; marker_share
calls run_cell_proportion_model(method="marker_share") with no QP-only kwargs."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402

PROFILE = {"A": {"Major": ["m0"], "Minor": []}, "B": {"Major": ["m1"], "Minor": []}}
NUC = pd.Series({"spot_1": 3, "spot_2": 4})


class FakeModel:
    last = None

    def __init__(self, **kw):
        FakeModel.last = self
        self.run_kwargs = None
        self.results = {}
        self.gene_expression_adata = kw["gene_expression_adata"]

    def filter_gex(self, **kw):
        pass

    def preprocess_gex(self, **kw):
        pass

    def preprocess_antibody(self, **kw):
        pass

    def load_cell_profile_dict(self, p):
        pass

    def run_cell_proportion_model(self, **kw):
        self.run_kwargs = kw
        df = pd.DataFrame({t: [0.5, 0.5] for t in PROFILE}, index=["spot_1", "spot_2"])
        self.results["cell_prop"] = df
        return df, df.copy()


@pytest.fixture
def fake(monkeypatch):
    import CITEgeist.model.citegeist_model as cm

    monkeypatch.setattr(cm, "CitegeistModel", FakeModel)
    monkeypatch.setattr(B, "load_nuclei_counts", lambda region: NUC)
    monkeypatch.delenv("PHASE1", raising=False)
    FakeModel.last = None


def _adatas():
    g = ad.AnnData(X=np.ones((2, 2)), obs=pd.DataFrame(index=["spot_1", "spot_2"]))
    return g, g.copy()


@pytest.mark.unit
def test_xenium_qp_default_is_pinned_qp(fake):
    g, c = _adatas()
    B.run_xenium_qp(g, c, PROFILE, 0)
    kw = FakeModel.last.run_kwargs
    assert kw["method"] == "qp" and kw["lambda_laplacian"] == 0.0 and kw["nuclei_counts"] is NUC


@pytest.mark.unit
def test_xenium_qp_marker_share(fake):
    g, c = _adatas()
    out = B.run_xenium_qp(g, c, PROFILE, 0, phase1="marker_share")
    assert FakeModel.last.run_kwargs == {"method": "marker_share"}
    assert list(out.columns) == list(PROFILE)


@pytest.mark.unit
def test_xenium_qp_env_var(fake, monkeypatch):
    monkeypatch.setenv("PHASE1", "marker_share")
    g, c = _adatas()
    B.run_xenium_qp(g, c, PROFILE, 0)
    assert FakeModel.last.run_kwargs == {"method": "marker_share"}


@pytest.mark.unit
def test_unknown_phase1_is_loud(fake):
    g, c = _adatas()
    with pytest.raises(ValueError):
        B.run_xenium_qp(g, c, PROFILE, 0, phase1="bogus")


@pytest.mark.unit
def test_sim_default_and_marker_share(fake):
    g, c = _adatas()
    B.run_sim(g, c, PROFILE, sace=False)
    assert FakeModel.last.run_kwargs == {"method": "qp"}
    B.run_sim(g, c, PROFILE, sace=False, phase1="marker_share")
    assert FakeModel.last.run_kwargs == {"method": "marker_share"}


@pytest.mark.unit
def test_sim_v41p_turns_on_brightness_only(fake):
    g, c = _adatas()
    B.run_sim(g, c, PROFILE, sace=False, phase1="marker_share_v41p")
    assert FakeModel.last.run_kwargs == {"method": "marker_share", "brightness_correction": True}
    B.run_sim(g, c, PROFILE, sace=False, phase1="marker_share")
    assert FakeModel.last.run_kwargs == {"method": "marker_share"}


@pytest.mark.unit
def test_marker_share_bright_passes_env_method_and_panel(fake, monkeypatch):
    monkeypatch.setenv("BRIGHTNESS_METHOD", "rna_slope")
    g, c = _adatas()
    prof = {"Epithelial": {"Major": ["m0"], "Minor": []}, "B cells": {"Major": ["m1"], "Minor": []}}
    B.run_xenium_qp(g, c, prof, 0, phase1="marker_share_bright")
    kw = FakeModel.last.run_kwargs
    assert kw["method"] == "marker_share" and kw["brightness_method"] == "rna_slope"
    assert kw["rna_panel"]["Epithelial"] == ["KRT7", "MUC1", "CDH1", "EHF"] and "B cells" in kw["rna_panel"]
    B.run_sim(g, c, PROFILE, sace=False, phase1="marker_share_bright")
    kw = FakeModel.last.run_kwargs
    assert kw["brightness_method"] == "rna_slope" and "Cancer Epithelial" in kw["rna_panel"]


@pytest.mark.unit
def test_marker_share_bright_needs_env_method(fake, monkeypatch):
    monkeypatch.delenv("BRIGHTNESS_METHOD", raising=False)
    g, c = _adatas()
    with pytest.raises(ValueError, match="BRIGHTNESS_METHOD"):
        B.run_xenium_qp(g, c, PROFILE, 0, phase1="marker_share_bright")
