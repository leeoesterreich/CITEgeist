"""Tests for discrete cell assignment with morphology nudge."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest
import torch

from CITEgeist.model.assignment.cell_assignment import (
    discretize_proportions,
    hungarian_assign_spot,
    assign_cells_to_types,
    extract_embeddings,
    fit_morphology_scores,
    assign_cells,
    bayesian_assign_cells,
)
from CITEgeist.model.morphology.prototype_contrastive import PrototypeContrastiveModel


class TestPrototypeContrastiveEmbeddingsMode:
    """Test that PrototypeContrastiveModel accepts precomputed embeddings."""

    def test_from_embeddings_forward(self):
        """Model with from_embeddings=True accepts (B, 384) tensors."""
        model = PrototypeContrastiveModel(
            num_types=6,
            embed_dim=128,
            from_embeddings=True,
        )
        model.eval()
        embeddings = torch.randn(10, 384)
        result = model(embeddings)
        assert result["q_proto"].shape == (10, 6)
        assert result["q_class"].shape == (10, 6)
        assert result["z"].shape == (10, 128)

    def test_from_embeddings_no_encoder(self):
        """Model with from_embeddings=True should not have a ViT encoder."""
        model = PrototypeContrastiveModel(
            num_types=6,
            embed_dim=128,
            from_embeddings=True,
        )
        assert model.encoder is None

    def test_from_embeddings_flag_stored(self):
        """from_embeddings attribute is stored on the model."""
        model = PrototypeContrastiveModel(
            num_types=4,
            embed_dim=64,
            from_embeddings=True,
        )
        assert model.from_embeddings is True

    def test_default_mode_has_encoder(self):
        """Default model (from_embeddings=False) still has a ViT encoder."""
        model = PrototypeContrastiveModel(
            num_types=4,
            embed_dim=64,
            from_embeddings=False,
        )
        assert model.encoder is not None
        assert model.from_embeddings is False

    def test_freeze_encoder_noop_when_from_embeddings(self):
        """freeze_encoder() is a no-op (no AttributeError) when from_embeddings=True."""
        model = PrototypeContrastiveModel(
            num_types=4,
            embed_dim=64,
            from_embeddings=True,
        )
        # Should not raise
        model.freeze_encoder()

    def test_unfreeze_last_n_blocks_noop_when_from_embeddings(self):
        """unfreeze_last_n_blocks() is a no-op when from_embeddings=True."""
        model = PrototypeContrastiveModel(
            num_types=4,
            embed_dim=64,
            from_embeddings=True,
        )
        # Should not raise
        model.unfreeze_last_n_blocks(2)

    def test_output_probabilities_sum_to_one(self):
        """Softmax outputs should sum to 1 across types."""
        model = PrototypeContrastiveModel(
            num_types=5,
            embed_dim=64,
            from_embeddings=True,
        )
        model.eval()
        embeddings = torch.randn(8, 384)
        result = model(embeddings)
        # q_proto and q_class are softmax outputs — must sum to 1
        torch.testing.assert_close(
            result["q_proto"].sum(dim=-1),
            torch.ones(8),
            atol=1e-5,
            rtol=1e-5,
        )
        torch.testing.assert_close(
            result["q_class"].sum(dim=-1),
            torch.ones(8),
            atol=1e-5,
            rtol=1e-5,
        )


class TestDiscretizeProportions:
    """Test Stage 1: proportion -> integer count rounding."""

    def test_basic_rounding(self):
        """Proportions * nuclei_counts -> integer counts summing to N."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.5, 0.3], "T_Cell": [0.3, 0.5], "B_Cell": [0.2, 0.2]},
            index=["spot_0", "spot_1"],
        )
        nuclei_counts = pd.Series([10, 5], index=["spot_0", "spot_1"])

        int_counts = discretize_proportions(cell_prop, nuclei_counts)

        assert int_counts.shape == (2, 3)
        assert (int_counts.sum(axis=1).values == [10, 5]).all()
        assert (int_counts.values >= 0).all()

    def test_single_cell_spot(self):
        """Spot with 1 nucleus gets exactly one cell type assigned."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.6], "T_Cell": [0.4]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([1], index=["spot_0"])

        int_counts = discretize_proportions(cell_prop, nuclei_counts)

        assert int_counts.sum(axis=1).values[0] == 1
        assert int_counts.values[0, 0] == 1  # Cancer gets the 1 cell (higher proportion)

    def test_zero_nuclei_spot(self):
        """Spot with 0 nuclei returns all zeros."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.5], "T_Cell": [0.5]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([0], index=["spot_0"])

        int_counts = discretize_proportions(cell_prop, nuclei_counts)

        assert int_counts.sum(axis=1).values[0] == 0

    def test_reuses_existing_integer_counts(self):
        """When existing_integer_counts provided, returns them directly."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.5], "T_Cell": [0.5]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([10], index=["spot_0"])
        existing = pd.DataFrame(
            {"Cancer": [7], "T_Cell": [3]},
            index=["spot_0"],
        )

        int_counts = discretize_proportions(
            cell_prop, nuclei_counts, existing_integer_counts=existing
        )

        pd.testing.assert_frame_equal(int_counts, existing)


class TestHungarianAssignSpot:
    """Test per-spot Hungarian assignment."""

    def test_no_morphology_respects_counts(self):
        """Without morphology, assigns correct number of each type."""
        counts = np.array([3, 2, 1])  # 3 Cancer, 2 T_Cell, 1 B_Cell
        cell_types = ["Cancer", "T_Cell", "B_Cell"]
        n_cells = 6

        assignments = hungarian_assign_spot(
            n_cells=n_cells,
            integer_counts=counts,
            cell_types=cell_types,
            morph_scores=None,
            morphology_weight=0.0,
        )

        assert len(assignments) == 6
        type_counts = pd.Series(assignments).value_counts()
        assert type_counts.get("Cancer", 0) == 3
        assert type_counts.get("T_Cell", 0) == 2
        assert type_counts.get("B_Cell", 0) == 1

    def test_morphology_nudge_changes_assignment(self):
        """Morphology scores change which specific cells get which type."""
        counts = np.array([1, 1])
        cell_types = ["Cancer", "T_Cell"]
        morph_scores = np.array([
            [0.9, 0.1],  # cell 0: strong Cancer signal
            [0.1, 0.9],  # cell 1: strong T_Cell signal
        ])

        assignments = hungarian_assign_spot(
            n_cells=2,
            integer_counts=counts,
            cell_types=cell_types,
            morph_scores=morph_scores,
            morphology_weight=1.0,
        )

        assert assignments[0] == "Cancer"
        assert assignments[1] == "T_Cell"

    def test_single_cell_spot(self):
        counts = np.array([1, 0])
        cell_types = ["Cancer", "T_Cell"]
        assignments = hungarian_assign_spot(1, counts, cell_types, None, 0.0)
        assert len(assignments) == 1
        assert assignments[0] == "Cancer"

    def test_deterministic_no_morphology(self):
        counts = np.array([2, 2])
        cell_types = ["Cancer", "T_Cell"]
        a1 = hungarian_assign_spot(4, counts, cell_types, None, 0.0)
        a2 = hungarian_assign_spot(4, counts, cell_types, None, 0.0)
        assert a1 == a2


class TestAssignCellsToTypes:
    """Test full Stage 3: multi-spot assignment."""

    def test_basic_assignment(self):
        integer_counts = pd.DataFrame(
            {"Cancer": [2, 1], "T_Cell": [1, 2]},
            index=["spot_0", "spot_1"],
        )
        cell_to_spot = np.array([0, 0, 0, 1, 1, 1])

        result = assign_cells_to_types(
            integer_counts=integer_counts,
            cell_to_spot=cell_to_spot,
            cell_types=["Cancer", "T_Cell"],
            morph_scores=None,
            morphology_weight=0.0,
            spot_proportions=pd.DataFrame(
                {"Cancer": [0.67, 0.33], "T_Cell": [0.33, 0.67]},
                index=["spot_0", "spot_1"],
            ),
        )

        assert len(result) == 6
        assert "assigned_type" in result.columns
        assert "confidence" in result.columns
        assert "spot_id" in result.columns
        spot0 = result[result["spot_id"] == "spot_0"]
        assert (spot0["assigned_type"] == "Cancer").sum() == 2
        assert (spot0["assigned_type"] == "T_Cell").sum() == 1


class TestEdgeCases:
    """Edge case tests for robustness."""

    def test_all_one_type(self):
        """Spot where one type has 100% proportion."""
        cell_prop = pd.DataFrame(
            {"Cancer": [1.0], "T_Cell": [0.0]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([5], index=["spot_0"])
        cell_to_spot = np.array([0, 0, 0, 0, 0])

        result = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        assert (result["assigned_type"] == "Cancer").all()

    def test_many_types_few_cells(self):
        """More types than cells in a spot — some types get 0."""
        cell_prop = pd.DataFrame(
            {"A": [0.3], "B": [0.3], "C": [0.2], "D": [0.1], "E": [0.1]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([2], index=["spot_0"])
        cell_to_spot = np.array([0, 0])

        result = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        assert len(result) == 2
        assert set(result["assigned_type"]).issubset({"A", "B", "C", "D", "E"})

    def test_empty_spot_no_cells(self):
        """Spot exists in cell_prop but no cells map to it."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.5, 0.5], "T_Cell": [0.5, 0.5]},
            index=["spot_0", "spot_1"],
        )
        nuclei_counts = pd.Series([3, 3], index=["spot_0", "spot_1"])
        cell_to_spot = np.array([0, 0, 0])  # all in spot 0, none in spot 1

        result = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        assert len(result) == 3
        assert (result["spot_id"] == "spot_0").all()

    def test_reconcile_cell_count_mismatch(self):
        """When cells_in_spot != nuclei_counts, re-rounds to actual cell count."""
        integer_counts = pd.DataFrame(
            {"Cancer": [5], "T_Cell": [5]},  # 10 total
            index=["spot_0"],
        )
        cell_to_spot = np.array([0, 0, 0, 0, 0, 0, 0, 0])  # 8 actual cells

        result = assign_cells_to_types(
            integer_counts=integer_counts,
            cell_to_spot=cell_to_spot,
            cell_types=["Cancer", "T_Cell"],
            morph_scores=None,
            morphology_weight=0.0,
            spot_proportions=pd.DataFrame(
                {"Cancer": [0.5], "T_Cell": [0.5]},
                index=["spot_0"],
            ),
        )

        assert len(result) == 8
        assert (result["assigned_type"] == "Cancer").sum() == 4
        assert (result["assigned_type"] == "T_Cell").sum() == 4


class TestFitMorphologyScores:
    """Test Stage 2: morphology scoring via prototype-contrastive LLP."""

    def test_returns_correct_shape(self):
        """Morphology scores are (C, T) probabilities."""
        np.random.seed(42)
        n_cells = 100
        n_types = 4
        embeddings = np.random.randn(n_cells, 384).astype(np.float32)
        cell_to_spot = np.repeat(np.arange(10), 10)
        oracle_props = np.random.dirichlet(np.ones(n_types), size=10).astype(np.float32)

        scores = fit_morphology_scores(
            embeddings=embeddings,
            cell_to_spot=cell_to_spot,
            oracle_props=oracle_props,
            num_types=n_types,
            n_epochs=5,  # minimal for test speed
            device="cpu",
        )

        assert scores.shape == (n_cells, n_types)
        # Probabilities sum to ~1 per cell
        row_sums = scores.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=0.01)
        # All non-negative
        assert (scores >= 0).all()


class TestExtractEmbeddings:
    """Test embedding extraction from patches with caching."""

    def test_extract_from_patches(self):
        """Extracts (C, 384) embeddings from (C, 3, 96, 96) patches."""
        np.random.seed(42)
        patches = np.random.randn(10, 3, 96, 96).astype(np.float32)

        with tempfile.TemporaryDirectory() as tmpdir:
            embeddings = extract_embeddings(
                patches=patches,
                encoder_checkpoint=None,
                output_folder=tmpdir,
                sample_name="test",
                device="cpu",
            )

            assert embeddings.shape == (10, 384)
            assert embeddings.dtype == np.float32

    def test_caching(self):
        """Second call loads from cache instead of re-extracting."""
        np.random.seed(42)
        patches = np.random.randn(5, 2, 96, 96).astype(np.float32)

        with tempfile.TemporaryDirectory() as tmpdir:
            e1 = extract_embeddings(
                patches=patches,
                encoder_checkpoint=None,
                output_folder=tmpdir,
                sample_name="test",
                device="cpu",
            )
            cache_files = [f for f in os.listdir(tmpdir) if f.startswith("morphology_embeddings_")]
            assert len(cache_files) == 1

            e2 = extract_embeddings(
                patches=patches,
                encoder_checkpoint=None,
                output_folder=tmpdir,
                sample_name="test",
                device="cpu",
            )
            np.testing.assert_array_equal(e1, e2)

    def test_rejects_invalid_patch_shape(self):
        """Rejects patches that aren't (C, ch, 96, 96) with ch in {2, 3}."""
        patches_bad_spatial = np.random.randn(5, 3, 224, 224).astype(np.float32)

        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match="96x96"):
                extract_embeddings(
                    patches=patches_bad_spatial,
                    encoder_checkpoint=None,
                    output_folder=tmpdir,
                    sample_name="test",
                    device="cpu",
                )


class TestAssignCellsFullPipeline:
    """Test the top-level assign_cells() orchestrator."""

    def test_no_morphology_pipeline(self):
        """Full pipeline without morphology: discretize -> assign."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.6, 0.2], "T_Cell": [0.4, 0.8]},
            index=["spot_0", "spot_1"],
        )
        nuclei_counts = pd.Series([5, 5], index=["spot_0", "spot_1"])
        cell_to_spot = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

        result = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        assert len(result) == 10
        assert "spot_id" in result.columns
        assert "cell_id" in result.columns
        assert "assigned_type" in result.columns
        assert "confidence" in result.columns
        assert "Cancer" in result.columns
        assert "T_Cell" in result.columns

        # Spot 0: 3 Cancer + 2 T_Cell (0.6 * 5 = 3, 0.4 * 5 = 2)
        spot0 = result[result["spot_id"] == "spot_0"]
        assert (spot0["assigned_type"] == "Cancer").sum() == 3
        assert (spot0["assigned_type"] == "T_Cell").sum() == 2

    def test_with_morphology_embeddings(self):
        """Full pipeline with morphology: discretize -> score -> assign."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.5], "T_Cell": [0.5]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([4], index=["spot_0"])
        cell_to_spot = np.array([0, 0, 0, 0])
        embeddings = np.random.randn(4, 384).astype(np.float32)

        result = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
            morphology_embeddings=embeddings,
            morphology_weight=0.5,
            n_morph_epochs=5,
            device="cpu",
        )

        assert len(result) == 4
        assert "assigned_type" in result.columns
        # 2 Cancer + 2 T_Cell (0.5 * 4 = 2 each)
        assert (result["assigned_type"] == "Cancer").sum() == 2
        assert (result["assigned_type"] == "T_Cell").sum() == 2

    def test_morphology_weight_zero_matches_no_morphology(self):
        """morphology_weight=0 produces same result as no morphology."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.6], "T_Cell": [0.4]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([5], index=["spot_0"])
        cell_to_spot = np.array([0, 0, 0, 0, 0])

        r_no_morph = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        r_zero_weight = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
            morphology_embeddings=np.random.randn(5, 384).astype(np.float32),
            morphology_weight=0.0,
            n_morph_epochs=5,
            device="cpu",
        )

        assert list(r_no_morph["assigned_type"]) == list(r_zero_weight["assigned_type"])

    def test_validates_cell_to_spot_alignment(self):
        """Raises ValueError if cell_to_spot indices exceed cell_prop."""
        cell_prop = pd.DataFrame(
            {"Cancer": [0.5], "T_Cell": [0.5]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([2], index=["spot_0"])
        cell_to_spot = np.array([0, 5])  # index 5 is out of range

        with pytest.raises(ValueError, match="cell_to_spot"):
            assign_cells(
                cell_prop=cell_prop,
                nuclei_counts=nuclei_counts,
                cell_to_spot=cell_to_spot,
            )


class TestCitegeistModelAssignCells:
    """Test CitegeistModel.assign_cells() wrapper."""

    def _make_model_with_results(self):
        """Create a model-like object with results already populated."""
        from CITEgeist.model.citegeist_model import CitegeistModel
        model = CitegeistModel.__new__(CitegeistModel)
        model.results = {
            "cell_prop": pd.DataFrame(
                {"Cancer": [0.6, 0.3], "T_Cell": [0.4, 0.7]},
                index=["spot_0", "spot_1"],
            ),
        }
        model.output_folder = "/tmp/test_citegeist"
        model.sample_name = "test_sample"
        return model

    def test_wrapper_stores_results(self):
        """Wrapper stores cell_assignments, cell_types_hard in self.results."""
        model = self._make_model_with_results()
        nuclei_counts = pd.Series([3, 3], index=["spot_0", "spot_1"])
        cell_to_spot = np.array([0, 0, 0, 1, 1, 1])

        result = model.assign_cells(
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        assert isinstance(result, pd.DataFrame)
        assert "cell_assignments" in model.results
        assert "cell_types_hard" in model.results
        assert len(model.results["cell_types_hard"]) == 6

    def test_raises_without_proportions(self):
        """Raises if run_cell_proportion_model() wasn't called."""
        from CITEgeist.model.citegeist_model import CitegeistModel
        model = CitegeistModel.__new__(CitegeistModel)
        model.results = {}

        with pytest.raises(RuntimeError, match="run_cell_proportion_model"):
            model.assign_cells(
                nuclei_counts=pd.Series([5]),
                cell_to_spot=np.array([0, 0, 0, 0, 0]),
            )

    def test_reuses_existing_integer_counts(self):
        """Wrapper passes existing cell_counts_integer if available."""
        model = self._make_model_with_results()
        model.results["cell_counts_integer"] = pd.DataFrame(
            {"Cancer": [2, 1], "T_Cell": [1, 2]},
            index=["spot_0", "spot_1"],
        )
        nuclei_counts = pd.Series([3, 3], index=["spot_0", "spot_1"])
        cell_to_spot = np.array([0, 0, 0, 1, 1, 1])

        result = model.assign_cells(
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
        )

        spot0 = result[result["spot_id"] == "spot_0"]
        assert (spot0["assigned_type"] == "Cancer").sum() == 2


class TestDetectionMaskPersistence:
    """Test that detection_mask_bool types and shapes are correct."""

    def test_detection_mask_bool_dtype_and_shape(self):
        """Boolean detection mask has correct dtype and shape."""
        detection = np.array([[True, False, True], [True, True, False]])
        assert detection.dtype == bool
        assert detection.shape == (2, 3)
        # Cast to float works correctly
        detection_float = detection.astype(np.float64)
        assert detection_float[0, 1] == 0.0
        assert detection_float[1, 0] == 1.0

    def test_rescue_all_false_row(self):
        """When all types masked for a spot, rescue should set all True."""
        detection = np.array([[False, False, False], [True, False, True]])
        # Simulate rescue logic
        for i in range(detection.shape[0]):
            if detection[i].sum() < 2:
                detection[i] = True
        assert detection[0].all()  # rescued
        assert detection[1, 0] and not detection[1, 1] and detection[1, 2]  # unchanged


class TestBayesianAssignCells:
    """Test Bayesian posterior cell assignment."""

    def test_basic_posterior(self):
        """Cells assigned to highest posterior type."""
        morph = np.array([
            [0.8, 0.1, 0.1],
            [0.1, 0.8, 0.1],
            [0.1, 0.1, 0.8],
        ])
        cell_to_spot = np.array([0, 0, 0])
        prior = np.array([[0.4, 0.3, 0.3]])
        mask = np.array([[True, True, True]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B", "C"],
            spot_ids=["spot_0"],
        )

        assert len(result) == 3
        assert result.iloc[0]["assigned_type"] == "A"
        assert result.iloc[1]["assigned_type"] == "B"
        assert result.iloc[2]["assigned_type"] == "C"

    def test_mask_eliminates_types(self):
        """Masked types get zero posterior regardless of morphology."""
        morph = np.array([
            [0.8, 0.1, 0.1],
        ])
        cell_to_spot = np.array([0])
        prior = np.array([[0.4, 0.3, 0.3]])
        mask = np.array([[False, True, True]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B", "C"],
            spot_ids=["spot_0"],
        )

        assert result.iloc[0]["assigned_type"] in ("B", "C")
        assert result.iloc[0]["A"] == 0.0

    def test_prior_weights_assignment(self):
        """When morphology is uniform, prior determines assignment."""
        morph = np.array([
            [0.33, 0.33, 0.34],
        ])
        cell_to_spot = np.array([0])
        prior = np.array([[0.9, 0.05, 0.05]])
        mask = np.array([[True, True, True]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B", "C"],
            spot_ids=["spot_0"],
        )

        assert result.iloc[0]["assigned_type"] == "A"

    def test_all_masked_fallback(self):
        """When all types masked, falls back to uniform."""
        morph = np.array([[0.5, 0.5]])
        cell_to_spot = np.array([0])
        prior = np.array([[0.5, 0.5]])
        mask = np.array([[False, False]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B"],
            spot_ids=["spot_0"],
        )

        assert len(result) == 1
        assert np.isfinite(result.iloc[0]["confidence"])

    def test_single_type_detected(self):
        """Only 1 type detected -> one-hot assignment."""
        morph = np.array([[0.3, 0.7]])
        cell_to_spot = np.array([0])
        prior = np.array([[0.5, 0.5]])
        mask = np.array([[True, False]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B"],
            spot_ids=["spot_0"],
        )

        assert result.iloc[0]["assigned_type"] == "A"
        assert result.iloc[0]["confidence"] == pytest.approx(1.0, abs=0.01)

    def test_no_count_constraints(self):
        """Bayesian does NOT enforce count budgets — all cells can be same type."""
        morph = np.array([
            [0.9, 0.1],
            [0.8, 0.2],
            [0.7, 0.3],
        ])
        cell_to_spot = np.array([0, 0, 0])
        prior = np.array([[0.5, 0.5]])
        mask = np.array([[True, True]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B"],
            spot_ids=["spot_0"],
        )

        assert (result["assigned_type"] == "A").sum() == 3

    def test_output_columns(self):
        """Output has correct columns."""
        morph = np.array([[0.5, 0.5]])
        cell_to_spot = np.array([0])
        prior = np.array([[0.5, 0.5]])
        mask = np.array([[True, True]])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B"],
            spot_ids=["spot_0"],
            cell_ids=np.array(["cell_0"]),
        )

        assert "spot_id" in result.columns
        assert "cell_id" in result.columns
        assert "assigned_type" in result.columns
        assert "confidence" in result.columns
        assert "A" in result.columns
        assert "B" in result.columns

    def test_multi_spot(self):
        """Works across multiple spots with different masks."""
        morph = np.array([
            [0.8, 0.1, 0.1],
            [0.1, 0.8, 0.1],
        ])
        cell_to_spot = np.array([0, 1])
        prior = np.array([[0.5, 0.3, 0.2], [0.1, 0.6, 0.3]])
        mask = np.array([
            [True, False, True],
            [False, True, True],
        ])

        result = bayesian_assign_cells(
            morph_scores=morph,
            cell_to_spot=cell_to_spot,
            proportion_prior=prior,
            detection_mask=mask,
            cell_types=["A", "B", "C"],
            spot_ids=["spot_0", "spot_1"],
        )

        assert result.iloc[0]["assigned_type"] == "A"
        assert result.iloc[1]["assigned_type"] == "B"


class TestAssignCellsBayesianDispatch:
    """Test assignment_method='bayesian' dispatch in assign_cells()."""

    def test_bayesian_dispatch(self):
        """assign_cells with assignment_method='bayesian' uses Bayesian posterior."""
        cell_prop = pd.DataFrame(
            {"A": [0.6], "B": [0.4]},
            index=["spot_0"],
        )
        nuclei_counts = pd.Series([3], index=["spot_0"])
        cell_to_spot = np.array([0, 0, 0])
        morph = np.array([
            [0.9, 0.1],
            [0.9, 0.1],
            [0.1, 0.9],
        ])
        mask = np.array([[True, True]])

        result = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
            assignment_method="bayesian",
            detection_mask=mask,
            proportion_prior=cell_prop.values,
            morph_scores_precomputed=morph,
        )

        assert len(result) == 3
        assert (result["assigned_type"] == "A").sum() == 2
        assert (result["assigned_type"] == "B").sum() == 1

    def test_bayesian_requires_morph(self):
        """Bayesian assignment raises without morphology scores."""
        cell_prop = pd.DataFrame({"A": [0.5]}, index=["spot_0"])
        mask = np.array([[True]])

        with pytest.raises(ValueError, match="[Mm]orphology"):
            assign_cells(
                cell_prop=cell_prop,
                nuclei_counts=pd.Series([1], index=["spot_0"]),
                cell_to_spot=np.array([0]),
                assignment_method="bayesian",
                detection_mask=mask,
            )

    def test_bayesian_requires_detection_mask(self):
        """Bayesian assignment raises without detection mask."""
        cell_prop = pd.DataFrame({"A": [0.5]}, index=["spot_0"])
        morph = np.array([[0.5]])

        with pytest.raises(ValueError, match="detection_mask"):
            assign_cells(
                cell_prop=cell_prop,
                nuclei_counts=pd.Series([1], index=["spot_0"]),
                cell_to_spot=np.array([0]),
                assignment_method="bayesian",
                morph_scores_precomputed=morph,
            )
