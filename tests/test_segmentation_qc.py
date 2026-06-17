"""Tests for segmentation QC PDF generation."""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


class TestSegmentationQcPdf:
    def test_generates_pdf_file(self):
        from CITEgeist.model.morphology.segmentation import SegmentationQC
        from CITEgeist.model.morphology.segmentation_qc import (
            generate_segmentation_qc_pdf,
        )

        N_spots = 20
        image = np.random.randint(0, 255, (200, 200, 3), dtype=np.uint8)
        masks = np.zeros((200, 200), dtype=np.int32)
        for i in range(1, 31):
            r, c = np.random.randint(10, 190, 2)
            masks[r - 3 : r + 3, c - 3 : c + 3] = i
        centroids_xy = np.random.rand(30, 2) * 200
        spot_centers_xy = np.random.rand(N_spots, 2) * 200
        spot_radius_px = 15.0
        nuclei_count_raw = pd.Series(
            np.random.randint(0, 5, N_spots),
            index=[f"spot_{i}" for i in range(N_spots)],
        )
        qc = SegmentationQC(
            n_raw=35,
            n_after_area_filter=30,
            n_removed_small=3,
            n_removed_large=2,
            pixel_size_um=0.45,
            min_area_um2=10.0,
            max_area_um2=500.0,
            area_median_um2=80.0,
            area_p25_um2=60.0,
            area_p75_um2=100.0,
            density_per_mm2=2500.0,
            density_flag="ok",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            pdf_path = generate_segmentation_qc_pdf(
                output_folder=tmpdir,
                sample_name="test_sample",
                image=image,
                masks=masks,
                centroids_xy=centroids_xy,
                spot_centers_xy=spot_centers_xy,
                spot_radius_px=spot_radius_px,
                nuclei_count_raw=nuclei_count_raw,
                qc=qc,
            )

            assert Path(pdf_path).exists()
            assert pdf_path.endswith(".pdf")
            assert Path(pdf_path).stat().st_size > 0

    def test_handles_dummy_masks(self):
        """Patchwise segmentation returns empty masks — should use centroid dots."""
        from CITEgeist.model.morphology.segmentation import SegmentationQC
        from CITEgeist.model.morphology.segmentation_qc import (
            generate_segmentation_qc_pdf,
        )

        image = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
        masks = np.zeros((100, 100), dtype=np.int32)
        centroids_xy = np.random.rand(10, 2) * 100
        spot_centers_xy = np.random.rand(5, 2) * 100
        nuclei_count_raw = pd.Series([2, 3, 1, 0, 4], index=[f"s{i}" for i in range(5)])
        qc = SegmentationQC(
            n_raw=12,
            n_after_area_filter=10,
            n_removed_small=1,
            n_removed_large=1,
            pixel_size_um=0.5,
            min_area_um2=10.0,
            max_area_um2=500.0,
            area_median_um2=80.0,
            area_p25_um2=60.0,
            area_p75_um2=100.0,
            density_per_mm2=2000.0,
            density_flag="ok",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            pdf_path = generate_segmentation_qc_pdf(
                output_folder=tmpdir,
                sample_name="patchwise_test",
                image=image,
                masks=masks,
                centroids_xy=centroids_xy,
                spot_centers_xy=spot_centers_xy,
                spot_radius_px=10.0,
                nuclei_count_raw=nuclei_count_raw,
                qc=qc,
            )
            assert Path(pdf_path).exists()
