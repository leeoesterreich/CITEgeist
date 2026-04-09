"""Tests for Module 3.5 projection helpers."""

from CITEgeist.model.annotation.module3_5_projection import should_enrich_module3_5_output


def test_should_enrich_module3_5_output_requires_module3_5_benchmark_pass():
    summary = {"benchmark_passed": True, "validated_pairs": ["CD8_T_Cells__PDCD1"]}
    assert should_enrich_module3_5_output(summary) is True

    failing = {"benchmark_passed": False, "validated_pairs": ["CD8_T_Cells__PDCD1"]}
    assert should_enrich_module3_5_output(failing) is False
