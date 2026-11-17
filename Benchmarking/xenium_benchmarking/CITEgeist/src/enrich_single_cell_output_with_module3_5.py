#!/usr/bin/env python
"""Benchmark-gated enrichment helpers for Module 3.5 single-cell outputs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import anndata as ad

from CITEgeist.model.module3_5_projection import should_enrich_module3_5_output



def load_summary(summary_path: Path) -> dict[str, Any]:
    return json.loads(summary_path.read_text())


def write_enriched_output(
    source_h5ad: Path,
    enriched_h5ad: Path,
    *,
    functional_annotations: Any = None,
    functional_metadata: Mapping[str, Any] | None = None,
) -> Path:
    adata = ad.read_h5ad(source_h5ad)
    if functional_annotations is not None and not getattr(functional_annotations, "empty", False):
        aligned = functional_annotations.reindex(adata.obs_names)
        for column in aligned.columns:
            adata.obs[column] = aligned[column]
    if functional_metadata is not None:
        adata.uns["functional_annotation_meta"] = dict(functional_metadata)
    adata.write_h5ad(enriched_h5ad)
    return enriched_h5ad


def should_write_enriched_output(summary: Mapping[str, Any]) -> bool:
    return should_enrich_module3_5_output(summary)


if __name__ == "__main__":
    raise SystemExit("This module is intended to be imported, not executed directly.")
