"""Tests for bulk RNA chaperone interaction analysis."""

import json
import numpy as np
import pandas as pd
import pytest
from pathlib import Path


def _make_synthetic_tpm(tmp_path, gene_id_map):
    """Create a synthetic TPM matrix with known interaction pattern."""
    rng = np.random.RandomState(42)
    n_genes = len(gene_id_map) + 5  # 12 targets + 5 noise

    gene_ids = list(gene_id_map.values()) + [str(99990 + i) for i in range(5)]

    groups = {
        "MCF7_WT": ["GSM_MW1", "GSM_MW2", "GSM_MW3", "GSM_MW4"],
        "MCF7_D538G": ["GSM_MD1", "GSM_MD2", "GSM_MD3", "GSM_MD4"],
        "T47D_WT": ["GSM_TW1", "GSM_TW2", "GSM_TW3", "GSM_TW4"],
        "T47D_D538G": ["GSM_TD1", "GSM_TD2", "GSM_TD3", "GSM_TD4"],
    }
    all_samples = [s for grp in groups.values() for s in grp]

    data = {}
    for sample in all_samples:
        data[sample] = rng.uniform(10, 100, n_genes)

    # Inject known interaction for MDK (id=4192): MCF7 D538G DOWN, T47D D538G UP
    mdk_idx = gene_ids.index("4192")
    for s in groups["MCF7_WT"]:
        data[s][mdk_idx] = 150 + rng.normal(0, 5)
    for s in groups["MCF7_D538G"]:
        data[s][mdk_idx] = 110 + rng.normal(0, 5)
    for s in groups["T47D_WT"]:
        data[s][mdk_idx] = 25 + rng.normal(0, 3)
    for s in groups["T47D_D538G"]:
        data[s][mdk_idx] = 35 + rng.normal(0, 3)

    # Inject opposite interaction for HSP90B1 (id=7184): MCF7 D538G UP, T47D DOWN
    hsp_idx = gene_ids.index("7184")
    for s in groups["MCF7_WT"]:
        data[s][hsp_idx] = 300 + rng.normal(0, 10)
    for s in groups["MCF7_D538G"]:
        data[s][hsp_idx] = 500 + rng.normal(0, 10)
    for s in groups["T47D_WT"]:
        data[s][hsp_idx] = 900 + rng.normal(0, 20)
    for s in groups["T47D_D538G"]:
        data[s][hsp_idx] = 650 + rng.normal(0, 20)

    df = pd.DataFrame(data, index=gene_ids)
    df.index.name = "GeneID"
    tpm_path = tmp_path / "test_tpm.tsv.gz"
    df.to_csv(tpm_path, sep="\t", compression="gzip")

    mapping_path = tmp_path / "gene_id_mapping.json"
    with open(mapping_path, "w") as f:
        json.dump(gene_id_map, f)

    return tpm_path, mapping_path, groups


GENE_ID_MAP = {
    "MDK": "4192", "HSP90B1": "7184", "HSPA5": "3309", "CALR": "811",
    "CANX": "821", "PDIA4": "9601", "PDIA6": "10130", "XBP1": "7494",
    "ATF6": "22926", "SEC61B": "10952", "MAN1A1": "4121", "SEC23A": "10484",
}


class TestBulkRnaInteraction:
    def test_loads_and_maps_genes(self, tmp_path):
        from midkine.scripts.bulk_rna_interaction import load_target_genes
        tpm_path, mapping_path, groups = _make_synthetic_tpm(tmp_path, GENE_ID_MAP)
        sub = load_target_genes(str(tpm_path), str(mapping_path))
        assert "MDK" in sub.index
        assert "HSP90B1" in sub.index
        assert len(sub) == 12

    def test_interaction_table_columns(self, tmp_path):
        from midkine.scripts.bulk_rna_interaction import load_target_genes, compute_interaction_table
        tpm_path, mapping_path, groups = _make_synthetic_tpm(tmp_path, GENE_ID_MAP)
        sub = load_target_genes(str(tpm_path), str(mapping_path))
        table = compute_interaction_table(sub, groups)
        assert "gene" in table.columns
        assert "mcf7_wt_mean" in table.columns
        assert "interaction_F" in table.columns
        assert "interaction_p" in table.columns
        assert "fdr_q" in table.columns
        assert len(table) == 12

    def test_mdk_interaction_significant(self, tmp_path):
        from midkine.scripts.bulk_rna_interaction import load_target_genes, compute_interaction_table
        tpm_path, mapping_path, groups = _make_synthetic_tpm(tmp_path, GENE_ID_MAP)
        sub = load_target_genes(str(tpm_path), str(mapping_path))
        table = compute_interaction_table(sub, groups)
        mdk_row = table[table["gene"] == "MDK"].iloc[0]
        assert mdk_row["interaction_p"] < 0.05

    def test_composite_interaction(self, tmp_path):
        from midkine.scripts.bulk_rna_interaction import load_target_genes, compute_composite_interaction
        tpm_path, mapping_path, groups = _make_synthetic_tpm(tmp_path, GENE_ID_MAP)
        sub = load_target_genes(str(tpm_path), str(mapping_path))
        result = compute_composite_interaction(sub, groups)
        assert "interaction_F" in result
        assert "interaction_p" in result
        assert "mcf7_wt_mean" in result
