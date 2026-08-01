"""Pure helpers for Fig 4A directional E2 scoring.

Empirical signed-contrast scoring of four MSigDB Li 2023 gene sets:
  LI_ESTROGENE_EARLY_E2_RESPONSE_{UP,DN}
  LI_ESTROGENE_LATE_E2_RESPONSE_{UP,DN}

Per-cancer-cell `(UP_score − DN_score)` for each phase, then spot-mean.
Per-call `gene_pool` keeps the current signature in the pool and excludes
the OTHER three signatures from control draws (scanpy requires
`gene_list ⊆ gene_pool`, so excluding the current signature would crash).

NaN-return (not zero-return) for missing-data cases so the diverging cmap
in the plotter can render those spots gray via na_color=RGBA(0.83,0.83,0.83,0.4).

No matplotlib, no I/O — testable in isolation.
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd

# Hardcoded per Codex review (control-pool leakage prevention requires
# the four sets to be known at import time). Fetched from MSigDB v2024.1.Hs
# on 2026-06-03; refresh requires re-pulling all four lists from gsea-msigdb.org.
# Source: Li et al. 2023, PMID 37272757 (EstroGene database).
SIGNATURE_GENES: dict[str, list[str]] = {
    "early_up": [
        "ADCY9",
        "ADORA1",
        "ADRB1",
        "AMZ1",
        "ASB13",
        "C1orf226",
        "CA12",
        "CASTOR1",
        "CCDC88C",
        "CCND1",
        "CELSR2",
        "CISH",
        "CXCL12",
        "CXXC5",
        "DDX10",
        "DEPTOR",
        "DOK7",
        "EEIG1",
        "EGR3",
        "ELOVL2",
        "FCMR",
        "FHL2",
        "FKBP4",
        "FMN1",
        "FOS",
        "FOXC1",
        "GREB1",
        "HCK",
        "HEY2",
        "HSPB8",
        "IGFBP4",
        "IL17RB",
        "JAK2",
        "KCNK5",
        "KCNK6",
        "KRT13",
        "LONRF2",
        "MANEAL",
        "MBOAT1",
        "MICB",
        "MREG",
        "MYB",
        "MYBL1",
        "MYC",
        "NEIL2",
        "NHERF1",
        "NKAIN1",
        "NOP16",
        "NPY1R",
        "NR5A2",
        "NRIP1",
        "OSGIN1",
        "PDZK1",
        "PFKFB3",
        "PGR",
        "PKIB",
        "PLIN5",
        "PMAIP1",
        "PPIF",
        "PRAG1",
        "PTGES",
        "PXK",
        "RAB31",
        "RAPGEFL1",
        "RARA",
        "RASGRP1",
        "RBM24",
        "RBP7",
        "RET",
        "SEC14L2",
        "SEMA3G",
        "SFXN2",
        "SGK1",
        "SGK3",
        "SIAH2",
        "SLC22A5",
        "SLC47A1",
        "SLC7A5",
        "SLITRK4",
        "STC2",
        "SYBU",
        "THBS1",
        "TIPARP",
        "TMPRSS3",
        "TPD52L1",
        "USP31",
        "WWC1",
        "ZNF703",
    ],
    "early_dn": [
        "AMIGO2",
        "ARHGEF37",
        "ARID5B",
        "ARL4D",
        "BAMBI",
        "BIK",
        "BLNK",
        "BMF",
        "BTG2",
        "C9orf152",
        "CCN2",
        "CCNG2",
        "CDYL2",
        "CFAP206",
        "CMYA5",
        "CREBRF",
        "DDIT4",
        "DRAM1",
        "DRD1",
        "EDN2",
        "EFNA1",
        "EGLN3",
        "ELF5",
        "ENC1",
        "FAM107B",
        "FAM110C",
        "FAM171B",
        "GRAMD2B",
        "GRB7",
        "GRHL3",
        "HCAR1",
        "HCAR2",
        "HILPDA",
        "ID1",
        "ID2",
        "ID3",
        "IL1R1",
        "KCNJ8",
        "LIMA1",
        "LIPH",
        "LNX1",
        "LYPD3",
        "NUAK1",
        "PDP1",
        "PIK3R3",
        "PLEKHF2",
        "PLK2",
        "PPFIBP2",
        "PTGER4",
        "RAB27B",
        "REL",
        "RGS2",
        "RIN2",
        "RND3",
        "RNF144B",
        "RNF43",
        "RPRM",
        "RUNDC3B",
        "S1PR3",
        "SGCG",
        "SH3BP4",
        "SMAD6",
        "SOCS2",
        "SOX2",
        "SP6",
        "SPRY1",
        "ST8SIA4",
        "STON1",
        "TENT5B",
        "TGFB2",
        "TGFB3",
        "TM4SF1",
        "TNFRSF11B",
        "TP53INP1",
        "TUFT1",
        "VEGFC",
        "ZNF467",
    ],
    "late_up": [
        "A4GALT",
        "ACOX2",
        "ADRA2A",
        "ARTN",
        "ATRNL1",
        "BFSP2",
        "C1QTNF6",
        "C5AR2",
        "CA8",
        "CCN5",
        "CD44",
        "CHAC1",
        "CHST8",
        "CTPS1",
        "CTSD",
        "CYP26B1",
        "DSCC1",
        "EXO1",
        "FHL1",
        "GAL",
        "GATA4",
        "GJB3",
        "GPR68",
        "GRIK3",
        "H19",
        "KCNH1",
        "KLK11",
        "LRP8",
        "MAT1A",
        "MCM10",
        "MCM7",
        "MGAT3",
        "MGP",
        "MPPED2",
        "NAV2",
        "NCAPG",
        "NEIL3",
        "NOD2",
        "NPY5R",
        "NXPH3",
        "PDLIM3",
        "PLAT",
        "PLCL1",
        "PPM1K",
        "PRSS23",
        "RRM2",
        "SERPINA1",
        "SERPINA5",
        "SOX3",
        "SPINK4",
        "TFCP2L1",
        "TH",
        "THSD4",
        "TMEM229B",
        "TRPC6",
        "UGT2B15",
    ],
    "late_dn": [
        "ABCA1",
        "ABCC5",
        "ACHE",
        "ACKR3",
        "ANTXR2",
        "AQP3",
        "C1orf115",
        "CA5B",
        "CAMK2N1",
        "CEMIP",
        "CLIP4",
        "CLMN",
        "COL28A1",
        "CRYM",
        "CSTA",
        "CTNND2",
        "CXCR4",
        "CYP4B1",
        "DAB2",
        "DKK1",
        "DLX2",
        "EFNB2",
        "EMP1",
        "ENTPD3",
        "EPAS1",
        "FAM110B",
        "FAM83A",
        "FBN2",
        "FBXO32",
        "FPR3",
        "GABBR2",
        "GRB14",
        "GRM4",
        "H2AC6",
        "H2BC8",
        "H4C8",
        "HDAC9",
        "HOPX",
        "HPX",
        "IGDCC3",
        "IGFBP3",
        "ITGB6",
        "KLHDC7A",
        "KRT24",
        "KRT4",
        "KRT81",
        "LDLRAD4",
        "LXN",
        "MAB21L4",
        "MMP16",
        "MYO1B",
        "NBEA",
        "NTN4",
        "P2RX2",
        "PALMD",
        "PCDH7",
        "PGM5",
        "PHEX",
        "PIK3IP1",
        "PLEKHF1",
        "PMP22",
        "PRKG1",
        "PVALB",
        "RFTN1",
        "RND1",
        "RTL5",
        "SEMA6D",
        "SERHL2",
        "SLIT2",
        "SSPOP",
        "ST3GAL1",
        "ST6GAL1",
        "TCP11L2",
        "TMEM86A",
        "TSPAN1",
        "UNC5C",
        "UPK1A",
        "UPK2",
        "VGLL1",
        "WNT11",
    ],
}


def signature_union() -> set[str]:
    """All genes appearing in any of the four directional signatures."""
    out: set[str] = set()
    for v in SIGNATURE_GENES.values():
        out |= set(v)
    return out


def build_gene_pool_for_call(
    var_names: Iterable[str],
    current_signature_key: str,
) -> list[str]:
    """Per-call control gene pool.

    Returns var_names minus the union of the OTHER three signatures.
    Scanpy requires gene_list \u2286 gene_pool, so the *current* signature
    must remain in the pool; only the other three are excluded.
    """
    if current_signature_key not in SIGNATURE_GENES:
        raise KeyError(
            f"unknown signature key: {current_signature_key!r}; " f"expected one of {sorted(SIGNATURE_GENES.keys())}"
        )
    excluded: set[str] = set()
    for key, genes in SIGNATURE_GENES.items():
        if key != current_signature_key:
            excluded |= set(genes)
    return [g for g in var_names if g not in excluded]


def compute_directional_score(
    cancer,  # AnnData of cancer-only cells
    phase: str,  # "early" or "late"
    min_overlap: int = 5,
    random_state: int = 0,
) -> Optional[pd.Series]:
    """Per-spot directional score = mean over cancer cells of (UP - DN).

    Returns None if EITHER constituent set is unscorable. Caller treats
    None as "this phase is invalid for this sample" and NaN-fills.
    """
    if phase not in ("early", "late"):
        raise ValueError(f"phase must be 'early' or 'late', got {phase!r}")

    # Guard 1: no cancer cells.
    if cancer.n_obs == 0:
        return None

    import scanpy as sc

    up_key, dn_key = f"{phase}_up", f"{phase}_dn"
    up_set = SIGNATURE_GENES[up_key]
    dn_set = SIGNATURE_GENES[dn_key]

    up_present = [g for g in up_set if g in cancer.var_names]
    dn_present = [g for g in dn_set if g in cancer.var_names]

    # Guard 2: low overlap in either constituent — invalidate whole phase.
    if len(up_present) < min_overlap or len(dn_present) < min_overlap:
        return None

    # Symmetric ctrl_size per phase.
    ctrl_size = min(len(up_present), len(dn_present))

    cancer = cancer.copy()  # don't mutate caller's AnnData
    var_names_list = list(cancer.var_names)
    up_pool = build_gene_pool_for_call(var_names_list, up_key)
    dn_pool = build_gene_pool_for_call(var_names_list, dn_key)

    try:
        sc.tl.score_genes(
            cancer,
            gene_list=up_present,
            score_name=f"{up_key}_score",
            ctrl_size=ctrl_size,
            gene_pool=up_pool,
            use_raw=False,
            random_state=random_state,
        )
        sc.tl.score_genes(
            cancer,
            gene_list=dn_present,
            score_name=f"{dn_key}_score",
            ctrl_size=ctrl_size,
            gene_pool=dn_pool,
            use_raw=False,
            random_state=random_state,
        )
    except (ValueError, RuntimeError):
        # Guard 3: scanpy bailed (insufficient pool, NaN expression, etc).
        return None

    per_cell = cancer.obs[f"{up_key}_score"].values - cancer.obs[f"{dn_key}_score"].values

    df = pd.DataFrame(
        {
            "spot_barcode": cancer.obs["spot_barcode"].values,
            "score": per_cell,
        }
    )
    spot_mean = df.groupby("spot_barcode")["score"].mean()
    spot_mean.name = f"{phase}_e2"
    return spot_mean


def compute_directional_score_layer(
    cancer_layer: "ad.AnnData",
    phase: str,
    min_overlap: int = 5,
    random_state: int = 0,
) -> Optional[pd.Series]:
    """Spot-level directional score on SACE cancer-layer AnnData.

    Each obs is a Visium spot (not a cell). Expression is the SACE-allocated
    cancer compartment. Returns per-spot (UP - DN) scores directly — no
    cell-to-spot aggregation needed.
    """
    if phase not in ("early", "late"):
        raise ValueError(f"phase must be 'early' or 'late', got {phase!r}")
    if cancer_layer.n_obs == 0:
        return None

    import scanpy as sc

    up_key, dn_key = f"{phase}_up", f"{phase}_dn"
    up_set = SIGNATURE_GENES[up_key]
    dn_set = SIGNATURE_GENES[dn_key]

    up_present = [g for g in up_set if g in cancer_layer.var_names]
    dn_present = [g for g in dn_set if g in cancer_layer.var_names]

    if len(up_present) < min_overlap or len(dn_present) < min_overlap:
        return None

    ctrl_size = min(len(up_present), len(dn_present))
    adata = cancer_layer.copy()
    var_names_list = list(adata.var_names)
    up_pool = build_gene_pool_for_call(var_names_list, up_key)
    dn_pool = build_gene_pool_for_call(var_names_list, dn_key)

    try:
        sc.tl.score_genes(
            adata,
            gene_list=up_present,
            score_name=f"{up_key}_score",
            ctrl_size=ctrl_size,
            gene_pool=up_pool,
            use_raw=False,
            random_state=random_state,
        )
        sc.tl.score_genes(
            adata,
            gene_list=dn_present,
            score_name=f"{dn_key}_score",
            ctrl_size=ctrl_size,
            gene_pool=dn_pool,
            use_raw=False,
            random_state=random_state,
        )
    except (ValueError, RuntimeError):
        return None

    scores = adata.obs[f"{up_key}_score"].values - adata.obs[f"{dn_key}_score"].values
    return pd.Series(scores, index=adata.obs_names, name=f"{phase}_e2")


def decide_percentile(combined_abs_finite: np.ndarray) -> tuple[float, str]:
    """Pick e2_vmax via objective 99/95 ratio trigger.

    Spec: ratio strictly > 1.5 drops to p95; <=1.5 keeps p99.
    """
    if combined_abs_finite.size == 0:
        return 1.0, "p99"
    p99 = float(np.nanpercentile(combined_abs_finite, 99))
    p95 = float(np.nanpercentile(combined_abs_finite, 95))
    if p95 > 0 and (p99 / p95) > 1.5:  # strict
        return p95, "p95"
    return p99, "p99"


def decision_gate_preconditions(
    per_sample_n_evaluable: dict[str, int],
    early_std: float,
    late_std: float,
    min_samples: int = 2,
    min_per_sample_n: int = 10,
    max_single_sample_share: float = 0.7,
) -> dict[str, bool]:
    """Check whether pooled Pearson r is informative enough to drive the
    collapse decision. Returns a dict of named booleans; caller treats
    all-True as "use r"; any False as "keep 3 rows + flag for human review".
    """
    counts = [n for n in per_sample_n_evaluable.values() if n >= min_per_sample_n]
    total = sum(per_sample_n_evaluable.values())
    max_share = max(per_sample_n_evaluable.values()) / total if total > 0 else 1.0
    return {
        "min_samples_evaluable": len(counts) >= min_samples,
        "balanced_contribution": max_share <= max_single_sample_share,
        "nonzero_variance": early_std > 0.0 and late_std > 0.0,
    }
