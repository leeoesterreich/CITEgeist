"""A12.2 input inventory: pooled depth reference M and per-specimen nuclei-floor parameters f_spec.

Observed raw counts only (filtered_feature_bc_matrix.h5 column sums) and the legacy StarDist nuclei counts
(read-only, used for n_cells_est only). No SACE/COMMOT/hurdle/Moran outcome is read.
M = median observed raw spot UMI over all spots of all 12 specimens present in BOTH the QP csv and the GEX object.
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))

import compartments as C  # noqa: E402
import provenance  # noqa: E402
import spatial_utils as S  # noqa: E402

SPECIMENS = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2", "HCC22-088-P2-S1", "HCC22-088-P2-S2", "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep", "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep", "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]


def spot_umi(h5_path: Path) -> pd.Series:
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        bc = np.array([b.decode() for b in g["barcodes"][:]])
        data, indptr = g["data"][:], g["indptr"][:]
    sums = np.add.reduceat(data.astype(np.float64), indptr[:-1]) if len(data) else np.zeros(len(bc))
    sums[np.diff(indptr) == 0] = 0.0          # reduceat returns data[start] for empty segments
    return pd.Series(sums, index=bc)


def main():
    umi, qps, per = {}, {}, {}
    for s in SPECIMENS:
        u = spot_umi(S.DATA_ROOT / s / "outs" / "filtered_feature_bc_matrix.h5")
        qp = pd.read_csv(S.QP_DIR / f"{s}_cell_prop_finetuned_results.csv", index_col=0)
        umi[s], qps[s] = u, qp
    M = C.pooled_depth_reference(umi, qps)
    for s in SPECIMENS:
        common = umi[s].index.intersection(qps[s].index)
        nuc = pd.read_csv(Path(S.REPO_ROOT) / "output" / "morphology_assignment_v3" / s / f"{s}_nuclei_counts_fullres_patchwise.csv", index_col=0)["nuclei_count_raw"]
        n = nuc.reindex(common)
        per[s] = {
            "n_spots_qp": int(len(qps[s])), "n_spots_gex": int(len(umi[s])), "n_common": int(len(common)),
            "n_common_missing_nuclei": int(n.isna().sum()),
            "median_spot_umi_common": float(umi[s].loc[common].median()),
            "median_nuclei_ge1": C.median_nuclei_ge1(n.to_numpy()),
            "f_spec": C.f_spec(n.to_numpy()),
            "frac_common_spots_zero_nuclei": float((n.fillna(0) < 1).mean()),
        }
    payload = {
        "M_pooled_median_spot_umi": M,
        "n_spots_in_pool": int(sum(v["n_common"] for v in per.values())),
        "definition": "median observed raw spot UMI over all spots of the 12 specimens present in both the QP csv and the GEX object",
        "nuclei_source": "output/morphology_assignment_v3/{specimen}/{specimen}_nuclei_counts_fullres_patchwise.csv (legacy StarDist counts, read-only)",
        "per_specimen": per,
    }
    out = HERE.parent / "inventory_M_fspec.json"
    provenance.write_stamped_json(out, payload, "Benchmarking/midkine_remediation/src/inventory_depth_nuclei.py")
    print(out, M, {k: v["f_spec"] for k, v in per.items()})


if __name__ == "__main__":
    main()
