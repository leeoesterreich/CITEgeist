"""rnaprot Xenium GT: signed-off RNA-cluster labels, minus cells whose protein confidently says a different lineage.

Spec: docs/superpowers/specs/2026-10-24-xenium-gt-rna-protein-concordant.md (Part A; Z = 3 fixed before scoring).
Per channel: z = (log1p x - negative mean) / negative SD, the negative component being the low mode of a 2-component
GMM on the nonzero cells. A lineage is positive when z > Z (strict). The call is the positive lineage with the highest
z, else "none". A cell is dropped only when its call is a different lineage ("none" keeps it). Plasma guard: an RNA
B cell with CD138 z > Z is never dropped. 7-type: also drop an RNA CD4/CD8 cell when protein is CD8A+/CD4- (resp.
CD4+/CD8A-) and the call is T. E-Cadherin is not used.

Usage: python -m gt_rebuild.rnaprot build [--out-root DIR] [--z 3.0]
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from . import gt_config, manifest

Z_DEFAULT = 3.0
SPEC_REL = "docs/superpowers/specs/2026-10-24-xenium-gt-rna-protein-concordant.md"
REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))
import provenance  # noqa: E402

DEFAULT_OUT = gt_config.BASE / "ground_truth_xenium_rnaprot"
RNA_ROOT = gt_config.BASE / "ground_truth_xenium_clusters"

# lineage -> protein channels (max over channels); order is the tie-break order of the reference analysis
LINEAGE_CHANNELS = {
    "B cells": ["CD20", "CD138"],
    "T cells": ["CD3E"],
    "Macrophages": ["CD68", "CD163"],
    "Endothelial": ["CD31"],
    "Epithelial": ["PanCK"],
    "Fibroblasts": ["alphaSMA"],
}
SUBTYPE_CHANNELS = ["CD4", "CD8A"]
NEEDED_CHANNELS = sorted({c for v in LINEAGE_CHANNELS.values() for c in v} | set(SUBTYPE_CHANNELS))
_LIN_OF = {"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"}
GMM_SUBSAMPLE = 200000


def negative_fit(x_log1p, seed=0, rng=None, include_zeros=False):
    """(mean, sd) of the negative component of one channel (log1p values).

    A 2-component GMM (random_state=42) is fitted on a <=200k subsample of the nonzero cells; the negative
    component is the lower-mean one, assigned to every nonzero cell. As in the reference analysis (concord.py /
    background.py) the statistics are taken over the low NONZERO component; ``include_zeros=True`` also pools the
    exact zeros (diagnostic only, not used by the build). ``rng`` lets protein_z share one RandomState across
    channels (the reference consumed it in channel order).
    """
    from sklearn.mixture import GaussianMixture

    v = np.asarray(x_log1p, dtype=np.float64)
    rng = np.random.RandomState(seed) if rng is None else rng
    nz = v > 0
    vz = v[nz]
    sub = rng.choice(len(vz), min(GMM_SUBSAMPLE, len(vz)), replace=False)
    g = GaussianMixture(2, random_state=42).fit(vz[sub].reshape(-1, 1))
    low = int(np.argmin(g.means_.ravel()))
    neg = np.zeros(len(v), bool)
    neg[nz] = g.predict(vz.reshape(-1, 1)) == low
    if include_zeros:
        neg |= ~nz
    return float(v[neg].mean()), float(v[neg].std())


def protein_z(protein: pd.DataFrame, keep: Optional[Sequence[str]] = None, seed=0) -> pd.DataFrame:
    """dsb-style z per channel from RAW protein counts. All channels are fitted in column order from one
    RandomState(seed) (so the subsample draws match the reference); only ``keep`` columns are returned.
    The per-channel (mean, sd) fits are stored in ``result.attrs["neg_fit"]``."""
    rng = np.random.RandomState(seed)
    fits, out = {}, {}
    for c in protein.columns:
        v = np.log1p(protein[c].to_numpy(dtype=np.float64))
        m, sd = negative_fit(v, rng=rng)
        fits[c] = {"neg_mean": m, "neg_sd": sd}
        if keep is None or c in keep:
            out[c] = (v - m) / sd
    z = pd.DataFrame(out, index=protein.index)
    if keep is not None:
        z = z[list(keep)]
    z.attrs["neg_fit"] = fits
    return z


def lineage_call(z: pd.DataFrame, Z: float = Z_DEFAULT) -> pd.Series:
    """Positive lineage (z > Z, strict) with the highest z, else "none"."""
    names = list(LINEAGE_CHANNELS)
    M = np.column_stack([z[chs].max(axis=1).to_numpy() for chs in LINEAGE_CHANNELS.values()])
    pos = M > Z
    best = np.where(pos, M, -np.inf).argmax(1)
    call = np.where(pos.any(1), np.array(names, dtype=object)[best], "none").astype(object)
    return pd.Series(call, index=z.index)


def apply_no_contradiction(
    labels6: pd.Series, labels7: pd.Series, z: pd.DataFrame, Z: float = Z_DEFAULT
) -> pd.DataFrame:
    """keep6 / keep7 / reason per RNA-labelled cell (index = labels6.index; labels7 must cover the same cells)."""
    if set(labels6.index) != set(labels7.index) or labels6.index.has_duplicates:
        raise ValueError("labels6 and labels7 must cover the same unique cells")
    ids = labels6.index
    labels7 = labels7.reindex(ids)
    zz = z.reindex(ids)
    if zz[NEEDED_CHANNELS].isna().any().any():
        raise ValueError("protein z missing for some labelled cells")
    call = lineage_call(zz, Z)
    guard = ((labels6 == "B cells") & (zz["CD138"] > Z)).to_numpy()
    lin7 = labels7.map(lambda v: _LIN_OF.get(v, v))
    if not (lin7.to_numpy() == labels6.to_numpy()).all():
        raise ValueError("labels7 does not collapse to labels6")
    contra = ((call != "none") & (call != labels6)).to_numpy() & ~guard
    cd4, cd8 = (zz["CD4"] > Z).to_numpy(), (zz["CD8A"] > Z).to_numpy()
    is_t_call = (call == "T cells").to_numpy()
    sub = is_t_call & (
        ((labels7 == "CD4+ T cells").to_numpy() & cd8 & ~cd4) | ((labels7 == "CD8+ T cells").to_numpy() & cd4 & ~cd8)
    )
    keep6 = ~contra
    keep7 = keep6 & ~sub
    reason = np.where(contra, "lineage:protein_calls_" + call.astype(str), np.where(sub, "subtype:cd4_cd8_swap", ""))
    return pd.DataFrame({"keep6": keep6, "keep7": keep7, "reason": reason}, index=ids)


# ----------------------------------------------------------------------------- build
def _git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def _spot_universe(root: Path) -> pd.MultiIndex:
    parts = []
    for r in range(5):
        p = root / "ground_truth_6type" / "ground_truth" / f"Xenium_region_{r}_prop.csv"
        parts.append(pd.DataFrame({"r": r, "s": pd.read_csv(p, index_col=0, usecols=[0]).index}))
    df = pd.concat(parts)
    return pd.MultiIndex.from_arrays([df.r.values, df.s.values], names=["r", "s"])


def keep_stats(res: pd.DataFrame, lab6: pd.Series, lab7: pd.Series, mapping: pd.DataFrame, universe) -> dict:
    """Kept fractions (all RNA-labelled cells) and spot coverage over the RNA-GT spot universe."""
    out = {"n_labelled": int(len(res))}
    for n, lab, col in ((6, lab6, "keep6"), (7, lab7, "keep7")):
        k = res[col]
        out[f"kept_{n}type"] = float(k.mean())
        out[f"kept_per_type_{n}type"] = {t: float(k[lab == t].mean()) for t in sorted(lab.unique())}
        out[f"n_per_type_{n}type"] = {t: int((lab == t).sum()) for t in sorted(lab.unique())}
        m = mapping.reindex(res.index)
        ok = ((m["spot_idx"] >= 0) & m["region_id"].notna()).to_numpy() & k.to_numpy()
        cnt = pd.DataFrame({"r": m["region_id"][ok].astype(int).to_numpy(), "s": m["spot_id"][ok].to_numpy()})
        nn = cnt.groupby(["r", "s"]).size().reindex(universe, fill_value=0)
        out[f"spots_{n}type"] = {
            "median_kept_per_spot": float(nn.median()),
            "p10": float(nn.quantile(0.1)),
            "spots_0_kept": int((nn == 0).sum()),
            "spots_lt3_kept": int((nn < 3).sum()),
            "n_spots_universe": int(len(universe)),
        }
    return out


def build(out_root: Path, Z: float = Z_DEFAULT, seed: int = 0) -> dict:
    out_root = Path(out_root)
    if out_root.exists() and any(out_root.iterdir()):
        raise FileExistsError(f"refusing to write into non-empty {out_root}")
    rna_manifest_path = RNA_ROOT / "manifest.json"
    manifest.require_rnagt_signoff(rna_manifest_path, labelset="rnagt")  # RNA GT must be signed off
    labels_dir = out_root / "labels"
    labels_dir.mkdir(parents=True, exist_ok=True)

    lab = {
        n: pd.read_csv(RNA_ROOT / "labels" / f"cell_labels_{n}type.csv", index_col=0, dtype={"cell_id": str})
        for n in (6, 7)
    }
    lab6, cluster = lab[6]["singler_label"], lab[7]["cluster"]
    lab7 = lab[7]["singler_label"].reindex(lab6.index)
    prot_csv = gt_config.OLD_LABELS_DIR / "xenium_cell_proteins.csv"
    ph = list(pd.read_csv(prot_csv, nrows=0).columns)
    chs = ph[ph.index("PD-1") :]
    prot = pd.read_csv(prot_csv, index_col=0, usecols=[ph[0]] + chs).astype("float32")
    print(f"protein {prot.shape}; RNA-labelled cells {len(lab6)}", flush=True)
    z = protein_z(prot, keep=NEEDED_CHANNELS, seed=seed)
    fits = z.attrs["neg_fit"]
    missing = lab6.index.difference(z.index)
    if len(missing):
        raise ValueError(f"{len(missing)} RNA-labelled cells have no protein row, e.g. {list(missing[:3])}")
    zl = z.reindex(lab6.index)
    res = apply_no_contradiction(lab6, lab7, zl, Z)
    call = lineage_call(zl, Z)

    mapping = pd.read_csv(gt_config.MAPPING_CSV, index_col=0)
    mapping.index = mapping.index.astype(str)
    stats = keep_stats(res, lab6, lab7, mapping, _spot_universe(RNA_ROOT))

    # kept labels (same format as the RNA GT)
    for n, lb, col in ((6, lab6, "keep6"), (7, lab7, "keep7")):
        k = res[col].to_numpy()
        pd.DataFrame(
            {"cell_id": lb.index[k], "singler_label": lb.to_numpy()[k], "cluster": cluster.reindex(lb.index).to_numpy()[k]}
        ).to_csv(labels_dir / f"cell_labels_{n}type.csv", index=False)
    drop = res[~(res.keep6 & res.keep7)]
    pd.DataFrame(
        {
            "cell_id": drop.index,
            "label6": lab6[drop.index].to_numpy(),
            "label7": lab7[drop.index].to_numpy(),
            "protein_call": call[drop.index].to_numpy(),
            "keep6": drop.keep6.to_numpy(),
            "keep7": drop.keep7.to_numpy(),
            "reason": drop.reason.to_numpy(),
        }
    ).to_csv(labels_dir / "dropped_cells.csv", index=False)
    rna_exc = pd.read_csv(RNA_ROOT / "labels" / "excluded_cells.csv", dtype={"cell_id": str})
    pd.concat(
        [rna_exc, pd.DataFrame({"cell_id": res.index[~res.keep7.to_numpy()], "cluster": "protein"})]
    ).to_csv(labels_dir / "excluded_cells.csv", index=False)

    # spot proportions + GEX layers: the same scripts the RNA GT uses (GEX is zero-filled over the props' spots)
    xb = gt_config.BASE
    for n in (6, 7):
        subprocess.run(
            [
                sys.executable, "-u", str(xb / "ground_truth_singler/src/aggregate_to_spots.py"),
                "--labels-csv", str(labels_dir / f"cell_labels_{n}type.csv"),
                "--mapping-csv", str(gt_config.MAPPING_CSV),
                "--output-dir", str(out_root / f"ground_truth_{n}type"),
                "--n-types", str(n), "--prefix", "Xenium",
            ],
            check=True,
        )
    subprocess.run(
        [
            sys.executable, "-u", str(xb / "evaluation/src/generate_singler_gex_gt.py"),
            "--labels-csv", str(labels_dir / "cell_labels_7type.csv"),
            "--out-root", str(out_root / "ground_truth_7type/ground_truth_gex"),
            "--props-dir", str(out_root / "ground_truth_7type/ground_truth"),
        ],
        check=True,
    )
    from .assemble_labels import excluded_cells_per_spot

    for n, col in ((6, "keep6"), (7, "keep7")):
        kept_ids = res.index[res[col].to_numpy()]
        exc_ids = set(rna_exc["cell_id"]) | set(res.index[~res[col].to_numpy()])
        excluded_cells_per_spot(mapping, kept_ids, exc_ids).to_csv(
            out_root / f"ground_truth_{n}type" / "excluded_cells_per_spot.csv", index=False
        )

    # the 9 integrity gates, on the 7-type label file (the RNA GT's own gate setup)
    from . import integrity_gates as ig
    from .cluster_names import read_cluster_names

    names = read_cluster_names(RNA_ROOT / "work" / "cluster_names" / "cluster_names.csv")
    dropped_ids = json.loads((RNA_ROOT / "work" / "level1_qc.json").read_text())["dropped_ids"]
    gates = ig.run_gates(
        out_root, labels_dir / "cell_labels_7type.csv", gt_config.MAPPING_CSV, names, dropped_ids,
        labels_dir / "excluded_cells.csv",
    )
    ig.write_report(gates, out_root / "reports" / "integrity_report.json")
    for g in gates:
        print(f"gate {g.name}: {'ok' if g.ok else 'FAIL'} {g.detail}", flush=True)

    uni, new = _spot_universe(RNA_ROOT), _spot_universe(out_root)
    stats["n_spots_excluded_zero_kept_6type"] = int(len(uni.difference(new)))
    for n in (6, 7):
        stats[f"regions_present_{n}type"] = [
            r for r in range(5) if (out_root / f"ground_truth_{n}type/ground_truth/Xenium_region_{r}_prop.csv").exists()
        ]
    stats["protein_call_counts"] = {k: int(v) for k, v in call.value_counts().items()}
    stats["drop_reasons"] = {k: int(v) for k, v in drop.reason.value_counts().items()}

    payload = {
        "created_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "kind": "rnaprot (RNA labels not contradicted by protein)",
        "Z": Z,
        "seed": seed,
        "spec": SPEC_REL if (REPO / SPEC_REL).is_file() else None,
        "spec_sha256": manifest.sha256_file(REPO / SPEC_REL) if (REPO / SPEC_REL).is_file() else None,
        "rna_gt_manifest": str(rna_manifest_path),
        "rna_gt_manifest_sha256": manifest.sha256_file(rna_manifest_path),
        "protein_csv_sha256": manifest.sha256_file(prot_csv),
        "mapping_csv_sha256": manifest.sha256_file(gt_config.MAPPING_CSV),
        "lineage_channels": LINEAGE_CHANNELS,
        "subtype_channels": SUBTYPE_CHANNELS,
        "negative_fit_log1p": {c: fits[c] for c in NEEDED_CHANNELS},
        "negative_fit_all_channels": fits,
        "keep_stats": stats,
        "git_head": provenance.code_commit(REPO),
        # None = no git checkout, so uncommitted changes could not be checked (never a faked False)
        "rnaprot_py_uncommitted_changes": (
            bool(_git("status", "--porcelain", "--", "Benchmarking/xenium_benchmarking/gt_rebuild/rnaprot.py"))
            if provenance.in_git_checkout(REPO)
            else None
        ),
        "gates": {g.name: {"ok": g.ok, "detail": g.detail} for g in gates},
        "versions": manifest.collect_versions(),
    }
    manifest.write_manifest(out_root / "manifest.json", payload)
    (out_root / "reports" / "rnaprot_stats.json").write_text(json.dumps(stats, indent=2))
    return payload


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    b.add_argument("--z", type=float, default=Z_DEFAULT)
    a = ap.parse_args(argv)
    p = build(a.out_root, a.z)
    manifest.write_complete_marker(
        a.out_root, [os.environ.get("SLURM_JOB_ID", "login")], manifest.sha256_file(a.out_root / "manifest.json")
    )
    print(json.dumps(p["keep_stats"], indent=2))
    if not all(g["ok"] for g in p["gates"].values()):
        sys.exit(2)


if __name__ == "__main__":
    main()
