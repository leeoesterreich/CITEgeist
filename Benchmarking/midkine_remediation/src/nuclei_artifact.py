"""The v2 rebuilt-nuclei input (stage B2): guards shared by every consumer entry point (B1 refit, B3 rebuild, B4 SACE layers).

B2 writes ``<nuclei_dir>/<specimen>/<specimen>_nuclei_counts_fullres_patchwise.csv`` (index = Visium barcode, column
``nuclei_count_raw``, one row per p0 spot) as a REAL, stamped, non-symlink file. A consumer never opens a nuclei file directly;
it calls ``load_nuclei``, which:

  * refuses ``output/morphology_assignment_v3`` (the v4 stand-in) by lexical path, by resolved path, and by resolved target
    of any stand-in entry, and refuses a resolved path that names another v41 artifact (``FORBIDDEN_TOKENS``);
  * refuses a symlinked file or a symlinked parent below the repo root;
  * checks the stamp (stage B2, this protocol) and the v41 denylist (``scripts/v42_stamp.verify_stamp``);
  * requires the barcodes to be unique and EXACTLY equal to the specimen's full p0 barcode set (no missing, no extra).

Zero is a legal count (an edge or blank patch). It is kept as 0: no floor is invented here.

ONE nuclei source (v42 fix round, author decision): B2 also writes the per-nucleus table ``<specimen>_nucleus_spot_map.csv``
(``nucleus_id, spot_id, x_pixel, y_pixel``: the nucleus-to-spot map and the centroids, fullres pixels) from the SAME StarDist run. The counts
file is DERIVED from that table (a spot's count = the number of its rows) and its stamp lists the table among its inputs. Counts and positions
are loaded only together (``load_nuclei_artifact``; ``load_nuclei`` is its counts view): a counts-only B2 (no table, a table the counts stamp
does not bind, or counts that are not the table's row counts) is refused with ``NucleiArtifactOld`` / ``NucleiCountMismatch``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import check_v42_lineage as L  # noqa: E402
import v42_stamp as W  # noqa: E402

STANDIN = "morphology_assignment_v3"
STANDIN_DIR = REPO / "output" / STANDIN
COUNT_COL = "nuclei_count_raw"
TABLE_COLS = ["nucleus_id", "spot_id", "x_pixel", "y_pixel"]


class StandInRefused(ValueError):
    """The path is (or resolves into) output/morphology_assignment_v3 or another v41 artifact."""


class NucleiCoverageError(ValueError):
    """Nuclei barcodes are not exactly the specimen's p0 barcodes, or a barcode is duplicated."""


def nuclei_csv(nuclei_dir, specimen: str) -> Path:
    return Path(nuclei_dir) / specimen / f"{specimen}_nuclei_counts_fullres_patchwise.csv"


def table_csv(nuclei_dir, specimen: str) -> Path:
    """The per-nucleus table (same name and columns as the segmentation's ``*_nucleus_spot_map.csv``)."""
    return Path(nuclei_dir) / specimen / f"{specimen}_nucleus_spot_map.csv"


def qc_json(nuclei_dir, specimen: str) -> Path:
    """Segmentation QC of the one StarDist run (``SegmentationQC.to_dict()``: raw/filtered/removed counts, area stats, density, density flag)."""
    return Path(nuclei_dir) / specimen / f"{specimen}_segmentation_qc.json"


class NucleiArtifactOld(NucleiCoverageError):
    """A counts-only (or unbound) B2 artifact: positions are missing, so the counts are not from the one nuclei source."""


class NucleiCountMismatch(NucleiCoverageError):
    """The counts are not exactly the per-spot row counts of the nucleus table."""


def refuse_standin(path, standin_dir=None) -> None:
    """Raise StandInRefused if ``path`` is, or resolves into, the v3 stand-in (or names another v41 artifact)."""
    p = Path(path)
    resolved = p.resolve()
    for what, parts in (("path", p.parts), ("resolved path", resolved.parts)):
        if STANDIN in parts:
            raise StandInRefused(f"{what} {p} is under {STANDIN} (a v4 stand-in): no v2 stage may read it")
    sd = Path(standin_dir if standin_dir is not None else STANDIN_DIR)
    if sd.is_dir():
        for dp, dns, fns in os.walk(sd):  # entries of the stand-in are symlinks into a v4 rebuild: refuse those targets too
            for name in dns + fns:
                entry = Path(dp) / name
                if entry.is_symlink():
                    tgt = entry.resolve()
                    if resolved == tgt or tgt in resolved.parents:
                        raise StandInRefused(f"{p} resolves into {tgt}, the target of stand-in entry {entry}")
    hit = L.forbidden_hit(str(resolved).replace(str(REPO), ""))
    if hit:
        raise StandInRefused(f"{p} resolves to {resolved}, which names a v41 artifact ({hit})")


def check_coverage(nuclei_index, p0_index, specimen: str = "?") -> None:
    """Exact set equality between the nuclei barcodes and the specimen's full p0 barcodes, both unique."""
    n, p = pd.Index(nuclei_index), pd.Index(p0_index)
    if len(p) == 0 or len(n) == 0:
        raise NucleiCoverageError(f"{specimen}: empty p0 ({len(p)}) or nuclei ({len(n)}) barcode set (vacuous coverage)")
    for name, idx in (("p0", p), ("nuclei", n)):
        if idx.has_duplicates:
            dup = list(idx[idx.duplicated()].unique()[:3])
            raise NucleiCoverageError(f"{specimen}: duplicate {name} barcodes ({int(idx.duplicated().sum())}), e.g. {dup}")
    missing, extra = p.difference(n), n.difference(p)
    if len(missing) or len(extra):
        raise NucleiCoverageError(
            f"{specimen}: nuclei barcodes are not the p0 barcodes: {len(missing)} p0 spots without nuclei "
            f"(e.g. {list(missing[:3])}), {len(extra)} nuclei barcodes not in p0 (e.g. {list(extra[:3])})"
        )


def read_counts(csv: Path) -> pd.Series:
    df = pd.read_csv(csv, index_col=0)
    if COUNT_COL not in df.columns:
        raise ValueError(f"{csv}: no {COUNT_COL} column (columns: {list(df.columns)})")
    s = df[COUNT_COL]
    if s.isna().any() or (s < 0).any() or not np.allclose(s, np.round(s)):
        raise ValueError(f"{csv}: {COUNT_COL} must be finite non-negative integers")
    return s.astype(np.int64)


def read_table(csv: Path) -> pd.DataFrame:
    df = pd.read_csv(csv, float_precision="round_trip")  # exact float round trip: positions are written and re-read bit for bit
    if list(df.columns) != TABLE_COLS:
        raise ValueError(f"{csv}: columns must be {TABLE_COLS}, found {list(df.columns)}")
    if len(df) and not np.isfinite(df[["x_pixel", "y_pixel"]].to_numpy(dtype=float)).all():
        raise ValueError(f"{csv}: non-finite nucleus positions")
    if not np.array_equal(df["nucleus_id"].to_numpy(), np.arange(1, len(df) + 1)):
        raise ValueError(f"{csv}: nucleus_id must be 1..N in order (the cell ids of the Hungarian step)")
    df["spot_id"] = df["spot_id"].astype(str)
    return df


def counts_from_table(table: pd.DataFrame, spot_index) -> pd.Series:
    """Per-spot nucleus counts of the table over ``spot_index`` (spots without a nucleus = 0). A nucleus in a spot outside the index raises."""
    idx = pd.Index(spot_index)
    stray = pd.Index(table["spot_id"].unique()).difference(idx)
    if len(stray):
        raise NucleiCountMismatch(f"nucleus table has {len(stray)} spots that are not in the spot set, e.g. {list(stray[:3])}")
    return table.groupby("spot_id").size().reindex(idx).fillna(0).astype(np.int64).rename(COUNT_COL)


def check_counts_are_table(counts: pd.Series, table: pd.DataFrame, specimen: str = "?") -> None:
    """Exact: same spot set, and every count equals the number of table rows of that spot."""
    derived = counts_from_table(table, counts.index)
    bad = (derived.to_numpy() != counts.to_numpy())
    if bad.any():
        i = np.where(bad)[0][:3]
        raise NucleiCountMismatch(
            f"{specimen}: counts are not the nucleus table's row counts for {int(bad.sum())} spots, e.g. "
            + ", ".join(f"{counts.index[j]}: counts {int(counts.iloc[j])} vs table {int(derived.iloc[j])}" for j in i))


def load_nuclei(*a, **k) -> pd.Series:
    """Counts view of ``load_nuclei_artifact``: the counts load only together with the positions table (see there)."""
    return load_nuclei_artifact(*a, **k)[0]


def load_nuclei_artifact(
    nuclei_dir,
    specimen: str,
    p0_index,
    *,
    protocol_sha12: str | None = None,
    run_id: str | None = None,
    repo_root=None,
    runs: dict | None = None,
    verify_stamp: bool = True,
):
    """The only way a v2 stage reads nuclei: ``(counts, table)`` from ONE stamped artifact. ``run_id`` defaults to the protocol's B2 run id.

    Refuses: the stand-in; a symlink; a missing table (``NucleiArtifactOld``: a counts-only B2 from before the one-source rerun); counts whose
    stamp does not list this very table among its inputs; counts that are not the table's per-spot row counts; barcodes that are not exactly
    the specimen's p0 barcodes. With ``verify_stamp`` it also checks both stamps (stage B2, this protocol) and the v41 denylist."""
    refuse_standin(nuclei_dir)
    f = nuclei_csv(nuclei_dir, specimen)
    t = table_csv(nuclei_dir, specimen)
    refuse_standin(f)
    refuse_standin(t)
    if not f.is_file():
        raise FileNotFoundError(f"B2 nuclei artifact missing: {f}")
    if not t.is_file():
        raise NucleiArtifactOld(
            f"{specimen}: {t.name} is missing in {t.parent}: this is a counts-only B2 artifact (no per-nucleus positions). v42 reads ONE nuclei source "
            f"(counts + positions from the same segmentation run); rerun B2 under the current V2_TS (v42_run_ids.yaml current_v2_ts) instead of using {nuclei_dir}")
    root = W.resolve_roots(repo_root)[0]  # the namespace root: the tree the stamped paths are relative to
    rels = {}
    for what, path in (("counts", f), ("table", t)):
        try:
            rels[what] = str(path.relative_to(root))
        except ValueError:
            rels[what] = None
        if rels[what] is not None:
            bad = L.path_problems(root, rels[what])
            if bad:
                raise StandInRefused("; ".join(bad))
        elif path.is_symlink() or any(q.is_symlink() for q in path.parents if q != Path("/")):
            raise StandInRefused(f"{path} is, or sits below, a symlink: the B2 artifact must be a real file")
    if verify_stamp:
        if protocol_sha12 is None:
            raise ValueError("protocol_sha12 is required to verify the B2 stamp")
        rid = run_id or L.resolve_run(protocol_sha12, "B2", runs if runs is not None else L.load_run_ids())
        stamps = {}
        for what, path in (("counts", f), ("table", t)):
            try:
                stamps[what] = W.verify_stamp(path, run_id=rid, protocol_sha12=protocol_sha12, repo_root=root)
            except W.StampError as e:
                raise NucleiCoverageError(f"{specimen}: nuclei {what} file is not a valid stamped B2 output: {e}") from e
        want = L.sha256_file(t)
        if rels["table"] is None or stamps["counts"]["input_sha256"].get(rels["table"]) != want:
            raise NucleiArtifactOld(
                f"{specimen}: the counts stamp does not list {t.name} (same hash) among its inputs: the counts are not derived from this positions table "
                f"(a counts-only B2 with a table added later is refused)")
    s = read_counts(f)
    check_coverage(s.index, p0_index, specimen)
    table = read_table(t)
    check_counts_are_table(s, table, specimen)
    return s, table
