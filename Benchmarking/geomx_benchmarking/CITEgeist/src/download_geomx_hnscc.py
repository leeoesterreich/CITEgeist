"""
download_geomx_hnscc.py
-----------------------
Download and format GeoMx HNSCC data from GEO into paired AnnData objects.

GEO accessions:
  - GSE288406 — GeoMx IPA protein (138 AOIs × 580 proteins)
  - GSE290057 — GeoMx WTA RNA (122 AOIs × 18,815 transcripts)

Usage:
  python download_geomx_hnscc.py --output-dir /path/to/output
  python download_geomx_hnscc.py --output-dir /path/to/output --skip-download
"""

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

try:
    import GEOparse

    HAS_GEOPARSE = True
except ImportError:
    HAS_GEOPARSE = False
    logging.warning("GEOparse not installed — series_matrix fallback disabled. " "Install with: pip install GEOparse")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# GEO accessions
IPA_ACCESSION = "GSE288406"
WTA_ACCESSION = "GSE290057"


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def download_geo_series(accession: str, output_dir: Path) -> Path:
    """Download supplementary files from GEO FTP for a given accession.

    URL pattern:
        https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/suppl/

    Parameters
    ----------
    accession:
        GEO series accession (e.g. "GSE288406").
    output_dir:
        Root directory under which a subdirectory named ``accession`` will
        be created and files mirrored into.

    Returns
    -------
    Path
        Path to the local directory containing the downloaded files.
    """
    prefix = accession[:6]  # e.g. "GSE288"
    ftp_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/suppl/"
    series_dir = output_dir / accession
    series_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading %s from %s", accession, ftp_url)
    logger.info("Local target: %s", series_dir)

    cmd = [
        "wget",
        "--mirror",
        "--no-host-directories",
        "--cut-dirs=7",
        "--no-parent",
        "--directory-prefix",
        str(series_dir),
        ftp_url,
    ]
    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"wget failed for {accession} (exit {result.returncode}). " "Check network access and FTP URL."
        )

    return series_dir


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def _read_possibly_gzipped(path: Path) -> pd.DataFrame:
    """Read a TSV/CSV file, handling optional gzip compression."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(str(path), "rt") as fh:
        sample = fh.read(2048)
    sep = "\t" if "\t" in sample else ","
    logger.debug("Reading %s (sep=%r)", path.name, sep)
    try:
        return pd.read_csv(str(path), sep=sep, index_col=0)
    except Exception:
        # Try without index_col in case the first column is not an index
        return pd.read_csv(str(path), sep=sep)


def _try_parse_supplementary_files(series_dir: Path):
    """Strategy 1: Look for count/expression/matrix/processed files."""
    patterns = [
        "*count*",
        "*Count*",
        "*expression*",
        "*Expression*",
        "*matrix*",
        "*Matrix*",
        "*processed*",
        "*Processed*",
        "*norm*",
        "*Norm*",
    ]
    candidates = []
    for pattern in patterns:
        candidates.extend(series_dir.glob(pattern))
        candidates.extend(series_dir.glob(f"**/{pattern}"))  # one level deep

    # Deduplicate, prefer non-gz first for readability during debugging
    seen = set()
    unique = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            unique.append(c)
    unique.sort(key=lambda p: (p.suffix == ".gz", p.name))

    errors = []
    for path in unique:
        try:
            df = _read_possibly_gzipped(path)
            logger.info("Parsed supplementary file: %s  shape=%s", path.name, df.shape)
            return df, None  # (counts_df, metadata_df)
        except Exception as exc:
            errors.append(f"{path.name}: {exc}")

    if errors:
        logger.debug("Supplementary file parse attempts: %s", "; ".join(errors))
    return None, None


def _try_parse_series_matrix(series_dir: Path):
    """Strategy 2: Fall back to GEOparse series_matrix files."""
    if not HAS_GEOPARSE:
        return None, None

    matrix_files = list(series_dir.glob("*series_matrix*"))
    if not matrix_files:
        return None, None

    try:
        gse = GEOparse.get_GEO(filepath=str(matrix_files[0]))
        pivot = gse.pivot_samples("VALUE")
        meta_rows = []
        for gsm_name, gsm in gse.gsms.items():
            row = {"sample_id": gsm_name}
            row.update(gsm.metadata)
            meta_rows.append(row)
        meta_df = pd.DataFrame(meta_rows).set_index("sample_id")
        logger.info(
            "Parsed series_matrix: counts shape=%s, metadata shape=%s",
            pivot.shape,
            meta_df.shape,
        )
        return pivot, meta_df
    except Exception as exc:
        logger.warning("GEOparse series_matrix failed: %s", exc)
        return None, None



# ---------------------------------------------------------------------------
# HNSCC-specific parsing
# ---------------------------------------------------------------------------


def parse_hnscc_ipa(series_dir: Path):
    """Parse HNSCC IPA protein data (GSE288406).

    Reads raw counts (proteins × AOIs with BOM) and metadata, transposes
    to AOIs × proteins, and extracts segment labels from AOI name suffixes.

    Returns
    -------
    counts_df : pd.DataFrame
        AOIs × proteins.
    metadata_df : pd.DataFrame
        AOI-level metadata with 'segment' column, indexed by AOI name.
    """
    series_dir = Path(series_dir)

    counts_file = series_dir / "GSE288406_HNSCC_IPA_raw_counts.csv.gz"
    meta_file = series_dir / "GSE288406_HNSCC_IPA_raw_metadata.csv.gz"

    if not counts_file.exists():
        raise FileNotFoundError(f"IPA counts file not found: {counts_file}")

    logger.info("Reading IPA raw counts: %s", counts_file)
    counts_df = pd.read_csv(counts_file, index_col=0, encoding="utf-8-sig")
    logger.info("  Raw shape (proteins × AOIs): %d × %d", *counts_df.shape)
    counts_df = counts_df.T
    logger.info("  Transposed to AOIs × proteins: %d × %d", *counts_df.shape)

    metadata_df = None
    if meta_file.exists():
        logger.info("Reading IPA metadata: %s", meta_file)
        metadata_df = pd.read_csv(meta_file, index_col="Sample_ID")
        metadata_df["segment"] = metadata_df.index.map(lambda x: x.rsplit("_", 1)[-1] if "_" in x else "Unknown")
        metadata_df["core_id"] = metadata_df.index.map(lambda x: x.rsplit("_", 1)[0] if "_" in x else x)
        logger.info("  Metadata: %d AOIs, columns: %s", len(metadata_df), list(metadata_df.columns))
    else:
        logger.warning("IPA metadata file not found: %s", meta_file)

    return counts_df, metadata_df


def parse_hnscc_wta(series_dir: Path):
    """Parse HNSCC WTA RNA data (GSE290057).

    Reads QC-filtered counts (genes × AOIs), transposes to AOIs × genes,
    and derives segment labels from AOI name suffixes.

    Returns
    -------
    counts_df : pd.DataFrame
        AOIs × genes.
    metadata_df : pd.DataFrame
        Minimal metadata with 'segment' and 'core_id' columns.
    """
    series_dir = Path(series_dir)

    qc_file = series_dir / "GSE290057_QC.csv.gz"
    if not qc_file.exists():
        raise FileNotFoundError(f"WTA QC file not found: {qc_file}")

    logger.info("Reading WTA QC counts: %s", qc_file)
    counts_df = pd.read_csv(qc_file, index_col=0)
    logger.info("  Raw shape (genes × AOIs): %d × %d", *counts_df.shape)
    counts_df = counts_df.T
    logger.info("  Transposed to AOIs × genes: %d × %d", *counts_df.shape)

    metadata_df = pd.DataFrame(index=counts_df.index)
    metadata_df["segment"] = metadata_df.index.map(lambda x: x.rsplit("_", 1)[-1] if "_" in x else "Unknown")
    metadata_df["core_id"] = metadata_df.index.map(lambda x: x.rsplit("_", 1)[0] if "_" in x else x)
    logger.info("  Segments: %s", dict(metadata_df["segment"].value_counts()))

    return counts_df, metadata_df


# ---------------------------------------------------------------------------
# AnnData construction
# ---------------------------------------------------------------------------


def _extract_spatial_coords(meta_df: pd.DataFrame) -> pd.DataFrame | None:
    """Extract x/y spatial coordinates from a metadata DataFrame.

    Looks for columns named (case-insensitive):
      x, y, xcoord, ycoord, x_coord, y_coord, roi_x, roi_y

    Returns a DataFrame with columns ["x", "y"] indexed by AOI, or None.
    """
    if meta_df is None:
        return None

    col_lower = {c.lower(): c for c in meta_df.columns}
    x_candidates = ["x", "xcoord", "x_coord", "roi_x", "loc_x", "centroid_x"]
    y_candidates = ["y", "ycoord", "y_coord", "roi_y", "loc_y", "centroid_y"]

    x_col = next((col_lower[k] for k in x_candidates if k in col_lower), None)
    y_col = next((col_lower[k] for k in y_candidates if k in col_lower), None)

    if x_col is None or y_col is None:
        return None

    coords = meta_df[[x_col, y_col]].copy()
    coords.columns = ["x", "y"]
    coords = coords.apply(pd.to_numeric, errors="coerce")
    n_valid = coords.dropna().shape[0]
    logger.info("Found spatial coords columns ('%s', '%s'); %d/%d valid", x_col, y_col, n_valid, len(coords))
    return coords


def _extract_segment_annotations(meta_df: pd.DataFrame) -> pd.Series | None:
    """Extract segment / ROI type annotations from metadata."""
    if meta_df is None:
        return None

    col_lower = {c.lower(): c for c in meta_df.columns}
    candidates = ["segment", "roi_type", "region_type", "annotation", "tissue_type", "cell_type", "pathology"]
    for key in candidates:
        if key in col_lower:
            col = col_lower[key]
            logger.info("Using '%s' column for segment annotations", col)
            return meta_df[col].rename("segment")
    return None



def build_anndata_pair(
    ipa_dir: Path,
    wta_dir: Path,
    output_dir: Path,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Parse IPA and WTA series directories and build a matched AnnData pair.

    Parameters
    ----------
    ipa_dir:
        Local directory for the IPA protein GEO series (GSE288406).
    wta_dir:
        Local directory for the WTA RNA GEO series (GSE290057).
    output_dir:
        Directory where ``geomx_protein.h5ad`` and ``geomx_gex.h5ad`` are saved.

    Returns
    -------
    protein_adata : AnnData
        AOIs × proteins, obs contains metadata + spatial coords if available.
    gex_adata : AnnData
        AOIs × transcripts, obs contains metadata + spatial coords if available.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Parse both datasets (HNSCC-specific) ---
    logger.info("Parsing IPA protein matrix from %s", ipa_dir)
    ipa_counts, ipa_meta = parse_hnscc_ipa(ipa_dir)

    logger.info("Parsing WTA RNA matrix from %s", wta_dir)
    wta_counts, wta_meta = parse_hnscc_wta(wta_dir)

    logger.info("IPA counts: %d AOIs × %d proteins", *ipa_counts.shape)
    logger.info("WTA counts: %d AOIs × %d transcripts", *wta_counts.shape)

    # --- Intersect to common AOIs ---
    common_aois = ipa_counts.index.intersection(wta_counts.index)
    if len(common_aois) == 0:
        logger.warning(
            "No overlapping AOI IDs found between IPA (%d) and WTA (%d). "
            "Check that row indices represent the same AOI identifiers. "
            "Proceeding with all AOIs in each modality (unpaired).",
            len(ipa_counts),
            len(wta_counts),
        )
        ipa_common = ipa_counts
        wta_common = wta_counts
    else:
        logger.info(
            "Intersected to %d common AOIs (IPA had %d, WTA had %d)",
            len(common_aois),
            len(ipa_counts),
            len(wta_counts),
        )
        ipa_common = ipa_counts.loc[common_aois]
        wta_common = wta_counts.loc[common_aois]

    # --- Build AnnData objects ---
    protein_adata = ad.AnnData(
        X=ipa_common.values.astype(np.float32),
        obs=pd.DataFrame(index=ipa_common.index),
        var=pd.DataFrame(index=ipa_common.columns),
    )
    protein_adata.var["feature_types"] = "Protein"

    gex_adata = ad.AnnData(
        X=wta_common.values.astype(np.float32),
        obs=pd.DataFrame(index=wta_common.index),
        var=pd.DataFrame(index=wta_common.columns),
    )
    gex_adata.var["feature_types"] = "Gene Expression"

    # --- Attach metadata ---
    for adata, meta_df, label in [
        (protein_adata, ipa_meta, "IPA"),
        (gex_adata, wta_meta, "WTA"),
    ]:
        if meta_df is not None:
            shared_index = adata.obs.index.intersection(meta_df.index)
            if len(shared_index) > 0:
                adata.obs = adata.obs.join(meta_df.loc[shared_index], how="left")
                logger.info("%s: attached %d metadata columns", label, meta_df.shape[1])
            else:
                logger.warning(
                    "%s: metadata index does not overlap with AOI index — skipping",
                    label,
                )

    # --- Spatial coordinates ---
    for adata, meta_df, label in [
        (protein_adata, ipa_meta, "IPA"),
        (gex_adata, wta_meta, "WTA"),
    ]:
        coords = _extract_spatial_coords(meta_df)
        if coords is not None:
            shared = adata.obs.index.intersection(coords.index)
            if len(shared) > 0:
                coord_arr = coords.loc[shared].values.astype(np.float32)
                adata.obsm["spatial"] = coord_arr
                logger.info("%s: stored spatial coords in obsm['spatial']", label)
        else:
            logger.warning(
                "%s: no spatial coordinate columns found in metadata. "
                "Using sequential dummy coordinates (0, 0), (1, 0), (2, 0), ...",
                label,
            )
            n = adata.n_obs
            dummy = np.column_stack([np.arange(n), np.zeros(n)]).astype(np.float32)
            adata.obsm["spatial"] = dummy
            adata.uns["spatial_coords_are_dummy"] = True

    # --- Segment annotations ---
    for adata, meta_df, label in [
        (protein_adata, ipa_meta, "IPA"),
        (gex_adata, wta_meta, "WTA"),
    ]:
        seg = _extract_segment_annotations(meta_df)
        if seg is not None:
            shared = adata.obs.index.intersection(seg.index)
            if len(shared) > 0:
                adata.obs.loc[shared, "segment"] = seg.loc[shared].values
                logger.info("%s: attached segment annotations", label)

    # --- Global metadata ---
    protein_adata.uns["geo_accession"] = IPA_ACCESSION
    protein_adata.uns["platform"] = "GeoMx_IPA"
    protein_adata.uns["disease"] = "HNSCC"

    gex_adata.uns["geo_accession"] = WTA_ACCESSION
    gex_adata.uns["platform"] = "GeoMx_WTA"
    gex_adata.uns["disease"] = "HNSCC"

    # --- Save ---
    protein_path = output_dir / "geomx_protein.h5ad"
    gex_path = output_dir / "geomx_gex.h5ad"

    protein_adata.write_h5ad(str(protein_path))
    logger.info("Saved protein AnnData: %s  (%d AOIs × %d proteins)", protein_path, *protein_adata.shape)

    gex_adata.write_h5ad(str(gex_path))
    logger.info("Saved GEX AnnData: %s  (%d AOIs × %d transcripts)", gex_path, *gex_adata.shape)

    return protein_adata, gex_adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Download and format GeoMx HNSCC data from GEO (GSE288406 + GSE290057) " "into paired AnnData objects."
        )
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Root output directory. GEO files go into <output-dir>/{accession}/; " "h5ad files into <output-dir>/.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        default=False,
        help="Skip wget download and parse already-downloaded files.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    ipa_dir = output_dir / IPA_ACCESSION
    wta_dir = output_dir / WTA_ACCESSION

    if not args.skip_download:
        logger.info("=== Downloading IPA protein series (%s) ===", IPA_ACCESSION)
        download_geo_series(IPA_ACCESSION, output_dir)

        logger.info("=== Downloading WTA RNA series (%s) ===", WTA_ACCESSION)
        download_geo_series(WTA_ACCESSION, output_dir)
    else:
        logger.info("--skip-download set; using existing files in %s", output_dir)
        for d in (ipa_dir, wta_dir):
            if not d.exists():
                raise FileNotFoundError(
                    f"Expected directory {d} does not exist. " "Run without --skip-download to fetch the data first."
                )

    logger.info("=== Building AnnData pair ===")
    protein_adata, gex_adata = build_anndata_pair(ipa_dir, wta_dir, output_dir)

    logger.info("Done.")
    logger.info("  Protein: %s", output_dir / "geomx_protein.h5ad")
    logger.info("  GEX:     %s", output_dir / "geomx_gex.h5ad")
    return protein_adata, gex_adata


if __name__ == "__main__":
    main()
