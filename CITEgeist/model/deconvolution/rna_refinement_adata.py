"""AnnData-level wrapper and artifact writer for reference-free RNA refinement of Phase 1."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..alignment import align_to, densify
from .cd3e_anchor import anchored_pi, anchored_z, fit_cd3e_anchor
from .rna_refinement import refine_proportions, select_genes

SCHEMA_VERSION = 1


def run_rna_refinement(
    gex_adata,
    p0: pd.DataFrame,
    nuclei: pd.Series,
    *,
    c: float = 0.1,
    n_passes: int = 3,
    layer: str = "raw_counts",
    n_hvg: int = 2000,
    subtype_groups: dict[str, tuple[str, ...]] | None = None,
    protein_counts: pd.DataFrame | None = None,
    subtype_markers: dict[str, list[str]] | None = None,
    subtype_scale: str | None = None,
    anchor_marker: str | None = None,
) -> tuple:
    """Refine protein proportions ``p0`` with raw spot RNA counts and nuclei counts.

    This is the default Phase 1 method (frozen v2: multinomial, ``c=0.1``, ``n_passes=3``, ``anchor_thr=0``). Its
    output is an RNA compartment (an RNA-mass fraction exported as a cell fraction). The ``cd3e_anchor`` split is
    opt-in and experimental: it failed its validation (no region passed V1-V4 jointly). See
    ``docs/superpowers/handoffs/2026-11-01-proportion-engine-v2-decision.md``.

    ``p0`` and ``nuclei`` are aligned to ``gex_adata.obs_names`` by barcode (strict: missing, extra or
    duplicate barcodes raise; never positional). ``n_hvg`` is the number of highly variable genes kept
    (default 2000; changed only for sensitivity analyses). Returns ``(p, info)`` where ``p`` has ``p0``'s index and
    columns and ``info`` holds the refinement diagnostics plus ``selected_genes`` and ``params``.

    Optional hierarchy (Amendment A, off unless ``subtype_groups`` is given): ``subtype_groups`` names exactly one
    parent with exactly two child types, e.g. ``{"T": ("CD4pos_T_cells", "CD8pos_T_cells")}``. RNA refinement runs
    on the merged parent (its p0 is the sum of the two child p0 columns); the parent share of each spot is then split
    ``pi*T`` / ``(1-pi)*T`` with ``pi = (x_a + 0.5) / (x_a + x_b + 1)`` (Jeffreys prior). ``x`` is the spot's raw
    antibody count summed over that child's ``subtype_markers`` columns of ``protein_counts`` (barcode-indexed,
    raw integer counts, aligned strictly to the GEX spots). ``info["hier_split"]`` holds pi and the zero-evidence
    spot count; ``info["type_names"]`` and the per-type diagnostic arrays refer to the merged type set.

    ``subtype_scale="cd3e_anchor"`` (Amendment B, opt-in; needs ``anchor_marker``, the CD3E column of
    ``protein_counts``) replaces the raw-count pi: CD3E is regressed on the two child evidences over all spots of
    the sample, background (median over T-low spots) is removed, slopes scale the evidence, and
    ``pi = (z_a + 0.5) / (z_a + z_b + 1)`` (see ``cd3e_anchor``). A sample that cannot be fitted raises
    ``AnchorNotEstimable`` before any refinement; there is no fallback. ``None`` keeps the raw split above exactly.
    """
    hier = None
    if subtype_scale not in (None, "cd3e_anchor"):
        raise ValueError(f"subtype_scale must be None or 'cd3e_anchor', got {subtype_scale!r}")
    if subtype_scale is None and anchor_marker is not None:
        raise ValueError("anchor_marker was given without subtype_scale='cd3e_anchor'")
    if subtype_scale == "cd3e_anchor" and not anchor_marker:
        raise ValueError("subtype_scale='cd3e_anchor' needs anchor_marker (the CD3E column of protein_counts)")
    if subtype_groups is None:
        if protein_counts is not None or subtype_markers is not None:
            raise ValueError("protein_counts / subtype_markers were given without subtype_groups")
        if subtype_scale is not None:
            raise ValueError("subtype_scale needs subtype_groups (the CD4/CD8 hierarchy)")
    else:
        hier = _validate_hierarchy(p0, subtype_groups, protein_counts, subtype_markers, anchor_marker)
    if layer not in gex_adata.layers:
        raise ValueError(f"GEX AnnData has no layer {layer!r} with raw counts; pass layer= naming the raw-count layer")
    if not isinstance(nuclei, pd.Series):
        raise ValueError("nuclei must be a barcode-indexed pd.Series; pass a Series, not a bare array")
    spots = pd.Index(gex_adata.obs_names)
    p0_al = align_to(spots, p0, what="proportions (p0) vs GEX spots")
    if hier is not None:
        x = _subtype_evidence(spots, protein_counts, hier, anchor_marker)
        anchored = None
        if subtype_scale == "cd3e_anchor":
            a, b = hier["children"]
            fit = fit_cd3e_anchor(x[ANCHOR_COL].to_numpy(), x[a].to_numpy(), x[b].to_numpy())
            anchored_z(x[a].to_numpy(), x[b].to_numpy(), fit)  # raises AnchorNotEstimable before any refinement
            anchored = {"fit": fit, "marker": anchor_marker, "x_anchor": x[ANCHOR_COL]}
        x = x.drop(columns=ANCHOR_COL, errors="ignore")
        p0_al = _merge_children(p0_al, hier)
    n_al = align_to(spots, nuclei.to_frame("n"), what="nuclei counts vs GEX spots")["n"]
    if n_al.isna().any():
        raise ValueError("nuclei counts contain NaN")
    Y = densify(gex_adata.layers[layer])
    genes = list(gex_adata.var_names)
    idx = select_genes(Y, genes, n_top=n_hvg)
    mask = np.zeros(Y.shape[1], dtype=bool)
    mask[idx] = True
    p, diag = refine_proportions(
        Y, n_al.to_numpy(dtype=float), p0_al.to_numpy(dtype=float), c=c, n_passes=n_passes, gene_mask=mask
    )
    out = pd.DataFrame(p, index=spots, columns=p0_al.columns).loc[p0.index]
    if hier is not None:
        out, split = _split_parent(out, x.loc[p0.index], p0.columns, hier, anchored)
    info = dict(diag)
    info["selected_genes"] = [genes[i] for i in idx]
    info["params"] = {
        "c": float(c),
        "n_passes": int(n_passes),
        "layer": layer,
        "n_hvg": int(n_hvg),
        "model": diag["model"],
        "anchor_thr": diag["anchor_thr"],
    }
    info["type_names"] = list(p0_al.columns)
    if hier is not None:
        info["params"]["hierarchy"] = {
            "parent": hier["parent"],
            "children": list(hier["children"]),
            "markers": {k: list(v) for k, v in hier["markers"].items()},
            "prior": JEFFREYS_PRIOR,
        }
        if anchored is not None:
            info["params"]["hierarchy"]["subtype_scale"] = "cd3e_anchor"
            info["params"]["hierarchy"]["anchor_marker"] = anchor_marker
        info["hier_split"] = split
    return out, info


JEFFREYS_PRIOR = 0.5
ANCHOR_COL = "__anchor__"  # internal name of the CD3E evidence column


def _validate_hierarchy(p0, subtype_groups, protein_counts, subtype_markers, anchor_marker=None) -> dict:
    """Check the scope of the hierarchy (one group, two children) and its inputs; return a normalized spec."""
    if not isinstance(subtype_groups, dict) or len(subtype_groups) != 1:
        raise ValueError("subtype_groups must hold exactly one parent group (Amendment A authorizes only CD4/CD8 T)")
    ((parent, children),) = subtype_groups.items()
    children = tuple(children)
    if len(children) != 2 or children[0] == children[1]:
        raise ValueError(f"subtype group {parent!r} must have exactly 2 distinct children, got {list(children)}")
    unknown = [c for c in children if c not in p0.columns]
    if unknown:
        raise ValueError(f"subtype children {unknown} are not columns of p0 {list(p0.columns)}")
    if parent in p0.columns:
        raise ValueError(f"parent name {parent!r} collides with an existing p0 column")
    if protein_counts is None or subtype_markers is None:
        raise ValueError("subtype_groups needs both protein_counts and subtype_markers")
    if set(subtype_markers) != set(children):
        raise ValueError(
            f"subtype_markers keys {sorted(subtype_markers)} must be exactly the children {sorted(children)}"
        )
    markers = {c: list(subtype_markers[c]) for c in children}
    for c, m in markers.items():
        if not m:
            raise ValueError(f"no subtype-defining markers given for {c!r}")
    needed = {m for ms in markers.values() for m in ms} | ({anchor_marker} if anchor_marker else set())
    missing = sorted(needed - set(protein_counts.columns))
    if missing:
        raise ValueError(f"subtype markers {missing} are not columns of protein_counts")
    return {"parent": parent, "children": children, "markers": markers}


def _merge_children(p0_al: pd.DataFrame, hier: dict) -> pd.DataFrame:
    """Replace the two child columns by their sum, placed where the first child was; other columns unchanged."""
    a, b = hier["children"]
    cols = list(p0_al.columns)
    first = min(cols.index(a), cols.index(b))
    merged = p0_al[a] + p0_al[b]
    out = p0_al.drop(columns=[a, b])
    out.insert(first, hier["parent"], merged)
    return out


def _subtype_evidence(spots: pd.Index, protein_counts: pd.DataFrame, hier: dict, anchor_marker=None) -> pd.DataFrame:
    """Per-spot raw marker-count sums for each child, aligned strictly to the GEX spots.

    With ``anchor_marker`` an extra ``ANCHOR_COL`` column holds that marker's raw count (same integer validation).
    """
    cols = sorted({m for ms in hier["markers"].values() for m in ms} | ({anchor_marker} if anchor_marker else set()))
    sub = protein_counts[cols]
    if any(isinstance(dt, pd.SparseDtype) for dt in sub.dtypes):
        sub = sub.sparse.to_dense()
    if sub.index.has_duplicates:
        raise ValueError("protein counts vs GEX spots: has duplicate labels")
    sub = align_to(spots, sub.loc[sub.index.isin(spots)], what="protein counts vs GEX spots")  # extra barcodes dropped
    vals = densify(sub.to_numpy())
    if not np.all(np.isfinite(vals)) or (vals < 0).any() or not np.all(vals == np.round(vals)):
        raise ValueError(
            "protein_counts must be raw non-negative integer counts (got normalized or non-integer values)"
        )
    sub = pd.DataFrame(vals, index=sub.index, columns=cols)
    ev = pd.DataFrame({c: sub[ms].sum(axis=1) for c, ms in hier["markers"].items()}, index=spots)
    if anchor_marker:
        ev[ANCHOR_COL] = sub[anchor_marker]
    return ev


def _split_parent(p_merged: pd.DataFrame, x: pd.DataFrame, p0_cols, hier: dict, anchored: dict | None = None):
    """Split the refined parent share into the two children; restore p0's columns.

    Default: Jeffreys ratio of raw child counts. With ``anchored`` (Amendment B): the CD3E-anchored pi.
    """
    a, b = hier["children"]
    xa, xb = x[a].to_numpy(), x[b].to_numpy()
    if anchored is None:
        pi = (xa + JEFFREYS_PRIOR) / (xa + xb + 2 * JEFFREYS_PRIOR)
        evidence = xa + xb
    else:
        pi = anchored_pi(xa, xb, anchored["fit"])
        za, zb = anchored_z(xa, xb, anchored["fit"])  # from x already aligned to p0's order
        evidence = za + zb
    parent = p_merged[hier["parent"]].to_numpy()
    out = p_merged.drop(columns=[hier["parent"]])
    out[a] = pi * parent
    out[b] = (1.0 - pi) * parent
    out = out[list(p0_cols)]
    zero = evidence == 0
    split = {
        "parent": hier["parent"],
        "children": [a, b],
        "markers": {k: list(v) for k, v in hier["markers"].items()},
        "prior": JEFFREYS_PRIOR,
        "pi": pd.Series(pi, index=p_merged.index),
        "x": x.copy(),
        "n_spots": int(len(pi)),
        "zero_evidence_spots": int(zero.sum()),
        "zero_evidence_fraction": float(zero.mean()) if len(zero) else float("nan"),
    }
    if anchored is not None:
        idx = p_merged.index
        split["subtype_scale"] = "cd3e_anchor"
        split["anchor_marker"] = anchored["marker"]
        split["anchor_fit"] = anchored["fit"].to_dict()
        split["z"] = pd.DataFrame({a: za, b: zb}, index=idx)
        split["x_anchor"] = anchored["x_anchor"].loc[idx]  # spot-indexed Series, re-ordered to p0's order
    return out, split


def _jsonable(x):
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def save_refinement_artifacts(results: dict, out_dir) -> None:
    """Write artifact schema v1: p0.csv, cell_prop.csv, selected_genes.txt, rna_refinement_diagnostics.json."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    info = results["rna_refinement"]
    results["cell_prop_base"].to_csv(out_dir / "p0.csv")
    results["cell_prop"].to_csv(out_dir / "cell_prop.csv")
    (out_dir / "selected_genes.txt").write_text("\n".join(info["selected_genes"]) + "\n")
    names = info["type_names"]
    sup = np.asarray(info["supported_types"], dtype=bool)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "params": info["params"],
        "type_names": names,
        "supported_types": [names[i] for i in np.flatnonzero(sup)],
        "pass_movement": info["pass_movement"],
        "em_iters": info["em_iters"],
        "em_converged": info["em_converged"],
        "fallback": info["fallback"],
        "profile_rank": info["profile_rank"],
        "profile_condition": info["profile_condition"],
        "profile_cosine": info["profile_cosine"],
        "high_cosine_pairs": [[names[a], names[b]] for a, b in info["high_cosine_pairs"]],
        "profile_bootstrap_cv": info.get("profile_bootstrap_cv"),
        "zero_library_spots": int(np.asarray(info["zero_library_spots"]).sum()),
        "zero_nuclei_spots": int(np.asarray(info["zero_nuclei_spots"]).sum()),
        "n_selected_genes": len(info["selected_genes"]),
    }
    hs = info.get("hier_split")
    if hs is not None:
        payload["hier_split"] = {
            k: hs[k]
            for k in (
                "parent",
                "children",
                "markers",
                "prior",
                "n_spots",
                "zero_evidence_spots",
                "zero_evidence_fraction",
            )
        }
        table = hs["x"].copy()
        table.columns = [f"x_{c}" for c in table.columns]
        table["pi"] = hs["pi"]
        if "anchor_fit" in hs:
            payload["hier_split"].update(
                subtype_scale=hs["subtype_scale"], anchor_marker=hs["anchor_marker"], anchor_fit=hs["anchor_fit"]
            )
            for c in hs["z"].columns:
                table[f"z_{c}"] = hs["z"][c]
            table["x_anchor"] = hs["x_anchor"]
        table.to_csv(out_dir / "hier_split.csv")
    (out_dir / "rna_refinement_diagnostics.json").write_text(json.dumps(_jsonable(payload), indent=2))
