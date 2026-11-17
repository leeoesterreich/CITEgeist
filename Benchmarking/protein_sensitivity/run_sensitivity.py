"""Run one dataset unit across a list of stress-test conditions.

--save-layers DIR (baseline condition only): writes DIR/layers/<type>_layer.csv from the QP+SACE run and
DIR/props.csv from the first QP run (the same `pred` the SACE step consumes inputs for; QP and SACE are
deterministic, so props equal the QP+SACE run's proportions; verified by the A11 gate equality check).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import perturb as P  # noqa: E402
from protein_sensitivity import variant as V  # noqa: E402

sys.path.insert(0, str(B.REPO / "scripts"))
import provenance  # noqa: E402

ROW_KEYS = [
    "dataset",
    "unit",
    "n_types",
    "axis",
    "param",
    "level",
    "seed",
    "panel",
    "sentinel",
    "status",
    "pearson_r",
    "rmse",
    "mae",
    "jsd",
    "gex_flat_r",
    "gex_nrmse",
    "per_type_json",
    "failure_message",
]
_NAN_METRICS = {k: float("nan") for k in ["pearson_r", "rmse", "mae", "jsd", "gex_flat_r", "gex_nrmse"]}


def load_conditions(spec: str) -> list:
    """Inline JSON list (starts with '[') or a path to a JSON file. Never stats an inline JSON string."""
    spec = spec.strip()
    return json.loads(spec if spec.startswith("[") else Path(spec).read_text())


def condition_id(cond) -> str:
    """Deterministic slug of axis/level/seed/panel, e.g. thinning_0.5_s1, ablation_one_per_type, baseline."""
    parts = [cond["axis"]]
    if cond["axis"] == "joint":  # round 2 joint factorial: level is a combo slug (settings_conditions.joint_level)
        return f"joint_{cond['level']}"
    if cond["axis"] == "setting":  # string levels (e.g. nuclei_source) allowed; param is part of the identity
        return f"setting_{cond['param']}_{cond['level']}".replace(" ", "").replace("/", "-")
    if cond.get("level") is not None and not pd.isna(cond.get("level")):
        parts.append(str(float(cond["level"])))
    if cond.get("panel"):
        parts.append(str(cond["panel"]))
    if "seed" in cond and cond["seed"] is not None and cond["seed"] != -1:
        parts.append(f"s{cond['seed']}")
    return "_".join(parts).replace(" ", "").replace("/", "-")


def apply_condition(adata_cite, profile, cond, rng):
    axis = cond["axis"]
    meta = {"status": "ok", "dropped": []}
    if axis in ("baseline", "setting", "joint"):  # setting = solver/nuclei/marker-set override, no protein change
        return adata_cite, profile, meta
    names = list(adata_cite.var_names)
    X = np.asarray(adata_cite.X.todense() if hasattr(adata_cite.X, "todense") else adata_cite.X)
    if axis == "thinning":
        adata_cite.X = P.thin_counts(X, cond["level"], rng).astype(float)
    elif axis == "dropout":
        Y, dropped = P.dropout_channels(X, names, P.qp_active_markers(profile), cond["level"], rng)
        adata_cite.X = Y.astype(float)
        meta["dropped"] = dropped
    elif axis == "noise":
        bg = P.estimate_background(X)
        adata_cite.X = P.add_background_noise(X, cond["level"], bg, rng).astype(float)
    elif axis in ("ablation", "type_removal"):
        full_active = set(P.qp_active_markers(profile))
        profile = P.active_profile(cond["profile"])  # JSON turns Soft tuples into lists
        if not profile:
            meta["status"] = "failed_no_markers"
            return adata_cite, profile, meta
        removed = sorted(full_active - set(P.qp_active_markers(profile)))
        adata_cite = adata_cite[:, [n for n in names if n not in removed]].copy()
        meta["dropped"] = removed
    else:
        raise ValueError(f"unknown axis {axis}")
    return adata_cite, profile, meta


def evaluate_unit(
    dataset, unit, n_types, conditions, *, loader, qp_fn, sace_fn, prop_metrics_fn, gex_metrics_fn, profile, on_row=None, save_pred=None, save_layers=None, save_diag=None, on_row_meta=None, on_perturbed=None
):
    gex, cite = loader(unit)
    rows = []
    for cond in conditions:
        rng = np.random.default_rng(cond.get("seed", 0))
        row = {
            "dataset": dataset,
            "unit": unit,
            "n_types": n_types,
            "axis": cond["axis"],
            "param": cond.get("param", ""),
            "level": cond.get("level", float("nan")),
            "seed": cond.get("seed", -1),
            "panel": cond.get("panel", ""),
            "sentinel": cond.get("sentinel", ""),
            **_NAN_METRICS,
            "per_type_json": "{}",
            "failure_message": "",
        }
        c2, prof, meta = apply_condition(cite.copy(), profile, cond, rng)
        row["status"] = meta["status"]
        if on_perturbed is not None:
            on_perturbed(condition_id(cond), c2 if meta["status"] == "ok" else None)
        if meta["status"] == "ok":
            try:
                extra = (cond,) if (cond["axis"] in ("setting", "joint") or cond.get("diagnostics")) else ()
                pred = qp_fn(gex.copy(), c2, prof, unit, *extra)
                pm = prop_metrics_fn(pred, unit)
                row.update({k: pm[k] for k in ["pearson_r", "rmse", "mae", "jsd"]})
                row["per_type_json"] = json.dumps(pm.get("per_type", {}))
                if cond.get("sace"):
                    layers = sace_fn(gex.copy(), c2, prof, unit)
                    row.update(gex_metrics_fn(layers, unit))
                    if save_layers is not None:
                        try:  # side output: a save failure must never wipe computed metrics
                            save_layers(condition_id(cond), layers, pred)
                        except Exception as e:
                            print(f"WARNING: save_layers failed for {condition_id(cond)}: {e!r}", file=sys.stderr)
                if save_diag is not None and getattr(pred, "attrs", {}).get("diagnostics") is not None:
                    try:  # side output
                        save_diag(condition_id(cond), pred.attrs["diagnostics"])
                    except Exception as e:
                        print(f"WARNING: save_diag failed for {condition_id(cond)}: {e!r}", file=sys.stderr)
                if save_pred is not None:
                    try:  # side output: an I/O error must never wipe computed metrics
                        save_pred(condition_id(cond), pred)
                    except Exception as e:
                        print(f"WARNING: save_pred failed for {condition_id(cond)}: {e!r}", file=sys.stderr)
            except Exception as e:  # any failure on a perturbed input is a recorded result
                msg = str(e).split(".")[0][:60].replace(" ", "_")
                row["failure_message"] = f"{type(e).__name__}: {e}"  # full, untruncated; status above stays 60-char slug
                row["status"] = "failed_" + (msg if isinstance(e, ValueError) else f"{type(e).__name__}_{msg}")
                row.update(_NAN_METRICS)
                row["per_type_json"] = "{}"
        out_row = {k: row[k] for k in ROW_KEYS}
        rows.append(out_row)
        if on_row_meta is not None:
            on_row_meta(
                {
                    "condition_id": condition_id(cond),
                    "status": row["status"],
                    "dropped": list(meta["dropped"]),
                    "active_profile": {t: {"Major": list(sp.get("Major", []))} for t, sp in prof.items()},
                }
            )
        if on_row is not None:
            on_row(out_row)
    return rows


def _build_wiring(dataset: str, n_types: int, variant=None) -> dict:
    """Wire the real baseline.* functions; unit strings are converted here."""
    if dataset == "xenium" and variant is not None:
        import benchmark_paths as bp

        gt_root = bp.gt_dir(n_types, "rnaprot") / "ground_truth"
        gex_root = bp.gt_gex_dir("rnaprot")
        return dict(
            loader=lambda u: B.load_xenium(int(u)),
            qp_fn=lambda g, c, prof, u, *extra: B.run_xenium_qp(
                g, c, prof, int(u), setting=extra[0] if extra else None, phase1="marker_share", variant=variant
            ),
            sace_fn=lambda g, c, prof, u: B.run_xenium_sace(g, c, prof, int(u), phase1="marker_share", variant=variant),
            prop_metrics_fn=lambda pred, u: B.xenium_prop_metrics(pred, int(u), n_types, gt_root=gt_root),
            gex_metrics_fn=lambda layers, u: B.xenium_gex_metrics(layers, int(u), gex_gt_dir=gex_root),
            profile=V.profile_for(n_types, variant),
        )
    if dataset == "xenium":
        return dict(
            loader=lambda u: B.load_xenium(int(u)),
            qp_fn=lambda g, c, prof, u, *extra: B.run_xenium_qp(
                g, c, prof, int(u), setting=extra[0] if extra else None
            ),
            sace_fn=lambda g, c, prof, u: B.run_xenium_sace(g, c, prof, int(u)),
            prop_metrics_fn=lambda pred, u: B.xenium_prop_metrics(pred, int(u), n_types),
            gex_metrics_fn=lambda layers, u: B.xenium_gex_metrics(layers, int(u)),
            profile=B.PROFILE[n_types],
        )

    def _split(u):
        cond, rep = u.split(":")
        return cond, int(rep)

    return dict(
        loader=lambda u: B.load_sim(*_split(u)),
        qp_fn=lambda g, c, prof, u, *extra: B.run_sim(g, c, prof, sace=False)[0],
        sace_fn=lambda g, c, prof, u: B.run_sim(g, c, prof, sace=True)[1],
        prop_metrics_fn=lambda pred, u: B.sim_prop_metrics(pred, *_split(u)),
        gex_metrics_fn=lambda layers, u: B.sim_gex_metrics(layers, *_split(u)),
        profile=B.sim_profile(),
    )


def _git_hash() -> str:
    return provenance.code_commit(B.REPO)


def _code_dirty():
    """True/False from git; None when there is no git checkout (a clean tree cannot be verified)."""
    if not provenance.in_git_checkout(B.REPO):
        return None
    out = subprocess.run(
        ["git", "--no-optional-locks", "status", "--porcelain", "-uno"],
        capture_output=True,
        text=True,
        cwd=str(B.REPO),
    ).stdout
    return bool(out.strip())


_LEGACY_ARGS = ("dataset", "unit", "n_types", "conditions", "out", "save_layers")


def build_run_meta(args, conditions_text: str, variant) -> dict:
    """run_meta.json payload. Legacy keys unchanged; with a variant adds run_variant + code_dirty."""
    HERE_P = Path(__file__).resolve().parent
    protocol = HERE_P / "protocol.json"
    settings_protocol = HERE_P / "settings_protocol.json"
    v2 = [p for p in (HERE_P / "settings_protocol_v2.json", HERE_P / "settings_protocol_v2.DRAFT.json") if p.exists()]
    meta = {
        "git_hash": _git_hash(),
        "conditions_sha256": hashlib.sha256(conditions_text.encode()).hexdigest(),
        "settings_protocol_v2_md5": hashlib.md5(v2[0].read_bytes()).hexdigest() if v2 else None,
        "settings_protocol_v2_file": v2[0].name if v2 else None,
        "settings_protocol_md5": hashlib.md5(settings_protocol.read_bytes()).hexdigest()
        if settings_protocol.exists()
        else None,
        "protocol_md5": hashlib.md5(protocol.read_bytes()).hexdigest() if protocol.exists() else None,
        "args": {k: getattr(args, k) for k in _LEGACY_ARGS} if variant is None else vars(args),
    }
    if variant is not None:
        meta["run_variant"] = V.meta(variant, args.n_types)
        meta["code_dirty"] = _code_dirty()
        meta["code_commit_source"] = provenance.code_commit_info(B.REPO)[1]
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True, choices=["xenium", "sim"])
    ap.add_argument("--unit", required=True, help="region (xenium) or <condition>:<rep> (sim)")
    ap.add_argument("--n-types", type=int, default=6, choices=[6, 7])
    ap.add_argument("--conditions", required=True, help="JSON list, or path to a JSON file")
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--save-layers",
        default=None,
        help="dir; for sace conditions write <dir>/layers/<type>_layer.csv and <dir>/props.csv (baseline only)",
    )
    ap.add_argument("--profile", default="achievable", choices=list(V.PROFILES))
    ap.add_argument("--brightness", type=int, default=0, choices=[0, 1])
    args = ap.parse_args(argv)
    variant = V.from_args(args.profile, args.brightness)
    if variant is not None and args.dataset != "xenium":
        ap.error("--profile ms_v41 / --brightness 1 apply to the xenium dataset only")

    conditions = load_conditions(args.conditions)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    metrics_path = out / "metrics.csv"
    pd.DataFrame(columns=ROW_KEYS).to_csv(metrics_path, index=False)

    def _append(r):
        pd.DataFrame([r], columns=ROW_KEYS).to_csv(metrics_path, mode="a", header=False, index=False)

    pred_dir = out / "predictions"

    def _save_pred(cid, pred):
        pred_dir.mkdir(exist_ok=True)
        pd.DataFrame(pred).to_csv(pred_dir / f"{cid}.csv.gz")

    diag_dir = out / "diagnostics"

    def _save_diag(cid, diag):
        diag_dir.mkdir(exist_ok=True)
        (diag_dir / f"{cid}.json").write_text(json.dumps(diag))

    layers_dir = Path(args.save_layers) if args.save_layers else None

    def _save_layers(cid, layers, pred):
        if cid != "baseline":
            return
        (layers_dir / "layers").mkdir(parents=True, exist_ok=True)
        for ct, df in layers.items():
            df.to_csv(layers_dir / "layers" / f"{ct}_layer.csv")
        pd.DataFrame(pred).to_csv(layers_dir / "props.csv")

    hooks = {}
    if variant is not None:
        row_meta_path = out / "row_meta.jsonl"
        row_meta_path.write_text("")
        perturbed = {}

        def _on_perturbed(cid, c2):
            perturbed[cid] = V.perturbed_sha256(c2) if c2 is not None else None

        def _on_row_meta(rm):
            rm = {**rm, "perturbed_sha256": perturbed.get(rm["condition_id"]), "variant": V.meta(variant, args.n_types)}
            with open(row_meta_path, "a") as fh:
                fh.write(json.dumps(rm) + "\n")

        hooks = {"on_row_meta": _on_row_meta, "on_perturbed": _on_perturbed}

    evaluate_unit(
        args.dataset,
        args.unit,
        args.n_types,
        conditions,
        on_row=_append,
        save_pred=_save_pred,
        save_layers=_save_layers if layers_dir else None,
        save_diag=_save_diag,
        **hooks,
        **_build_wiring(args.dataset, args.n_types, variant),
    )

    cond_text = (
        args.conditions.strip() if args.conditions.strip().startswith("[") else Path(args.conditions).read_text()
    )
    meta = build_run_meta(args, cond_text, variant)
    (out / "run_meta.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
