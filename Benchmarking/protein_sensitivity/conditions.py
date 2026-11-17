"""Expand protocol.json into the list of stress-test conditions (JSON string)."""
from __future__ import annotations

import json
from pathlib import Path

PROTOCOL = Path(__file__).resolve().parent / "protocol.json"


def profile_key(dataset: str, n_types: int) -> str:
    return "sim" if dataset == "sim" else f"xenium{n_types}"


def get_profile(dataset: str, n_types: int, variant=None) -> dict:
    from protein_sensitivity import baseline as B

    if dataset == "sim":
        return B.sim_profile()
    if variant is None:
        return B.PROFILE[n_types]
    from protein_sensitivity import variant as V

    return V.profile_for(n_types, variant)


def _sace_allowed(dataset: str, n_types: int, sentinel: str) -> bool:
    # SACE runs only on sentinel conditions, and never on Xenium 6-type.
    return sentinel != "" and not (dataset == "xenium" and n_types == 6)


def build(mode: str, dataset: str, n_types: int = 6, variant=None) -> str:
    from protein_sensitivity import perturb as P

    if mode == "gate":
        return json.dumps(
            [{"axis": "baseline", "sentinel": "baseline", "sace": _sace_allowed(dataset, n_types, "baseline")}]
        )
    if mode != "sweep":
        raise ValueError(f"unknown mode {mode}")

    proto = json.loads(PROTOCOL.read_text())
    key = profile_key(dataset, n_types)
    profile = get_profile(dataset, n_types, variant)

    # sentinel lookup: (axis, level-or-panel) -> tag
    tags: dict[tuple, str] = {}
    for tag, items in proto["sentinels"].items():
        for it in items:
            ax = it["axis"]
            if "level" in it:
                tags[(ax, float(it["level"]))] = tag
            elif "panel" in it:
                panel = it["panel"]
                if isinstance(panel, dict):
                    panel = panel[key]
                tags[(ax, panel)] = tag
            else:
                tags[(ax, None)] = tag

    def cond(**kw):
        sent = kw.pop("sentinel", "")
        return {**kw, "sentinel": sent, "sace": _sace_allowed(dataset, n_types, sent)}

    out = [cond(axis="baseline", sentinel=tags[("baseline", None)])]
    for ax in ("thinning", "dropout", "noise"):
        for lvl in proto["axes"][ax]:
            for seed in proto["seeds"]:
                out.append(cond(axis=ax, level=lvl, seed=seed, sentinel=tags.get((ax, float(lvl)), "")))
    for ax, panels in (
        ("ablation", P.ablation_panels(profile) + P.leave_one_type_out(profile)),
        ("type_removal", P.type_removal_series(profile)),
    ):
        for name, prof in panels:
            out.append(cond(axis=ax, panel=name, profile=prof, sentinel=tags.get((ax, name), "")))
    return json.dumps(out)
