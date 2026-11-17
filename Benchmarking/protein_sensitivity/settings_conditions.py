"""Expand settings_protocol.json into settings-sensitivity conditions and resolve them to solver inputs.

Pure module: heavy imports are avoided so it is cheap on CPU and unit-testable. Data paths are
constants; loaders take explicit frames/paths where tests need to inject synthetic data.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SETTINGS_PROTOCOL = HERE / "settings_protocol.json"
XENIUM_DATA = REPO / "Benchmarking/xenium_pseudovisium/data_protein_gt"
STARDIST_DIR = REPO / "Benchmarking/xenium_benchmarking/CITEgeist/output/stardist_counts"

ON_OFF = {"on": True, "off": False}
PROTOCOL_V2_DRAFT = HERE / "settings_protocol_v2.DRAFT.json"
PROTOCOL_V2 = HERE / "settings_protocol_v2.json"  # name after freeze
JOINT_KEYS = ("alpha", "lambda_reg", "lambda_laplacian", "detection_gate_ub")
DETECTION_VARIANT_PARAMS = ("gmm_n_components", "gmm_base_threshold", "gmm_adaptive_threshold")
JOINT_DEFAULTS = {"alpha": 0.5, "lambda_reg": 1, "lambda_laplacian": 0.0, "detection_gate_ub": 0.01}  # production values (see v2 protocol)


def load_protocol(path: Path | str = SETTINGS_PROTOCOL) -> dict:
    return json.loads(Path(path).read_text())


def protocol_md5(path: Path | str = SETTINGS_PROTOCOL) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def norm_level(level) -> str:
    """Canonical string for a level so 2, 2.0 and "2" compare equal (CSV round trips change the type)."""
    try:
        return format(float(level), "g")
    except (TypeError, ValueError):
        return str(level)


def build_list(n_types: int = 6, proto: dict | None = None) -> list[dict]:
    """Baseline + one-at-a-time non-baseline setting conditions (order = protocol order)."""
    proto = proto or load_protocol()
    base = proto["baseline"]["levels"]
    out = [{"axis": "baseline", "sentinel": "", "sace": False}]
    for param, spec in proto["settings"].items():
        for lvl in spec["levels"]:
            if lvl == base[param] and type(lvl) is type(base[param]):
                continue
            out.append({"axis": "setting", "param": param, "level": lvl, "sentinel": "", "sace": False})
    return out


def build(n_types: int = 6) -> str:
    return json.dumps(build_list(n_types))


def qp_overrides(param: str, level) -> dict:
    """Keyword overrides for run_cell_proportion_model for one (param, level).

    `markers` and `nuclei_source` do not map to a solver kwarg here (profile / nuclei_counts are
    resolved by `resolve`), so they return {}.
    """
    if param in ("detection_gate_ub", "lambda_laplacian", "lambda_reg"):
        return {param: level}
    if param in ("use_detection_gating", "use_gex_detection"):
        if level not in ON_OFF:
            raise ValueError(f"{param}: level must be on/off, got {level!r}")
        return {param: ON_OFF[level]}
    if param in ("fusion_mode", "gex_detection_k", "gex_detection_min_corr"):  # round 2: genuine production kwargs
        return {param: level}
    if param in ("markers", "nuclei_source") or param in DETECTION_VARIANT_PARAMS:
        return {}
    raise ValueError(f"unknown setting param {param}")


def promote_minor(profile: dict) -> dict:
    """Return a copy with every Minor marker appended to Major (QP ignores Minor). Order kept, no dups."""
    out = copy.deepcopy(profile)
    for info in out.values():
        major = list(info.get("Major", []))
        for m in info.get("Minor", []):
            if m not in major:
                major.append(m)
        info["Major"] = major
        info["Minor"] = []
    return out


def _mapped(mapping: pd.DataFrame, region: int) -> pd.DataFrame:
    return mapping[(mapping["region_id"] == region) & (mapping["spot_idx"] != -1)]


def label_derived_counts(mapping: pd.DataFrame, cell_types: pd.DataFrame, region: int) -> pd.Series:
    """Per-spot count of mapped cells whose type is not 'Unknown' (the published nuclei_counts)."""
    known = _mapped(mapping, region).merge(cell_types, left_index=True, right_index=True)
    return known[known["cell_type"] != "Unknown"].groupby("spot_id").size()


def all_segmented_counts(mapping: pd.DataFrame, region: int) -> pd.Series:
    """Per-spot count of ALL mapped cells, Unknown included."""
    return _mapped(mapping, region).groupby("spot_id").size()


def stardist_counts(region: int, directory: Path | str = STARDIST_DIR) -> pd.Series:
    """Precomputed StarDist-DAPI nuclei counts, indexed by spot_<spot_idx> (same ids as spot_id)."""
    csv = Path(directory) / f"stardist_counts_region_{region}.csv"
    if not csv.exists():
        raise FileNotFoundError(f"StarDist counts not found: {csv}")
    return pd.read_csv(csv, index_col=0)["nuclei_count_raw"]


def nuclei_counts_for(source: str, region: int, loader_published) -> pd.Series | None:
    """Resolve a nuclei_source level to a per-spot Series (or None). `loader_published(region)` is the
    published label-derived loader, injected so baseline behaviour is bit-for-bit what it was."""
    if source == "label_derived":
        return loader_published(region)
    if source == "all_segmented":
        mapping = pd.read_csv(XENIUM_DATA / "cell_to_spot_mapping.csv", index_col=0)
        return all_segmented_counts(mapping, region)
    if source == "stardist_dapi":
        return stardist_counts(region)
    if source == "none":
        return None
    raise ValueError(f"unknown nuclei_source {source}")


def resolve(setting: dict | None, region: int, profile: dict, loader_published) -> tuple[dict, dict]:
    """(profile, run_cell_proportion_model kwargs) for a condition. setting=None -> published call."""
    kw = {"lambda_laplacian": 0.0}
    nuclei = None
    use_source = "label_derived"
    if setting and setting.get("axis") == "joint":
        kw.update({k: setting[k] for k in JOINT_KEYS})
    if setting and setting.get("axis") == "setting":
        param, level = setting["param"], setting["level"]
        if param == "markers":
            if level == "major_plus_minor":
                profile = promote_minor(profile)
            elif level != "major":
                raise ValueError(f"unknown markers level {level}")
        elif param == "nuclei_source":
            use_source = level
        else:
            kw.update(qp_overrides(param, level))
    nuclei = nuclei_counts_for(use_source, region, loader_published)
    kw["nuclei_counts"] = nuclei
    return profile, kw


# ---------------------------------------------------------------------------------------------
# Round 2 (settings_protocol_v2): condition expansion, benchmark-only detection variants, diagnostics
# ---------------------------------------------------------------------------------------------
def load_protocol_v2(path: Path | str | None = None) -> dict:
    p = Path(path) if path else (PROTOCOL_V2 if PROTOCOL_V2.exists() else PROTOCOL_V2_DRAFT)
    return json.loads(p.read_text())


def joint_level(alpha, lambda_reg, lambda_laplacian, detection_gate_ub) -> str:
    return f"a{norm_level(alpha)}_r{norm_level(lambda_reg)}_l{norm_level(lambda_laplacian)}_g{norm_level(detection_gate_ub)}"


def build_list_v2(proto: dict | None = None) -> list[dict]:
    """Baseline (+diagnostics) + OAT detection settings + joint factorial minus the baseline cell."""
    proto = proto or load_protocol_v2()
    out = [{"axis": "baseline", "diagnostics": True, "sentinel": "", "sace": False}]
    for param, spec in proto["settings"].items():
        for lvl in spec["levels_non_baseline"]:
            out.append(
                {"axis": "setting", "param": param, "level": lvl, "family": spec["family"], "sentinel": "", "sace": False}
            )
    jb = proto["joint_block"]
    base = tuple(norm_level(jb["baseline_cell"][k]) for k in JOINT_KEYS)
    import itertools

    for a, r, l, g in itertools.product(*(jb["levels"][k] for k in JOINT_KEYS)):
        if (norm_level(a), norm_level(r), norm_level(l), norm_level(g)) == base:
            continue  # the baseline condition stands in for this cell
        out.append(
            {
                "axis": "joint", "param": "joint", "level": joint_level(a, r, l, g),
                "alpha": a, "lambda_reg": r, "lambda_laplacian": l, "detection_gate_ub": g,
                "sentinel": "", "sace": False,
            }  # fmt: skip
        )
    return out


def build_v2() -> str:
    return json.dumps(build_list_v2())


def detection_variant(setting: dict | None) -> dict | None:
    """Benchmark-only protein-GMM variant for a condition, or None. Keys: n_components, threshold, adaptive_threshold."""
    if not setting or setting.get("axis") != "setting":
        return None
    p, lvl = setting["param"], setting["level"]
    if p == "gmm_n_components":
        return {"n_components": int(lvl)}
    if p == "gmm_base_threshold":
        return {"threshold": float(lvl)}
    if p == "gmm_adaptive_threshold":
        if lvl != "off":
            raise ValueError("gmm_adaptive_threshold only supports 'off'")
        return {"adaptive_threshold": False}
    return None


class ThresholdRecorder:
    """Passive record of every _compute_adaptive_threshold call: (signal weight, base threshold, effective threshold)."""

    def __init__(self):
        self.calls: list[tuple[float, float, float]] = []


def make_detect_wrapper(orig, variant: dict):
    """Wrapper for detection.detect_cell_types that REPLACES the existing threshold=/adaptive_threshold= keywords
    (production passes threshold=0.5 and adaptive_threshold=True by keyword, citegeist_model.py:989-995)."""

    def wrapped(*args, **kwargs):
        if "threshold" in variant:
            kwargs["threshold"] = variant["threshold"]
        if "adaptive_threshold" in variant:
            kwargs["adaptive_threshold"] = variant["adaptive_threshold"]
        return orig(*args, **kwargs)

    return wrapped


def make_gmm_factory(orig, n_components: int):
    """GaussianMixture replacement forcing n_components (production passes n_components=2 by keyword)."""

    def factory(*args, **kwargs):
        kwargs["n_components"] = n_components
        return orig(*args, **kwargs)

    return factory


def detection_context(variant: dict | None, recorder: ThresholdRecorder | None):
    """Context manager applying the benchmark-only protein-GMM variant and the passive recorder to
    CITEgeist.model.deconvolution.detection for the duration of ONE condition (restored on exit, also on error).
    Does not touch production source; works because citegeist_model imports these names inside the function and
    detection.py resolves GaussianMixture / _compute_adaptive_threshold at call time."""
    import contextlib

    @contextlib.contextmanager
    def _cm():
        from CITEgeist.model.deconvolution import detection as D

        saved = (D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold)
        try:
            if variant and ("threshold" in variant or "adaptive_threshold" in variant):
                D.detect_cell_types = make_detect_wrapper(saved[0], variant)
            if variant and "n_components" in variant:
                D.GaussianMixture = make_gmm_factory(saved[1], variant["n_components"])
            if recorder is not None:
                orig_thr = saved[2]

                def rec(gmm, signal_cluster, base_threshold=0.5):
                    t = orig_thr(gmm, signal_cluster, base_threshold=base_threshold)
                    recorder.calls.append((float(gmm.weights_[signal_cluster]), float(base_threshold), float(t)))
                    return t

                D._compute_adaptive_threshold = rec
            yield
        finally:
            D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold = saved

    return _cm()


def _mask_hash(mask) -> str | None:
    import numpy as np

    if mask is None:
        return None
    m = np.ascontiguousarray(np.asarray(mask, dtype=bool))
    return hashlib.sha256(m.shape.__repr__().encode() + np.packbits(m).tobytes()).hexdigest()


def gex_candidate_counts(H, k: int, min_corr: float) -> list[int]:
    """Per-type number of qualifying GEX genes, replicating detect_cell_types_gex gene selection
    (detection_refinement.py:141-157). < 3 means the all-True fallback (no GEX gate) for that type."""
    import numpy as np

    H = np.asarray(H)
    best = np.argmax(H, axis=0)
    counts = []
    for t in range(H.shape[0]):
        n = 0
        for g in np.argsort(-H[t]):
            if H[t][g] < min_corr:
                break
            if best[g] != t:
                continue
            n += 1
            if n >= k:
                break
        counts.append(n)
    return counts


def collect_diagnostics(results: dict, qp_kwargs: dict, recorder: ThresholdRecorder | None) -> dict:
    """Detection-mask hashes, effective thresholds and GEX candidate counts for one condition."""
    k = qp_kwargs.get("gex_detection_k", 10)
    mc = qp_kwargs.get("gex_detection_min_corr", 0.15)
    diag = {
        "mask_hash_protein": _mask_hash(results.get("detection_protein")),
        "mask_hash_gex": _mask_hash(results.get("detection_gex")),
        "mask_hash_fused": _mask_hash(results.get("detection_fused")),
        "effective_thresholds": [list(c) for c in (recorder.calls if recorder else [])],
    }
    H = results.get("gex_detection_corr")
    if H is not None:
        counts = gex_candidate_counts(H, k, mc)
        diag["gex_candidate_counts"] = counts
        diag["gex_allTrue_fallback_types"] = [i for i, c in enumerate(counts) if c < 3]
    fused = results.get("detection_fused")
    if fused is not None:
        import numpy as np

        diag["fused_detection_rate_per_type"] = [float(x) for x in np.asarray(fused, dtype=bool).mean(axis=0)]
    return diag
