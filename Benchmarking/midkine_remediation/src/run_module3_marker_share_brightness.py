"""Module 3 marker-share on one specimen with brightness correction ON (v41p) or OFF (gate).

Same as run_module3_marker_share.py (same profile, preprocessing, spot set); the only change is that
``run_cell_proportion_model`` is called with ``brightness_correction=<on/off>``. The GMM brightness is fit on
this specimen's own spots (same convention as the Xenium v41p producer). Per-marker GMM diagnostics
(components, brightness, fallback) are written to <sample>_brightness_diag.json.

    python run_module3_marker_share_brightness.py --sample S --output-dir DIR --brightness {on,off}
"""
from __future__ import annotations

import argparse
import functools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # import CITEgeist from THIS checkout, not the editable install
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_module3_marker_share as base  # noqa: E402

from CITEgeist.model.citegeist_model import CitegeistModel  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, choices=sorted(base.SPECIMENS))
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--brightness", required=True, choices=("on", "off"))
    a = ap.parse_args()
    on = a.brightness == "on"
    captured = {}
    orig_run = CitegeistModel.run_cell_proportion_model

    @functools.wraps(orig_run)  # keep the real signature visible to inspect.signature
    def run(self, *args, **kw):
        kw["brightness_correction"] = on
        out = orig_run(self, *args, **kw)
        captured["brightness"] = self.results.get("marker_brightness")
        captured["diag"] = self.results.get("marker_brightness_diagnostics")
        captured["n_zero"] = self.results.get("n_all_zero_marker_spots")
        return out

    CitegeistModel.run_cell_proportion_model = run

    from CITEgeist.model.deconvolution import marker_share as ms
    from CITEgeist.model.discovery.marker_gmm import fit_marker_gmm

    orig_est = ms.estimate_marker_brightness_with_diagnostics

    def est(X, var_names, markers, *args, **kw):
        out, diag = orig_est(X, var_names, markers, *args, **kw)
        Xc = ms.cosine_normalize(X)
        for mk in markers:  # same seed/max_components as the estimator -> identical fit; records n_components + means
            fit = fit_marker_gmm(Xc[:, list(var_names).index(mk)])
            diag[mk]["n_components"] = int(fit.n_components)
            diag[mk]["component_means"] = [float(v) for v in fit.means]
            diag[mk]["component_weights"] = [float(v) for v in fit.weights]
        return out, diag

    ms.estimate_marker_brightness_with_diagnostics = est
    m = base.load_unified()
    m.DATA_ROOT = base.DATA_ROOT
    m.SAMPLE_METADATA = {k: {"min_counts": v} for k, v in base.SPECIMENS.items()}
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    res = m.run_module3(a.sample, out, profile=m.UNIFIED_PROFILE, skip_gex=True, phase1_method="marker_share")
    res["profile"] = {ct: list(v.get("Major", [])) for ct, v in m.UNIFIED_PROFILE.items()}
    res["brightness_correction"] = on
    (out / f"{a.sample}_module3_meta.json").write_text(json.dumps(res, indent=2))
    (out / f"{a.sample}_brightness_diag.json").write_text(json.dumps(captured, indent=2))
    print(f"done {a.sample} brightness={a.brightness}: {res['n_spots']} spots; diag={captured.get('diag')}")


if __name__ == "__main__":
    main()
