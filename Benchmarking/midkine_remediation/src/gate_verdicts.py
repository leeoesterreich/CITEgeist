"""Verdict functions over the stamped G6/G7 synthetic result files (fail closed: missing/empty/tampered -> raise)."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
import provenance  # noqa: E402

MAX_TYPE1 = 0.07
SYNTH = HERE.parent / "synthetic"
G6_NULL_CASES = ("null_autocorrelated", "null_noisy_nst", "null_depth_confounded")     # protocol-required cases
G6_STRESS_CASES = ("null_autocorrelated_stress",)                                      # reported, not gating


def load_stamped(path: Path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    ok, msg = provenance.verify(path)
    if not ok:
        raise ValueError(f"{path.name}: {msg}")
    return json.loads(path.read_text())["result"]


def _rates(res: dict) -> Dict[str, float]:
    out = {}
    for part in ("magnitude", "participation"):
        if res[part]["n_ok"] <= 0:
            raise ValueError("no estimable datasets: gate cannot pass vacuously")
        out[part] = res[part]["type1"]
    return out


def g6_type1(prefix: str, cases=G6_NULL_CASES, synth: Path = SYNTH) -> Dict[str, dict]:
    """{case: {magnitude, participation, pass}} for result files g6_<prefix>_<case>.json."""
    if not cases:
        raise ValueError("empty case list")
    out = {}
    for c in cases:
        r = load_stamped(synth / f"g6_{prefix}_{c}.json")
        rates = _rates(r)
        out[c] = {**rates, "n_datasets": r["n_datasets"], "pass": all(v <= MAX_TYPE1 for v in rates.values())}
    return out


def g6_injection(name: str = "inj", synth: Path = SYNTH) -> dict:
    r = load_stamped(synth / f"g6_{name}.json")
    return {**r, "pass": bool(r["recovery_rate"] >= 0.90)}


def g7_definition_c(synth: Path = SYNTH) -> dict:
    r = load_stamped(synth / "g7_definition_c.json")
    rows = r["settings"]
    if not rows:
        raise ValueError("no settings")
    msr = [x["msr_raw"]["type1_two_sided"] for x in rows]
    hi = [x["naive_raw"]["type1_two_sided"] for x in rows if x["rho"] == max(y["rho"] for y in rows)]
    return {"n_settings": len(rows), "msr_max": max(msr), "msr_pass": max(msr) <= MAX_TYPE1,
            "naive_mean_high_rho": sum(hi) / len(hi), "naive_inflated": sum(hi) / len(hi) > MAX_TYPE1}


FIX_VARIANTS = ("v_b3_nocov", "v_b3_cov", "v_b6_nocov", "v_b6_cov")
FIX_CASES = ("null_autocorrelated", "null_autocorrelated_stress", "null_noisy_nst_kshare", "null_depth_confounded")


def g6fix_rule(synth: Path = SYNTH / "g6_fix") -> Dict[str, dict]:
    """The pre-declared selection rule (g6_fix_round_selection_rule.md): per variant, pass iff every null case has
    magnitude AND participation type-I <= 0.07 and injection recovery >= 0.80 with BH false positives <= 0.10."""
    out = {}
    for v in FIX_VARIANTS:
        fails = []
        for c in FIX_CASES:
            r = load_stamped(synth / f"t1_{v}_{c}.json")
            for part in ("magnitude", "participation"):
                if r[part]["n_ok"] <= 0:
                    raise ValueError("no estimable datasets")
                if r[part]["type1"] > MAX_TYPE1:
                    fails.append(f"{c}:{part}={r[part]['type1']:.3f}")
        inj = load_stamped(synth / f"inj_{v}.json")
        if inj["recovery_rate"] < 0.80:
            fails.append(f"injection recovery {inj['recovery_rate']:.2f}")
        if inj["bh_false_positive_rate_other_pairs"] > 0.10:
            fails.append(f"BH false positives {inj['bh_false_positive_rate_other_pairs']:.2f}")
        out[v] = {"pass": not fails, "fails": fails}
    return out
