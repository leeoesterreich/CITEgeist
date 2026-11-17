"""H1 claim tree, DESCRIPTIVE regime (protocol rev 6; decision_rules.hurdle_branches_for_H1_per_part; spec 9.4.8, G10).

SCOPE: the HEADLINE specimen only (HCC22-088-P4-S2_1i_rep). The inputs ``c`` (agreement across N-P, N-A, N-S; runs R02
exist only for P4-S2) and ``f`` (nuclei-free R11 agreement; P4-S2 is below the segmentation gate, so ``f`` is always
required) are defined only there. The other 11 specimens get plain descriptive reporting and no claim tree.

No p-value, q-value, interval property or control is an input: the wording does not depend on whether the nominal
interval excludes 1. Each of the 5 conditions is a standalone Boolean function of the inputs. ``evaluate`` returns EVERY
condition that holds; ``decide`` requires exactly one (no first-match ordering). ``enumerate_inputs`` yields all
2^4 * 3 = 48 combinations of (est, s, c, d in {agree, disagree, not_estimable}, f). Applied per part (magnitude GMR,
participation odds ratio).
"""
from __future__ import annotations

import itertools
from collections import namedtuple
from typing import Callable, Dict, Iterator, List, Optional

Inputs = namedtuple("Inputs", "est s c d f")
D_VALUES = ("agree", "disagree", "not_estimable")
SCOPE_SPECIMEN = "HCC22-088-P4-S2_1i_rep"
BRANCH_IDS = (
    "not_estimable",
    "not_higher",
    "direction_not_robust",
    "positive_between_spot_agrees",
    "positive_native_only",
)


def _robust(i: Inputs) -> bool:
    """P4-S2 is below the segmentation gate (fixed by preflight): robustness needs the normalisations (c) AND the nuclei-free run (f)."""
    return i.c and i.f


CONDITIONS: Dict[str, Callable[[Inputs], bool]] = {
    "not_estimable": lambda i: not i.est,
    "not_higher": lambda i: i.est and not i.s,
    "direction_not_robust": lambda i: i.est and i.s and not _robust(i),
    "positive_between_spot_agrees": lambda i: i.est and i.s and _robust(i) and i.d == "agree",
    "positive_native_only": lambda i: i.est and i.s and _robust(i) and i.d != "agree",
}


# The wording each branch allows (single code source; protocol decision_rules.hurdle_branches_for_H1_per_part must carry the same text). {estimand} is filled by the caller.
ALLOWED_CLAIM: Dict[str, str] = {
    "not_estimable": "could not be estimated for {estimand} (state group sizes and the typed reason_code); remove the claim",
    "not_higher": "the point estimate was not above 1 for {estimand} (report the estimate, and the nominal interval with its measured coverage for the magnitude); no supportive wording",
    "direction_not_robust": "positive point estimate for {estimand}, but the sign is not stable across normalisation or segmentation; no claim of higher signalling",
    "positive_between_spot_agrees": "higher in the native COMMOT analysis for {estimand}, and the between-spot-only rerun agreed in direction (point estimate above 1); descriptive; interval (magnitude) is nominal with its measured coverage stated",
    "positive_native_only": "higher in the native COMMOT analysis for {estimand}; the between-spot-only rerun did not agree or was not estimable; no claim of signalling between spots; descriptive",
}
# The condition of each branch as the Boolean text used in the protocol JSON and the spec (env: est, s, c, f booleans; d one of agree/disagree/not_estimable; robust = c and f).
CONDITION_TEXT: Dict[str, str] = {
    "not_estimable": "not est",
    "not_higher": "est and not s",
    "direction_not_robust": "est and s and not robust",
    "positive_between_spot_agrees": "est and s and robust and d==agree",
    "positive_native_only": "est and s and robust and d!=agree",
}


def enumerate_inputs() -> Iterator[Inputs]:
    for est, s, c, f in itertools.product((False, True), repeat=4):
        for d in D_VALUES:
            yield Inputs(est, s, c, d, f)


def evaluate(i: Inputs, conditions: Optional[Dict[str, Callable]] = None) -> List[str]:
    """All branch ids whose condition is true (independent evaluation)."""
    conds = CONDITIONS if conditions is None else conditions
    return [k for k, fn in conds.items() if fn(i)]


def decide(i: Inputs, conditions: Optional[Dict[str, Callable]] = None) -> dict:
    """The unique branch (plus the recorded reason for the robustness branch). Raises unless exactly one holds."""
    hits = evaluate(i, conditions)
    if len(hits) != 1:
        raise AssertionError(f"{len(hits)} branches true for {i}: {hits}")
    branch = hits[0]
    reason = ""
    if branch == "direction_not_robust":
        reason = "segmentation-dependent" if (i.c and not i.f) else "not robust to normalisation"
    return {"branch": branch, "reason": reason}


def check_exhaustive(conditions: Optional[Dict[str, Callable]] = None) -> dict:
    """Run all 48 combinations; return counts of combinations with 0, 1 and >1 true conditions."""
    n = zero = multi = 0
    for i in enumerate_inputs():
        k = len(evaluate(i, conditions))
        n += 1
        zero += k == 0
        multi += k > 1
    return {"n": n, "zero": zero, "multi": multi, "exactly_one": n - zero - multi}
