"""A12.2: independent Hungarian allocation gate (gate 1, spec 9.3), synthetic data + mutation checks."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking" / "midkine_remediation" / "src"))
import hungarian_gate as G  # noqa: E402

pytestmark = pytest.mark.unit

TYPES = ["T0", "T1", "T2"]


def _n(prop, n):
    return pd.Series(n, index=prop.index)


def _labels_from_counts(counts: pd.DataFrame) -> pd.DataFrame:
    rows = [(spot, t) for spot, r in counts.iterrows() for t in counts.columns for _ in range(int(r[t]))]
    return pd.DataFrame(rows, columns=["spot_barcode", "cell_type"])


# ---- literal (hand-computed) cases -------------------------------------------------
def test_largest_remainder_literals():
    assert G.largest_remainder([5.0, 3.0, 2.0], 10).tolist() == [5, 3, 2]
    assert G.largest_remainder([0.5, 0.5], 3).tolist() == [2, 1]  # exact tie -> lower index
    assert G.largest_remainder([0.0, 0.0, 0.0], 4).tolist() == [2, 1, 1]  # all-zero -> uniform, remainder to first
    assert G.largest_remainder([1.0, 2.0], 0).tolist() == [0, 0]


def test_two_stage_literal_differs_from_single_stage():
    """p=[.67,.16,.17], nuclei_count_raw=11, 8 mapped cells.
    Stage 1: [7,2,2]. Stage 2 re-round to 8: [5.09,1.45,1.45] -> floors [5,1,1], tie on remainder -> type1 -> [5,2,1].
    A single-stage 8*p rounding would give [6,1,1] and must NOT be the expectation."""
    prop = pd.DataFrame([[0.67, 0.16, 0.17]], index=["s"], columns=TYPES)
    exp = G.expected_counts(prop, pd.Series({"s": 11}), pd.Series({"s": 8}))
    assert exp.loc["s"].tolist() == [5, 2, 1]
    assert G.largest_remainder(8 * prop.loc["s"].values, 8).tolist() == [6, 1, 1]
    ok = G.verify_allocation(prop, pd.Series({"s": 11}), _n(prop, 8), _labels_from_counts(pd.DataFrame([[5, 2, 1]], index=["s"], columns=TYPES)))
    assert ok.passed
    single_stage = G.verify_allocation(prop, pd.Series({"s": 11}), _n(prop, 8), _labels_from_counts(pd.DataFrame([[6, 1, 1]], index=["s"], columns=TYPES)))
    assert not single_stage.passed and single_stage.n_mismatch_spots == 1


def test_zero_nuclei_but_mapped_cells_uses_uniform_rule():
    # raw nuclei 0 -> stage-1 zeros -> production re-round of zeros to 4 mapped cells is uniform [2,1,1]
    prop = pd.DataFrame([[0.9, 0.05, 0.05]], index=["s"], columns=TYPES)
    exp = G.expected_counts(prop, pd.Series({"s": 0}), pd.Series({"s": 4}))
    assert exp.loc["s"].tolist() == [2, 1, 1]


def test_spot_without_mapped_cells_expects_zero():
    prop = pd.DataFrame([[0.5, 0.5, 0.0], [0.2, 0.3, 0.5]], index=["a", "b"], columns=TYPES)
    exp = G.expected_counts(prop, pd.Series({"a": 4, "b": 5}), pd.Series({"b": 5}))
    assert exp.loc["a"].tolist() == [0, 0, 0] and exp.loc["b"].tolist() == [1, 2, 2]


# ---- randomised cross-check against the production code ------------------------------
def _production():
    try:
        from CITEgeist.model.assignment import cell_assignment as CA
        from CITEgeist.model.assignment.cellularity_utils import round_counts_largest_remainder as prod_lr
    except Exception:  # cuOPT / GPU import problems on CPU nodes must not error collection
        pytest.skip("CITEgeist assignment modules not importable here")
    return CA, prod_lr


def test_largest_remainder_matches_production_randomised():
    _, prod_lr = _production()
    rng = np.random.default_rng(1)
    for _ in range(3000):
        t = int(rng.integers(2, 8))
        c = rng.dirichlet(np.ones(t)) * rng.integers(1, 40)
        if rng.random() < 0.3:
            c = np.round(c)  # induce ties and zeros
        n = int(rng.integers(0, 40))
        np.testing.assert_array_equal(G.largest_remainder(c, n), prod_lr(c, n))


def _synthetic_world(seed=0, n_spots=60):
    rng = np.random.default_rng(seed)
    idx = [f"S{i:03d}" for i in range(n_spots)]
    prop = pd.DataFrame(rng.dirichlet(np.ones(4), n_spots) * 0.92, index=idx, columns=["A", "B", "C", "D"])  # rows ~0.92, like the QP files
    raw = pd.Series(rng.integers(0, 15, n_spots), index=idx)
    mapped = np.clip(raw.values + rng.integers(-3, 4, n_spots), 0, None)  # nuclei_count_raw != mapped in many spots
    mapped[-1] = raw.iloc[-1] + 3  # guarantee a mismatch in the LAST spot
    return idx, prop, raw, mapped


def _run_production(idx, prop, raw, mapped):
    CA, _ = _production()
    cell_to_spot = np.repeat(np.arange(len(idx)), mapped)
    res = CA.assign_cells(prop, raw, cell_to_spot, assignment_method="hungarian", morphology_weight=0.0)
    return pd.DataFrame({"spot_barcode": res["spot_id"].values, "cell_type": res["assigned_type"].values})


def test_gate_accepts_production_hungarian_output_with_unequal_counts():
    idx, prop, raw, mapped = _synthetic_world()
    assert (raw.values != mapped).sum() >= 10 and raw.iloc[-1] != mapped[-1]
    labels = _run_production(idx, prop, raw, mapped)
    assert len(labels) == mapped.sum()
    res = G.verify_allocation(prop, raw, _n(prop, mapped), labels)
    assert res.passed, res.to_dict()
    assert res.n_cells == int(mapped.sum()) and res.n_spots_checked == len(idx)


def test_mutation_wrong_allocation_in_last_spot_fails():
    idx, prop, raw, mapped = _synthetic_world()
    labels = _run_production(idx, prop, raw, mapped)
    last = labels["spot_barcode"] == idx[-1]
    assert last.sum() >= 2
    bad = labels.copy()
    pos = np.where(last)[0]
    # move one cell of the last spot from its label to a different type
    cur = bad.loc[bad.index[pos[0]], "cell_type"]
    bad.loc[bad.index[pos[0]], "cell_type"] = "A" if cur != "A" else "B"
    res = G.verify_allocation(prop, raw, _n(prop, mapped), bad)
    assert not res.passed and res.n_mismatch_spots == 1 and res.mismatches[0]["spot"] == idx[-1]


def test_mutation_drop_last_type_label_swap_and_scramble_fail():
    idx, prop, raw, mapped = _synthetic_world()
    labels = _run_production(idx, prop, raw, mapped)
    # swap two type names everywhere
    swapped = labels.replace({"A": "B", "B": "A"})
    assert not G.verify_allocation(prop, raw, _n(prop, mapped), swapped).passed
    # scramble spot assignment
    scr = labels.copy()
    scr["spot_barcode"] = np.roll(scr["spot_barcode"].values, 7)
    assert not G.verify_allocation(prop, raw, _n(prop, mapped), scr).passed
    # drop the last type column from the proportions: label set no longer inside type set
    assert not G.verify_allocation(prop.iloc[:, :-1], raw, _n(prop, mapped), labels).passed
    # delete the LAST cell
    assert not G.verify_allocation(prop, raw, _n(prop, mapped), labels.iloc[:-1]).passed


def test_empty_inputs_fail_closed():
    idx, prop, raw, mapped = _synthetic_world()
    labels = _run_production(idx, prop, raw, mapped)
    assert not G.verify_allocation(prop, raw, _n(prop, mapped), labels.iloc[0:0]).passed
    assert not G.verify_allocation(prop.iloc[0:0], raw, _n(prop, mapped), labels).passed
    # broken matcher: labels all reference unknown spots -> failure, not a vacuous pass
    ghost = labels.copy()
    ghost["spot_barcode"] = "nope"
    assert not G.verify_allocation(prop, raw, _n(prop, mapped), ghost).passed


def test_missing_nuclei_count_fails_closed():
    idx, prop, raw, mapped = _synthetic_world()
    labels = _run_production(idx, prop, raw, mapped)
    res = G.verify_allocation(prop, raw.drop(idx[-1]), _n(prop, mapped), labels)
    assert not res.passed and "nuclei_count_raw" in " ".join(res.reasons)
