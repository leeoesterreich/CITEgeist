"""A12.3 G9/G10/G11: hash loader, typed statuses, output namespace, v3 check, R00 distance assertion."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import protocol_io as P

pytestmark = pytest.mark.unit

DRAFT = Path(__file__).resolve().parents[2] / "Benchmarking" / "midkine_remediation" / "protocol.DRAFT.json"


def _frozen(tmp_path, raw=None):
    raw = raw if raw is not None else DRAFT.read_bytes()
    out = tmp_path / "protocol.json"
    out.write_bytes(P.freeze_text(raw))
    return out


def test_draft_refused_unless_allowed():
    with pytest.raises(P.ProtocolHashError):
        P.load_protocol(DRAFT)
    assert P.load_protocol(DRAFT, allow_draft=True)["status"] == "DRAFT"


def test_frozen_roundtrip_and_hash_is_of_null_field(tmp_path):
    f = _frozen(tmp_path)
    proto = P.load_protocol(f)
    assert proto["status"] == "FROZEN"
    raw = f.read_bytes()
    # independent recomputation: null the field by plain string replacement
    nulled = raw.replace(f'"frozen_sha256": "{proto["frozen_sha256"]}"'.encode(), b'"frozen_sha256": null')
    assert hashlib.sha256(nulled).hexdigest() == proto["frozen_sha256"] == proto["_protocol_sha256"]


def test_one_byte_edit_anywhere_aborts(tmp_path):
    f = _frozen(tmp_path)
    raw = bytearray(f.read_bytes())
    # edit one byte inside the body (a digit in the draft date), then inside the recorded hash itself
    i = raw.index(b"2026-10-14")
    for pos in (i, raw.index(b'"frozen_sha256": "') + len('"frozen_sha256": "') + 3):
        mut = bytearray(raw)
        mut[pos] = ord("0") if mut[pos] != ord("0") else ord("1")
        g = tmp_path / "mut.json"
        g.write_bytes(bytes(mut))
        with pytest.raises(P.ProtocolHashError):
            P.load_protocol(g)


def test_hash_matcher_must_find_field():
    with pytest.raises(P.ProtocolHashError):
        P.protocol_hash(b'{"status": "FROZEN"}')


def test_find_to_fill_values_only():
    assert P.find_to_fill({"a": "x TO_FILL_PRE_FREEZE y", "b": [{"c": "TO_FILL_PRE_FREEZE"}], "TO_FILL_PRE_FREEZE": "k"}) == ["/a", "/b[0]/c"]
    assert P.find_to_fill({"a": "clean"}) == []


def test_typed_status_contract():
    P.validate_typed({"status": "ok", "value": 0.3})
    P.validate_typed({"status": "not_estimable", "reason": "n<30"})
    for bad in ({}, {"status": "ok"}, {"status": "ok", "value": float("nan")}, {"status": "not_estimable"}, {"status": "weird", "reason": "x"}):
        with pytest.raises((P.GateFailure, ValueError)):
            P.validate_typed(bad)


def test_output_target_refusals(tmp_path):
    ok = tmp_path / "abc123def456" / "run1" / "P1"
    assert P.check_output_target(ok, "abc123def456", "run1") == ok
    with pytest.raises(P.OutputTargetError):
        P.check_output_target(tmp_path / "morphology_assignment_v3" / "abc123def456" / "run1", "abc123def456", "run1")
    with pytest.raises(P.OutputTargetError):
        P.check_output_target(tmp_path / "other" / "run1", "abc123def456", "run1")
    ok.mkdir(parents=True)
    P.check_output_target(ok, "abc123def456", "run1")  # empty existing dir is fine
    (ok / "x.txt").write_text("1")
    with pytest.raises(P.OutputTargetError):
        P.check_output_target(ok, "abc123def456", "run1")  # non-empty is refused


def test_atomic_write_leaves_no_temp(tmp_path):
    P.atomic_write_json(tmp_path / "d" / "o.json", {"a": 1})
    assert json.loads((tmp_path / "d" / "o.json").read_text()) == {"a": 1}
    assert [p.name for p in (tmp_path / "d").iterdir()] == ["o.json"]
    with pytest.raises(TypeError):
        P.atomic_write_json(tmp_path / "d" / "bad.json", {"a": object()})
    assert [p.name for p in (tmp_path / "d").iterdir()] == ["o.json"]


def test_v3_hash_check_detects_change_and_empty(tmp_path):
    f = tmp_path / "a.h5ad"
    f.write_bytes(b"abc")
    man = tmp_path / "m.sha256"
    man.write_text(f"{hashlib.sha256(b'abc').hexdigest()}  a.h5ad\n")
    assert P.v3_hashes_ok(man, tmp_path) == []
    f.write_bytes(b"abd")
    assert P.v3_hashes_ok(man, tmp_path) == [str(f)]
    man.write_text("")
    with pytest.raises(P.GateFailure):
        P.v3_hashes_ok(man, tmp_path)


def test_real_v3_manifest_has_12_entries():
    man = DRAFT.parent / "v3_manifest.sha256"
    if not man.exists():
        pytest.skip("v3 manifest is not part of the public tree")
    assert len([ln for ln in man.read_text().splitlines() if ln.strip()]) == 12


class _Fake:
    def __init__(self, keys):
        self.obsp = {k: None for k in keys}


def test_r00_has_no_spatial_distance():
    P.assert_no_spatial_distance(_Fake([]), "R00")
    P.assert_no_spatial_distance(_Fake(["spatial_distance"]), "R10")
    with pytest.raises(P.GateFailure):
        P.assert_no_spatial_distance(_Fake(["spatial_distance"]), "R00")
    with pytest.raises(P.GateFailure):
        P.assert_no_spatial_distance(_Fake([]), "R10")


def test_unfilled_ignores_only_the_two_meta_fields():
    doc = {"freeze_rule": "needs zero TO_FILL_PRE_FREEZE", "pre_freeze_checklist": ["no TO_FILL_PRE_FREEZE left"], "a": {"b": "TO_FILL_PRE_FREEZE x"}}
    assert P.unfilled(doc) == ["/a/b"]
    assert P.unfilled({**doc, "a": {"b": "filled"}}) == []


def test_draft_protocol_has_no_unfilled_placeholder():
    proto = json.loads(DRAFT.read_text())
    assert P.unfilled(proto) == [], P.unfilled(proto)


def test_g9_outputs_carry_the_verified_protocol_hash(tmp_path):
    proto = P.load_protocol(_frozen(tmp_path))
    out = P.with_protocol_hash({"x": 1}, proto)
    assert out["protocol_sha256"] == proto["frozen_sha256"] and out["protocol_sha12"] == proto["frozen_sha256"][:12] and out["x"] == 1
    draft = P.load_protocol(DRAFT, allow_draft=True)
    with pytest.raises(P.ProtocolHashError):
        P.with_protocol_hash({"x": 1}, draft)             # a draft can never be stamped into outputs
