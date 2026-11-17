"""Tests for gt_rebuild/record_inputs.py."""

import hashlib
import os
import sys
from datetime import datetime
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import record_inputs as ri  # noqa: E402

pytestmark = pytest.mark.unit


def _file(tmp_path, name, data=b"abc", mtime="2026-10-01T00:00:00"):
    p = tmp_path / name
    p.write_bytes(data)
    t = datetime.fromisoformat(mtime).timestamp()
    os.utime(p, (t, t))
    return p


def test_inventory_hashes_match_sha256sum(tmp_path):
    p = _file(tmp_path, "a.txt", b"hello")
    (e,) = ri.inventory([p], base=tmp_path)
    assert e["sha256"] == hashlib.sha256(b"hello").hexdigest()
    assert e["relpath"] == "a.txt" and e["size"] == 5
    assert e["mtime_iso"] == "2026-10-01T00:00:00"


def test_assert_predates_passes_when_older(tmp_path):
    es = ri.inventory([_file(tmp_path, "a")], base=tmp_path)
    ri.assert_predates(es, "2026-10-24T14:43:33")


def test_assert_predates_fails_when_one_file_newer(tmp_path):
    es = ri.inventory([_file(tmp_path, "a"), _file(tmp_path, "b", mtime="2026-10-25T00:00:00")], base=tmp_path)
    with pytest.raises(AssertionError, match="b"):
        ri.assert_predates(es, "2026-10-24T14:43:33")


def test_empty_inventory_raises():
    with pytest.raises(ValueError):
        ri.inventory([])
    with pytest.raises(ValueError):
        ri.assert_predates([], "2026-10-24T14:43:33")


def test_earliest_gt_built_at_takes_min():
    m = {
        "gt_built": {"gt_built_at": "2026-10-24T15:42:21"},
        "superseded_A": {"previous_sections": {"gt_built": {"gt_built_at": "2026-10-24T14:43:33"}}},
        "superseded_B": {"previous_sections": {"gt_built": {"gt_built_at": "2026-10-24T16:00:00"}}},
        "other": {},
    }
    assert ri.earliest_gt_built_at(m) == "2026-10-24T14:43:33"
    with pytest.raises(ValueError):
        ri.earliest_gt_built_at({})


def test_zero_match_group_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(ri, "XB", tmp_path)
    monkeypatch.setattr(ri, "group_globs", lambda: {"g": ["nothing/*.csv"]})
    with pytest.raises(FileNotFoundError, match="g"):
        ri.resolve_groups()


def test_missing_reference_dir_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(ri, "XB", tmp_path)
    with pytest.raises(FileNotFoundError, match="reference"):
        ri.group_globs()
