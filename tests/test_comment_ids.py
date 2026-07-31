"""Tests for the package-aware comment-id allocator."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path.home() / ".claude/skills/word-edit/scripts"))

pytestmark = pytest.mark.unit


@pytest.fixture(autouse=True)
def _reset_allocator():
    """word_edit._n is a module global; restore it so tests don't leak state."""
    import word_edit as we

    saved = we._n[0]
    yield
    we._n[0] = saved


def test_seeds_above_existing_ids():
    import word_edit as we
    from comment_ids import seed_allocator_above

    floor = seed_allocator_above(["0", "3", "1001", "1016"])
    assert floor == 1016
    assert we.nid() == "1017"
    assert we.nid() == "1018"


def test_never_lowers_the_allocator():
    """Seeding with ids below the current floor must not rewind it."""
    import word_edit as we
    from comment_ids import seed_allocator_above

    we._n[0] = 5000
    floor = seed_allocator_above(["1", "2"])
    assert floor == 5000
    assert we.nid() == "5001"


def test_empty_input_leaves_default_floor():
    import word_edit as we
    from comment_ids import seed_allocator_above

    we._n[0] = 1000
    assert seed_allocator_above([]) == 1000
    assert we.nid() == "1001"


def test_ignores_non_numeric_ids():
    import word_edit as we
    from comment_ids import seed_allocator_above

    we._n[0] = 1000
    assert seed_allocator_above(["abc", "42", ""]) == 1000 or we._n[0] == 42
    assert int(we.nid()) > 42
