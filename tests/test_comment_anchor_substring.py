"""Unit tests for substring-level comment anchoring.

migrate_comments_to_v33.inject_comment_anchor wraps the WHOLE paragraph, so
several comments on one paragraph get identical ranges -> invalid package.
This injector wraps only the named clause.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
NSMAP = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def make_para(*run_texts):
    p = etree.SubElement(etree.Element(W + "body", nsmap=NSMAP), W + "p")
    for txt in run_texts:
        r = etree.SubElement(p, W + "r")
        t = etree.SubElement(r, W + "t")
        t.text = txt
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    return p


def para_text(p):
    return "".join(t.text or "" for t in p.iter(W + "t"))


def anchored_text(p, cid):
    """Text lying between commentRangeStart and commentRangeEnd for cid."""
    out, inside = [], False
    for el in p.iter():
        if el.tag == W + "commentRangeStart" and el.get(W + "id") == cid:
            inside = True
        elif el.tag == W + "commentRangeEnd" and el.get(W + "id") == cid:
            inside = False
        elif inside and el.tag == W + "t":
            out.append(el.text or "")
    return "".join(out)


def test_anchors_only_the_named_substring():
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("Hello world. ", "The quick brown fox. ", "Goodbye.")
    inject_comment_anchor_substring(p, "7", "The quick brown fox.")
    assert anchored_text(p, "7") == "The quick brown fox."


def test_paragraph_text_is_preserved_exactly():
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("alpha beta gamma delta")
    before = para_text(p)
    inject_comment_anchor_substring(p, "1", "beta gamma")
    assert para_text(p) == before


def test_substring_spanning_two_runs():
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("one two thr", "ee four five")
    inject_comment_anchor_substring(p, "2", "two three four")
    assert anchored_text(p, "2") == "two three four"
    assert para_text(p) == "one two three four five"


def test_two_disjoint_comments_do_not_overlap():
    """The whole point: two threads on one paragraph get separate ranges."""
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("First clause here. Second clause there.")
    inject_comment_anchor_substring(p, "10", "First clause here.")
    inject_comment_anchor_substring(p, "20", "Second clause there.")
    assert anchored_text(p, "10") == "First clause here."
    assert anchored_text(p, "20") == "Second clause there."


def test_missing_needle_raises():
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("nothing to see")
    with pytest.raises(ValueError, match="not found"):
        inject_comment_anchor_substring(p, "3", "absent clause")


def test_ambiguous_needle_raises():
    """A silent wrong-occurrence anchor is worse than a loud failure."""
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("repeat me and repeat me")
    with pytest.raises(ValueError, match="ambiguous|occurs"):
        inject_comment_anchor_substring(p, "4", "repeat me")


def test_comment_reference_run_is_emitted():
    from comment_anchor_substring import inject_comment_anchor_substring

    p = make_para("anchor this clause please")
    inject_comment_anchor_substring(p, "9", "this clause")
    refs = [e for e in p.iter(W + "commentReference") if e.get(W + "id") == "9"]
    assert len(refs) == 1
