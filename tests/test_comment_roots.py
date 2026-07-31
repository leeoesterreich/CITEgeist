"""Tests for add_root_comment."""

import sys
from pathlib import Path

import pytest
from lxml import etree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path.home() / ".claude/skills/shared"))
sys.path.insert(0, str(Path.home() / ".claude/skills/word-edit/scripts"))

pytestmark = pytest.mark.unit

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
W15 = "{http://schemas.microsoft.com/office/word/2012/wordml}"
W16CID = "{http://schemas.microsoft.com/office/word/2016/wordml/cid}"
W16CEX = "{http://schemas.microsoft.com/office/word/2018/wordml/cex}"

NSMAP = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "w14": "http://schemas.microsoft.com/office/word/2010/wordml",
    "w15": "http://schemas.microsoft.com/office/word/2012/wordml",
    "w16cid": "http://schemas.microsoft.com/office/word/2016/wordml/cid",
    "w16cex": "http://schemas.microsoft.com/office/word/2018/wordml/cex",
}


@pytest.fixture(autouse=True)
def _reset_allocator():
    import word_edit as we

    saved = we._n[0]
    yield
    we._n[0] = saved


@pytest.fixture
def parts():
    """Minimal but structurally valid empty comment parts plus a one-paragraph body."""
    com = ET.Element(W + "comments", nsmap=NSMAP)
    ext = ET.Element(W15 + "commentsEx", nsmap=NSMAP)
    cids = ET.Element(W16CID + "commentsIds", nsmap=NSMAP)
    cex = ET.Element(W16CEX + "commentsExtensible", nsmap=NSMAP)
    people = ET.Element(W15 + "people", nsmap=NSMAP)
    body = ET.Element(W + "body", nsmap=NSMAP)
    p = ET.SubElement(body, W + "p")
    r = ET.SubElement(p, W + "r")
    t = ET.SubElement(r, W + "t")
    t.text = "The quick brown fox jumps over the lazy dog."
    return dict(com=com, ext=ext, cids=cids, cex=cex, people=people, body=body, para=p)


def _call(parts, **over):
    from comment_roots import add_root_comment

    kwargs = dict(
        com=parts["com"],
        ext=parts["ext"],
        cids=parts["cids"],
        cex=parts["cex"],
        people=parts["people"],
        body=parts["body"],
        author="Russell Schwartz",
        initials="RS",
        date="2026-07-24T18:35:00Z",
        text="I might temper this a bit.",
        anchor_para=parts["para"],
        anchor_needle="quick brown fox",
    )
    kwargs.update(over)
    return add_root_comment(**kwargs)


def test_writes_all_five_parts(parts):
    cid = _call(parts)

    comments = parts["com"].findall(W + "comment")
    assert len(comments) == 1
    assert comments[0].get(W + "id") == cid
    assert comments[0].get(W + "author") == "Russell Schwartz"
    assert comments[0].get(W + "date") == "2026-07-24T18:35:00Z"
    assert "I might temper this a bit." in "".join(comments[0].itertext())

    assert len(parts["ext"].findall(W15 + "commentEx")) == 1
    assert len(parts["cids"].findall(W16CID + "commentId")) == 1
    assert len(parts["cex"].findall(W16CEX + "commentExtensible")) == 1
    assert len(parts["people"].findall(W15 + "person")) == 1


def test_root_has_no_parent_link(parts):
    _call(parts)
    ce = parts["ext"].find(W15 + "commentEx")
    assert ce.get(W15 + "paraIdParent") is None


def test_durable_id_is_shared_between_cids_and_extensible(parts):
    _call(parts)
    did = parts["cids"].find(W16CID + "commentId").get(W16CID + "durableId")
    assert parts["cex"].find(W16CEX + "commentExtensible").get(W16CEX + "durableId") == did


def test_anchor_is_narrow_not_whole_paragraph(parts):
    cid = _call(parts)
    p = parts["para"]
    starts = [e for e in p.iter(W + "commentRangeStart") if e.get(W + "id") == cid]
    ends = [e for e in p.iter(W + "commentRangeEnd") if e.get(W + "id") == cid]
    refs = [e for e in p.iter(W + "commentReference") if e.get(W + "id") == cid]
    assert len(starts) == 1 and len(ends) == 1 and len(refs) == 1
    # Text between the markers must be the needle, not the whole sentence.
    covered = []
    seen_start = False
    for node in p.iter():
        if node.tag == W + "commentRangeStart" and node.get(W + "id") == cid:
            seen_start = True
        elif node.tag == W + "commentRangeEnd" and node.get(W + "id") == cid:
            break
        elif seen_start and node.tag == W + "t":
            covered.append(node.text or "")
    assert "".join(covered) == "quick brown fox"


def test_anchor_markers_are_direct_children_of_the_paragraph(parts):
    """add_threaded_reply searches only direct children; a nested ref breaks replies."""
    cid = _call(parts)
    p = parts["para"]
    assert any(e.get(W + "id") == cid for e in p.findall(W + "commentRangeStart"))
    assert any(e.get(W + "id") == cid for e in p.findall(W + "commentRangeEnd"))
    assert any(
        r.find(W + "commentReference") is not None and r.find(W + "commentReference").get(W + "id") == cid
        for r in p.findall(W + "r")
    )


def test_ids_do_not_collide_with_replies(parts):
    import word_edit as we

    first = _call(parts)
    second = _call(parts, text="second thread", anchor_needle="lazy dog")
    assert first != second
    assert int(second) > int(first)
    # A reply minted afterwards must also be distinct.
    reply = we.add_threaded_reply(
        second,
        "our answer",
        com=parts["com"],
        ext=parts["ext"],
        cids=parts["cids"],
        cex=parts["cex"],
        body=parts["body"],
    )
    assert len({first, second, reply}) == 3


def test_author_registration_is_idempotent(parts):
    _call(parts)
    _call(parts, text="second", anchor_needle="lazy dog")
    assert len(parts["people"].findall(W15 + "person")) == 1


def test_missing_needle_raises(parts):
    with pytest.raises(ValueError):
        _call(parts, anchor_needle="not present in this paragraph")
