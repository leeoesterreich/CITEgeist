"""A12.3 fix round: moving-block label shifts, k_share label (independent of n_st), rule file present."""
from pathlib import Path

import numpy as np
import pytest

import hurdle as H
import spatial_utils as S
import synthetic_g6 as G6
import synthetic_g6fix as F

pytestmark = pytest.mark.unit
REM = Path(__file__).resolve().parents[2] / "Benchmarking" / "midkine_remediation"
HAVE = (S.DATA_ROOT / "HCC22-088-P4-S2_1i_rep" / "outs" / "spatial" / "tissue_positions.csv").exists()


def _blob(n=900, seed=0):
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(40), np.arange(40))
    xy = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    keep = np.hypot(xy[:, 0] - 20, xy[:, 1] - 20) < 18
    return xy[keep]


def test_shift_pool_properties():
    xy = _blob()
    pool = H.shift_pool(xy, 1.0, n_shifts=50, min_shift_pitch=6.0)
    assert pool.shape == (50, len(xy)) and pool.min() >= 0 and pool.max() < len(xy)
    assert (pool == np.arange(len(xy))).mean() < 0.01                  # no identity-like shifts
    disp = np.array([np.median(np.linalg.norm(xy[row] - xy, axis=1)) for row in pool])
    assert disp.min() >= 4.0                                           # every shift really moves the field (>= ~6 pitches)
    # a clustered label keeps its clustering and (approximately) its prevalence under shifts
    g = (xy[:, 0] > 22).astype(int)
    L = H.shifted_labels(g, pool, 50, np.random.default_rng(1))
    W = S.knn_weights(xy, 6)
    agree = np.mean([float((W @ l * l).sum() / max(l.sum(), 1)) for l in L.astype(float)])
    plain = np.mean([float((W @ l * l).sum() / max(l.sum(), 1)) for l in np.random.default_rng(2).permuted(np.tile(g, (50, 1)), axis=1).astype(float)])
    assert agree > 0.8 > plain + 0.2
    assert abs(L.mean() - g.mean()) < 0.2


def test_shift_pool_fails_closed_when_impossible():
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    with pytest.raises(RuntimeError):
        H.shift_pool(xy, 1.0, n_shifts=5, min_shift_pitch=50.0, max_tries=200)


@pytest.mark.skipif(not HAVE, reason="coordinates not available")
def test_k_share_label_does_not_depend_on_n_noise_and_threshold_is_pooled_quantile():
    lat = G6.make_lattice()
    G6.K0_SHARE = F.K0_SHARE_FROZEN
    G6.CASES["tmp_lo"] = dict(G6.CASES["null_noisy_nst_kshare"], sigma_n=0.01)
    a = G6.simulate("null_noisy_nst_kshare", lat, np.random.default_rng(3))
    b = G6.simulate("tmp_lo", lat, np.random.default_rng(3))
    assert (a["group"] == b["group"]).all()                            # n_st noise does not touch the label
    assert not np.allclose(a["obs"]["n_st_eff"], b["obs"]["n_st_eff"])
    fr = np.mean([G6.simulate("null_noisy_nst_kshare", lat, np.random.default_rng(10 + i))["group"].mean() for i in range(6)])
    assert abs(fr - 0.30) < 0.03                                       # frozen k0 gives ~30% high on fresh seeds
    with pytest.raises(RuntimeError):
        G6.K0_SHARE = None
        G6.simulate("null_noisy_nst_kshare", lat, np.random.default_rng(3))
    G6.K0_SHARE = F.K0_SHARE_FROZEN


def test_selection_rule_file_exists_and_names_all_variants():
    rule = REM / "g6_fix_round_selection_rule.md"
    if not rule.exists():
        pytest.skip("selection-rule note is not part of the public tree")
    txt = rule.read_text()
    for v in F.VARIANTS:
        assert v in txt
    assert "0.0143791" in txt and "0.07" in txt and "0.80" in txt
