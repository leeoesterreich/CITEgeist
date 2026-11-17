"""A12.3 hurdle inference: statistics checked against independent algebra; gates; typed statuses."""
import numpy as np
import pytest
import statsmodels.api as sm

import hurdle as H
import hurdle_g6_historical as HG

pytestmark = pytest.mark.unit


def _grid(nx=30, ny=30):
    xs, ys = np.meshgrid(np.arange(nx), np.arange(ny))
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(float)


def test_blocks_tiling_merge_and_median_pitch():
    xy = _grid(24, 12)                       # 288 points, pitch 1
    b = H.make_blocks(xy, pitch=1.0, side=6.0, min_n=8)
    assert len(np.unique(b)) == 8            # 4 x 2 tiles of 6x6
    assert np.bincount(b).tolist() == [36] * 8
    assert H.median_pitch(xy) == 1.0
    # a 2-point island merges into its nearest tile
    xy2 = np.vstack([xy, [[40.0, 0.0], [40.0, 1.0]]])
    b2 = H.make_blocks(xy2, pitch=1.0, side=6.0, min_n=8)
    assert np.bincount(b2).min() >= 8 and len(b2) == 290
    assert b2[-1] == b2[-2] and b2[-1] == b2[np.argmin(np.linalg.norm(xy - [40, 0], axis=1))] 


def test_benjamini_hochberg_literal():
    q = HG.historical_benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
    # sorted p: .01 .03 .04 .2 -> p*m/rank = .04 .06 .0533 .2 -> monotone from the top: .04 .0533 .0533 .2
    assert q == pytest.approx([0.04, 0.05333333, 0.05333333, 0.2])


def test_estimability_gate_branches():
    g = np.tile([0, 1], 150)                                  # 300 units
    blk = np.repeat(np.arange(30), 10)                        # 30 blocks of 10
    assert H.estimability(g, blk, "participation")["ok"]
    assert "compartments per group" in H.estimability(g[:40], blk[:40], "participation")["reason"]
    # one block only -> 0 or 1 mixed blocks
    r = H.estimability(g, np.zeros(300, int), "participation")
    assert not r["ok"] and "mixed-label blocks" in r["reason"]
    # 9 mixed blocks (rest pure) is below the gate, 10 passes
    lab = np.zeros(300, int)
    blk2 = np.repeat(np.arange(30), 10)
    lab[:] = (blk2 < 9).astype(int) * (np.arange(300) % 2) + (blk2 >= 9) * (blk2 % 2)
    assert not H.estimability(lab, blk2, "participation")["ok"]
    lab[blk2 == 9] = np.tile([0, 1], 5)
    assert H.estimability(lab, blk2, "participation")["ok"]
    # distinct-permutation arithmetic: 10 mixed size-2 blocks -> 2^10 = 1024 >= 999
    _, logd = HG._mixed_and_distinct(np.repeat(np.arange(10), 2), np.tile([0, 1], 10), "participation")
    assert np.exp(logd) == pytest.approx(1024.0)
    _, logd9 = HG._mixed_and_distinct(np.repeat(np.arange(9), 2), np.tile([0, 1], 9), "participation")
    assert np.exp(logd9) == pytest.approx(512.0)
    # magnitude: too few active senders in a group
    sc = np.zeros(300)
    sc[:12] = 1.0
    r = H.estimability(g, blk, "magnitude", sc > 0)
    assert not r["ok"] and "active senders" in r["reason"]


def test_small_cases_return_typed_not_estimable():
    rng = np.random.default_rng(1)
    n = 200
    g = rng.integers(0, 2, n)
    blk = np.repeat(np.arange(4), 50)                         # 4 blocks only
    sc = rng.lognormal(size=n) * (rng.random(n) < 0.6)
    out = HG.historical_h1_permutation_test(sc, g, np.zeros((n, 0)), blk, B=99)
    for part in ("magnitude", "participation"):
        assert out[part]["status"] == "not_estimable" and "mixed-label blocks" in out[part]["reason"]


def test_magnitude_t_matches_ols_with_block_fixed_effects():
    rng = np.random.default_rng(2)
    n = 400
    g = rng.integers(0, 2, n)
    blk = np.repeat(np.arange(20), 20)
    Z = rng.normal(size=(n, 2))
    y = 0.5 * g + 0.3 * Z[:, 0] + rng.normal(size=20)[blk] + rng.normal(size=n)
    score = np.exp(y)
    res = HG.historical_magnitude_permutation_test(score, g, Z, blk, B=199)
    D = np.eye(20)[blk][:, 1:]
    ols = sm.OLS(y, np.column_stack([np.ones(n), g, Z, D])).fit()
    assert res["coef"] == pytest.approx(ols.params[1], abs=1e-10)
    assert res["t"] == pytest.approx(ols.tvalues[1], abs=1e-9)
    assert res["gmr"] == pytest.approx(np.exp(ols.params[1]))
    assert res["p"] < 0.01 and res["sign_positive"]


def test_fl_p_value_agrees_with_classical_when_unclustered():
    rng = np.random.default_rng(3)
    n = 300
    g = rng.integers(0, 2, n)
    blk = np.repeat(np.arange(15), 20)
    y = 0.15 * g + rng.normal(size=n)
    res = HG.historical_magnitude_permutation_test(np.exp(y), g, np.zeros((n, 0)), blk, B=4999)
    D = np.eye(15)[blk][:, 1:]
    ols = sm.OLS(y, np.column_stack([np.ones(n), g, D])).fit()
    assert res["p"] == pytest.approx(ols.pvalues[1], abs=0.04)


def test_participation_score_matches_full_information_matrix_algebra():
    rng = np.random.default_rng(4)
    n = 300
    g = rng.integers(0, 2, n)
    Z = rng.normal(size=(n, 3))
    a = (rng.random(n) < 1 / (1 + np.exp(-(0.2 + 0.5 * Z[:, 0])))).astype(float)
    blk = np.repeat(np.arange(15), 20)
    res = HG.historical_participation_permutation_test(a, g, Z, blk, B=199)
    # independent route: score vector / information of the FULL model evaluated at the null fit
    X0 = np.column_stack([np.ones(n), Z])
    null = sm.GLM(a, X0, family=sm.families.Binomial()).fit()
    mu = null.fittedvalues
    Xf = np.column_stack([X0, g])
    s = Xf.T @ (a - mu)
    I = Xf.T @ ((mu * (1 - mu))[:, None] * Xf)
    stat = float(s @ np.linalg.solve(I, s))
    assert res["score"] == pytest.approx(stat, rel=1e-6)
    full = sm.GLM(a, Xf, family=sm.families.Binomial()).fit()
    assert res["odds_ratio"] == pytest.approx(np.exp(full.params[-1]), rel=1e-5)


def test_participation_detects_effect_and_not_null():
    rng = np.random.default_rng(5)
    n = 600
    blk = np.repeat(np.arange(30), 20)
    g = (rng.random(n) < 0.4).astype(int)
    Z = rng.normal(size=(n, 1))
    a_eff = (rng.random(n) < 1 / (1 + np.exp(-(-0.5 + 1.5 * g)))).astype(float)
    a_null = (rng.random(n) < 0.5).astype(float)
    assert HG.historical_participation_permutation_test(a_eff, g, Z, blk, B=999)["p"] < 0.01
    assert HG.historical_participation_permutation_test(a_null, g, Z, blk, B=999)["p"] > 0.01


def test_few_distinct_permutations_are_sampled_without_replacement():
    n = 20
    blk = np.repeat(np.arange(10), 2)
    g = np.tile([0, 1], 10)
    rng = np.random.default_rng(6)
    a = (rng.random(n) < 0.5).astype(float)
    # participation needs >= MIN_PER_GROUP, so call the permutation helper directly
    P = HG._perm_matrix(blk, 4999, np.random.default_rng(0), "block", np.log(1024.0))
    assert len(P) <= 1023 and len({r.tobytes() for r in P}) == len(P)
    assert all((blk[r] == blk).all() for r in P)               # within-block only


def test_block_perms_stay_in_block_and_plain_do_not():
    blk = np.repeat(np.arange(6), 7)
    P = HG.within_block_perms(blk, 50, np.random.default_rng(0))
    assert all((blk[r] == blk).all() for r in P)
    Q = HG.plain_perms(len(blk), 50, np.random.default_rng(0))
    assert not all((blk[r] == blk).all() for r in Q)
    assert all(sorted(r) == list(range(len(blk))) for r in P[:5])


def test_swap_nonzero_filter_mutation_changes_result():
    rng = np.random.default_rng(7)
    n = 600
    blk = np.repeat(np.arange(30), 20)
    g = (rng.random(n) < 0.5).astype(int)
    act = rng.random(n) < 0.7
    score = np.where(act, np.exp(1.0 * g + rng.normal(size=n)), 0.0)
    good = HG.historical_magnitude_permutation_test(score, g, np.zeros((n, 0)), blk, B=199)
    assert good["status"] == "ok" and good["p"] < 0.05
    # mutation (swap_nonzero_filter): the inverted filter selects zero scores -> the guard refuses (the gate can fail)
    with pytest.raises(ValueError):
        HG.historical_magnitude_permutation_test(score, g, np.zeros((n, 0)), blk, B=199, active=~act)
    # mutation: hand the zero-inflated vector to the magnitude model without filtering -> refused too
    with pytest.raises(ValueError):
        HG.historical_magnitude_permutation_test(score, g, np.zeros((n, 0)), blk, B=199, active=np.ones(n, bool))


def _obs4(**over):
    import pandas as pd

    d = {"spot_umi": [100.0, 200, 400, 800], "rna_per_cell_equiv": [1.0, 1, 2, 4], "qp_prop": [.1, .2, .3, .4], "n_st_eff": [1.0, 2, 4, 8],
         "rna_per_cell_equiv_nf": [2.0, 2, 3, 5]}
    d.update(over)
    return pd.DataFrame(d)


def test_covariate_matrix_defaults_to_the_three_primary_covariates():
    Z = H.covariate_matrix(_obs4())
    assert Z.shape == (4, 3) and np.isfinite(Z).all()
    assert H.PRIMARY_COVARIATES == ("log_spot_umi", "log_rna_per_cell_equiv", "qp_prop")
    assert np.allclose(Z.mean(axis=0), 0) and np.allclose(Z.std(axis=0), 1)          # frozen standardisation: population SD (ddof=0)
    # log(n_st_eff) only by explicit opt-in
    assert H.covariate_matrix(_obs4(), include_log_n_st_eff=True).shape == (4, 4)
    assert H.covariate_matrix(_obs4(), drop=("log_spot_umi", "log_rna_per_cell_equiv", "qp_prop")).shape == (4, 0)


def test_covariate_matrix_nuclei_free_swaps_only_the_rna_column():
    a, b = H.covariate_matrix(_obs4()), H.covariate_matrix(_obs4(), nuclei_free=True)
    assert np.allclose(a[:, [0, 2]], b[:, [0, 2]]) and not np.allclose(a[:, 1], b[:, 1])


def test_covariate_domain_is_gated_not_floored():
    with pytest.raises(H.CovariateDomainError):
        H.covariate_matrix(_obs4(rna_per_cell_equiv=[0.0, 1, 2, 4]))
    with pytest.raises(H.CovariateDomainError):
        H.covariate_matrix(_obs4(spot_umi=[100.0, np.nan, 400, 800]))
    with pytest.raises(H.CovariateDomainError):
        H.covariate_matrix(_obs4(qp_prop=[.1, np.inf, .3, .4]))
    # the G6-era floor survives only in the historical module
    assert np.isfinite(HG.historical_covariate_matrix(_obs4(rna_per_cell_equiv=[0.0, 1, 2, 4]))).all()
    assert HG.historical_covariate_matrix(_obs4()).shape == (4, 4)


def test_protocol_covariates_match_primary_covariate_set():
    import json
    from pathlib import Path

    proto = json.loads((Path(__file__).resolve().parents[2] / "Benchmarking/midkine_remediation/protocol.DRAFT.json").read_text())
    assert proto["hurdle"]["models"]["covariates"] == ["log(spot_umi)", "log(rna_per_cell_equiv)", "qp_prop"]
    assert H.covariate_matrix(_obs4()).shape == (4, len(proto["hurdle"]["models"]["covariates"]))


def test_primary_api_has_no_pvalue_functions():
    for name in ("magnitude_test", "participation_test", "h1_test", "benjamini_hochberg", "surrogate_labels", "within_block_perms", "plain_perms", "safe_log"):
        assert not hasattr(H, name), name


def test_distinct_permutation_gate_is_implied_by_the_mixed_block_gate():
    # protocol rev 5 drops the permutation-count requirement: every mixed block has >= 2 units and k in [1, m-1],
    # so it contributes a factor >= 2 (participation: C(m,k) >= 2; magnitude: m! >= 2), and 10 mixed blocks give >= 2^10 = 1024 >= 999.
    rng = np.random.default_rng(0)
    for _ in range(50):
        sizes = rng.integers(2, 6, size=10)
        blk = np.repeat(np.arange(10), sizes)
        g = np.concatenate([rng.permutation(np.r_[1, np.zeros(m - 1, int), rng.integers(0, 2, 0)]) for m in sizes])
        for kind in ("participation", "magnitude"):
            mixed, logd = HG._mixed_and_distinct(blk, g, kind)
            assert mixed == 10 and np.exp(logd) >= 1024 - 1e-6
