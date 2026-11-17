"""Synthetic tests for Benchmarking/midkine_remediation/keratin_secretory_association.py (no real data)."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "midkine_remediation"))
import keratin_secretory_association as KS  # noqa: E402
import morans as M  # noqa: E402
from spatial_utils import gaussian_smoother, knn_weights  # noqa: E402

pytestmark = pytest.mark.unit


def lattice(side=14, seed=0):
    rng = np.random.default_rng(seed)
    g = np.array([(i, j) for i in range(side) for j in range(side)], float)
    return g + rng.uniform(-0.2, 0.2, g.shape)


XY = lattice()
SN = M.SpectralNull(knn_weights(XY, 6))
KERN = gaussian_smoother(XY, 3.0)


def smooth(rng):
    return KERN @ rng.standard_normal(len(XY))


def test_spearman_matches_scipy():
    rng = np.random.default_rng(1)
    x, y = rng.normal(size=81), rng.normal(size=81)
    assert KS.spearman(x, y) == pytest.approx(spearmanr(x, y)[0])
    assert KS.msr_spearman(x, y, M.SpectralNull(knn_weights(lattice(9), 6)), B=99)["rho"] == pytest.approx(spearmanr(x, y)[0])


def test_planted_association_detected():
    rng = np.random.default_rng(2)
    k = smooth(rng) + rng.standard_normal(len(XY))
    sec = 0.8 * k + 0.4 * rng.standard_normal(len(XY))
    r = KS.msr_spearman(np.exp(k), sec, SN, B=999, seed=1)
    assert r["rho"] > 0.5 and r["p_two_sided"] < 0.01
    r_neg = KS.msr_spearman(np.exp(k), -sec, SN, B=999, seed=1)
    assert r_neg["rho"] < -0.5 and r_neg["p_two_sided"] < 0.01


def test_independent_autocorrelated_fields_are_calibrated():
    rng = np.random.default_rng(3)
    n_sim = 300
    p_msr, p_naive = [], []
    for _ in range(n_sim):
        x, y = smooth(rng), smooth(rng)
        p_msr.append(KS.msr_spearman(x, y, SN, B=199, seed=7)["p_two_sided"])
        rho = KS.spearman(x, y)
        perm = np.array([KS.spearman(x, rng.permutation(y)) for _ in range(99)])
        p_naive.append((1 + np.sum(np.abs(perm) >= abs(rho))) / 100)
    p_msr, p_naive = np.array(p_msr), np.array(p_naive)
    assert (p_msr < 0.05).mean() < 0.09                  # nominal 0.05 (binomial sd ~0.013 at 300 sims)
    assert 0.40 < p_msr.mean() < 0.60                    # roughly uniform
    assert (p_naive < 0.05).mean() > 0.20                # the free permutation is anti-conservative here: the MSR test is doing real work


def _inputs(n=60, seed=5, kshare=None, mdk_nonzero=40):
    rng = np.random.default_rng(seed)
    genes = KS.SECRETORY + ["MDK", "ACTB"]
    X = rng.normal(size=(n, len(genes)))
    X[:, genes.index("MDK")] = 0.0
    X[rng.choice(n, min(n, mdk_nonzero), replace=False), genes.index("MDK")] = 1.0
    xy = lattice(10, seed)[:n]
    k = rng.gamma(2.0, 0.01, n) if kshare is None else kshare
    return k, X, genes, xy


def test_estimability_rules():
    k, X, g, xy = _inputs(n=29)
    r = KS.analyze_specimen("s", k, X, g, xy, B=19)
    assert all(r[t]["status"] == "not_estimable" and r[t]["reason_code"] == "eligibility_lt30" for t in KS.TESTS)
    k, X, g, xy = _inputs(n=60, kshare=np.zeros(60))
    r = KS.analyze_specimen("s", k, X, g, xy, B=19)
    assert r["primary"]["reason_code"] == "constant_keratin"
    k, X, g, xy = _inputs(n=60, mdk_nonzero=5)
    r = KS.analyze_specimen("s", k, X, g, xy, covariates=np.random.default_rng(0).normal(size=(60, 3)), B=19)
    assert r["primary"]["status"] == "ok" and r["partial_primary"]["status"] == "ok"
    assert r["secondary"]["reason_code"] == "mdk_nonzero_lt10" and r["partial_secondary"]["reason_code"] == "mdk_nonzero_lt10"
    k, X, g, xy = _inputs(n=60, mdk_nonzero=40)
    r = KS.analyze_specimen("s", k, X, g, xy, covariates=None, cov_error="bad covariate", B=19)
    assert r["primary"]["status"] == r["secondary"]["status"] == "ok"
    assert r["partial_primary"]["reason_code"] == "covariate_domain"


def test_deterministic_with_seed():
    rng = np.random.default_rng(9)
    x, y = smooth(rng), smooth(rng)
    a = KS.msr_spearman(x, y, SN, B=600, seed=11)
    b = KS.msr_spearman(x, y, SN, B=600, seed=11)
    c = KS.msr_spearman(x, y, SN, B=600, seed=12)
    assert a == b
    assert a["rho"] == c["rho"]
    k, X, g, xy = _inputs()
    assert KS.analyze_specimen("s", k, X, g, xy, B=99, seed=3) == KS.analyze_specimen("s", k, X, g, xy, B=99, seed=3)


def test_assemble_bh_and_cohort():
    def ok(rho, p):
        return {"status": "ok", "rho": rho, "p_two_sided": p}

    ne = {"status": "not_estimable", "reason_code": "eligibility_lt30"}
    per = {"a": {"n_eligible": 50, **{t: ok(0.3, 0.001) for t in KS.TESTS}},
           "HCC22-088-P4-S2_1i_rep": {"n_eligible": 50, **{t: ok(-0.2, 0.5) for t in KS.TESTS}},
           "c": {"n_eligible": 5, **{t: dict(ne) for t in KS.TESTS}}}
    out = KS.assemble(per, {"a": "ESR1_WT", "HCC22-088-P4-S2_1i_rep": "ESR1_D538G", "c": "ESR1_WT"})
    co = out["cohort"]["primary"]
    assert co["n_estimable"] == 2 and co["n_not_estimable"] == 1 and co["n_positive"] == 1
    assert co["n_significant_q_lt_0.05"] == 1 and co["significant_specimens"] == ["a"]
    assert per["a"]["primary"]["q"] == pytest.approx(0.002)
    assert out["focal_P4_S2"]["primary"]["rho"] == -0.2
    assert "DESCRIPTIVE" in out["descriptive_mutant_vs_wt"]["note"]


def test_refuses_to_overwrite_existing_ks1(tmp_path):
    sha12 = (REPO / "Benchmarking" / "midkine_remediation" / "ACTIVE_PROTOCOL").read_text().strip()  # the protocol.json hash
    out = tmp_path / sha12 / "ks1"
    out.mkdir(parents=True)
    (out / "x.csv").write_text("a")
    with pytest.raises(Exception, match="not an empty directory|never overwritten"):
        KS.run(tmp_path, out, REPO / "Benchmarking" / "midkine_remediation" / "protocol.json")
