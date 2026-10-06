"""Line-by-line verification of the official metric.

Every assertion here corresponds to a sentence on the competition problem-description page
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gems49 import metric as M  # noqa: E402
from gems49.spec import ALPHA, BETA, R_M  # noqa: E402


def test_official_constants():
    """The page states alpha = 0.2, beta = 0.8 and a 300 m triangular kernel."""
    assert (ALPHA, BETA, R_M) == (0.2, 0.8, 300.0)


def test_kernel_is_triangular_with_300m_support():
    d = np.array([0.0, 100.0, 150.0, 200.0, 300.0, 400.0])
    k = M.triangle_kernel(d)
    assert np.allclose(k, [1.0, 2 / 3, 0.5, 1 / 3, 0.0, 0.0])


def test_worked_example_from_the_page():
    """The page's own numbers: TP_w 3.00, FP_w 1.89, FN_w 2.00 -> 0.60."""
    v = M.dti_from_components(tp=3.00, fp=1.89, n_truth=3 + 2)
    assert abs(v - 0.602652) < 1e-5, v  # 3.00/(3.00+0.378+1.60)


def test_identity_tp_plus_fn_equals_n_truth():
    rng = np.random.default_rng(3)
    pred = (rng.random((40, 40)) < 0.1).astype(float)
    truth = (rng.random((40, 40)) < 0.1)
    c = M.dti_components(pred, truth)
    assert abs(c["TPw"] + c["FNw"] - c["n_truth"]) < 1e-9


def test_perfect_prediction_scores_one():
    truth = np.zeros((60, 60), dtype=bool)
    truth[20:40, 30] = True
    c = M.dti_components(truth.astype(float), truth)
    assert abs(c["DTI"] - 1.0) < 1e-9
    assert c["FPw"] == 0.0


def test_zero_prediction_is_zero():
    truth = np.zeros((60, 60), dtype=bool)
    truth[20:40, 30] = True
    c = M.dti_components(np.zeros((60, 60)), truth)
    assert c["TPw"] == 0.0 and c["DTI"] == 0.0


def _brute_force(pred, truth, res=100.0, r_m=300.0, alpha=0.2, beta=0.8):
    """Literal transcription of the published equations, one nested loop per term."""
    h, w = pred.shape
    tp = 0.0
    for gy in range(h):
        for gx in range(w):
            if not truth[gy, gx]:
                continue
            best = 0.0
            for yy in range(h):
                for xx in range(w):
                    d = np.hypot((yy - gy) * res, (xx - gx) * res)
                    if d <= r_m:
                        best = max(best, pred[yy, xx] * max(1.0 - d / r_m, 0.0))
            tp += best
    fp = 0.0
    for yy in range(h):
        for xx in range(w):
            p = pred[yy, xx]
            if p <= 0:
                continue
            bestk = 0.0
            for gy in range(h):
                for gx in range(w):
                    if truth[gy, gx]:
                        d = np.hypot((yy - gy) * res, (xx - gx) * res)
                        bestk = max(bestk, max(1.0 - d / r_m, 0.0))
            fp += p * (1.0 - bestk)
    n_g = int(truth.sum())
    fn = n_g - tp
    return tp, fp, fn, tp / (tp + alpha * fp + beta * fn + 1e-9)


def test_matches_literal_brute_force():
    rng = np.random.default_rng(11)
    pred = np.zeros((14, 14))
    pred[5, 6] = 0.9
    pred[7, 7] = 0.4
    pred[3, 9] = 1.0
    pred[10, 4] = 0.55
    truth = np.zeros((14, 14), dtype=bool)
    truth[5, 7] = truth[9, 5] = truth[1, 2] = True

    fast = M.dti_components(pred, truth)
    tp, fp, fn, dti = _brute_force(pred, truth)
    assert abs(fast["TPw"] - tp) < 1e-9
    assert abs(fast["FPw"] - fp) < 1e-9
    assert abs(fast["DTI"] - dti) < 1e-9


def test_masking_removes_credit_and_penalty():
    """A prediction pixel on a masked pixel must neither earn nor radiate credit."""
    truth = np.zeros((20, 20), dtype=bool)
    truth[10, 10] = True
    mask = np.zeros((20, 20), dtype=bool)
    mask[10, 9] = True  # one pixel from truth: would otherwise earn 2/3 credit
    pred = np.zeros((20, 20))
    pred[10, 9] = 1.0
    assert M.dti_components(pred, truth)["TPw"] > 0.6
    assert M.dti_components(M.scoring_domain(pred, mask), truth)["TPw"] == 0.0


def test_marginal_rule_algebra():
    """Adding a unit pixel raises D by exactly alpha; the bar is k > alpha*s."""
    assert abs(M.marginal_bar(0.30) - 0.06) < 1e-12
    assert abs(M.marginal_max_distance(0.30) - 282.0) < 1e-9


def test_marginal_bar_is_exact():
    """Direct computation of the marginal rule, including the shadowing caveat.

    A newly added pixel pays for itself iff BOTH hold:
      (a) it becomes the best-weighted prediction for at least one scored truth pixel
          (otherwise its credit is 0 -- it is *shadowed*), and
      (b) that credit exceeds ``alpha * DTI``.
    """
    # two isolated truth pixels 20 px apart, so neither can shadow the other
    truth = np.zeros((41, 41), dtype=bool)
    truth[20, 5] = True
    truth[20, 25] = True
    pred = np.zeros((41, 41))
    pred[20, 5] = 1.0
    s = M.dti_components(pred, truth)["DTI"]
    bar = M.marginal_bar(s)
    assert abs(s - 1.0 / 1.8) < 1e-9       # D = 0.2*(1+0) + 0.8*2 = 1.8 (+ the 1e-9 guard)
    assert abs(bar - 1.0 / 9.0) < 1e-9     # 0.2*0.5556

    for d_px in (1, 2, 3):
        credit = float(M.triangle_kernel(d_px * 100.0))
        p = pred.copy()
        p[20 - d_px, 25] = 1.0
        after = M.dti_components(p, truth)["DTI"]
        if credit > bar:
            assert after > s, (d_px, credit, bar, after, s)
        else:
            assert after < s, (d_px, credit, bar, after, s)

    # and the closed form: the added pixel moves D by exactly alpha
    p = pred.copy()
    p[18, 25] = 1.0
    c0 = M.dti_components(pred, truth)
    c1 = M.dti_components(p, truth)
    d0 = 0.2 * (c0["TPw"] + c0["FPw"]) + 0.8 * c0["n_truth"]
    d1 = 0.2 * (c1["TPw"] + c1["FPw"]) + 0.8 * c1["n_truth"]
    assert abs((d1 - d0) - 0.2) < 1e-12


def test_shadowed_pixel_never_pays():
    """A pixel whose credit is below the incumbent best for every truth pixel adds FP only."""
    truth = np.zeros((41, 41), dtype=bool)
    truth[20, 20:29] = True
    pred = np.zeros((41, 41))
    pred[20, 20] = pred[20, 23] = pred[20, 26] = 1.0
    s = M.dti_components(pred, truth)["DTI"]
    shadowed = pred.copy()
    shadowed[19, 23] = 1.0     # credit 2/3 to truth 23 but the incumbent there is 1.0 -> shadowed
    assert M.dti_components(shadowed, truth)["DTI"] < s


def test_dti_reduces_to_published_denominator_identity():
    """With alpha+beta=1 and FN=|G|-TP, the denominator equals alpha*(TP+FP)+beta*|G|."""
    rng = np.random.default_rng(5)
    pred = (rng.random((30, 30)) < 0.2) * rng.random((30, 30))
    truth = rng.random((30, 30)) < 0.15
    c = M.dti_components(pred, truth)
    lhs = c["TPw"] + ALPHA * c["FPw"] + BETA * c["FNw"]
    rhs = ALPHA * (c["TPw"] + c["FPw"]) + BETA * c["n_truth"]
    assert abs(lhs - rhs) < 1e-9


def test_graded_prediction_is_dominated_by_its_support():
    """Metric-optimality of binary emission, proved by direct computation."""
    rng = np.random.default_rng(7)
    truth = np.zeros((80, 90), dtype=bool)
    truth[10:60, 40] = True
    truth[30, 20:70] = True
    base = np.zeros((80, 90))
    base[9:61, 40] = 1.0
    base[30, 19:71] = 1.0
    scores = [M.dti_components(base * lam, truth)["DTI"] for lam in (0.2, 0.4, 0.6, 0.8, 1.0)]
    assert scores == sorted(scores)
    assert scores[-1] > scores[0]


def test_discrete_lattice_offsets_within_support():
    offs = M.kernel_offsets(100.0, 300.0)
    d = np.array([np.hypot(dy, dx) * 100.0 for dy, dx, _ in offs])
    assert d.max() <= 300.0 + 1e-9
    assert len(offs) == 29  # all integer offsets (dy,dx) with 0 <= hypot*100 <= 300


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([__file__, "-q"]))
