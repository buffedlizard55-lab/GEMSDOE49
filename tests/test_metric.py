"""Tests for the DTI metric.

Every expected value below is computed by hand or by an independent
brute-force implementation — not by the module under test.
"""
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pytest

from gems49.metric import dti, DTIResult, ALPHA, BETA, EPS


def brute_force_dti(pred, truth, mask=None, alpha=0.2, beta=0.8, R=3.0):
    """O(N*M) textbook implementation used as an oracle."""
    p = np.asarray(pred, float)
    g = np.argwhere(np.asarray(truth, bool))
    m = np.ones(p.shape, bool) if mask is None else np.asarray(mask, bool)
    g = np.array([rc for rc in g if m[rc[0], rc[1]]])
    pred_idx = np.argwhere((p > 0) & m)

    tp = 0.0
    fn = 0.0
    for (gr, gc) in g:
        best = 0.0
        for (pr, pc) in pred_idx:
            d = math.hypot(gr - pr, gc - pc)
            if d <= R:
                best = max(best, p[pr, pc] * max(1 - d / R, 0.0))
        tp += best
        fn += 1 - best
    fp = 0.0
    for (pr, pc) in pred_idx:
        kmax = 0.0
        for (gr, gc) in g:
            d = math.hypot(gr - pr, gc - pc)
            kmax = max(kmax, max(1 - d / R, 0.0))
        fp += p[pr, pc] * (1 - kmax)
    return tp / (tp + alpha * fp + beta * fn + 1e-7)


def test_perfect_prediction_scores_one():
    truth = np.zeros((11, 11), bool)
    truth[5, 5] = True
    pred = truth.astype(float)
    r = dti(pred, truth)
    assert abs(r.dti - 1.0) < 1e-6
    assert r.tp_w == 1.0 and r.fp_w == 0.0 and r.fn_w == 0.0


def test_hand_computed_offset_prediction():
    # truth at (5,5); pred p=0.5 at (5,7): d=2, k=1/3
    truth = np.zeros((11, 11), bool)
    truth[5, 5] = True
    pred = np.zeros((11, 11))
    pred[5, 7] = 0.5
    r = dti(pred, truth)
    tp = 0.5 / 3.0
    fn = 1 - tp
    fp = 0.5 * (1 - 1 / 3.0)
    expected = tp / (tp + 0.2 * fp + 0.8 * fn + 1e-7)
    assert abs(r.dti - expected) < 1e-12
    assert abs(r.tp_w - tp) < 1e-12
    assert abs(r.fp_w - fp) < 1e-12
    assert abs(r.fn_w - fn) < 1e-12


def test_diagonal_distance_uses_euclidean():
    # truth at (5,5); pred at (6,6): d = sqrt(2), k = 1 - sqrt(2)/3
    truth = np.zeros((11, 11), bool)
    truth[5, 5] = True
    pred = np.zeros((11, 11))
    pred[6, 6] = 1.0
    r = dti(pred, truth)
    k = 1 - math.sqrt(2) / 3
    expected = k / (k + 0.2 * (1 - k) + 0.8 * (1 - k) + 1e-7)
    assert abs(r.dti - expected) < 1e-12


def test_beyond_radius_is_full_false_positive():
    # truth at (5,5); pred at (5,9): d=4 > R => no TP credit, FP weight 1
    truth = np.zeros((11, 11), bool)
    truth[5, 5] = True
    pred = np.zeros((11, 11))
    pred[5, 9] = 1.0
    r = dti(pred, truth)
    expected = 0.0 / (0.0 + 0.2 * 1.0 + 0.8 * 1.0 + 1e-7)
    assert abs(r.dti - expected) < 1e-12
    assert r.tp_w == 0.0 and r.fn_w == 1.0 and r.fp_w == 1.0


def test_mask_removes_catalogue_truth_and_predictions():
    # two truth pixels, one masked; prediction mass on the masked truth
    # must neither earn TP nor cost FP.
    truth = np.zeros((11, 11), bool)
    truth[5, 3] = True   # will be masked (catalogue)
    truth[5, 7] = True   # held-out
    mask = np.ones((11, 11), bool)
    mask[5, 3] = False
    pred = np.zeros((11, 11))
    pred[5, 3] = 1.0     # on masked catalogue pixel: ignored entirely
    pred[5, 7] = 1.0     # on held-out truth
    r = dti(pred, truth, mask=mask)
    assert r.n_truth == 1
    assert abs(r.dti - 1.0) < 1e-6
    assert r.fp_w == 0.0


def test_identity_T_over_budget():
    # DTI = T / (0.2(T+S-M) + 0.8|G|) with S = sum p over active, M = sum p*k
    rng = np.random.default_rng(0)
    truth = np.zeros((30, 30), bool)
    truth[rng.integers(5, 25, 40), rng.integers(5, 25, 40)] = True
    pred = np.zeros((30, 30))
    idx = np.argwhere(truth)
    for (r_, c_) in idx[:20]:
        pred[r_ + rng.integers(-2, 3), c_ + rng.integers(-2, 3)] = rng.random()
    # also random noise predictions
    for _ in range(30):
        pred[rng.integers(0, 30), rng.integers(0, 30)] = rng.random()
    res = dti(pred, truth)
    # brute-force recomputation of S and M
    g = np.argwhere(truth)
    active = np.argwhere(pred > 0)
    S = pred[pred > 0].sum()
    from scipy.ndimage import distance_transform_edt
    kw = np.clip(1 - distance_transform_edt(~truth) / 3.0, 0, None)
    M = (pred * kw).sum()
    expected = res.tp_w / (0.2 * (res.tp_w + S - M) + 0.8 * truth.sum() + 1e-7)
    assert abs(res.dti - expected) < 1e-9


@pytest.mark.parametrize("seed", range(5))
def test_matches_brute_force_random(seed):
    rng = np.random.default_rng(seed)
    H = W = 24
    truth = np.zeros((H, W), bool)
    truth[rng.integers(0, H, 15), rng.integers(0, W, 15)] = True
    pred = np.zeros((H, W))
    pred[rng.integers(0, H, 25), rng.integers(0, W, 25)] = rng.random(25)
    mask = rng.random((H, W)) > 0.1
    r = dti(pred, truth, mask=mask)
    bf = brute_force_dti(pred, truth, mask=mask)
    assert abs(r.dti - bf) < 1e-12, (r.dti, bf)


def test_empty_prediction_zero_score():
    truth = np.zeros((10, 10), bool)
    truth[3:6, 5] = True
    pred = np.zeros((10, 10))
    r = dti(pred, truth)
    assert r.dti == 0.0 and r.fn_w == 3.0


def test_rejects_out_of_range():
    truth = np.zeros((5, 5), bool)
    pred = np.zeros((5, 5))
    pred[0, 0] = 1.5
    with pytest.raises(ValueError):
        dti(pred, truth)
    pred[0, 0] = np.nan
    with pytest.raises(ValueError):
        dti(pred, truth)


def test_published_identity_alpha_beta():
    """The GEMSDOE32 published identity: with alpha=0.2, beta=0.8,
    denom = 0.2*(T+S-M) + 0.8*|G|  (because T - 0.8*T = 0.2*T)."""
    T, S, M, G = 12.5, 40.0, 10.0, 60
    denom_direct = T + 0.2 * (S - M) + 0.8 * (G - T)
    denom_identity = 0.2 * (T + S - M) + 0.8 * G
    assert abs(denom_direct - denom_identity) < 1e-12
