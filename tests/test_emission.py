import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np

from gems49.emission import nms_topk, dots_to_prediction


def test_nms_emits_exactly_n_dots_not_suppression_discs():
    rng = np.random.default_rng(0)
    h = w = 400
    score = rng.random((h, w)).astype(np.float32)
    footprint = np.ones((h, w), bool)
    catalogue = np.zeros((h, w), bool)
    chosen = nms_topk(score, catalogue, footprint, n_dots=1000, nms_radius_px=2.0)
    assert int(chosen.sum()) == 1000, int(chosen.sum())


def test_nms_minimum_spacing():
    rng = np.random.default_rng(1)
    h = w = 300
    score = rng.random((h, w)).astype(np.float32)
    footprint = np.ones((h, w), bool)
    catalogue = np.zeros((h, w), bool)
    chosen = nms_topk(score, catalogue, footprint, n_dots=500, nms_radius_px=2.0)
    rc = np.argwhere(chosen)
    from scipy.spatial import cKDTree
    tree = cKDTree(rc)
    d, _ = tree.query(rc, k=2)
    assert d[:, 1].min() > 2.0 - 1e-9


def test_nms_never_emits_on_catalogue():
    rng = np.random.default_rng(2)
    h = w = 200
    score = rng.random((h, w)).astype(np.float32)
    footprint = np.ones((h, w), bool)
    catalogue = np.zeros((h, w), bool)
    catalogue[50:150, 50:150] = True
    score[50:150, 50:150] += 5.0  # make catalogue the global maximum
    chosen = nms_topk(score, catalogue, footprint, n_dots=100, nms_radius_px=2.0)
    assert not (chosen & catalogue).any()
    assert int(chosen.sum()) == 100


def test_dots_to_prediction_values():
    footprint = np.ones((10, 10), bool)
    chosen = np.zeros((10, 10), bool)
    chosen[3, 4] = True
    p = dots_to_prediction(chosen, footprint, value=1.0)
    assert p.dtype == np.float32 and p.sum() == 1.0 and p.max() == 1.0
