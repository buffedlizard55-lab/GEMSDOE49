"""Hide-and-recover instrument: segmentation, folds, leakage guards and scorer parity."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import holdout  # noqa: E402
from gems49 import metric as M  # noqa: E402


def _synthetic_two_lines(shape=(60, 60)):
    cat = np.zeros(shape, dtype=bool)
    cat[10, 5:40] = True
    cat[40, 5:40] = True
    return cat


def test_segmentation_finds_two_components():
    cat = _synthetic_two_lines()
    seg, sizes = holdout.segment_catalogue(cat)
    assert seg.max() == 2
    assert sorted(sizes.tolist()) == [35, 35]
    assert (seg[cat] > 0).all()
    assert (seg[~cat] == 0).all()


def test_speck_filter_drops_tiny_segments():
    cat = _synthetic_two_lines()
    cat[5, 5] = True  # a single-pixel speck
    seg, sizes = holdout.segment_catalogue(cat, min_px=3)
    assert seg.max() == 2


def test_fold_assignment_is_complete_and_balanced():
    cat = np.zeros((120, 120), dtype=bool)
    for i in range(0, 100, 4):
        cat[10 + i % 80, 10:90] = True
    seg, sizes = holdout.segment_catalogue(cat)
    fa = holdout.fold_assignment(seg, seg.max(), 4, seed=1)
    assert fa.shape[0] == seg.max()
    assert set(np.unique(fa)) == {0, 1, 2, 3}
    load = np.zeros(4)
    for s, f in enumerate(fa, start=1):
        load[f] += sizes[s - 1]
    assert load.max() / max(load.min(), 1) < 2.0


def test_leakage_probe_is_exactly_zero():
    """Scoring the retained catalogue itself must earn nothing; otherwise the fold leaks."""
    cat = _synthetic_two_lines()
    fp = np.ones(cat.shape, dtype=bool)
    seg, sizes = holdout.segment_catalogue(cat)
    fa = holdout.fold_assignment(seg, seg.max(), 2, seed=0)
    for k in (0, 1):
        fold = holdout.build_fold(cat, fp, fa, seg, k)
        assert holdout.leakage_probe(fold)["TPw"] == 0.0
        assert holdout.null_probe(fold)["DTI"] == 0.0


def test_dots_on_holdout_truth_score_high():
    cat = _synthetic_two_lines()
    fp = np.ones(cat.shape, dtype=bool)
    seg, _ = holdout.segment_catalogue(cat)
    fa = holdout.fold_assignment(seg, seg.max(), 2, seed=0)
    fold = holdout.build_fold(cat, fp, fa, seg, 0, sep_m=0.0)
    truth = fold["truth_unrestricted"]
    assert truth.any()
    pred = truth.astype(np.float64)
    r = holdout.score_fold(pred, fold, truth_key="truth_unrestricted")
    assert r["DTI"] > 0.99


def test_collared_truth_is_a_subset_of_unrestricted():
    cat = _synthetic_two_lines()
    fp = np.ones(cat.shape, dtype=bool)
    seg, _ = holdout.segment_catalogue(cat)
    fa = holdout.fold_assignment(seg, seg.max(), 2, seed=0)
    fold = holdout.build_fold(cat, fp, fa, seg, 0, sep_m=600.0)
    assert (fold["truth_collared"] & ~fold["truth_unrestricted"]).sum() == 0


def test_fold_scorer_matches_the_reference_metric():
    rng = np.random.default_rng(0)
    cat = np.zeros((80, 80), dtype=bool)
    cat[20, 10:70] = True
    cat[60, 10:70] = True
    fp = np.ones(cat.shape, dtype=bool)
    seg, _ = holdout.segment_catalogue(cat)
    fa = holdout.fold_assignment(seg, seg.max(), 2, seed=0)
    fold = holdout.build_fold(cat, fp, fa, seg, 0, sep_m=0.0)

    pred = np.zeros(cat.shape, dtype=np.float32)
    idx = rng.choice(np.flatnonzero(fp), size=50, replace=False)
    pred.ravel()[idx] = 1.0
    pred[fold["train_mask"]] = 0.0

    sc = holdout.FoldScorer(fold, "truth_unrestricted")
    got = sc.score(pred)
    ref = M.dti_components(np.where(fold["valid"], pred, 0.0), fold["truth_unrestricted"])
    assert abs(got["DTI"] - ref["DTI"]) < 1e-9
    assert abs(got["TPw"] - ref["TPw"]) < 1e-9
    assert abs(got["FPw"] - ref["FPw"]) < 1e-6


def test_scorable_domain_excludes_the_retained_map():
    cat = _synthetic_two_lines()
    fp = np.ones(cat.shape, dtype=bool)
    seg, _ = holdout.segment_catalogue(cat)
    fa = holdout.fold_assignment(seg, seg.max(), 2, seed=0)
    fold = holdout.build_fold(cat, fp, fa, seg, 0)
    dom = holdout.scorable_domain(fold)
    assert not (dom & fold["train_mask"]).any()
    assert dom.sum() == int(fp.sum() - fold["train_mask"].sum())
