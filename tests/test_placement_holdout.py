import numpy as np

from gemsdoe49.holdout import apply_analysis_role, evaluate_spatial_holdout, spatial_column_folds
from gemsdoe49.placement import greedy_kernel_suppress


def test_greedy_suppression_is_deterministic_and_respects_distance():
    yy, xx = np.indices((40, 40))
    score = (yy * 101 + xx).astype(np.float32)
    valid = np.ones(score.shape, dtype=bool)
    a, receipt = greedy_kernel_suppress(score, valid, 30, min_distance_px=3.0)
    b, _ = greedy_kernel_suppress(score, valid, 30, min_distance_px=3.0)
    np.testing.assert_array_equal(a, b)
    assert int((a > 0).sum()) == 30
    coords = np.argwhere(a > 0)
    for i in range(len(coords)):
        if i + 1 < len(coords):
            d = np.sqrt(((coords[i + 1:] - coords[i]) ** 2).sum(axis=1))
            assert np.all(d > 3.0)
    assert receipt["selected"] == 30


def test_spatial_blocks_are_disjoint_and_hold_guard_gaps():
    foot = np.ones((50, 100), dtype=bool)
    folds = spatial_column_folds(foot.shape, foot, n_folds=5, boundary_guard_px=3)
    domains = [fold["domain"] for fold in folds]
    for i in range(len(domains)):
        for j in range(i + 1, len(domains)):
            assert not np.any(domains[i] & domains[j])
    assert sum(int(x.sum()) for x in domains) < foot.size
    for left, right in zip(folds, folds[1:]):
        assert right["col_start"] - left["col_end"] >= 6


def test_spatial_holdout_emits_equal_mass_and_reports_every_fold():
    h, w = 100, 250
    rows, cols = np.indices((h, w))
    footprint = np.ones((h, w), dtype=bool)
    support = footprint.copy()
    labels = np.zeros((h, w), dtype=bool)
    for x in [20, 60, 100, 150, 200, 230]:
        labels[20:80, x] = True
    candidate = np.exp(-((cols % 50) - 20) ** 2 / 300.0).astype(np.float32)
    baseline = np.exp(-((cols % 50) - 30) ** 2 / 300.0).astype(np.float32)
    random = np.random.default_rng(7).random((h, w), dtype=np.float32)
    result = evaluate_spatial_holdout(
        {
            "cross_gradient": candidate,
            "single_gravity": baseline,
            "control_random": random,
            "control_shift_test": np.roll(candidate, 8, axis=1),
        },
        labels,
        footprint,
        support,
        total_budget=300,
        min_distance_px=2.0,
        n_folds=5,
        boundary_guard_px=3,
    )
    assert len(result["results"]["cross_gradient"]["folds"]) == 5
    assert sum(result["fold_budgets"]) == 300
    assert result["results"]["cross_gradient"]["pooled"]["prediction_mass"] == 300


def test_sensitivity_thresholds_cannot_promote_candidate():
    raw = {"promotion_pass": True, "results": {"fixture": "unchanged"}}
    sensitivity = apply_analysis_role(raw, "sensitivity")
    assert sensitivity["thresholds_met_on_this_run"] is True
    assert sensitivity["promotion_pass"] is False
    assert sensitivity["results"] is raw["results"]
    primary = apply_analysis_role(raw, "primary")
    assert primary["promotion_pass"] is True


def test_unknown_analysis_role_fails_closed():
    with np.testing.assert_raises(ValueError):
        apply_analysis_role({"promotion_pass": True}, "tuned")
