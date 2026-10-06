import numpy as np

from gemsdoe49.metric import dti_components, dti_from_components


def brute_force(pred, truth, radius=3.0):
    truth_rc = np.argwhere(truth)
    pred_rc = np.argwhere(pred > 0)
    tp = 0.0
    fn = 0.0
    for tr, tc in truth_rc:
        if pred_rc.size:
            dist = np.sqrt(((pred_rc - np.array([tr, tc])) ** 2).sum(axis=1))
            credit = np.max(pred[pred_rc[:, 0], pred_rc[:, 1]] * np.maximum(1 - dist / radius, 0))
        else:
            credit = 0.0
        tp += float(credit)
        fn += 1.0 - float(credit)
    fp = 0.0
    for pr, pc in pred_rc:
        if truth_rc.size:
            dist = np.sqrt(((truth_rc - np.array([pr, pc])) ** 2).sum(axis=1))
            nearest = np.max(np.maximum(1 - dist / radius, 0))
        else:
            nearest = 0.0
        fp += float(pred[pr, pc]) * (1.0 - float(nearest))
    return tp, fp, fn


def test_official_metric_worked_example_equation():
    # The competition page reports 0.60 to two decimals; the unrounded value is 0.60265.
    assert np.isclose(dti_from_components(3.0, 1.89, 2.0), 0.6026516673, atol=1e-9)


def test_fast_metric_matches_bruteforce_random_cases():
    rng = np.random.default_rng(49)
    for _ in range(15):
        truth = rng.random((11, 13)) < 0.08
        pred = np.zeros((11, 13), dtype=np.float32)
        pred[rng.random((11, 13)) < 0.06] = 1.0
        # Add non-binary confidence values without changing the support.
        pred[pred > 0] = rng.uniform(0.1, 1.0, size=int((pred > 0).sum()))
        expected = brute_force(pred, truth)
        actual = dti_components(pred, truth)
        np.testing.assert_allclose(
            [actual.tp_weighted, actual.fp_weighted, actual.fn_weighted],
            expected,
            rtol=1e-11,
            atol=1e-11,
        )


def test_metric_rejects_invalid_values_and_nan():
    truth = np.zeros((4, 4), dtype=bool)
    truth[1, 1] = True
    with np.testing.assert_raises(ValueError):
        dti_components(np.full((4, 4), np.nan), truth)
    with np.testing.assert_raises(ValueError):
        dti_components(np.full((4, 4), 1.01), truth)
