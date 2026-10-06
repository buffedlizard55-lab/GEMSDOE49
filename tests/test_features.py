import numpy as np

from gemsdoe49.features import cross_gradient_alignment, shifted_alignment_control


def grid_fields(shape=(61, 61)):
    rows, cols = np.indices(shape, dtype=np.float32)
    return cols * 100.0, cols * -250.0, cols * 30.0


def test_parallel_and_antiparallel_gradients_have_high_alignment():
    a, b, c = grid_fields()
    valid = np.ones(a.shape, dtype=bool)
    result = cross_gradient_alignment(
        {"gravity": a, "rtp": b, "topography": c}, valid, sigma_px=1.0
    )
    center = (30, 30)
    assert result.score[center] > 0.95
    assert result.orthogonality_flag[center] < 0.05
    assert result.pair_support[center] > 0.95


def test_orthogonality_is_a_flag_not_positive_alignment():
    rows, cols = np.indices((61, 61), dtype=np.float32)
    result = cross_gradient_alignment(
        {"gravity": cols, "rtp": rows, "topography": rows * 2.0},
        np.ones((61, 61), dtype=bool),
        sigma_px=1.0,
    )
    center = (30, 30)
    assert np.isclose(result.score[center], 1.0 / 3.0, atol=0.03)
    assert result.orthogonality_flag[center] > 0.60
    # The operator does not turn a high orthogonality flag into extra evidence.
    assert result.score[center] < result.orthogonality_flag[center]


def test_zero_gradient_pair_contributes_nothing_but_other_pairs_remain():
    rows, cols = np.indices((61, 61), dtype=np.float32)
    zeros = np.zeros((61, 61), dtype=np.float32)
    result = cross_gradient_alignment(
        {"gravity": zeros, "rtp": cols, "topography": cols * 2.0},
        np.ones((61, 61), dtype=bool),
        sigma_px=1.0,
    )
    # Two zero-gradient pairs contribute zero; the remaining parallel rtp/topography pair contributes 1/3.
    assert np.isclose(result.score[30, 30], 1.0 / 3.0, atol=0.03)
    assert result.orthogonality_flag[30, 30] == 0.0


def test_operator_rejects_degenerate_scale_and_invalid_grid():
    a, b, c = grid_fields(shape=(8, 9))
    valid = np.ones(a.shape, dtype=bool)
    with np.testing.assert_raises(ValueError):
        cross_gradient_alignment({"a": a, "b": b, "c": c}, valid, sigma_px=0.0)
    with np.testing.assert_raises(ValueError):
        cross_gradient_alignment({"a": a, "b": b, "c": c}, valid, pixel_size_m=0.0)
    with np.testing.assert_raises(ValueError):
        cross_gradient_alignment({"a": a.ravel(), "b": b.ravel(), "c": c.ravel()}, valid.ravel())


def test_shift_control_excludes_wrap_and_preserves_range():
    a, b, c = grid_fields()
    result = cross_gradient_alignment(
        {"gravity": a, "rtp": b, "topography": c},
        np.ones(a.shape, dtype=bool),
        sigma_px=1.0,
    )
    score, support = shifted_alignment_control(result, shifted_field="rtp", shift_px=(0, 7))
    assert np.all((score >= 0) & (score <= 1))
    assert not support[:, :7].any()
    assert support[:, 7:].any()
