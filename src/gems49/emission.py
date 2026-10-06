"""Metric-aware placement of prediction support.

The published metric has an exact marginal rule.  Because ``alpha + beta = 1`` and
``FN_w = |G| - TP_w``, the denominator is ``D = alpha*(TP_w + FP_w) + beta*|G|``.  Adding one
prediction pixel of unit weight whose kernel credit is ``k`` raises ``TP_w`` by ``k`` and the
false-positive term by ``1 - k``, so ``TP_w + FP_w`` rises by exactly 1 and ``D`` rises by
exactly ``alpha``, whatever ``k`` is.  Hence

    the move pays  <=>  (TP_w + k)/(D + alpha) > s  <=>  k > alpha * s.

At ``s = 0.30`` the bar is 0.06, i.e. a pixel pays whenever it lands within
``300*(1-0.06) = 282 m`` of a *scored* truth pixel.  The metric is therefore a distance
filter, not a probability-calibration filter, and the whole placement problem reduces to:
put as much support as possible inside 300 m of the hidden truth, and keep total support
small.  This module implements that selection directly.
"""

from __future__ import annotations

import numpy as np

from .metric import marginal_bar, marginal_max_distance  # re-exported for callers
from .spec import GRID

__all__ = [
    "marginal_bar",
    "marginal_max_distance",
    "select_support",
    "spacing_from_kernel",
    "coverage_per_dot",
]


def spacing_from_kernel(px: int = 3) -> float:
    """Distance in pixels at which a chain of dots maximises credit delivered per dot.

    For a dot train of spacing ``s`` pixels laid along a 1-px-wide truth trace, each dot owns
    about ``s`` truth pixels and delivers an average kernel credit of ``1 - s/(2R_px)``.
    Credit per dot is ``s*(1 - s/(2R_px))``, maximised at ``s = R_px`` (3 px at 100 m).
    """
    r_px = GRID["res_m"] and 3
    return float(px if px <= r_px else r_px)


def coverage_per_dot(spacing_px: float, r_px: float = 3.0) -> float:
    """Total truth credit delivered per dot for a perfectly-placed dot train."""
    if spacing_px <= r_px:
        return spacing_px * (1.0 - spacing_px / (2.0 * r_px))
    return float(r_px * r_px / (2.0 * spacing_px)) * 2.0  # decoupled tails


def select_support(
    score: np.ndarray,
    valid: np.ndarray,
    *,
    n_target: int,
    min_sep_px: int = 0,
    score_floor: float | None = None,
) -> np.ndarray:
    """Greedy top-score selection with a minimum-separation (Poisson-disk) constraint.

    ``score`` is the continuous field in [0,1]; ``valid`` the allowed domain.  Returns a
    float32 field that is 0 where nothing was selected and 1.0 where a dot was placed --
    a *binary* support, because the metric strictly prefers it: for a support scaled by
    ``lambda``, ``DTI(lambda) = lambda*T / (a*lambda*(T+F) + b*|G|)`` is increasing in
    ``lambda``, so any graded map is dominated by its own thresholded support.
    """
    score = np.asarray(score, dtype=np.float64)
    allow = np.asarray(valid, dtype=bool) & np.isfinite(score)
    s = np.where(allow, score, -np.inf)
    if score_floor is not None:
        s = np.where(s >= score_floor, s, -np.inf)

    flat = s.ravel()
    order = np.argsort(-flat, kind="stable")
    h, w = s.shape
    blocked = np.zeros(flat.shape, dtype=bool)
    radii = _disc_offsets(min_sep_px)

    out = np.zeros(flat.shape, dtype=np.float32)
    taken = 0
    for i in order:
        if not np.isfinite(flat[i]) or flat[i] <= 0.0:
            break
        if blocked[i]:
            continue
        y, x = divmod(int(i), w)
        out[i] = 1.0
        taken += 1
        for dy, dx in radii:
            yy, xx = y + dy, x + dx
            if 0 <= yy < h and 0 <= xx < w:
                blocked[yy * w + xx] = True
        if taken >= n_target:
            break
    return out.reshape(h, w)


def _disc_offsets(radius: int) -> np.ndarray:
    if radius <= 0:
        return np.zeros((1, 2), dtype=np.int32)
    r = int(radius)
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    m = (yy * yy + xx * xx) <= r * r
    return np.stack([yy[m], xx[m]], axis=1).astype(np.int32)


def selection_order(score: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Flat indices of valid, finite pixels ordered by descending score (computed once)."""
    s = np.where(np.asarray(valid, dtype=bool) & np.isfinite(score),
                 np.asarray(score, dtype=np.float64), -np.inf).ravel()
    return np.argsort(-s, kind="stable")


def select_from_order(
    order: np.ndarray,
    shape: tuple[int, int],
    valid: np.ndarray,
    *,
    n_target: int,
    min_sep_px: int = 0,
    score: np.ndarray | None = None,
) -> np.ndarray:
    """Greedy Poisson-disk selection from a precomputed descending-score order."""
    h, w = shape
    flat_valid = np.asarray(valid, dtype=bool).ravel()
    blocked = np.zeros(order.shape, dtype=bool)
    out = np.zeros(order.shape, dtype=np.float32)
    radii = _disc_offsets(min_sep_px)
    taken = 0
    for i in order:
        if not flat_valid[i] or blocked[i]:
            continue
        if score is not None and not (score.ravel()[i] > 0.0):
            break
        y, x = divmod(int(i), w)
        out[i] = 1.0
        taken += 1
        for dy, dx in radii:
            yy, xx = y + dy, x + dx
            if 0 <= yy < h and 0 <= xx < w:
                blocked[yy * w + xx] = True
        if taken >= n_target:
            break
    return out.reshape(shape)


def candidate_pool(score: np.ndarray, valid: np.ndarray, n_candidates: int) -> np.ndarray:
    """Indices of the ``n_candidates`` highest-scoring valid pixels (flat, row-major)."""
    s = np.where(valid & np.isfinite(score), np.asarray(score, dtype=np.float64), -np.inf).ravel()
    k = min(int(n_candidates), int(np.isfinite(s).sum()))
    pool = np.argpartition(-s, k - 1)[:k]
    return pool[np.argsort(-s[pool])]
