"""Deterministic, kernel-aware sparse point placement."""
from __future__ import annotations

import numpy as np


def _disk_offsets(radius_px: float) -> tuple[tuple[int, int], ...]:
    r = int(np.ceil(radius_px))
    return tuple(
        (dy, dx)
        for dy in range(-r, r + 1)
        for dx in range(-r, r + 1)
        if (dy * dy + dx * dx) <= radius_px * radius_px + 1.0e-12
    )


def greedy_kernel_suppress(
    score: np.ndarray,
    valid: np.ndarray,
    count: int,
    *,
    min_distance_px: float = 3.0,
) -> tuple[np.ndarray, dict]:
    """Keep highest-evidence points while suppressing redundant centers within R=3 px.

    This is deterministic greedy non-maximum suppression aligned with the metric's
    300 m support. It does not inspect labels. Returned values are binary point
    emissions; they are valid [0,1] probabilities/confidence scores.
    """
    x = np.asarray(score, dtype=np.float32)
    m = np.asarray(valid, dtype=bool)
    if x.ndim != 2 or x.shape != m.shape:
        raise ValueError("score and valid must be matching 2D arrays")
    if count < 0:
        raise ValueError("count must be nonnegative")
    if min_distance_px < 0:
        raise ValueError("min_distance_px must be nonnegative")
    good = m & np.isfinite(x) & (x > 0.0)
    candidates = np.flatnonzero(good)
    output = np.zeros(x.shape, dtype=np.float32)
    if count == 0:
        return output, {"requested": 0, "selected": 0, "candidate_pixels": int(candidates.size), "min_distance_px": min_distance_px}
    if candidates.size < count:
        raise ValueError(f"Only {candidates.size} positive-evidence cells for requested budget {count}")

    flat = x.ravel()
    # Stable sorting means exact ties resolve by row-major pixel index.
    order = candidates[np.argsort(-flat[candidates], kind="mergesort")]
    h, w = x.shape
    suppressed = np.zeros(x.shape, dtype=bool)
    chosen: list[int] = []
    offsets = _disk_offsets(min_distance_px)
    for index in order:
        row, col = divmod(int(index), w)
        if suppressed[row, col]:
            continue
        chosen.append(int(index))
        for dy, dx in offsets:
            rr, cc = row + dy, col + dx
            if 0 <= rr < h and 0 <= cc < w:
                suppressed[rr, cc] = True
        if len(chosen) == count:
            break
    if len(chosen) != count:
        raise ValueError(f"Metric-aware suppression selected {len(chosen)} of {count} requested points")
    output.ravel()[np.asarray(chosen, dtype=np.int64)] = 1.0
    return output, {
        "requested": int(count),
        "selected": int(len(chosen)),
        "candidate_pixels": int(candidates.size),
        "min_distance_px": float(min_distance_px),
        "tie_break": "stable row-major flat index",
    }
