"""Distance-weighted Tversky metric transcribed from the official competition page."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
from scipy.ndimage import distance_transform_edt

ALPHA = 0.2
BETA = 0.8
RADIUS_PX = 3.0  # 300 m at the official 100 m pixel size
EPS = 1.0e-12


@dataclass(frozen=True)
class MetricComponents:
    tp_weighted: float
    fp_weighted: float
    fn_weighted: float
    truth_pixels: int
    prediction_mass: float
    dti: float

    def to_dict(self) -> dict:
        return asdict(self)


def triangular_kernel(distance_px: np.ndarray | float, radius_px: float = RADIUS_PX):
    return np.maximum(1.0 - np.asarray(distance_px) / radius_px, 0.0)


def _offsets(radius_px: float = RADIUS_PX):
    limit = int(math.ceil(radius_px))
    for dy in range(-limit, limit + 1):
        for dx in range(-limit, limit + 1):
            distance = math.hypot(dy, dx)
            if distance < radius_px:
                yield dy, dx, 1.0 - distance / radius_px


def dti_components(
    prediction: np.ndarray,
    truth: np.ndarray,
    *,
    radius_px: float = RADIUS_PX,
    alpha: float = ALPHA,
    beta: float = BETA,
) -> MetricComponents:
    """Exact metric on an in-memory 2D raster, using the nearest 300 m kernel.

    Predictions are interpreted as p(x) in [0,1]. NaNs are never accepted here;
    callers must crop/mask to a valid evaluation domain first.
    """
    p = np.asarray(prediction, dtype=np.float64)
    g = np.asarray(truth, dtype=bool)
    if p.ndim != 2 or p.shape != g.shape:
        raise ValueError("prediction and truth must be 2D arrays with the same shape")
    if not np.all(np.isfinite(p)):
        raise ValueError("prediction contains NaN/Inf; mask to the scored domain first")
    if np.any(p < 0.0) or np.any(p > 1.0):
        raise ValueError("prediction values must be within [0, 1]")
    if radius_px <= 0:
        raise ValueError("radius_px must be positive")

    n_truth = int(g.sum())
    mass = float(p.sum(dtype=np.float64))
    if n_truth == 0:
        fp = mass
        return MetricComponents(0.0, fp, 0.0, 0, mass, 0.0)

    # For each predicted pixel, the largest truth-kernel value is based on its
    # Euclidean distance to the nearest truth cell.
    dist_to_truth = distance_transform_edt(~g)
    nearest_truth_weight = triangular_kernel(dist_to_truth, radius_px)
    fp = float(np.sum(p * (1.0 - nearest_truth_weight), dtype=np.float64))

    # For each truth pixel, find max_x p(x) k(d(x,g)); offsets outside support
    # are omitted because their kernel weight is exactly zero.
    best_at_truth = np.zeros(p.shape, dtype=np.float64)
    h, w = p.shape
    for dy, dx, weight in _offsets(radius_px):
        tr0 = max(0, -dy)
        tr1 = min(h, h - dy)
        tc0 = max(0, -dx)
        tc1 = min(w, w - dx)
        pr0, pr1 = tr0 + dy, tr1 + dy
        pc0, pc1 = tc0 + dx, tc1 + dx
        candidate = p[pr0:pr1, pc0:pc1] * weight
        target = best_at_truth[tr0:tr1, tc0:tc1]
        np.maximum(target, candidate, out=target)

    best = best_at_truth[g]
    tp = float(best.sum(dtype=np.float64))
    fn = float(np.sum(1.0 - best, dtype=np.float64))
    denom = tp + alpha * fp + beta * fn + EPS
    score = float(tp / denom) if denom > 0 else 0.0
    return MetricComponents(tp, fp, fn, n_truth, mass, score)


def dti_from_components(tp: float, fp: float, fn: float, *, alpha: float = ALPHA, beta: float = BETA) -> float:
    denom = float(tp) + alpha * float(fp) + beta * float(fn) + EPS
    return float(tp / denom) if denom > 0 else 0.0


def credit_bar(current_dti: float, *, alpha: float = ALPHA) -> float:
    """Expected kernel credit needed for a unit-mass addition to improve DTI."""
    return alpha * float(current_dti)
