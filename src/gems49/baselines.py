"""The single-layer gradient baseline family.

The baseline is deliberately strong: *every* one of the 19 supplied bands is turned into a
single-layer edge field, and the best of the 19 is chosen **on the training folds only**, at
a matched emission budget and with the same metric-aware placement as the candidate.  A
candidate that cannot beat the best single band has not earned a submission slot.
"""

from __future__ import annotations

import numpy as np

from .crossgrad import horizontal_gradient, magnitude, robust_unit, smooth
from .spec import BANDS


def single_layer_edge(field: np.ndarray, footprint: np.ndarray, sigma: float = 2.0) -> np.ndarray:
    """|grad(smoothed field)|, robustly normalised to [0, 1] over the footprint."""
    f = smooth(field, sigma)
    gx, gy = horizontal_gradient(f)
    m = magnitude(gx, gy)
    valid = footprint & np.isfinite(m)
    n, _ = robust_unit(m, valid, 95.0)
    out = np.where(valid, n, np.nan)
    return out


def all_single_layer_edges(
    bands: dict[str, np.ndarray], footprint: np.ndarray, sigma: float = 2.0
) -> dict[str, np.ndarray]:
    return {name: single_layer_edge(bands[name], footprint, sigma) for _i, name, _c, _d in BANDS
            if name in bands}


def curvature(field: np.ndarray, footprint: np.ndarray, sigma: float = 2.0) -> np.ndarray:
    """|Laplacian| of the smoothed field -- a second-order single-layer edge field."""
    from scipy.ndimage import laplace

    f = smooth(field, sigma)
    lap = np.abs(laplace(np.nan_to_num(f, nan=0.0), mode="nearest"))
    valid = footprint & np.isfinite(f)
    n, _ = robust_unit(lap, valid, 95.0)
    return np.where(valid, n, np.nan)
