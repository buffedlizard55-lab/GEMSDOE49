"""Shared field-construction helpers used by the screening scripts."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter

from .crossgrad import horizontal_gradient, magnitude, robust_unit, smooth


def rank_normalise(field: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Map a field to its within-domain rank scaled to [0,1] (heavy-tail immune)."""
    f = np.where(valid & np.isfinite(field), field, np.nan)
    idx = np.flatnonzero(np.isfinite(f.ravel()))
    vals = f.ravel()[idx]
    order = np.argsort(vals, kind="stable")
    r = np.empty(len(idx), dtype=np.float32)
    r[order] = np.arange(len(idx), dtype=np.float32) / max(len(idx) - 1, 1)
    out = np.full(f.size, np.nan, dtype=np.float32)
    out[idx] = r
    return out.reshape(field.shape)


def multiscale_edge(field: np.ndarray, footprint: np.ndarray, sigmas) -> np.ndarray:
    """Max over scales of |grad G_sigma(field)|, each robustly normalised then combined."""
    acc = np.zeros(field.shape, dtype=np.float32)
    for s in sigmas:
        gx, gy = horizontal_gradient(smooth(field, s))
        m = magnitude(gx, gy)
        ok = footprint & np.isfinite(m)
        n, _ = robust_unit(m, ok, 95.0)
        acc = np.maximum(acc, np.where(ok, n, 0.0).astype(np.float32))
    return np.where(footprint, acc, np.float32(np.nan))


def coherence(field: np.ndarray, footprint: np.ndarray, sigma: float = 2.0,
              win: float = 3.0) -> np.ndarray:
    """Structure-tensor coherence R = |sum w e^{2i theta}| / sum w: how line-like an edge is."""
    gx, gy = horizontal_gradient(smooth(field, sigma))
    ok = footprint & np.isfinite(gx)
    jxx = gaussian_filter(np.where(ok, gx * gx, 0.0), win, mode="nearest")
    jyy = gaussian_filter(np.where(ok, gy * gy, 0.0), win, mode="nearest")
    jxy = gaussian_filter(np.where(ok, gx * gy, 0.0), win, mode="nearest")
    num = np.hypot(jxx - jyy, 2.0 * jxy)
    den = jxx + jyy
    r = np.where(den > 0, num / np.maximum(den, 1e-30), 0.0).astype(np.float32)
    return np.where(ok, r, np.float32(np.nan))
