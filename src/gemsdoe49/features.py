"""Scale-aware horizontal-gradient and cross-gradient features.

The score is a weighted sum of pairwise alignment operators, not a product of
separately trained layer scores. Orthogonality is returned as a diagnostic only.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.ndimage import gaussian_filter


@dataclass(frozen=True)
class GradientSet:
    gx: np.ndarray
    gy: np.ndarray
    magnitude: np.ndarray
    normalized_magnitude: np.ndarray
    p99_magnitude: float


@dataclass(frozen=True)
class CrossGradientResult:
    score: np.ndarray
    orthogonality_flag: np.ndarray
    pair_support: np.ndarray
    single_layer_baselines: dict[str, np.ndarray]
    gradients: dict[str, GradientSet]
    valid_support: np.ndarray
    diagnostics: dict[str, float]


def _robust_scale(magnitude: np.ndarray, valid: np.ndarray, q: float = 99.0) -> tuple[np.ndarray, float]:
    if not np.isfinite(q) or not 0.0 < q <= 100.0:
        raise ValueError("robust percentile must be in (0, 100]")
    values = np.asarray(magnitude[valid], dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.zeros(magnitude.shape, dtype=np.float32), 0.0
    scale = float(np.percentile(values, q))
    if not np.isfinite(scale) or scale <= np.finfo(np.float32).eps:
        return np.zeros(magnitude.shape, dtype=np.float32), 0.0
    norm = np.clip(magnitude / scale, 0.0, 1.0).astype(np.float32, copy=False)
    norm[~valid] = 0.0
    return norm, scale


def _masked_gaussian(
    field: np.ndarray,
    valid: np.ndarray,
    sigma_px: float,
    *,
    truncate: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalized Gaussian convolution that never mixes nodata sentinels into values."""
    x = np.asarray(field, dtype=np.float32)
    m = np.asarray(valid, dtype=bool) & np.isfinite(x)
    if x.ndim != 2 or x.shape != m.shape:
        raise ValueError("field and valid mask must have identical 2D shapes")
    if not np.isfinite(sigma_px) or sigma_px <= 0:
        raise ValueError("sigma_px must be finite and positive")
    weights = m.astype(np.float32)
    denom = gaussian_filter(weights, sigma=sigma_px, mode="constant", cval=0.0, truncate=truncate)
    numer = gaussian_filter(np.where(m, x, 0.0), sigma=sigma_px, mode="constant", cval=0.0, truncate=truncate)
    smooth = np.zeros_like(numer, dtype=np.float32)
    np.divide(numer, denom, out=smooth, where=denom > 1.0e-6)
    # Values are normalized wherever the stencil has nonzero support. The >=99%
    # mask is applied to the score, not to the smoothed field before differencing;
    # zeroing unsupported cells first would create artificial edge gradients.
    support = denom >= 0.99
    return smooth, support


def horizontal_gradient(
    field: np.ndarray,
    valid: np.ndarray,
    *,
    sigma_px: float,
    pixel_size_m: float = 100.0,
    robust_percentile: float = 99.0,
) -> tuple[GradientSet, np.ndarray]:
    """Return central-difference gradient (east/north axes), robust amplitude and support.

    Raster row coordinates increase southward; that sign is consistent for every
    input field, so pairwise orientation comparisons are unaffected.
    """
    if not np.isfinite(pixel_size_m) or pixel_size_m <= 0:
        raise ValueError("pixel_size_m must be finite and positive")
    smooth, support = _masked_gaussian(field, valid, sigma_px)
    gy, gx = np.gradient(smooth, pixel_size_m, pixel_size_m, edge_order=1)
    gx = np.asarray(gx, dtype=np.float32)
    gy = np.asarray(gy, dtype=np.float32)
    magnitude = np.hypot(gx, gy).astype(np.float32)
    normalized, scale = _robust_scale(magnitude, support, robust_percentile)
    return GradientSet(gx, gy, magnitude, normalized, scale), support


def cross_gradient_alignment(
    fields: Mapping[str, np.ndarray],
    valid: np.ndarray,
    *,
    sigma_px: float = 5.0,
    pixel_size_m: float = 100.0,
    robust_percentile: float = 99.0,
) -> CrossGradientResult:
    """Compute magnitude-weighted, sign-invariant structural alignment.

    For each field pair i,j, the cosine term is |grad_i dot grad_j| / (|grad_i||grad_j|)
    and the cross-gradient term is |grad_i x grad_j| / (|grad_i||grad_j|). Both are
    computed directly from the vectors. Pair weight is sqrt(a_i*a_j), with a_i the
    robustly scaled gradient magnitude. The evidence score is the arithmetic mean of
    weighted alignment over the three pairs. The orthogonality flag uses the same
    weights and |sin(theta)|, but is deliberately not included in score.
    """
    if len(fields) != 3:
        raise ValueError("This preregistered H49-XG operator requires exactly three fields")
    if not np.isfinite(sigma_px) or sigma_px <= 0:
        raise ValueError("sigma_px must be finite and positive")
    if not np.isfinite(pixel_size_m) or pixel_size_m <= 0:
        raise ValueError("pixel_size_m must be finite and positive")
    if not np.isfinite(robust_percentile) or not 0.0 < robust_percentile <= 100.0:
        raise ValueError("robust_percentile must be in (0, 100]")
    keys = tuple(fields.keys())
    shape = np.asarray(fields[keys[0]]).shape
    if len(shape) != 2:
        raise ValueError("input fields must be 2D arrays")
    common = np.asarray(valid, dtype=bool).copy()
    if common.shape != shape:
        raise ValueError("valid mask shape does not match input fields")
    cleaned: dict[str, np.ndarray] = {}
    for name in keys:
        arr = np.asarray(fields[name], dtype=np.float32)
        if arr.shape != shape:
            raise ValueError(f"Field {name!r} has a different shape")
        common &= np.isfinite(arr)
        cleaned[name] = arr

    gradients: dict[str, GradientSet] = {}
    supports: list[np.ndarray] = []
    for name in keys:
        grad, support = horizontal_gradient(
            cleaned[name], common,
            sigma_px=sigma_px,
            pixel_size_m=pixel_size_m,
            robust_percentile=robust_percentile,
        )
        gradients[name] = grad
        supports.append(support)
    support = common.copy()
    for item in supports:
        support &= item

    score = np.zeros(shape, dtype=np.float32)
    orth = np.zeros(shape, dtype=np.float32)
    pair_support = np.zeros(shape, dtype=np.float32)
    pair_count = 0
    eps = np.float32(1.0e-12)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a = gradients[keys[i]]
            b = gradients[keys[j]]
            denom = np.maximum(a.magnitude * b.magnitude, eps)
            dot = np.abs(a.gx * b.gx + a.gy * b.gy) / denom
            cross = np.abs(a.gx * b.gy - a.gy * b.gx) / denom
            cos_abs = np.clip(dot, 0.0, 1.0)
            sin_abs = np.clip(cross, 0.0, 1.0)
            weight = np.sqrt(a.normalized_magnitude * b.normalized_magnitude)
            score += (weight * cos_abs).astype(np.float32)
            orth += (weight * sin_abs).astype(np.float32)
            pair_support += weight.astype(np.float32)
            pair_count += 1

    score /= float(pair_count)
    orth /= float(pair_count)
    pair_support /= float(pair_count)
    score[~support] = 0.0
    orth[~support] = 0.0
    pair_support[~support] = 0.0
    score = np.clip(score, 0.0, 1.0).astype(np.float32)
    orth = np.clip(orth, 0.0, 1.0).astype(np.float32)
    pair_support = np.clip(pair_support, 0.0, 1.0).astype(np.float32)
    baselines = {name: gradients[name].normalized_magnitude.copy() for name in keys}
    for baseline in baselines.values():
        baseline[~support] = 0.0

    diagnostics = {
        "sigma_px": float(sigma_px),
        "sigma_m": float(sigma_px * pixel_size_m),
        "fwhm_m": float(sigma_px * pixel_size_m * 2.354820045),
        "valid_support_pixels": int(support.sum()),
        "score_p50": float(np.median(score[support])) if support.any() else 0.0,
        "score_p90": float(np.percentile(score[support], 90)) if support.any() else 0.0,
        "score_p99": float(np.percentile(score[support], 99)) if support.any() else 0.0,
        "orthogonality_p99": float(np.percentile(orth[support], 99)) if support.any() else 0.0,
    }
    return CrossGradientResult(score, orth, pair_support, baselines, gradients, support, diagnostics)


def shifted_alignment_control(
    result: CrossGradientResult,
    *,
    shifted_field: str = "rtp",
    shift_px: tuple[int, int] = (0, 10),
) -> tuple[np.ndarray, np.ndarray]:
    """Build a spatial-shift null that preserves each gradient field's texture.

    One field's vector components and amplitudes are shifted together, destroying
    only local registration. Wrapped edge pixels are excluded, never scored.
    Returns (alignment score, valid support).
    """
    if shifted_field not in result.gradients:
        raise KeyError(f"Unknown field for shifted control: {shifted_field}")
    dy, dx = (int(shift_px[0]), int(shift_px[1]))
    if (dy, dx) == (0, 0):
        raise ValueError("shift_px must be nonzero")
    h, w = result.valid_support.shape
    if abs(dy) >= h or abs(dx) >= w:
        raise ValueError("shift is larger than the raster")

    grads = dict(result.gradients)
    original = grads[shifted_field]
    shifted = GradientSet(
        gx=np.roll(original.gx, shift=(dy, dx), axis=(0, 1)),
        gy=np.roll(original.gy, shift=(dy, dx), axis=(0, 1)),
        magnitude=np.roll(original.magnitude, shift=(dy, dx), axis=(0, 1)),
        normalized_magnitude=np.roll(original.normalized_magnitude, shift=(dy, dx), axis=(0, 1)),
        p99_magnitude=original.p99_magnitude,
    )
    grads[shifted_field] = shifted
    support = result.valid_support & np.roll(result.valid_support, shift=(dy, dx), axis=(0, 1))
    if dy > 0:
        support[:dy, :] = False
    elif dy < 0:
        support[dy:, :] = False
    if dx > 0:
        support[:, :dx] = False
    elif dx < 0:
        support[:, dx:] = False

    score = np.zeros(result.score.shape, dtype=np.float32)
    eps = np.float32(1.0e-12)
    pair_count = 0
    names = tuple(grads)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a = grads[names[i]]
            b = grads[names[j]]
            denom = np.maximum(a.magnitude * b.magnitude, eps)
            dot = np.abs(a.gx * b.gx + a.gy * b.gy) / denom
            weight = np.sqrt(a.normalized_magnitude * b.normalized_magnitude)
            score += (weight * np.clip(dot, 0.0, 1.0)).astype(np.float32)
            pair_count += 1
    score /= float(pair_count)
    score[~support] = 0.0
    return np.clip(score, 0.0, 1.0).astype(np.float32), support


def normalize_rank01(score: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Map a finite score field to robust [0,1] ranks; invalid/outside cells are zero."""
    x = np.asarray(score, dtype=np.float32)
    m = np.asarray(valid, bool) & np.isfinite(x)
    out = np.zeros(x.shape, dtype=np.float32)
    vals = x[m]
    if vals.size == 0:
        return out
    # Empirical rank normalisation is invariant to arbitrary input units.
    order = np.argsort(vals, kind="mergesort")
    ranks = np.empty(order.size, dtype=np.float32)
    if order.size == 1:
        ranks[order] = 1.0
    else:
        ranks[order] = np.linspace(0.0, 1.0, order.size, dtype=np.float32)
    out[m] = ranks
    return out
