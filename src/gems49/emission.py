
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

import hashlib
import os

import numpy as np
import rasterio
from rasterio.transform import from_origin

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


# ---------------------------------------------------------------------------
# Compatibility layer: placement and writing helpers from the earlier (merged)
# session, kept so its scripts and receipts remain reproducible.
# ---------------------------------------------------------------------------

def nms_topk(score: np.ndarray, catalogue: np.ndarray, footprint: np.ndarray,
             n_dots: int, nms_radius_px: float = 2.0) -> np.ndarray:
    """Greedy NMS selection: repeatedly take the highest remaining score,
    suppress a disc of radius nms_radius_px around it.

    Returns a boolean raster of emitted pixels (exactly n_dots True if the
    candidate pool is large enough).
    """
    s = np.where(footprint & ~catalogue & np.isfinite(score), score, -np.inf)
    flat = s.ravel()
    order = np.argpartition(-flat, min(n_dots * 8, flat.size - 1))[: n_dots * 8]
    order = order[np.argsort(-flat[order])]
    h, w = score.shape
    r = int(np.ceil(nms_radius_px))
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    disc = (xx * xx + yy * yy) <= nms_radius_px ** 2
    emitted = np.zeros(h * w, dtype=bool)   # the selected dots
    dead = np.zeros(h * w, dtype=bool)      # suppressed neighbourhoods
    picked = 0
    for idx in order:
        if flat[idx] == -np.inf:
            break
        if dead[idx]:
            continue
        emitted[idx] = True
        picked += 1
        if picked >= n_dots:
            break
        rc, cc = divmod(int(idx), w)
        rr = np.clip(rc + yy, 0, h - 1)
        ccc = np.clip(cc + xx, 0, w - 1)
        dead[(rr * w + ccc)[disc]] = True
    out = emitted.reshape(h, w)
    assert int(out.sum()) == picked
    return out


def dots_to_prediction(chosen: np.ndarray, footprint: np.ndarray,
                       value: float = 1.0) -> np.ndarray:
    p = np.zeros(footprint.shape, dtype=np.float32)
    p[chosen] = value
    return p


def write_submission(pred: np.ndarray, footprint: np.ndarray, out_path: str,
                     mode: str) -> dict:
    """Write the submission GeoTIFF.  mode='zeros': 0 outside footprint,
    no nodata tag.  mode='nan': NaN outside footprint, no nodata tag."""
    geom = D.submission_footprint()[1]
    transform = from_origin(geom["transform"][2], geom["transform"][5],
                            geom["transform"][0], abs(geom["transform"][4]))
    arr = pred.astype(np.float32).copy()
    if mode == "zeros":
        arr[~footprint] = 0.0
        if not np.isfinite(arr).all():
            raise ValueError("zeros mode requires all-finite predictions inside footprint")
        arr = np.nan_to_num(arr, nan=0.0)
    elif mode == "nan":
        arr[~footprint] = np.nan
    else:
        raise ValueError("mode must be 'zeros' or 'nan'")
    lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if lo < 0.0 or hi > 1.0:
        raise ValueError(f"values outside [0,1]: min={lo}, max={hi}")
    profile = dict(
        driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
        dtype="float32", crs=geom["crs"], transform=transform,
        compress="lzw", tiled=True, blockxsize=256, blockysize=256,
        BIGTIFF="IF_SAFER",
    )  # deliberately no 'nodata' key -> no nodata tag
    with rasterio.open(out_path, "w", **profile) as ds:
        ds.write(arr, 1)
    return verify_submission(out_path, footprint)


def verify_submission(path: str, footprint: np.ndarray | None = None) -> dict:
    """Independent re-read of a written file; mirrors the competition's
    validator requirements (problem description, 'Submission format')."""
    with rasterio.open(path) as ds:
        arr = ds.read(1)
        checks = {
            "single_band": ds.count == 1,
            "dtype_float32": ds.dtypes[0] == "float32",
            "dimensions_3730x3292": (ds.height, ds.width) == (3730, 3292),
            "crs_epsg_32611": ds.crs is not None and ds.crs.to_epsg() == 32611,
            "resolution_100m": ds.transform[0] == 100.0 and ds.transform[4] == -100.0,
            "origin_matches_template": (round(ds.transform[2], 6), round(ds.transform[5], 6)) == (243350.0, 4508550.0),
            "nodata_tag_absent_or_nan": ds.nodata is None or np.isnan(ds.nodata),
        }
    if footprint is None:
        footprint = D.submission_footprint()[0]
    inside = arr[footprint]
    outside = arr[~footprint]
    checks.update({
        "in_footprint_all_finite": bool(np.isfinite(inside).all()),
        "in_footprint_range_0_1": bool((inside.min() >= 0.0) and (inside.max() <= 1.0)),
        "outside_footprint_only_0_or_nan": bool(np.all((outside == 0) | np.isnan(outside))),
        "no_sentinel_values": bool(not np.any(np.abs(arr[np.isfinite(arr)]) > 1e30)),
    })
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    return {
        "path": path,
        "sha256": sha,
        "size_bytes": os.path.getsize(path),
        "emitted_positive_pixels": int((inside > 0).sum()),
        "in_footprint_min": float(inside.min()),
        "in_footprint_max": float(inside.max()),
        "checks": checks,
        "all_checks_passed": bool(np.all(list(checks.values()))),
    }
