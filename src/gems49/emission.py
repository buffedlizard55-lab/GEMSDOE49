"""Metric-aware placement and submission writing.

Metric-aware placement: the DTI is a mass budget — every unit of predicted
mass that does not cover truth costs alpha = 0.2, and credit per emitted
pixel saturates at the kernel scale (R = 3 px).  Placement therefore:

1. ranks candidate pixels by the evidence score,
2. applies non-maximum suppression with radius ~R/2 so emitted dots are not
   redundant copies of the same cover (one dot already covers truth within
   3 px),
3. zeroes every catalogue pixel (masked in live scoring; emission there is
   pure cost — staff forum 11516),
4. selects the top-N dots, where N is calibrated on the hide-and-recover
   holdout at the point where the marginal dot's credit equals the
   break-even bar (the family measured that bar at ~0.052-0.055 credit per
   dot; GEMSDOE32 derives 0.2*0.26 = 0.0520 from the published formula).

Submission writing follows the format the portal accepts (the failure mode
"Predicted values must be in range [0, 1]" is caused either by values
outside [0,1] or by a nodata sentinel such as -3.4e38 — GEMSDOE32 measured
both).  We therefore write an ALL-FINITE float32 raster with no nodata tag:
outside the footprint = 0.0 in the '-zeros' file (the safest variant), and
NaN outside the footprint in the '-nan' twin (the competition's own
convention).  Every written file is re-opened from disk and re-checked.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass

import numpy as np
import rasterio
from rasterio.transform import from_origin

from . import data as D


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
