#!/usr/bin/env python3
"""Audit observed raster smoothness without mistaking it for native survey spacing."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe49.raster import _find_band, sha256_file  # noqa: E402

FIELDS = ("iso_grav_anom", "rtp", "det_elev")
LAGS_PX = (1, 2, 3, 5, 10, 20, 30)
STRIDE = 2


def pearson(x: np.ndarray, y: np.ndarray, valid: np.ndarray) -> tuple[float | None, int]:
    xv = np.asarray(x[valid], dtype=np.float64)
    yv = np.asarray(y[valid], dtype=np.float64)
    if xv.size < 2:
        return None, int(xv.size)
    xv -= xv.mean()
    yv -= yv.mean()
    denom = float(np.sqrt(np.dot(xv, xv) * np.dot(yv, yv)))
    if denom == 0.0:
        return None, int(xv.size)
    return float(np.dot(xv, yv) / denom), int(xv.size)


def lag_correlation(array: np.ndarray, footprint: np.ndarray, lag: int, axis: int) -> dict:
    if axis == 1:
        x0 = array[::STRIDE, : array.shape[1] - lag : STRIDE]
        x1 = array[::STRIDE, lag::STRIDE]
        m0 = footprint[::STRIDE, : array.shape[1] - lag : STRIDE]
        m1 = footprint[::STRIDE, lag::STRIDE]
        name = "east_west"
    else:
        x0 = array[: array.shape[0] - lag : STRIDE, ::STRIDE]
        x1 = array[lag::STRIDE, ::STRIDE]
        m0 = footprint[: array.shape[0] - lag : STRIDE, ::STRIDE]
        m1 = footprint[lag::STRIDE, ::STRIDE]
        name = "north_south"
    valid = m0 & m1 & np.isfinite(x0) & np.isfinite(x1)
    value, n = pearson(x0, x1, valid)
    return {"direction": name, "lag_pixels": int(lag), "lag_m": int(lag * 100), "pearson_r": value, "sample_pairs": n}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / ".cache" / "gems_data")
    parser.add_argument("--out", type=Path, default=ROOT / "research" / "input_spacing_audit_20261006.json")
    args = parser.parse_args()
    feature_path = args.data_dir / "training_features.tif"
    sample_path = args.data_dir / "sample_submission.tif"

    with rasterio.open(sample_path) as template:
        footprint = np.isfinite(template.read(1, masked=False))
        grid = {
            "crs": template.crs.to_string() if template.crs else None,
            "width": template.width,
            "height": template.height,
            "pixel_size_m": [float(template.res[0]), float(template.res[1])],
            "footprint_pixels": int(footprint.sum()),
            "sample_sha256": sha256_file(sample_path),
        }
    field_receipts = {}
    band_context = {}
    with rasterio.open(feature_path) as ds:
        if ds.shape != footprint.shape:
            raise ValueError("feature raster does not match the sample grid")
        for name in FIELDS:
            index = _find_band(ds, name)
            band = ds.read(index, masked=False).astype(np.float32, copy=False)
            valid = footprint & np.isfinite(band)
            if ds.nodata is not None and np.isfinite(ds.nodata):
                valid &= band != np.float32(ds.nodata)
            band_receipt = {
                "band_index_1_based": index,
                "description": ds.descriptions[index - 1],
                "tags": ds.tags(index),
                "valid_pixels": int(valid.sum()),
                "nodata_pixels_in_footprint": int(footprint.sum() - valid.sum()),
                "min": float(np.min(band[valid])),
                "p01": float(np.percentile(band[valid], 1)),
                "median": float(np.percentile(band[valid], 50)),
                "p99": float(np.percentile(band[valid], 99)),
                "max": float(np.max(band[valid])),
                "lag_correlations": [],
            }
            for lag in LAGS_PX:
                for axis in (1, 0):
                    band_receipt["lag_correlations"].append(lag_correlation(band, valid, lag, axis))
            field_receipts[name] = band_receipt
        for index in (5, 11, 18, 19):
            band_context[str(index)] = {
                "description": ds.descriptions[index - 1],
                "tags": ds.tags(index),
            }
        feature_hash = sha256_file(feature_path)

    receipt = {
        "audit_id": "GEMSDOE49-input-spacing-audit",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "interpretation_limit": "Raster autocorrelation and observed lag structure measure the resampled competition grid's smoothness, not native station spacing or independent spatial resolution.",
        "provenance": "Pinned owner-hosted public mirror; byte hash establishes mirror consistency only and is not organizer authentication.",
        "feature_file_sha256": feature_hash,
        "sampling": {"systematic_stride_px": STRIDE, "lags_px": list(LAGS_PX)},
        "grid": grid,
        "bands": field_receipts,
        "related_derivative_band_descriptions": band_context,
        "external_source_context": {
            "official_usgs_isostatic_page": "https://pubs.usgs.gov/ds/2006/234/nv_iso.htm",
            "official_usgs_bouguer_page": "https://pubs.usgs.gov/ds/2006/234/nv_boug.htm",
            "official_usgs_data_directory": "https://pubs.usgs.gov/ds/2006/234/data/",
            "verified_statement": "USGS DS 234 says Nevada isostatic gravity was derived from complete Bouguer data; the companion complete-Bouguer page states 71,055 station measurements and a 1-km grid. The competition mirror does not identify the precise source/version of its iso_grav_anom band, so this is a cautionary comparator, not a provenance match.",
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "out": str(args.out),
        "feature_sha256": feature_hash,
        "grid": grid,
        "gravity_lags": [
            x for x in field_receipts["iso_grav_anom"]["lag_correlations"]
            if x["direction"] == "east_west"
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
