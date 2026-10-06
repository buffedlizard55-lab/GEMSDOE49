"""Raster loading, grid checks, submission writing and independent format receipt."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import rasterio


BAND_NAMES = {
    "rtp": "rtp",
    "det_elev": "det_elev",
    "iso_grav_anom": "iso_grav_anom",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _grid_signature(ds: rasterio.io.DatasetReader) -> dict:
    return {
        "width": int(ds.width),
        "height": int(ds.height),
        "crs": ds.crs.to_string() if ds.crs else None,
        "transform": tuple(float(v) for v in ds.transform)[:6],
        "bounds": tuple(float(v) for v in (ds.bounds.left, ds.bounds.bottom, ds.bounds.right, ds.bounds.top)),
    }


def assert_same_grid(reference: rasterio.io.DatasetReader, other: rasterio.io.DatasetReader) -> None:
    if (reference.width, reference.height) != (other.width, other.height):
        raise ValueError("Raster shape does not match the sample/template")
    if reference.crs != other.crs or reference.transform != other.transform:
        raise ValueError("Raster CRS/geotransform does not match the sample/template")


def _find_band(ds: rasterio.io.DatasetReader, requested: str) -> int:
    matches = []
    for index in range(1, ds.count + 1):
        tags = ds.tags(index)
        tag_name = tags.get("band_name", "").strip()
        description = (ds.descriptions[index - 1] or "").strip()
        if tag_name == requested or description.split(" - ", 1)[0].strip() == requested:
            matches.append(index)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {requested!r} band, found indices {matches}")
    return matches[0]


def load_inputs(data_dir: Path) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, dict]:
    """Load named scalar fields, labels and the sample-footprint mask; fail closed on grids."""
    data_dir = Path(data_dir)
    feature_path = data_dir / "training_features.tif"
    label_path = data_dir / "labels.tif"
    sample_path = data_dir / "sample_submission.tif"
    fields: dict[str, np.ndarray] = {}
    with rasterio.open(sample_path) as template, rasterio.open(feature_path) as features:
        assert_same_grid(template, features)
        if features.count < 1:
            raise ValueError("Feature stack has no bands")
        footprint_template = template.read(1, masked=False)
        footprint = np.isfinite(footprint_template)
        feature_grid = _grid_signature(features)
        band_map: dict[str, dict] = {}
        for name in BAND_NAMES:
            index = _find_band(features, name)
            arr = features.read(index, masked=False).astype(np.float32, copy=False)
            valid = np.isfinite(arr)
            if features.nodata is not None and np.isfinite(features.nodata):
                valid &= arr != np.float32(features.nodata)
            arr = arr.copy()
            arr[~valid] = np.nan
            fields[name] = arr
            band_map[name] = {
                "band_index_1_based": int(index),
                "description": features.descriptions[index - 1],
                "tags": features.tags(index),
                "valid_in_template_footprint": int((valid & footprint).sum()),
                "invalid_in_template_footprint": int((~valid & footprint).sum()),
            }
    with rasterio.open(label_path) as labels_ds, rasterio.open(sample_path) as template:
        assert_same_grid(template, labels_ds)
        raw_labels = labels_ds.read(1, masked=False)
        labels = np.isfinite(raw_labels) & (raw_labels > 0) & footprint
        label_nodata = labels_ds.nodata
    metadata = {
        "feature_grid": feature_grid,
        "footprint_pixels": int(footprint.sum()),
        "outside_footprint_pixels": int((~footprint).sum()),
        "label_positive_pixels": int(labels.sum()),
        "label_nodata": float(label_nodata) if label_nodata is not None else None,
        "band_map": band_map,
        "feature_sha256": sha256_file(feature_path),
        "labels_sha256": sha256_file(label_path),
        "sample_sha256": sha256_file(sample_path),
    }
    return fields, labels, footprint, metadata


def write_submission(
    output_path: Path,
    template_path: Path,
    footprint: np.ndarray,
    prediction: np.ndarray,
    *,
    tags: dict[str, str] | None = None,
) -> dict:
    """Write one float32 GeoTIFF, NaN only outside the official valid footprint."""
    output_path = Path(output_path)
    template_path = Path(template_path)
    foot = np.asarray(footprint, dtype=bool)
    pred = np.asarray(prediction, dtype=np.float32)
    if pred.shape != foot.shape:
        raise ValueError("prediction shape differs from footprint")
    if not np.all(np.isfinite(pred[foot])):
        raise ValueError("NaN/Inf inside the sample footprint is prohibited")
    if np.any(pred[foot] < 0.0) or np.any(pred[foot] > 1.0):
        raise ValueError("in-footprint predictions must be in [0,1]")
    output = pred.copy()
    output[~foot] = np.nan
    with rasterio.open(template_path) as template:
        template_foot = np.isfinite(template.read(1, masked=False))
        if not np.array_equal(template_foot, foot):
            raise ValueError("provided footprint differs from the official sample-template footprint")
        profile = template.profile.copy()
        profile.update(
            driver="GTiff", count=1, dtype="float32", nodata=np.nan,
            compress="deflate", predictor=3, zlevel=6,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(output, 1)
            dst.set_band_description(1, "Fault confidence / probability")
            dst.update_tags(AREA_OR_POINT="Area")
            if tags:
                dst.update_tags(**{str(k): str(v) for k, v in tags.items()})
    return validate_submission(output_path, template_path)


def validate_submission(path: Path, template_path: Path) -> dict:
    """Re-open the written file and check format, NaN placement, values and grid."""
    path = Path(path)
    template_path = Path(template_path)
    with rasterio.open(template_path) as template, rasterio.open(path) as source:
        assert_same_grid(template, source)
        if source.count != 1:
            raise ValueError(f"Expected a single band; found {source.count}")
        if source.dtypes != ("float32",):
            raise ValueError(f"Expected float32; found {source.dtypes}")
        if source.crs is None or source.crs.to_epsg() != 32611:
            raise ValueError(f"Expected EPSG:32611; found {source.crs}")
        if abs(source.transform.a - 100.0) > 1.0e-9 or abs(source.transform.e + 100.0) > 1.0e-9:
            raise ValueError(f"Expected 100 m pixels; got {source.transform}")
        footprint = np.isfinite(template.read(1, masked=False))
        array = source.read(1, masked=False)
        inside = array[footprint]
        outside = array[~footprint]
        if not np.isfinite(inside).all():
            raise ValueError("Submission contains NaN/Inf inside the official footprint")
        if np.any(inside < 0.0) or np.any(inside > 1.0):
            raise ValueError("Submission contains values outside [0,1] inside the footprint")
        if outside.size and not np.isnan(outside).all():
            raise ValueError("Only NaN is permitted outside the official footprint")
        if source.nodata is None or not np.isnan(source.nodata):
            raise ValueError(f"Expected NaN nodata metadata; got {source.nodata!r}")
        if not np.array_equal(source.read_masks(1) > 0, footprint):
            raise ValueError("GeoTIFF validity mask does not match the sample footprint")
        signature = _grid_signature(source)
    return {
        "status": "PASS",
        "file": path.name,
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "bands": 1,
        "dtype": "float32",
        "crs": signature["crs"],
        "shape": [signature["height"], signature["width"]],
        "transform": list(signature["transform"]),
        "bounds": list(signature["bounds"]),
        "nodata": "NaN",
        "footprint_pixels": int(footprint.sum()),
        "outside_pixels": int((~footprint).sum()),
        "nan_inside": 0,
        "non_nan_outside": 0,
        "min_inside": float(np.min(inside)),
        "max_inside": float(np.max(inside)),
        "positive_pixels": int((inside > 0).sum()),
        "unique_values_inside": int(np.unique(inside).size),
        "checks": [
            "single band", "float32", "EPSG:32611", "100 m", "sample shape/geotransform/bounds",
            "in-footprint finite and within [0,1]", "NaN only outside footprint", "NaN nodata tag",
            "dataset validity mask matches sample footprint",
        ],
    }
