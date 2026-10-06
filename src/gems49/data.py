"""Data access for the GEMS Prize rasters.

All files under data/ were obtained from the official competition release
(DrivenData #306 data tab) and verified by sha256 against the pins in
data/bridge/manifest.json (which mirrors the official Dropbox mirrors listed
on the competition data page).  See docs/sources.html for the full list.

Band numbering is 1-based as in the official band descriptions.
"""
from __future__ import annotations

import json
import os

import numpy as np
import rasterio

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FEATURES = os.path.join(ROOT, "data", "training_features.tif")
LABELS = os.path.join(ROOT, "data", "bridge", "existing_faults.tif")
SAMPLE_SUBMISSION = os.path.join(ROOT, "data", "bridge", "example_submission.tif")
LIDAR_SCARP = os.path.join(ROOT, "data", "external", "lidar_scarp_features_u8.tif")

FEATURE_SENTINEL = np.float32(-3.4028234663852886e38)

# 1-based band index -> short name (from training_features.tif descriptions,
# cross-checked against the reference-solution notebook band list).
BANDS = {
    1: "mag_anom",
    2: "rtp",
    3: "tmi_hg",
    4: "geod_2ndinv",
    5: "iso_grav_anom_slope",
    6: "tc",
    7: "geod_shearrate",
    8: "geod_dilaterate",
    9: "tmi_vg",
    10: "deq_n100a15",
    11: "iso_grav_anom_vg",
    12: "det_elev",
    13: "iso_grav_anom",
    14: "tmi",
    15: "depth_to_base_surf",
    16: "ieq_n100a15",
    17: "cond_surf",
    18: "iso_grav_anom_hg",
    19: "det_elev_slope",
}
NAME_TO_BAND = {v: k for k, v in BANDS.items()}


def band(name_or_index) -> np.ndarray:
    """Read one feature band as float32 with nodata -> NaN.

    float32 is deliberate: the sandbox has 4 GB RAM and the grid holds
    12.28 M pixels (49 MB per plane); the features are delivered in
    float32, so no precision is lost.
    """
    if isinstance(name_or_index, str):
        name_or_index = NAME_TO_BAND[name_or_index]
    with rasterio.open(FEATURES) as ds:
        b = ds.read(name_or_index).astype(np.float32)
    b[b <= np.float32(-1e38)] = np.nan  # reference notebook convention
    return b


def labels() -> np.ndarray:
    """Catalogue raster: True where a known fault pixel is present."""
    with rasterio.open(LABELS) as ds:
        lab = ds.read(1)
    return lab == 1


def submission_footprint() -> tuple[np.ndarray, dict]:
    """Finite mask of the official sample_submission (= scoring footprint).

    Returns (mask, geometry) where geometry holds crs/transform/shape.
    """
    with rasterio.open(SAMPLE_SUBMISSION) as ds:
        s = ds.read(1)
        geom = {
            "crs": str(ds.crs),
            "transform": tuple(ds.transform),
            "height": ds.height,
            "width": ds.width,
            "nodata": ds.nodata,
        }
    return np.isfinite(s), geom


def load_lidar_scarp() -> tuple[np.ndarray, dict]:
    """LiDAR scarp matched-filter field (u8 quantised, 0 = nodata).

    From GEMSDOE24 data/external/lidar_scarp_features_u8.tif (built by the
    7GEMSDOE pipeline from GeoDAWN 1 m lidar; band table in the sibling
    lidar_scarp_features.json).
    """
    with rasterio.open(LIDAR_SCARP) as ds:
        a = ds.read(1)
        meta = json.load(open(os.path.join(ROOT, "data", "external", "lidar_scarp_features.json")))
    return a, meta
