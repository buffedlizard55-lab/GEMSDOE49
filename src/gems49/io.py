"""Raster input/output for the GEMS grid.

Reading is fail-closed: every band is returned as float32 with the feature stack's
most-negative-float32 sentinel converted to NaN, and the official footprint is always
recomputed from the example-submission template rather than inferred.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio

from .spec import BAND_INDEX, GRID, PINS, sha256_file

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


def _path(name: str) -> Path:
    return DATA / name


def check_pins(names: list[str] | None = None) -> dict[str, bool]:
    """Re-hash pinned data files and report which match."""
    names = names or list(PINS)
    out = {}
    for n in names:
        p = _path(n)
        out[n] = p.exists() and sha256_file(p) == PINS[n]
    return out


def read_template_footprint() -> np.ndarray:
    """Boolean mask of pixels that are NOT NaN in the official example submission."""
    with rasterio.open(_path("example_submission.tif")) as src:
        arr = src.read(1)
    return np.isfinite(arr)


def read_catalogue() -> np.ndarray:
    """Boolean mask of the known (USGS Quaternary + INGENIOUS) fault pixels."""
    with rasterio.open(_path("existing_faults.tif")) as src:
        arr = src.read(1)
    return arr == 1


def band_names() -> list[str]:
    return list(BAND_INDEX)


def read_band(name: str) -> np.ndarray:
    """Read one named band of training_features.tif as float32 with sentinel -> NaN."""
    if name not in BAND_INDEX:
        raise KeyError(f"unknown band {name!r}; known: {sorted(BAND_INDEX)}")
    idx = BAND_INDEX[name]
    with rasterio.open(_path("training_features.tif")) as src:
        arr = src.read(idx).astype(np.float32)
    arr[arr <= GRID["nodata_sentinel"] / 2.0] = np.nan
    return arr


def read_bands(names: list[str]) -> dict[str, np.ndarray]:
    """Read several bands in one pass over the file."""
    wanted = {BAND_INDEX[n]: n for n in names}
    out: dict[str, np.ndarray] = {}
    with rasterio.open(_path("training_features.tif")) as src:
        for idx, nm in sorted(wanted.items()):
            arr = src.read(idx).astype(np.float32)
            arr[arr <= GRID["nodata_sentinel"] / 2.0] = np.nan
            out[nm] = arr
    return out


def read_aux(name: str, band: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Read an auxiliary uint8 bridge raster; returns (values, valid) with 0 == nodata."""
    with rasterio.open(_path(f"external/{name}")) as src:
        arr = src.read(band)
    ok = arr > 0
    return arr.astype(np.float32), ok


def read_grid_meta() -> dict:
    with rasterio.open(_path("example_submission.tif")) as src:
        return {
            "height": src.height,
            "width": src.width,
            "crs": str(src.crs),
            "res": src.res,
            "transform": tuple(src.transform)[:6],
            "dtype": src.dtypes[0],
        }


def write_submission(
    path: str | Path,
    values: np.ndarray,
    footprint: np.ndarray,
    *,
    nan_outside: bool = True,
) -> Path:
    """Write the required single-band float32 GeoTIFF.

    Contract enforced here (all four are re-checked by scripts/validate_submission.py):

    * EPSG:32611, 100 m, 3730 x 3292, geotransform identical to the organizer template;
    * single band, dtype float32;
    * every FINITE value inside the footprint lies in [0, 1] (no NaN, no inf inside);
    * outside the footprint: NaN when ``nan_outside`` (byte-identical mask to the official
      template), otherwise 0.0.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.shape != footprint.shape:
        raise ValueError(f"shape mismatch {values.shape} vs footprint {footprint.shape}")
    inside = footprint
    finite_inside = np.isfinite(values[inside])
    if not finite_inside.all():
        raise ValueError(f"{int((~finite_inside).sum())} non-finite values INSIDE the footprint")
    vmin = float(values[inside].min())
    vmax = float(values[inside].max())
    if vmin < 0.0 - 1e-6 or vmax > 1.0 + 1e-6:
        raise ValueError(f"values inside footprint out of [0,1]: min={vmin} max={vmax}")

    out = np.where(inside, np.clip(values, 0.0, 1.0), np.nan if nan_outside else 0.0).astype(np.float32)
    if not nan_outside:
        # belt and braces: nothing outside the footprint is ever anything but zero
        out[~inside] = 0.0

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=GRID["height"],
        width=GRID["width"],
        count=1,
        dtype="float32",
        crs=GRID["crs"],
        transform=rasterio.Affine(*GRID["transform"]),
        nodata=np.nan if nan_outside else None,
        compress="deflate",
        predictor=3,
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        dst.write(out, 1)
    return path
