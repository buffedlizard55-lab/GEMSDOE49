import numpy as np
import rasterio
from affine import Affine

from gemsdoe49.raster import validate_submission, write_submission


def make_template(path):
    arr = np.full((20, 30), np.nan, dtype=np.float32)
    arr[2:18, 3:27] = 0.0
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=arr.shape[0],
        width=arr.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:32611",
        transform=Affine(100.0, 0.0, 500000.0, 0.0, -100.0, 4500000.0),
        nodata=np.nan,
    ) as dst:
        dst.write(arr, 1)


def test_submission_is_float32_and_nan_only_outside(tmp_path):
    template = tmp_path / "sample.tif"
    output = tmp_path / "candidate.tif"
    make_template(template)
    with rasterio.open(template) as ds:
        footprint = np.isfinite(ds.read(1))
    pred = np.zeros(footprint.shape, dtype=np.float32)
    pred[6, 8] = 1.0
    receipt = write_submission(output, template, footprint, pred)
    assert receipt["status"] == "PASS"
    assert receipt["nan_inside"] == 0
    assert receipt["non_nan_outside"] == 0
    assert receipt["positive_pixels"] == 1
    verified = validate_submission(output, template)
    assert verified["dtype"] == "float32"
    assert verified["crs"] == "EPSG:32611"


def test_writer_rejects_nan_inside_and_out_of_range(tmp_path):
    template = tmp_path / "sample.tif"
    make_template(template)
    with rasterio.open(template) as ds:
        footprint = np.isfinite(ds.read(1))
    bad = np.zeros(footprint.shape, dtype=np.float32)
    bad[3, 3] = np.nan
    with np.testing.assert_raises(ValueError):
        write_submission(tmp_path / "bad_nan.tif", template, footprint, bad)
    bad = np.zeros(footprint.shape, dtype=np.float32)
    bad[3, 3] = 1.2
    with np.testing.assert_raises(ValueError):
        write_submission(tmp_path / "bad_range.tif", template, footprint, bad)
