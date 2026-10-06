"""The submission format gate, including the two live rejection modes.

Mode A and mode B below are the only two ways a GeoTIFF can produce the DrivenData form error
``Predicted values must be in range [0, 1]`` while looking plausible in a GIS viewer.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gems49 import io  # noqa: E402
from gems49.spec import GRID  # noqa: E402
from validate_submission import validate  # noqa: E402

TEMPLATE = ROOT / "data" / "example_submission.tif"
HAVE_DATA = TEMPLATE.exists()


def _write(tmp_path: Path, arr: np.ndarray) -> Path:
    p = tmp_path / "cand.tif"
    with rasterio.open(
        p, "w", driver="GTiff", height=GRID["height"], width=GRID["width"], count=1,
        dtype="float32", crs=GRID["crs"], transform=rasterio.Affine(*GRID["transform"]),
        nodata=np.nan,
    ) as dst:
        dst.write(arr.astype(np.float32), 1)
    return p


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_mode_a_nan_inside_footprint_is_rejected(tmp_path):
    with rasterio.open(TEMPLATE) as s:
        arr = s.read(1).copy()
    fp = np.isfinite(arr)
    ys, xs = np.nonzero(fp)
    arr[ys[0], xs[0]] = np.nan  # one NaN strictly inside the footprint
    p = _write(tmp_path, arr)
    ok, checks, _summary = validate(p, TEMPLATE)
    assert not ok
    bad = [c for c in checks if not c["ok"]]
    assert any("mode A" in c["check"] for c in bad)
    # and a naive whole-array min/max check fails too, which is the live failure
    assert any("naive form check" in c["check"] and not c["ok"] for c in checks)


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_mode_b_out_of_range_is_rejected(tmp_path):
    with rasterio.open(TEMPLATE) as s:
        arr = s.read(1).copy()
    fp = np.isfinite(arr)
    ys, xs = np.nonzero(fp)
    arr[ys[0], xs[0]] = 2.5
    p = _write(tmp_path, arr)
    ok, checks, _ = validate(p, TEMPLATE)
    assert not ok
    assert any("mode B" in c["check"] and not c["ok"] for c in checks)


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_nodata_sentinel_left_in_the_raster_is_rejected(tmp_path):
    with rasterio.open(TEMPLATE) as s:
        arr = s.read(1).copy()
    fp = np.isfinite(arr)
    ys, xs = np.nonzero(fp)
    arr[ys[0], xs[0]] = -3.4028235e38  # the feature stack's sentinel
    p = _write(tmp_path, arr)
    ok, checks, _ = validate(p, TEMPLATE)
    assert not ok


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_valid_file_passes_every_check(tmp_path):
    with rasterio.open(TEMPLATE) as s:
        arr = s.read(1)
    fp = np.isfinite(arr)
    out = np.where(fp, 0.5, np.nan).astype(np.float32)
    p = _write(tmp_path, out)
    ok, checks, summary = validate(p, TEMPLATE)
    assert ok, [c for c in checks if not c["ok"]]
    assert summary["n_positive"] == int(fp.sum())


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_zero_outside_variant_is_immune_to_the_naive_check(tmp_path):
    """The -zeros variant must pass the naive whole-array min/max check that NaN fails."""
    fp = io.read_template_footprint()
    arr = np.where(fp, 0.0, 0.0).astype(np.float32)
    p = _write(tmp_path, arr)
    ok, checks, _ = validate(p, TEMPLATE)
    naive = [c for c in checks if "naive form check" in c["check"]][0]
    assert naive["ok"]
    # but it does NOT echo the template's NaN mask, which is expected and acceptable
    echo = [c for c in checks if "echoes the official NaN mask" in c["check"]][0]
    assert not echo["ok"]


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_writer_refuses_nan_inside_footprint(tmp_path):
    fp = io.read_template_footprint()
    bad = np.where(fp, 0.5, np.nan)
    ys, xs = np.nonzero(fp)
    bad[ys[0], xs[0]] = np.nan
    with pytest.raises(ValueError, match="non-finite values INSIDE"):
        io.write_submission(tmp_path / "x.tif", bad, fp)


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_writer_refuses_out_of_range(tmp_path):
    fp = io.read_template_footprint()
    bad = np.where(fp, 1.5, np.nan)
    with pytest.raises(ValueError, match=r"out of \[0,1\]"):
        io.write_submission(tmp_path / "x.tif", bad, fp)


@pytest.mark.skipif(not HAVE_DATA, reason="official data not placed")
def test_writer_roundtrip_contract(tmp_path):
    fp = io.read_template_footprint()
    support = np.zeros(fp.shape, dtype=np.float32)
    ys, xs = np.nonzero(fp)
    support[ys[::500], xs[::500]] = 1.0
    for nan_outside in (True, False):
        p = io.write_submission(tmp_path / f"r{int(nan_outside)}.tif", support, fp,
                                nan_outside=nan_outside)
        with rasterio.open(p) as s:
            assert s.count == 1 and s.dtypes[0] == "float32"
            assert str(s.crs) == GRID["crs"]
            assert s.res == (100.0, 100.0)
            assert tuple(s.transform)[:6] == GRID["transform"]
            a = s.read(1)
        assert np.array_equal(np.isnan(a), ~fp) if nan_outside else not np.isnan(a).any()
        assert np.nanmax(a) <= 1.0 and np.nanmin(a) >= 0.0
