"""Audit the placed competition rasters and write data/prepared_manifest.json.

Fail-closed: exits non-zero if any grid, CRS, dtype or label statistic disagrees with the
values measured from the official bytes and pinned in src/gems49/spec.py.

Run:  python scripts/prepare_data.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import io  # noqa: E402
from gems49.spec import (  # noqa: E402
    BANDS,
    CATALOGUE_PIXELS,
    FOOTPRINT_PIXELS,
    GRID,
    OUTSIDE_PIXELS,
    PINS,
    sha256_file,
)

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {msg}")
    if not cond:
        FAILURES.append(msg)


def main() -> int:
    print("== hash pins ==")
    for name, want in PINS.items():
        p = ROOT / "data" / name
        if not p.exists():
            print(f"  SKIP  {name} (absent)")
            continue
        got = sha256_file(p)
        check(got == want, f"{name} sha256 {got[:16]}... == {want[:16]}...")

    print("== grid contract ==")
    meta = io.read_grid_meta()
    check(meta["height"] == GRID["height"] and meta["width"] == GRID["width"],
          f"shape {meta['height']}x{meta['width']} == {GRID['height']}x{GRID['width']}")
    check(meta["crs"] == GRID["crs"], f"CRS {meta['crs']} == {GRID['crs']}")
    check(meta["res"] == (GRID["res_m"], GRID["res_m"]), f"resolution {meta['res']} == 100 m")
    check(meta["transform"] == GRID["transform"], f"geotransform {meta['transform']}")
    check(meta["dtype"] == "float32", f"template dtype {meta['dtype']} == float32")

    print("== template, labels, footprint ==")
    fp = io.read_template_footprint()
    cat = io.read_catalogue()
    check(int(fp.sum()) == FOOTPRINT_PIXELS, f"footprint pixels {int(fp.sum())} == {FOOTPRINT_PIXELS}")
    check(int((~fp).sum()) == OUTSIDE_PIXELS, f"outside pixels {int((~fp).sum())} == {OUTSIDE_PIXELS}")
    check(int(cat.sum()) == CATALOGUE_PIXELS, f"catalogue pixels {int(cat.sum())} == {CATALOGUE_PIXELS}")
    check(bool((cat & ~fp).sum() == 0), "no catalogue pixel lies outside the template footprint")

    print("== feature stack ==")
    band_registry = []
    import rasterio

    with rasterio.open(ROOT / "data" / "training_features.tif") as src:
        check(src.count == len(BANDS), f"band count {src.count} == {len(BANDS)}")
        check(src.crs is not None and str(src.crs) == GRID["crs"], f"features CRS {src.crs}")
        check(src.res == (GRID["res_m"], GRID["res_m"]), f"features resolution {src.res}")
        for idx, name, cat_, _desc in BANDS:
            desc = src.tags(idx).get("description", "")
            band_registry.append({"index": idx, "band": name, "category": cat_, "tag": desc})
            check(bool(desc), f"band {idx:2d} {name:22s} carries a description tag")

    print("== measured label statistics ==")
    stats = {
        "catalogue_pixels": int(cat.sum()),
        "footprint_pixels": int(fp.sum()),
        "catalogue_density_in_footprint": round(float(cat.sum() / fp.sum()), 8),
        "valid_fraction": round(float(fp.sum() / fp.size), 6),
    }
    for k, v in stats.items():
        print(f"    {k} = {v}")

    out = ROOT / "data" / "prepared_manifest.json"
    out.write_text(json.dumps({
        "grid": {k: v for k, v in GRID.items() if k != "transform"} | {"transform": list(GRID["transform"])},
        "pins": PINS,
        "bands": band_registry,
        "stats": stats,
        "failures": FAILURES,
    }, indent=1))
    print(f"wrote {out}")
    if FAILURES:
        print(f"\n{len(FAILURES)} CHECK(S) FAILED")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
