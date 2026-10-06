#!/usr/bin/env python3
"""Offline verification of the bytes that are actually committed and linked from the site.

This is the check that runs unattended in CI, so it is deliberately narrow and cannot be satisfied
by a description: it hashes the file on disk, compares it with the hash recorded in
``docs/data/submission.json``, and re-derives the format contract from the raster itself.

Run: python3 scripts/verify_committed_submission.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
GRID = {"width": 3292, "height": 3730, "crs": "EPSG:32611", "res": (100.0, 100.0),
        "transform": (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    sub = json.loads((DATA / "submission.json").read_text())
    fails: list[str] = []
    paths = {}
    for key in ("filename", "filename_zeros"):
        rel = sub["download_url"] if key == "filename" else sub["download_url_zeros"]
        # download_url is a site-relative path such as downloads/<file>; map it back to docs/
        p = ROOT / "docs" / rel.lstrip("/")
        if not p.exists():
            p = ROOT / rel
        paths[key] = p
        if not p.exists():
            fails.append(f"{key}: {p} does not exist")
    if fails:
        print("\n".join(fails))
        return 1

    nan_p = paths["filename"]
    if sha256(nan_p) != sub["sha256"]:
        fails.append(f"primary artifact sha256 {sha256(nan_p)} != recorded {sub['sha256']}")

    with rasterio.open(nan_p) as s:
        if s.count != 1:
            fails.append(f"band count {s.count} != 1")
        if s.dtypes[0] != "float32":
            fails.append(f"dtype {s.dtypes[0]} != float32")
        if str(s.crs) != GRID["crs"]:
            fails.append(f"crs {s.crs} != {GRID['crs']}")
        if (s.height, s.width) != (GRID["height"], GRID["width"]):
            fails.append(f"shape {s.height}x{s.width} != {GRID['height']}x{GRID['width']}")
        if tuple(round(v, 6) for v in s.res) != GRID["res"]:
            fails.append(f"resolution {s.res} != {GRID['res']}")
        if tuple(round(v, 6) for v in tuple(s.transform)[:6]) != GRID["transform"]:
            fails.append(f"transform {tuple(s.transform)[:6]} != {GRID['transform']}")
        arr = s.read(1)
    footprint = np.isfinite(arr)
    inside = arr[footprint]
    if not np.isfinite(inside).all():
        fails.append("non-finite values inside the footprint")
    if inside.size and (float(inside.min()) < 0.0 or float(inside.max()) > 1.0):
        fails.append(f"values outside [0,1]: min {inside.min()} max {inside.max()}")
    if int(np.count_nonzero(inside)) != int(sub["n_positive"]):
        fails.append(f"positive pixels {int(np.count_nonzero(inside))} != recorded {sub['n_positive']}")

    z_p = paths["filename_zeros"]
    with rasterio.open(z_p) as s:
        z = s.read(1)
    if np.isnan(z).any():
        fails.append("the zero-filled variant contains NaN")
    if z.dtype != np.float32:
        fails.append(f"zero variant dtype {z.dtype} != float32")
    zp = z > 0
    if not np.array_equal(zp, arr > 0):
        fails.append("the two variants disagree about which pixels are positive")
    if (z < 0).any() or (z > 1).any():
        fails.append("zero variant has values outside [0,1]")

    reported = {
        "sha256": sha256(nan_p),
        "n_positive": int(np.count_nonzero(inside)),
        "n_footprint": int(footprint.sum()),
        "footprint_share": float(footprint.mean()),
        "min_inside": float(inside.min()) if inside.size else None,
        "max_inside": float(inside.max()) if inside.size else None,
        "size_bytes": int(nan_p.stat().st_size),
    }
    print(json.dumps(reported, indent=1))
    if fails:
        print("\nFAILED:")
        print("\n".join(f"  - {f}" for f in fails))
        return 1
    print("\nOK: committed artifact matches docs/data/submission.json and the published contract.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
