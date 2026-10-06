"""Hard format gate for a GEMS submission GeoTIFF.

Every check corresponds to a sentence in the published submission format:

  "same projected coordinate reference system as the training data (EPSG 32611)"
  "same resolution as the training data (100m)"
  "same bounds as the training data, and data outside the bounds is null or nan"
  "a single layer with datatype of 32-bit float (float32) with values between 0 and 1"

It also reproduces the two failure modes behind the DrivenData form rejection
``Predicted values must be in range [0, 1]`` as *constructed test cases*, so the gate can be
shown to catch them (see tests/test_gate.py):

  mode A  a single NaN INSIDE the footprint  -> np.nanmin/np.nanmax or a min/max check fails
  mode B  a finite value outside [0, 1]      -> a plain min/max check fails

Run:  python scripts/validate_submission.py docs/downloads/<file>.tif
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49.spec import GRID  # noqa: E402


def validate(path: Path, template: Path) -> tuple[bool, list[dict], dict]:
    checks: list[dict] = []

    def add(name: str, ok: bool, detail: str, required: bool = True) -> None:
        checks.append({"check": name, "ok": bool(ok), "detail": detail, "required": required})

    with rasterio.open(template) as t:
        tmpl = t.read(1)
        tmpl_meta = {"shape": (t.height, t.width), "crs": str(t.crs), "res": t.res,
                     "transform": tuple(t.transform)[:6]}

    try:
        src = rasterio.open(path)
    except Exception as exc:  # noqa: BLE001
        return False, [{"check": "openable", "ok": False, "detail": str(exc)}], {}

    with src:
        arr = src.read(1)
        meta = {"driver": src.driver, "count": src.count, "dtype": src.dtypes[0],
                "shape": (src.height, src.width), "crs": str(src.crs), "res": src.res,
                "transform": tuple(src.transform)[:6], "nodata": src.nodata}

        add("driver is GTiff", src.driver == "GTiff", src.driver)
        add("single band", src.count == 1, f"count={src.count}")
        add("dtype float32", src.dtypes[0] == "float32", src.dtypes[0])
        add("CRS EPSG:32611", str(src.crs) == GRID["crs"], str(src.crs))
        add("resolution 100 m", (src.res[0], src.res[1]) == (100.0, 100.0), str(src.res))
        add("shape 3730x3292", (src.height, src.width) == (GRID["height"], GRID["width"]),
            f"{src.height}x{src.width}")
        add("geotransform matches template", tuple(src.transform)[:6] == tmpl_meta["transform"],
            str(tuple(src.transform)[:6]))
        add("bounds match template", tuple(src.bounds) == tuple(rasterio.open(template).bounds),
            str(tuple(round(b, 1) for b in src.bounds)))

        inside = np.isfinite(tmpl)
        out_mask = ~inside
        finite_inside = np.isfinite(arr[inside])
        add("no NaN/Inf INSIDE the footprint (rejection mode A)",
            bool(finite_inside.all()),
            f"{int((~finite_inside).sum())} bad pixels of {int(inside.sum())}")

        if finite_inside.any():
            vmin = float(arr[inside][finite_inside].min())
            vmax = float(arr[inside][finite_inside].max())
        else:
            vmin = vmax = float("nan")
        add("values inside footprint within [0, 1] (rejection mode B)",
            bool(finite_inside.any() and vmin >= -1e-6 and vmax <= 1.0 + 1e-6),
            f"min={vmin:.6f} max={vmax:.6f}")

        # whole-array min/max, exactly what a naive submission-form check computes
        raw_min = float(np.min(arr))
        raw_max = float(np.max(arr))
        # INFORMATIONAL ONLY: a raw min/max over the whole array must not be a hard gate,
        # because NaN outside the footprint is explicitly allowed by the submission format
        # (the organizer's own example_submission.tif is NaN there).
        add("whole-array min/max are finite and in [0,1] (informational: a naive form check)",
            bool(np.isfinite(raw_min) and np.isfinite(raw_max) and raw_min >= -1e-6 and raw_max <= 1.0 + 1e-6),
            f"min={raw_min} max={raw_max}. A NaN outside the footprint fails this naive check; use the "
            f"0.0-outside variant if a form rejects the file on it.",
            required=False)

        add("no positive values outside the footprint",
            bool(np.nan_to_num(arr[out_mask], nan=0.0).max(initial=0.0) <= 0.0),
            f"max outside = {float(np.nan_to_num(arr[out_mask], nan=0.0).max(initial=0.0))}")

        n_pos = int((np.nan_to_num(arr, nan=0.0) > 0).sum())
        uniq = np.unique(arr[np.isfinite(arr)])
        add("echoes the official NaN mask exactly (informational)",
            bool(np.array_equal(np.isnan(arr), np.isnan(tmpl))),
            f"nan pixels {int(np.isnan(arr).sum())} vs template {int(np.isnan(tmpl).sum())}",
            required=False)

    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {
        "path": str(path), "bytes": path.stat().st_size, "sha256": sha,
        "n_positive": n_pos, "unique_finite_values": [float(v) for v in uniq[:12]],
        "n_unique_finite_values": int(uniq.size),
        "fraction_of_footprint_emitted": round(n_pos / float(int(inside.sum())), 6),
        "meta": meta,
    }
    return all(c["ok"] for c in checks if c.get("required", True)), checks, summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("tif")
    ap.add_argument("--template", default=str(ROOT / "data" / "example_submission.tif"))
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    ok, checks, summary = validate(Path(args.tif), Path(args.template))
    for c in checks:
        mark = "PASS" if c["ok"] else ("INFO" if not c.get("required", True) else "FAIL")
        print(f"  {mark}  {c['check']}: {c['detail']}")
    print()
    print(json.dumps(summary, indent=1))
    if args.json:
        Path(args.json).write_text(json.dumps({"ok": ok, "checks": checks, "summary": summary}, indent=1))
    print("\n" + ("ALL CHECKS PASSED" if ok else "SUBMISSION REJECTED BY LOCAL GATE"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
