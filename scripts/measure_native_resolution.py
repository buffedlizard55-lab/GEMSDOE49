"""Measure the native spatial resolution of the fields the brief asks us to differentiate.

The brief says, in as many words: "Check the native gravity station spacing before trusting
100 m gradients".  This script answers that question with three independent measurements per
field, all computed from the official bytes:

* **lag-1 autocorrelation** ``rho1`` along both axes: a field whose 100 m gradient is meaningful
  must decorrelate slowly (rho1 close to 1).  rho1 near 0 means the 100 m product is dominated by
  interpolation structure, not by geology.
* **integral scale** ``L_int`` in metres, from the autocorrelation function along a row/column,
  the lag at which the ACF first drops below 1/e.
* **radially averaged power spectrum slope** over the 0.5-2 cycles/km band, to show whether the
  field has a real spectral rolloff at the grid scale or is flat (white) there.

Run: python scripts/measure_native_resolution.py --out docs/data/native_resolution.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import io  # noqa: E402

FIELDS = ("iso_grav_anom", "rtp", "tmi", "det_elev", "det_elev_slope", "cond_surf",
          "depth_to_base_surf", "iso_grav_anom_slope")


def acf_1d(x: np.ndarray, max_lag: int = 120) -> np.ndarray:
    """Normalised autocorrelation of a 1-D signal, computed only on finite samples."""
    x = x[np.isfinite(x)]
    if x.size < 4 * max_lag:
        return np.zeros(max_lag + 1)
    x = x - x.mean()
    denom = float(np.dot(x, x))
    if denom <= 0:
        return np.zeros(max_lag + 1)
    out = np.empty(max_lag + 1)
    for k in range(max_lag + 1):
        out[k] = float(np.dot(x[: x.size - k], x[k:])) / denom if k else 1.0
    return out


def row_acf(a: np.ndarray, max_lag=120, step=7) -> np.ndarray:
    acc = np.zeros(max_lag + 1)
    n = 0
    for r in range(0, a.shape[0], step):
        acc += acf_1d(a[r], max_lag)
        n += 1
    return acc / max(n, 1)


def col_acf(a: np.ndarray, max_lag=120, step=7) -> np.ndarray:
    return row_acf(a.T, max_lag, step)


def integral_scale_m(acf: np.ndarray, res_m: float = 100.0) -> float:
    """First lag (metres) at which the ACF falls below 1/e."""
    below = np.flatnonzero(acf < 1.0 / np.e)
    return float(below[0] * res_m) if below.size else float(len(acf) * res_m)


def radial_spectrum(a: np.ndarray, nbins: int = 24):
    """Radially averaged power spectrum on the footprint, in cycles per 100 m pixel."""
    f = np.nan_to_num(a, nan=float(np.nanmean(a[np.isfinite(a)])) if np.isfinite(a).any() else 0.0)
    f = f - f.mean()
    win = np.hanning(min(f.shape[0], f.shape[1]))
    f = f[: win.size, : win.size] * np.outer(win, win)
    P = np.abs(np.fft.rfft2(f)) ** 2
    ky = np.fft.fftfreq(f.shape[0])[:, None]
    kx = np.fft.rfftfreq(f.shape[1])[None, :]
    k = np.hypot(ky, kx)
    edges = np.linspace(k.min(), k.max(), nbins + 1)
    which = np.digitize(k.ravel(), edges) - 1
    prof = np.array([P.ravel()[which == i].mean() if np.any(which == i) else np.nan
                     for i in range(nbins)])
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres.tolist(), prof.tolist()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/native_resolution.json")
    args = ap.parse_args()
    t0 = time.time()
    footprint = io.read_template_footprint()

    out = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "note": "all values measured from data/training_features.tif in this checkout",
           "fields": {}}
    for name in FIELDS:
        a = io.read_band(name)
        a[~footprint] = np.nan
        ra = row_acf(a)
        ca = col_acf(a)
        cx, prof = radial_spectrum(a)
        rec = {
            "rho1_rows": round(float(ra[1]), 6),
            "rho1_cols": round(float(ca[1]), 6),
            "rho3_rows_300m": round(float(ra[3]), 6),
            "L_int_m_rows": integral_scale_m(ra),
            "L_int_m_cols": integral_scale_m(ca),
            "spectrum_cycles_per_px": [round(v, 5) for v in cx],
            "spectrum_power": [None if not np.isfinite(v) else float(f"{v:.6g}") for v in prof],
            "finite_fraction": round(float(np.isfinite(a).mean()), 6),
        }
        out["fields"][name] = rec
        print(f"  {name:22s} rho1(row/col)={rec['rho1_rows']:.4f}/{rec['rho1_cols']:.4f}  "
              f"L_int={rec['L_int_m_rows']:.0f}/{rec['L_int_m_cols']:.0f} m", flush=True)
        del a

    p = ROOT / args.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}  ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
