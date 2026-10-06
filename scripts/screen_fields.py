"""Field-family screen at fixed placement geometry, with no-information controls.

Controls are mandatory: without a matched-geometry random/uniform baseline there is no way to
tell a detector's skill from the geometry of the placement itself.

Run: python -u scripts/screen_fields.py --out docs/data/field_screen.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import emission, holdout, io, methods  # noqa: E402
from gems49.baselines import single_layer_edge  # noqa: E402
from gems49.crossgrad import horizontal_gradient, magnitude, robust_unit, smooth  # noqa: E402
from gems49.spec import GRID  # noqa: E402

N_FOLDS = 4
SEED = 49
SEP = 3
N_TARGET = 60_000
GEO = ((60_000, 3), (90_000, 3), (60_000, 4), (40_000, 3), (90_000, 4), (140_000, 4))


def rank_normalise(field: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Map a field to its within-domain rank scaled to [0,1] (heavy-tail immune)."""
    f = np.where(valid & np.isfinite(field), field, np.nan)
    idx = np.flatnonzero(np.isfinite(f.ravel()))
    vals = f.ravel()[idx]
    order = np.argsort(vals, kind="stable")
    r = np.empty(len(idx), dtype=np.float32)
    r[order] = np.arange(len(idx), dtype=np.float32) / max(len(idx) - 1, 1)
    out = np.full(f.size, np.nan, dtype=np.float32)
    out[idx] = r
    return out.reshape(field.shape)


def multiscale_edge(field: np.ndarray, footprint: np.ndarray, sigmas) -> np.ndarray:
    """Max over scales of |grad G_sigma(field)|, each robustly normalised then combined."""
    acc = np.zeros(field.shape, dtype=np.float32)
    for s in sigmas:
        gx, gy = horizontal_gradient(smooth(field, s))
        m = magnitude(gx, gy)
        ok = footprint & np.isfinite(m)
        n, _ = robust_unit(m, ok, 95.0)
        acc = np.maximum(acc, np.where(ok, n, 0.0).astype(np.float32))
    return np.where(footprint, acc, np.float32(np.nan))


def coherence(field: np.ndarray, footprint: np.ndarray, sigma: float = 2.0, win: float = 3.0) -> np.ndarray:
    """Structure-tensor coherence R = |sum w e^{2i theta}| / sum w -- how line-like the field is."""
    gx, gy = horizontal_gradient(smooth(field, sigma))
    ok = footprint & np.isfinite(gx)
    jxx = gaussian_filter(np.where(ok, gx * gx, 0.0), win, mode="nearest")
    jyy = gaussian_filter(np.where(ok, gy * gy, 0.0), win, mode="nearest")
    jxy = gaussian_filter(np.where(ok, gx * gy, 0.0), win, mode="nearest")
    num = np.hypot(jxx - jyy, 2.0 * jxy)
    den = jxx + jyy
    r = np.where(den > 0, num / np.maximum(den, 1e-30), 0.0).astype(np.float32)
    return np.where(ok, r, np.float32(np.nan))


def build_fields(bands: dict, footprint: np.ndarray) -> dict[str, np.ndarray]:
    g, m = bands["iso_grav_anom"], bands["rtp"]
    s19 = bands["det_elev_slope"]
    out: dict[str, np.ndarray] = {}

    out["single:det_elev_slope"] = single_layer_edge(s19, footprint, 2.0)
    out["topo_ms"] = multiscale_edge(s19, footprint, (1.0, 1.5, 2.5, 4.0))
    coh = coherence(s19, footprint)
    out["topo_ms_coh"] = (np.nan_to_num(out["topo_ms"], nan=0.0) * np.nan_to_num(coh, nan=0.0)).astype(np.float32)

    out["triple"] = methods.family("triple", bands, footprint, sigma=2.0)
    out["tg_g"] = methods.family("tg_g", bands, footprint, sigma=2.0)
    out["tg_m"] = methods.family("tg_m", bands, footprint, sigma=2.0)
    out["triple_s1"] = methods.family("triple", bands, footprint, sigma=1.0)
    out["triple_s4"] = methods.family("triple_s4", bands, footprint, sigma=4.0)
    out["cg_gxm"] = methods.family("cg_gxm", bands, footprint, sigma=2.0)

    # rank-average fusion of the topographic carrier and the coupling operator
    r_t = rank_normalise(out["topo_ms"], footprint)
    r_3 = rank_normalise(out["triple"], footprint)
    out["rankavg_topo_triple"] = (0.5 * np.nan_to_num(r_t, nan=0.0) + 0.5 * np.nan_to_num(r_3, nan=0.0)).astype(np.float32)
    r_g = rank_normalise(out["cg_gxm"], footprint)
    out["rankavg_topo_triple_gxm"] = (
        0.4 * np.nan_to_num(r_t, nan=0.0)
        + 0.4 * np.nan_to_num(r_3, nan=0.0)
        + 0.2 * np.nan_to_num(r_g, nan=0.0)
    ).astype(np.float32)

    # ---- no-information controls --------------------------------------------
    h, w = footprint.shape
    yy, xx = np.mgrid[0:h, 0:w]
    out["control:uniform_lattice4"] = np.where(footprint & (yy % 4 == 0) & (xx % 4 == 0), 1.0, 0.0).astype(np.float32)
    rng = np.random.default_rng(SEED)
    out["control:random"] = np.where(footprint, rng.random(footprint.shape), np.float32(np.nan)).astype(np.float32)
    # constant field: selection then follows raster order (top-left first)
    out["control:raster_order"] = np.where(footprint, 1.0, np.float32(np.nan)).astype(np.float32)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/field_screen.json")
    ap.add_argument("--geo", default=",".join(f"{n}:{s}" for n, s in GEO))
    args = ap.parse_args()

    t0 = time.time()
    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    seg, sizes = holdout.segment_catalogue(catalogue)
    fold_of_seg = holdout.fold_assignment(seg, len(sizes), N_FOLDS, seed=SEED)
    folds = [holdout.build_fold(catalogue, footprint, fold_of_seg, seg, k) for k in range(N_FOLDS)]
    del seg
    scorers_c = [holdout.FoldScorer(f, "truth_collared") for f in folds]
    scorers_u = [holdout.FoldScorer(f, "truth_unrestricted") for f in folds]
    print(f"[{time.time()-t0:.1f}s] folds ready; collared truth per fold "
          f"{[s.n_truth for s in scorers_c]}", flush=True)

    bands = io.read_bands(sorted(set(methods.CORE_BANDS) | {"det_elev_slope"}))
    fields = build_fields(bands, footprint)
    del bands
    print(f"[{time.time()-t0:.1f}s] {len(fields)} fields built", flush=True)

    geo = [tuple(int(v) for v in tok.split(":")) for tok in args.geo.split(",")]
    results: dict[str, dict] = {}
    for name, field in fields.items():
        f0 = np.nan_to_num(field, nan=0.0).astype(np.float32)
        rows = []
        for n_target, sep in geo:
            sup = emission.select_support(f0, footprint, n_target=n_target, min_sep_px=sep)
            cs = [s.score(sup)["DTI"] for s in scorers_c]
            us = [s.score(sup)["DTI"] for s in scorers_u]
            rows.append({"n_target": n_target, "min_sep_px": sep, "n_emitted": int(sup.sum()),
                         "collared_mean": float(np.mean(cs)),
                         "unrestricted_mean": float(np.mean(us)),
                         "collared_folds": [round(float(v), 6) for v in cs]})
            del sup
        best = max(rows, key=lambda r: r["collared_mean"])
        results[name] = {"rows": rows, "best": best}
        print(f"  {name:28s} best collared={best['collared_mean']:.4f} "
              f"(n={best['n_target']} sep={best['min_sep_px']})  "
              f"unres@best={best['unrestricted_mean']:.4f}   [{time.time()-t0:.0f}s]", flush=True)
        del field, f0

    payload = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": SEED, "n_folds": N_FOLDS, "geometry_grid": [list(g) for g in geo],
        "results": results, "runtime_s": round(time.time() - t0, 1),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    print(f"wrote {out} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
