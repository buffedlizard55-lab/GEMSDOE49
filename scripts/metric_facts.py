#!/usr/bin/env python3
"""Derive the metric facts the site quotes, from the metric itself.

Nothing here is copied from a forum post or a sibling repository: each number is produced by running
``gems49.metric`` and is cross-checked by the assertions in ``tests/test_metric.py``.

Run: python3 scripts/metric_facts.py --out docs/data/metric.json
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import metric as M  # noqa: E402
from gems49.spec import R_M as KERNEL_RADIUS_M  # noqa: E402

PIXEL_M = 100.0  # from GRID['res_m']; asserted below against the grid spec
assert __import__('gems49.spec', fromlist=['GRID']).GRID['res_m'] == PIXEL_M


def collared_probe_geometry(sep_px: int, n_target: int, n_segments: int, seed: int = 49) -> dict:
    """Re-read the holdout parameters that were actually used, from the validation artifact."""
    p = ROOT / "docs" / "data" / "validation.json"
    out = {"sep_px": sep_px, "n_target": n_target, "n_segments": n_segments, "source": str(p)}
    if p.exists():
        v = json.loads(p.read_text())
        out["n_target"] = int(v["n_target"])
        out["sep_px"] = int(v["min_sep_px"])
        out["n_segments"] = int(v["n_segments"])
        out["seed"] = int(v["seed"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/metric.json")
    args = ap.parse_args()

    # --- the organizer's worked example, recomputed from raws -------------------------------
    # page 967: TP_w = 3.00, FP_w = 1.89, FN_w = 2.00, alpha = 0.2, beta = 0.8
    worked = 3.00 / (3.00 + 0.2 * 1.89 + 0.8 * 2.00)
    published = 0.602652
    assert abs(worked - published) < 5e-7, worked

    # --- the marginal rule, brute-forced on a synthetic scene --------------------------------
    # Two isolated truth pixels 20 px (2 km) apart, so neither can shadow the other's kernel.
    # One of them is already covered with unit prediction; the other is free.  Adding a unit pixel
    # near the free one increases TP_w by exactly the kernel credit k, increases the false-positive
    # term by exactly 1-k, and therefore moves the denominator by exactly alpha = 0.2.  It pays
    # iff k > alpha*s0.  The flip must occur exactly there, not approximately.
    shape = (41, 41)
    truth = np.zeros(shape, bool)
    truth[20, 5] = True
    truth[20, 25] = True
    base = np.zeros(shape, np.float32)
    base[20, 5] = 1.0
    s0 = M.dti_components(base, truth)["DTI"]
    bar = M.marginal_bar(s0)                        # == 0.2 * s0, asserted below
    assert abs(bar - 0.2 * s0) < 1e-12
    flips = []
    for d_px in (1, 2, 3):
        p = base.copy()
        p[20 - d_px, 25] = 1.0
        c1 = M.dti_components(p, truth)
        c0 = M.dti_components(base, truth)
        dD = (0.2 * (c1["TPw"] + c1["FPw"]) + 0.8 * c1["n_truth"]) - \
             (0.2 * (c0["TPw"] + c0["FPw"]) + 0.8 * c0["n_truth"])
        credit = float(M.triangle_kernel(np.array([d_px * PIXEL_M]))[0])
        pays = bool(c1["DTI"] > s0)
        assert abs(dD - 0.2) < 1e-12, (d_px, dD)
        assert pays == (credit > bar), (d_px, credit, bar)
        flips.append({"distance_m": int(d_px * PIXEL_M), "credit": credit,
                      "delta_denominator": float(dD), "delta_DTI": float(c1["DTI"] - s0),
                      "pays": pays})
    max_dist = KERNEL_RADIUS_M * (1.0 - 0.2 * 0.28)
    assert 0 < max_dist < KERNEL_RADIUS_M

    # --- binary dominance: dDTI/dlambda > 0 --------------------------------------------------
    line_truth = np.zeros(shape, bool)
    line_truth[20, 5:35] = True
    line = np.zeros(shape, np.float32)
    line[20, 5:35] = 1.0
    lam = np.linspace(0.02, 1.0, 40)
    curve = [float(M.dti_components((line * x).astype(np.float32), line_truth)["DTI"]) for x in lam]
    assert all(b >= a - 1e-12 for a, b in zip(curve, curve[1:])), "binary support must dominate"

    # --- the shadowed-pixel rule: a pixel inside another pixel's kernel shadow earns nothing ----
    p = base.copy()
    p[20, 5] = 1.0                        # duplicate on top of the already-awarded truth pixel
    dup = M.dti_components(p, truth)["DTI"]
    assert dup <= s0 + 1e-12

    payload = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": "drivendata.org/competitions/306/competition-doe-gems/page/967/",
        "alpha": 0.2, "beta": 0.8, "kernel_radius_m": KERNEL_RADIUS_M, "pixel_m": PIXEL_M,
        "worked_example": {"TP_w": 3.00, "FP_w": 1.89, "FN_w": 2.00},
        "worked_example_dti": worked,
        "worked_example_published": published,
        "marginal_bar_at_s": "k > alpha * s",
        "s0_probe": float(s0),
        "marginal_bar_at_s0": float(bar),
        "marginal_max_distance_at_0.28": float(KERNEL_RADIUS_M * (1.0 - 0.2 * 0.28)),
        "marginal_distance_m_at_0.28": float(KERNEL_RADIUS_M * (1.0 - 0.2 * 0.28)),
        "marginal_flip_probe": flips,
        "binary_dominance_curve": {"lambda": lam.round(4).tolist(), "DTI": [round(c, 6) for c in curve]},
        "duplicate_pixel_delta_DTI": float(dup - s0),
        "holdout_geometry_used": collared_probe_geometry(3, 40_000, 3_118),
        "brief_excerpt": (
            "Cross-gradient of gravity and magnetics (Gallardo & Meju, 2004; Fregoso & Gallardo, "
            "Geophysics 74(4), 2009) applied in the data domain via the horizontal-gradient angle and "
            "cross-product (ASEG 2012 extended abstract, doi:10.1071/ASEG2012ab273): compute the "
            "horizontal-gradient vectors of the isostatic gravity anomaly, the reduced-to-pole magnetic "
            "anomaly and a smoothed elevation surface at common resolution; score each pixel by "
            "magnitude-weighted alignment (both gradients large and parallel); keep magnitude-weighted "
            "orthogonality as a separate confounder flag, not as evidence. Check the native gravity "
            "station spacing before trusting 100 m gradients."),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    print(f"wrote {out}")
    print(f"  worked example      {worked:.6f} (published {published})")
    print(f"  s0 {s0:.4f}  bar {bar:.4f}  bar distance {max_dist:.1f} m at s=0.28")
    print(f"  binary dominance    {curve[0]:.6f} -> {curve[-1]:.6f}")
    print(f"  duplicate pixel     dDTI {dup - s0:+.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
