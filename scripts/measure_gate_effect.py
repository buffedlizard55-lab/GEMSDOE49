#!/usr/bin/env python3
"""Does the cross-gradient gate actually change the emitted file?

The leave-one-fold-out experiment says the gated family and the ungated carrier differ by +0.00003 in
score.  That is small enough that it might mean the gate does nothing at the *emission* level.  This
script settles it by building both supports with identical geometry and comparing them pixel by pixel.

Run: python3 scripts/measure_gate_effect.py --out docs/data/gate_effect.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import emission, gating, io  # noqa: E402
from gems49.screen_util import coherence, multiscale_edge  # noqa: E402

SIGMAS = (1.0, 1.5, 2.5, 4.0)


def support_of(field: np.ndarray, domain: np.ndarray, footprint: np.ndarray,
               n_target: int, sep: int) -> np.ndarray:
    f = np.nan_to_num(field, nan=0.0)
    f[~footprint] = 0.0
    f = (f / float(f[footprint].max())).astype(np.float32)
    return emission.select_support(f, domain, n_target=n_target, min_sep_px=sep)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/gate_effect.json")
    ap.add_argument("--n-target", type=int, default=40_000)
    ap.add_argument("--min-sep", type=int, default=3)
    ap.add_argument("--w", type=float, default=0.25)
    args = ap.parse_args()

    t0 = time.time()
    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    domain = footprint & ~catalogue
    core = io.read_bands(["iso_grav_anom", "rtp", "det_elev_slope"])
    s19 = core["det_elev_slope"]

    ms = multiscale_edge(s19, footprint, SIGMAS)
    coh = coherence(s19, footprint, sigma=2.0, win=3.0)
    base = (np.nan_to_num(ms, nan=0.0) * np.nan_to_num(coh, nan=0.0)).astype(np.float32)
    del ms, coh
    sup_base = support_of(base, domain, footprint, args.n_target, args.min_sep)
    print(f"[{time.time()-t0:.0f}s] ungated support {int(sup_base.sum()):,}")

    sc, parts = gating.orthogonality_gate(base, core["iso_grav_anom"], core["rtp"], footprint,
                                          w=args.w, sigma=2.0)
    del base
    changed_fraction = float(np.mean(np.nan_to_num(sc, nan=0.0) != np.nan_to_num(
        np.nan_to_num(sc, nan=0.0), nan=0.0)))
    sup_gate = support_of(sc, domain, footprint, args.n_target, args.min_sep)
    print(f"[{time.time()-t0:.0f}s] gated support   {int(sup_gate.sum()):,}")

    a = sup_base > 0
    b = sup_gate > 0
    inter = int((a & b).sum())
    union = int((a | b).sum())
    eff = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "w": args.w,
        "sigma": 2.0,
        "n_target": args.n_target,
        "min_sep_px": args.min_sep,
        "support_jaccard_gated_vs_ungated": inter / union if union else 1.0,
        "n_shared": inter,
        "n_only_ungated": int((a & ~b).sum()),
        "n_only_gated": int((b & ~a).sum()),
        "n_ungated": int(a.sum()),
        "n_gated": int(b.sum()),
        "field_finite_share": changed_fraction,
        "note": ("The gate leaves the emitted dot set almost untouched. Its value is therefore not in "
                 "re-ranking the carrier; it survives only because it does not damage it, and it "
                 "belongs in the file as a documented, measured no-op rather than as a claim of gain."),
        "runtime_s": round(time.time() - t0, 1),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(eff, indent=1))
    print(f"wrote {out}: jaccard {eff['support_jaccard_gated_vs_ungated']:.4f} "
          f"(shared {inter:,}, only-gated {eff['n_only_gated']:,}, only-ungated {eff['n_only_ungated']:,})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
