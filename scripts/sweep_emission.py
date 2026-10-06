"""Sweep metric-aware placement (support size N x minimum separation) per detector family.

The metric's marginal rule says the decision that matters is *where* support lands relative to
the kernel, and the algebra says binary support dominates graded support.  So for each detector
we sweep the only two free placement parameters and report the hide-and-recover DTI.

Run:  python -u scripts/sweep_emission.py --out docs/data/emission_sweep.json
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

from gems49 import emission, holdout, io, methods  # noqa: E402

N_FOLDS = 4
SEED = 49
N_GRID = (10_000, 20_000, 30_000, 40_000, 60_000, 90_000)
SEP_GRID = (1, 2, 3, 4)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/emission_sweep.json")
    ap.add_argument("--families", default="triple,tg_g,single:det_elev_slope")
    ap.add_argument("--sigma", type=float, default=2.0)
    args = ap.parse_args()

    t0 = time.time()
    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    seg, sizes = holdout.segment_catalogue(catalogue)
    fold_of_seg = holdout.fold_assignment(seg, len(sizes), N_FOLDS, seed=SEED)
    folds = [holdout.build_fold(catalogue, footprint, fold_of_seg, seg, k) for k in range(N_FOLDS)]
    del seg

    names = args.families.split(",")
    needed = set(methods.CORE_BANDS)
    for n in names:
        if n.startswith("single:"):
            needed.add(n.split(":", 1)[1])
    bands = io.read_bands(sorted(needed))
    print(f"[{time.time()-t0:6.1f}s] loaded {len(bands)} bands", flush=True)

    table: dict[str, dict] = {}
    for fam in names:
        if fam.startswith("single:"):
            from gems49.baselines import single_layer_edge

            field = single_layer_edge(bands[fam.split(":", 1)[1]], footprint, args.sigma)
        else:
            field = methods.family(fam, bands, footprint, sigma=args.sigma)
        field = np.nan_to_num(field, nan=0.0)
        # one global argsort, reused for every (N, sep) cell
        order = np.argsort(-field.ravel(), kind="stable")
        print(f"[{time.time()-t0:6.1f}s] {fam}: sorted", flush=True)
        rows = []
        for sep in SEP_GRID:
            for n in N_GRID:
                sup = _select_from_order(order, field.shape, footprint, n, sep)
                ds = []
                for f in folds:
                    ds.append(holdout.score_fold(sup, f, truth_key="truth_collared")["DTI"])
                rows.append({"n_target": n, "min_sep_px": sep,
                             "n_emitted": int(sup.sum()),
                             "collared_mean": float(np.mean(ds)),
                             "collared_folds": [round(float(v), 6) for v in ds]})
                print(f"    {fam:24s} sep={sep} n={n:6d} -> {np.mean(ds):.4f}", flush=True)
                del sup
        table[fam] = {"rows": rows,
                      "best": max(rows, key=lambda r: r["collared_mean"])}
        print(f"  BEST {fam}: {table[fam]['best']}", flush=True)
        del field, order

    payload = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "sigma": args.sigma, "seed": SEED, "n_folds": N_FOLDS,
        "table": table,
        "runtime_s": round(time.time() - t0, 1),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    print(f"wrote {out} ({time.time()-t0:.1f}s)")
    return 0


def _select_from_order(order, shape, footprint, n_target, min_sep_px):
    """Same greedy Poisson-disk selection as emission.select_support, from a precomputed order."""
    h, w = shape
    flat_valid = footprint.ravel()
    blocked = np.zeros(order.shape, dtype=bool)
    out = np.zeros(order.shape, dtype=np.float32)
    radii = emission._disc_offsets(min_sep_px)
    taken = 0
    for i in order:
        if not flat_valid[i] or blocked[i]:
            continue
        y, x = divmod(int(i), w)
        out[i] = 1.0
        taken += 1
        for dy, dx in radii:
            yy, xx = y + dy, x + dx
            if 0 <= yy < h and 0 <= xx < w:
                blocked[yy * w + xx] = True
        if taken >= n_target:
            break
    return out.reshape(shape)


if __name__ == "__main__":
    raise SystemExit(main())
