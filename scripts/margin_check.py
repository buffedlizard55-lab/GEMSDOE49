#!/usr/bin/env python3
"""Does the catalogue-buffer rule help the candidate more than it helps the baseline?

``prune_experiment`` showed that excluding a k-pixel buffer around the retained catalogue raises the
score.  But a placement rule handed to one field and not to its competitor is not a comparison.  This
script re-measures the headline margin with the rule applied to *both* arms, at each arm's own
best geometry, for k in {-1, 0, 2, 3, 4}.  Only the margin that survives this is reportable.

Run: python3 scripts/margin_check.py --out docs/data/margin_check.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import emission, gating, holdout, io  # noqa: E402
from gems49.baselines import single_layer_edge  # noqa: E402
from gems49.screen_util import coherence, multiscale_edge  # noqa: E402

N_FOLDS = 4
SEED = 49
KS = (-1, 0, 2, 3, 4)
ARMS = {
    "candidate:gate_ortho_w0.25": (40_000, 3),
    "baseline:single:det_elev_slope": (60_000, 3),
}


def build(name: str, footprint: np.ndarray) -> np.ndarray:
    core = io.read_bands(["iso_grav_anom", "rtp", "det_elev_slope"])
    if name.startswith("baseline:"):
        f = single_layer_edge(core["det_elev_slope"], footprint, 2.0)
    else:
        ms = multiscale_edge(core["det_elev_slope"], footprint, (1.0, 1.5, 2.5, 4.0))
        coh = coherence(core["det_elev_slope"], footprint, sigma=2.0, win=3.0)
        base = np.nan_to_num(ms, nan=0.0) * np.nan_to_num(coh, nan=0.0)
        del ms, coh
        f, _ = gating.orthogonality_gate(base, core["iso_grav_anom"], core["rtp"], footprint,
                                         w=0.25, sigma=2.0)
    f = np.nan_to_num(f, nan=0.0).astype(np.float32)
    f[~footprint] = 0.0
    return (f / float(f[footprint].max())).astype(np.float32)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/margin_check.json")
    args = ap.parse_args()
    t0 = time.time()

    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    seg, sizes = holdout.segment_catalogue(catalogue)
    fold_of_seg = holdout.fold_assignment(seg, len(sizes), N_FOLDS, seed=SEED)
    folds = [holdout.build_fold(catalogue, footprint, fold_of_seg, seg, k) for k in range(N_FOLDS)]
    sc = {(f["fold"], key): holdout.FoldScorer(f, key)
          for f in folds for key in ("truth_collared", "truth_unrestricted")}
    dist = {f["fold"]: distance_transform_edt(~f["train_mask"]) for f in folds}

    out = {}
    for arm, (n_target, sep) in ARMS.items():
        field = build(arm, footprint)
        order = emission.selection_order(field, footprint)
        del field
        rows = []
        for k in KS:
            col, unr = [], []
            for f in folds:
                domain = footprint if k < 0 else (footprint & (dist[f["fold"]] > k))
                sup = emission.select_from_order(order, footprint.shape, domain,
                                                 n_target=n_target, min_sep_px=sep)
                col.append(sc[(f["fold"], "truth_collared")].score(sup)["DTI"])
                unr.append(sc[(f["fold"], "truth_unrestricted")].score(sup)["DTI"])
                del sup, domain
            rows.append({"exclusion_px": k, "collared_mean": float(np.mean(col)),
                         "unrestricted_mean": float(np.mean(unr)),
                         "collared_folds": [round(v, 6) for v in col],
                         "unrestricted_folds": [round(v, 6) for v in unr]})
            print(f"  {arm:34s} k={k:>2}  collared {rows[-1]['collared_mean']:.4f}  "
                  f"unrestricted {rows[-1]['unrestricted_mean']:.4f}", flush=True)
        out[arm] = {"n_target": n_target, "min_sep_px": sep, "rows": rows}
        print(f"[{time.time()-t0:.0f}s] {arm} done", flush=True)

    margins = []
    for k in KS:
        a = next(r for r in out["candidate:gate_ortho_w0.25"]["rows"] if r["exclusion_px"] == k)
        b = next(r for r in out["baseline:single:det_elev_slope"]["rows"] if r["exclusion_px"] == k)
        margins.append({"exclusion_px": k,
                        "collared_margin": a["collared_mean"] - b["collared_mean"],
                        "unrestricted_margin": a["unrestricted_mean"] - b["unrestricted_mean"]})
    payload = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "protocol": ("both arms re-selected under identical exclusion radii, each at its own "
                     "best geometry; collared = held-out truth at least 600 m from the retained map, "
                     "unrestricted = all held-out truth"),
        "arms": out,
        "margins": margins,
        "runtime_s": round(time.time() - t0, 1),
    }
    p = Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=1))
    print(f"wrote {p}")
    for m in margins:
        print(f"  k={m['exclusion_px']:>2}  collared margin {m['collared_margin']:+.4f}   "
              f"unrestricted margin {m['unrestricted_margin']:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
