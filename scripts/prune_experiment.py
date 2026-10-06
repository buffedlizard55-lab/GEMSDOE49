#!/usr/bin/env python3
"""Does excluding a buffer around the known-fault catalogue help or hurt?

The highest-scoring artifact in our cohort prunes every dot within 2 px (200 m) of the catalogue and
reports a live-mirror gain of +0.004870 in 4/4 folds.  That is a claim about the live scorer, which
we cannot reach; but the mechanism is exactly reproducible on our holdout, where the *retained* part
of the catalogue plays the role of the masked known faults:

    for k in 0..4 px:  select the dots with a k-pixel exclusion zone around the retained map,
                       score on the *held-out* segments (which is where truth lives).

``truth_collared`` drops held-out truth pixels that sit within the exclusion distance of the retained
map, so it answers "what happens away from the interface"; ``truth_unrestricted`` keeps them, so it
answers "what happens when a new fault is mapped right against a known one" -- which the staff ruling
says is allowed.  Reporting both is the only honest way to use this instrument, because production
new faults are *adjacent* to their catalogue parents while our holdout truth is a disjoint segment.

Run: python3 scripts/prune_experiment.py --out docs/data/prune_experiment.json
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
from gems49.screen_util import coherence, multiscale_edge  # noqa: E402

N_FOLDS = 4
SEED = 49
SIGMAS = (1.0, 1.5, 2.5, 4.0)
KS = (-1, 0, 1, 2, 3, 4, 6)  # -1 = no exclusion (the convention used by the earlier sweep)
GEOMETRY = ((40_000, 3), (60_000, 3))


def build_field(footprint: np.ndarray) -> np.ndarray:
    core = io.read_bands(["iso_grav_anom", "rtp", "det_elev_slope"])
    ms = multiscale_edge(core["det_elev_slope"], footprint, SIGMAS)
    coh = coherence(core["det_elev_slope"], footprint, sigma=2.0, win=3.0)
    base = np.nan_to_num(ms, nan=0.0) * np.nan_to_num(coh, nan=0.0)
    del ms, coh
    sc, _ = gating.orthogonality_gate(base, core["iso_grav_anom"], core["rtp"], footprint,
                                      w=0.25, sigma=2.0)
    del base
    f = np.nan_to_num(sc, nan=0.0).astype(np.float32)
    f[~footprint] = 0.0
    return (f / float(f[footprint].max())).astype(np.float32)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/prune_experiment.json")
    args = ap.parse_args()
    t0 = time.time()

    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    field = build_field(footprint)
    order = emission.selection_order(field, footprint)
    del field
    print(f"[{time.time()-t0:.0f}s] field built and ranked")

    _print_catalogue = int(catalogue.sum())
    seg, sizes = holdout.segment_catalogue(catalogue)
    fold_of_seg = holdout.fold_assignment(seg, len(sizes), N_FOLDS, seed=SEED)
    folds = [holdout.build_fold(catalogue, footprint, fold_of_seg, seg, k) for k in range(N_FOLDS)]
    scorers = {(f["fold"], key): holdout.FoldScorer(f, key)
               for f in folds for key in ("truth_collared", "truth_unrestricted")}
    print(f"[{time.time()-t0:.0f}s] folds ready")

    rows = []
    for n_target, sep in GEOMETRY:
        for k in KS:
            per_fold = {"collared": [], "unrestricted": [], "kept": []}
            for f in folds:
                # exclusion zone around the retained map only -- the production analogue of the
                # masked known-fault catalogue
                # The production analogue of the masked known-fault catalogue is this fold's
                # retained map (``train_mask``); held-out pixels are the unknown truth and must
                # never take part in the exclusion.
                if k < 0:
                    domain = footprint                      # no exclusion at all (sweep convention)
                else:
                    d = distance_transform_edt(~f["train_mask"])
                    domain = footprint & (d > k)            # k = 0 excludes the exact pixels
                sup = emission.select_from_order(order, footprint.shape, domain,
                                                 n_target=n_target, min_sep_px=sep)
                for key in ("collared", "unrestricted"):
                    per_fold[key].append(
                        scorers[(f["fold"], f"truth_{key}")].score(sup)["DTI"])
                per_fold["kept"].append(int(sup.sum()))
                del sup, domain
            row = {
                "n_target": n_target, "min_sep_px": sep, "exclusion_px": k,
                "collared_mean": float(np.mean(per_fold["collared"])),
                "collared_folds": [round(v, 6) for v in per_fold["collared"]],
                "unrestricted_mean": float(np.mean(per_fold["unrestricted"])),
                "unrestricted_folds": [round(v, 6) for v in per_fold["unrestricted"]],
                "dots_emitted": per_fold["kept"],
            }
            rows.append(row)
            print(f"  n={n_target:>6,} sep={sep} excl={k}px  collared {row['collared_mean']:.4f}  "
                  f"unrestricted {row['unrestricted_mean']:.4f}  dots {row['dots_emitted']}",
                  flush=True)

    payload = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "field": "gate_ortho_w0.25 (the shipped family)",
        "protocol": ("dots selected with a k-pixel exclusion zone around the *retained* catalogue of "
                     "each fold; scored on held-out segments.  See the module docstring."),
        "rows": rows,
        "runtime_s": round(time.time() - t0, 1),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
