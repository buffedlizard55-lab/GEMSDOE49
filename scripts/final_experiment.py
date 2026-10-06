"""Final detector screen under a leave-one-fold-out (LOFO) selection protocol.

Why LOFO: the placement geometry (support size N, minimum separation) is a free parameter.
Choosing it on the same folds that are then reported inflates the winner by the maximum of
many noisy draws.  LOFO removes that bias -- for each held-out fold the geometry is chosen on
the *other three* folds only, and only the held-out fold's score is reported.

Run: python -u scripts/final_experiment.py --out docs/data/final_experiment.json
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

from gems49 import emission, gating, holdout, io, methods  # noqa: E402
from gems49.baselines import single_layer_edge  # noqa: E402
from gems49.crossgrad import horizontal_gradient, magnitude, robust_unit, smooth  # noqa: E402
from gems49.screen_util import coherence, multiscale_edge, rank_normalise  # noqa: E402
from gems49.spec import BANDS  # noqa: E402

N_FOLDS = 4
SEED = 49
N_GRID = (15_000, 20_000, 30_000, 40_000, 60_000)
SEP_GRID = (3, 4)


CORE = ("iso_grav_anom", "rtp", "det_elev", "det_elev_slope")


def iter_fields(footprint: np.ndarray):
    """Yield (name, field) one at a time: build -> score -> discard keeps peak RAM ~1 GB."""
    core = io.read_bands(list(CORE))
    g, m = core["iso_grav_anom"], core["rtp"]
    s19 = core["det_elev_slope"]

    # ---- single-layer gradient baselines (the reference family) --------------
    for _i, nm, _c, _d in BANDS:
        if nm in CORE:
            field = single_layer_edge(core[nm], footprint, 2.0)
        else:
            field = single_layer_edge(io.read_band(nm), footprint, 2.0)
        yield f"single:{nm}", field
        del field

    # ---- topographic structural carriers ------------------------------------
    ms = multiscale_edge(s19, footprint, (1.0, 1.5, 2.5, 4.0))
    yield "topo_ms", ms
    base = None
    for tag, sg, win in (("coh1", 1.0, 2.0), ("coh", 2.0, 3.0), ("coh4", 3.0, 5.0)):
        c = coherence(s19, footprint, sigma=sg, win=win)
        v = (np.nan_to_num(ms, nan=0.0) * np.nan_to_num(c, nan=0.0)).astype(np.float32)
        del c
        yield f"topo_ms_{tag}", v
        if tag == "coh":
            base = v
        else:
            del v
    del ms

    # ---- the brief's coupling, unmodified ------------------------------------
    cg = methods.family("cg_gxm", core, footprint, sigma=2.0)
    yield "cg_gxm", cg
    tr = methods.family("triple", core, footprint, sigma=2.0)
    yield "triple", tr
    yield "tg_g", methods.family("tg_g", core, footprint, sigma=2.0)

    # ---- cross-gradient used as a confounder GATE on the topographic carrier --
    for w in (0.25, 0.5, 0.75, 1.0):
        sc, _d = gating.orthogonality_gate(base, g, m, footprint, w=w, sigma=2.0)
        yield f"gate_ortho_w{w:g}", sc
        del sc
    for w in (0.25, 0.5, 0.75):
        sc, _d = gating.alignment_gate(base, g, m, footprint, w=w, sigma=2.0)
        yield f"gate_align_w{w:g}", sc
        del sc

    # ---- fusions -------------------------------------------------------------
    r_b = rank_normalise(base, footprint)
    r_t = rank_normalise(tr, footprint)
    r_x = rank_normalise(cg, footprint)
    nb = np.nan_to_num(r_b, nan=0.0)
    yield "rankavg_base_cg", (0.7 * nb + 0.3 * np.nan_to_num(r_x, nan=0.0)).astype(np.float32)
    yield "rankavg_base_triple", (0.7 * nb + 0.3 * np.nan_to_num(r_t, nan=0.0)).astype(np.float32)
    del r_b, r_t, r_x, nb, cg, tr
    del core

    # ---- no-information controls --------------------------------------------
    h, w_ = footprint.shape
    yy, xx = np.mgrid[0:h, 0:w_]
    yield "control:uniform_lattice4", np.where(footprint & (yy % 4 == 0) & (xx % 4 == 0), 1.0, 0.0).astype(np.float32)
    rng = np.random.default_rng(SEED)
    yield "control:random", np.where(footprint, rng.random(footprint.shape), np.float32(np.nan)).astype(np.float32)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/final_experiment.json")
    args = ap.parse_args()

    t0 = time.time()
    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    seg, sizes = holdout.segment_catalogue(catalogue)
    fold_of_seg = holdout.fold_assignment(seg, len(sizes), N_FOLDS, seed=SEED)
    folds = [holdout.build_fold(catalogue, footprint, fold_of_seg, seg, k) for k in range(N_FOLDS)]
    del seg
    sc_c = [holdout.FoldScorer(f, "truth_collared") for f in folds]
    sc_u = [holdout.FoldScorer(f, "truth_unrestricted") for f in folds]
    print(f"[{time.time()-t0:.0f}s] folds: collared truth {[s.n_truth for s in sc_c]} "
          f"unrestricted {[s.n_truth for s in sc_u]}", flush=True)

    geo = [(n, s) for s in SEP_GRID for n in N_GRID]
    table: dict[str, dict] = {}
    for name, field in iter_fields(footprint):
        f0 = np.nan_to_num(field, nan=0.0).astype(np.float32)
        del field
        order = emission.selection_order(f0, footprint)
        cells = []
        for n_target, sep in geo:
            sup = emission.select_from_order(order, f0.shape, footprint,
                                             n_target=n_target, min_sep_px=sep)
            cs = [s.score(sup)["DTI"] for s in sc_c]
            us = [s.score(sup)["DTI"] for s in sc_u]
            cells.append({"n_target": n_target, "min_sep_px": sep, "n_emitted": int(sup.sum()),
                          "collared": [round(float(v), 6) for v in cs],
                          "unrestricted": [round(float(v), 6) for v in us]})
            del sup
        # leave-one-fold-out geometry selection
        lofo = []
        chosen = []
        for k in range(N_FOLDS):
            others = [i for i in range(N_FOLDS) if i != k]
            best_j = max(range(len(cells)),
                         key=lambda j: np.mean([cells[j]["collared"][i] for i in others]))
            lofo.append(cells[best_j]["collared"][k])
            chosen.append({"held_out_fold": k, "n_target": cells[best_j]["n_target"],
                           "min_sep_px": cells[best_j]["min_sep_px"]})
        lofo_u = []
        for k in range(N_FOLDS):
            others = [i for i in range(N_FOLDS) if i != k]
            best_j = max(range(len(cells)),
                         key=lambda j: np.mean([cells[j]["unrestricted"][i] for i in others]))
            lofo_u.append(cells[best_j]["unrestricted"][k])
        best_insample = max(cells, key=lambda c: float(np.mean(c["collared"])))
        table[name] = {
            "cells": cells,
            "lofo_collared_mean": float(np.mean(lofo)),
            "lofo_collared_folds": [round(float(v), 6) for v in lofo],
            "lofo_unrestricted_mean": float(np.mean(lofo_u)),
            "lofo_unrestricted_folds": [round(float(v), 6) for v in lofo_u],
            "lofo_chosen": chosen,
            "insample_best": {"n_target": best_insample["n_target"],
                              "min_sep_px": best_insample["min_sep_px"],
                              "collared_mean": float(np.mean(best_insample["collared"]))},
        }
        print(f"  {name:26s} LOFO={table[name]['lofo_collared_mean']:.4f} "
              f"(unres {table[name]['lofo_unrestricted_mean']:.4f})  "
              f"in-sample-best={table[name]['insample_best']['collared_mean']:.4f}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        del f0, order

    def rank(key):
        return sorted(table.items(), key=lambda kv: -kv[1][key])

    singles = [(k, v) for k, v in table.items() if k.startswith("single:")]
    best_single = max(singles, key=lambda kv: kv[1]["lofo_collared_mean"])
    cands = [(k, v) for k, v in table.items() if not k.startswith("single:") and not k.startswith("control:")]
    best_cand = max(cands, key=lambda kv: kv[1]["lofo_collared_mean"])

    print("\n=== LOFO-validated leaderboard (collared instrument) ===")
    for k, v in rank("lofo_collared_mean"):
        print(f"  {k:26s} {v['lofo_collared_mean']:.4f}   (unres {v['lofo_unrestricted_mean']:.4f})")
    print(f"\n  best single-layer baseline : {best_single[0]} {best_single[1]['lofo_collared_mean']:.4f}")
    print(f"  best candidate             : {best_cand[0]} {best_cand[1]['lofo_collared_mean']:.4f}")
    d = best_cand[1]["lofo_collared_mean"] - best_single[1]["lofo_collared_mean"]
    print(f"  delta = {d:+.4f}  ->  {'PASS' if d > 0 else 'FAIL'}")

    payload = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "protocol": "leave-one-fold-out geometry selection",
        "seed": SEED, "n_folds": N_FOLDS,
        "n_grid": list(N_GRID), "sep_grid": list(SEP_GRID),
        "fold_truth_sizes": {"collared": [s.n_truth for s in sc_c],
                             "unrestricted": [s.n_truth for s in sc_u]},
        "table": table,
        "best_single": {"name": best_single[0], **best_single[1]},
        "best_candidate": {"name": best_cand[0], **best_cand[1]},
        "gate_pass": bool(d > 0),
        "runtime_s": round(time.time() - t0, 1),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    print(f"\nwrote {out} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
