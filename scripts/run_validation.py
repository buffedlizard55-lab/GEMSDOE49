"""Hide-and-recover experiment: does cross-gradient coupling beat a single-layer gradient?

Memory-lean by construction: this host has ~4 GB of RAM, so bands are read and released one
family at a time and only supports (one byte per pixel) are retained between steps.

Run:  python -u scripts/run_validation.py --out docs/data/validation.json
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
from gems49.baselines import single_layer_edge  # noqa: E402
from gems49.metric import marginal_max_distance  # noqa: E402
from gems49.spec import BANDS  # noqa: E402

N_FOLDS = 4
SEED = 49
from gems49.methods import CORE_BANDS  # noqa: E402


def score_support(support: np.ndarray, folds: list[dict]) -> dict:
    out = {}
    for key, label in (("truth_collared", "collared"), ("truth_unrestricted", "unrestricted")):
        vals = []
        for f in folds:
            r = holdout.score_fold(support, f, truth_key=key)
            vals.append({"fold": f["fold"], "DTI": r["DTI"], "TPw": r["TPw"],
                         "FPw": r["FPw"], "n_truth": r["n_truth"], "n_pred_px": r["n_pred_px"]})
        out[label] = {
            "mean": float(np.mean([v["DTI"] for v in vals])),
            "folds": vals,
        }
    return out


def record(support: np.ndarray, folds: list[dict]) -> dict:
    res = score_support(support, folds)
    return {
        "n_emitted": int(support.sum()),
        "collared_mean": res["collared"]["mean"],
        "unrestricted_mean": res["unrestricted"]["mean"],
        "collared_folds": [round(v["DTI"], 6) for v in res["collared"]["folds"]],
        "unrestricted_folds": [round(v["DTI"], 6) for v in res["unrestricted"]["folds"]],
        "collared_TPw": [round(v["TPw"], 3) for v in res["collared"]["folds"]],
        "collared_FPw": [round(v["FPw"], 3) for v in res["collared"]["folds"]],
        "collared_n_truth": [v["n_truth"] for v in res["collared"]["folds"]],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/data/validation.json")
    ap.add_argument("--n-target", type=int, default=40000)
    ap.add_argument("--min-sep", type=int, default=2)
    ap.add_argument("--sigma", type=float, default=2.0)
    ap.add_argument("--skip-singles", action="store_true")
    ap.add_argument("--families", default=None, help="comma-separated subset of CANDIDATE_FAMILIES")
    args = ap.parse_args()

    t0 = time.time()
    footprint = io.read_template_footprint()
    catalogue = io.read_catalogue()
    print(f"[{time.time()-t0:6.1f}s] footprint={int(footprint.sum())} catalogue={int(catalogue.sum())}", flush=True)

    seg, sizes = holdout.segment_catalogue(catalogue)
    n_seg = len(sizes)
    fold_of_seg = holdout.fold_assignment(seg, n_seg, N_FOLDS, seed=SEED)
    folds = [holdout.build_fold(catalogue, footprint, fold_of_seg, seg, k) for k in range(N_FOLDS)]
    del seg
    print(f"[{time.time()-t0:6.1f}s] segments={n_seg} (min/med/max "
          f"{sizes.min()}/{int(np.median(sizes))}/{sizes.max()})", flush=True)
    for f in folds:
        print(f"    fold {f['fold']}: train={f['n_train_px']:6d} held={f['n_held_px']:6d} "
              f"truth_un={f['n_truth_unrestricted']:6d} truth_col={f['n_truth_collared']:6d} "
              f"({100*f['collared_fraction']:.1f}% survive collar)", flush=True)

    controls = {"leakage_probe_TPw": [], "null_probe_DTI": []}
    for f in folds:
        controls["leakage_probe_TPw"].append(holdout.leakage_probe(f)["TPw"])
        controls["null_probe_DTI"].append(holdout.null_probe(f)["DTI"])
    print(f"    control leakage_probe max TPw = {max(controls['leakage_probe_TPw'])} "
          f"(must be 0) | null DTI max = {max(controls['null_probe_DTI'])} (must be 0)", flush=True)

    # oracle: dots placed on the collared truth with the same geometry -> instrument ceiling
    oracle = {}
    for f in folds:
        sup = emission.select_support(
            f["truth_collared"].astype(np.float32),
            footprint & ~f["train_mask"],
            n_target=args.n_target, min_sep_px=args.min_sep)
        oracle[str(f["fold"])] = holdout.score_fold(sup, f, truth_key="truth_collared")["DTI"]
        del sup
    print(f"    oracle (on-truth dots, n<={args.n_target}, sep={args.min_sep}): "
          f"mean DTI = {np.mean(list(oracle.values())):.4f}", flush=True)

    results: dict[str, dict] = {}

    if not args.skip_singles:
        for idx, name, _cat, _desc in BANDS:
            field = single_layer_edge(io.read_band(name), footprint, sigma=args.sigma)
            sup = emission.select_support(np.nan_to_num(field, nan=0.0), footprint,
                                          n_target=args.n_target, min_sep_px=args.min_sep)
            results[f"single:{name}"] = record(sup, folds)
            r = results[f"single:{name}"]
            print(f"    single:{name:22s} n={r['n_emitted']:6d} col={r['collared_mean']:.4f} "
                  f"unres={r['unrestricted_mean']:.4f}", flush=True)
            del field, sup

    bands = io.read_bands(list(CORE_BANDS))
    print(f"[{time.time()-t0:6.1f}s] loaded {len(bands)} core bands", flush=True)
    fam_list = args.families.split(",") if args.families else list(methods.CANDIDATE_FAMILIES)
    for fam in fam_list:
        t1 = time.time()
        field = methods.family(fam, bands, footprint, sigma=args.sigma)
        sup = emission.select_support(np.nan_to_num(field, nan=0.0), footprint,
                                      n_target=args.n_target, min_sep_px=args.min_sep)
        results[fam] = record(sup, folds)
        r = results[fam]
        print(f"    {fam:24s} n={r['n_emitted']:6d} col={r['collared_mean']:.4f} "
              f"unres={r['unrestricted_mean']:.4f}   [{time.time()-t1:.1f}s]", flush=True)
        del field, sup

    singles = {k: v for k, v in results.items() if k.startswith("single:")}
    cands = {k: v for k, v in results.items() if not k.startswith("single:")}
    best_single = max(singles.items(), key=lambda kv: kv[1]["collared_mean"]) if singles else None
    best_cand = max(cands.items(), key=lambda kv: kv[1]["collared_mean"])

    print("\n=== PROMOTION GATE: candidate must beat the best single-layer gradient ===")
    if best_single:
        print(f"  best single-layer baseline : {best_single[0]:34s} {best_single[1]['collared_mean']:.4f}")
    print(f"  best candidate             : {best_cand[0]:34s} {best_cand[1]['collared_mean']:.4f}")
    if best_single:
        d = best_cand[1]["collared_mean"] - best_single[1]["collared_mean"]
        print(f"  delta                      : {d:+.4f}   {'PASS' if d > 0 else 'FAIL'}")

    payload = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_folds": N_FOLDS, "seed": SEED, "n_segments": int(n_seg),
        "n_target": args.n_target, "min_sep_px": args.min_sep, "sigma": args.sigma,
        "folds": [{k: v for k, v in f.items() if not isinstance(v, np.ndarray)} for f in folds],
        "controls": {"leakage_probe_TPw_max": float(max(controls["leakage_probe_TPw"])),
                     "null_probe_DTI_max": float(max(controls["null_probe_DTI"]))},
        "oracle_collared_DTI": oracle,
        "oracle_collared_mean": float(np.mean(list(oracle.values()))),
        "results": results,
        "best_single": ({"name": best_single[0], **best_single[1]} if best_single else None),
        "best_candidate": {"name": best_cand[0], **best_cand[1]},
        "marginal_max_distance_at_0.28": marginal_max_distance(0.28),
        "runtime_s": round(time.time() - t0, 1),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1))
    print(f"\nwrote {out}  ({time.time()-t0:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
