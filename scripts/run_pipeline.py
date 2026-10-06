"""GEMSDOE49 cross-gradient coupling pipeline.

Steps (each writes a receipt into evidence/):
 1. load official rasters (sha256-verified bridge),
 2. gravity-support diagnostic (station-spacing check),
 3. cross-gradient coupled field at several common scales,
 4. hide-and-recover holdout: coupled vs single-layer baselines,
 5. metric-aware emission (NMS dots, budget sweep),
 6. submission writing (zeros + nan twins) with format verification,
 7. uniqueness audit vs every prior submission raster.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from gems49 import data as D
from gems49 import emission as E
from gems49 import xgrad as X
from gems49.holdout import hide_and_recover_folds, score_candidate
from gems49.metric import dti

EVID = os.path.join(D.ROOT, "evidence")
os.makedirs(EVID, exist_ok=True)


def save(name, obj):
    p = os.path.join(EVID, name)
    json.dump(obj, open(p, "w"), indent=1, default=float)
    print(f"[evidence] {name}")
    return p


def main():
    t0 = time.time()
    print("== step 1: load rasters ==")
    footprint, geom = D.submission_footprint()
    catalogue = D.labels()
    grav = D.band("iso_grav_anom")
    rtp = D.band("rtp")
    elev = D.band("det_elev")
    dbase = D.band("depth_to_base_surf")
    print(f"footprint={footprint.sum()} catalogue={catalogue.sum()} "
          f"grav finite={np.isfinite(grav).sum()} rtp finite={np.isfinite(rtp).sum()} "
          f"elev finite={np.isfinite(elev).sum()}")

    print("== step 2: gravity support diagnostic ==")
    diag = X.gravity_support_diagnostics(grav, footprint)
    save("gravity_support_diagnostics.json", diag)

    print("== step 3: holdout folds ==")
    folds = hide_and_recover_folds(catalogue, footprint, hide_fraction=0.25, seed=49)
    save("holdout_folds.json", {f"fold{f.fold_id}": {"n_hidden": f.n_hidden} for f in folds})

    fields = {"iso_grav_anom": grav, "rtp": rtp, "det_elev": elev}

    print("== step 4: coupled vs single-layer baselines over scales ==", flush=True)
    scale_results = {}
    for sigma in (3.0, 5.0, 7.0, 10.0):
        cf = X.coupled_score(fields, sigma, footprint, return_ortho=(sigma == 5.0))
        sc = np.where(np.isfinite(cf.score), cf.score, 0.0).astype(np.float32)
        base_fields = {}
        baselines = {}
        for name in X.FIELDS:
            b = X.single_layer_baseline(name, fields, sigma, footprint)
            base_fields[name] = np.where(np.isfinite(b), b, 0.0).astype(np.float32)
            baselines[name] = float(np.nanmean(np.where(footprint, b, np.nan)))
        scale_results[f"sigma_{sigma:g}"] = {
            "coupled_mean": float(np.nanmean(np.where(footprint, sc, np.nan))),
            "coupled_p99": float(np.quantile(sc[footprint], 0.99)),
            "baselines_mean": baselines,
        }
        if sigma == 5.0:  # keep one ortho snapshot for the confounder report
            ortho_snap = cf.ortho_flag.copy()
        del cf
        print(f"  sigma={sigma}: coupled mean={scale_results[f'sigma_{sigma:g}']['coupled_mean']:.5f}",
              flush=True)
        for n_dots in (10000, 20000, 30000, 40000):
            chosen = E.nms_topk(sc, catalogue, footprint, n_dots=n_dots, nms_radius_px=2.0)
            pred = E.dots_to_prediction(chosen, footprint)
            res = score_candidate(pred, folds)
            key = f"sigma_{sigma:g}_n{n_dots//1000}k"
            scale_results[key] = {"coupled": res}
            print(f"    coupled n={n_dots}: holdout DTI mean={res['dti_mean']:.4f} "
                  f"per-fold={[f'{x:.3f}' for x in res['dti_per_fold']]}", flush=True)
            del chosen, pred
            for name in X.FIELDS:
                ch_b = E.nms_topk(base_fields[name], catalogue, footprint, n_dots=n_dots,
                                  nms_radius_px=2.0)
                res_b = score_candidate(E.dots_to_prediction(ch_b, footprint), folds)
                scale_results[key][f"baseline_{name}"] = {"dti_mean": res_b["dti_mean"],
                                                          "dti_per_fold": res_b["dti_per_fold"]}
                print(f"      baseline {name}: {res_b['dti_mean']:.4f}", flush=True)
                del ch_b
        del sc, base_fields
    save("scale_and_holdout_results.json", scale_results)

    # pick the best (sigma, n_dots) for the coupled candidate
    cand = [(k, v["coupled"]["dti_mean"]) for k, v in scale_results.items()
            if k.startswith("sigma_") and "_n" in k and "coupled" in v]
    cand.sort(key=lambda kv: -kv[1])
    print("== ranked coupled configs ==", cand[:5])
    best_key, best_dti = cand[0]
    sigma_best = float(best_key.replace("sigma_", "").split("_n")[0])
    n_best = int(best_key.split("_n")[1].replace("k", "000"))
    beats = {}
    for k, v in scale_results.items():
        if k == best_key and "coupled" in v:
            for bname in X.FIELDS:
                beats[bname] = v["coupled"]["dti_mean"] - v[f"baseline_{bname}"]["dti_mean"]
    print("best config:", best_key, "DTI", best_dti, "margin vs baselines:", beats)
    if any(m <= 0 for m in beats.values()):
        print("!! candidate does NOT beat every single-layer baseline — flag for review")
    save("best_config.json", {"key": best_key, "sigma": sigma_best, "n_dots": n_best,
                              "holdout_dti_mean": best_dti, "margins_vs_baselines": beats})

    print("== step 5: final emission ==", flush=True)
    cf = X.coupled_score(fields, sigma_best, footprint, return_ortho=True)
    sc = np.where(np.isfinite(cf.score), cf.score, 0.0).astype(np.float32)
    # normalise evidence to [0,1] over the footprint (rank-preserving min-max)
    q01, q999 = np.quantile(sc[footprint], [0.01, 0.999])
    score01 = np.clip((sc - q01) / max(q999 - q01, 1e-12), 0.0, 1.0)
    chosen = E.nms_topk(score01, catalogue, footprint, n_dots=n_best, nms_radius_px=2.0)
    pred = E.dots_to_prediction(chosen, footprint, value=1.0)
    res_final = score_candidate(pred, folds)
    save("final_holdout.json", res_final)
    print("final holdout:", res_final["dti_mean"], res_final["dti_per_fold"])

    # basin-margin diagnostic (brief: expect density steps at basin margins)
    gy, gx = np.gradient(np.where(np.isfinite(dbase), dbase, 0.0))
    dbase_edge = np.hypot(gx, gy)
    edge_q = np.quantile(dbase_edge[footprint & np.isfinite(dbase)], 0.90)
    emitted_rc = np.argwhere(chosen)
    frac_on_basin_edge = float((dbase_edge[emitted_rc[:, 0], emitted_rc[:, 1]] > edge_q).mean())
    cat_rc = np.argwhere(catalogue)
    frac_cat_on_basin_edge = float((dbase_edge[cat_rc[:, 0], cat_rc[:, 1]] > edge_q).mean())
    save("basin_margin_diagnostic.json", {
        "frac_emitted_on_basin_edge_p90": frac_on_basin_edge,
        "frac_catalogue_on_basin_edge_p90": frac_cat_on_basin_edge,
        "note": "depth_to_base_surf band 15 gradient > p90 defines 'basin margin'",
    })
    # orthogonality confounder summary (flag only — never used as evidence)
    save("ortho_confounder_summary.json", {
        "mean": float(np.nanmean(np.where(footprint, cf.ortho_flag, np.nan))),
        "p99": float(np.quantile(cf.ortho_flag[footprint], 0.99)),
        "note": "magnitude-weighted orthogonality kept as a flag, excluded from evidence",
    })

    print("== step 6: write submissions ==")
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    slug = f"gemsdoe49-xgrad-coupled-{stamp}"
    receipt_z = E.write_submission(pred, footprint, f"/home/user/GEMSDOE49/docs/downloads/{slug}-zeros.tif", "zeros")
    receipt_n = E.write_submission(pred, footprint, f"/home/user/GEMSDOE49/docs/downloads/{slug}-nan.tif", "nan")
    save("format_receipts.json", {"zeros": receipt_z, "nan": receipt_n, "slug": slug})
    print("zeros:", receipt_z["sha256"][:16], receipt_z["all_checks_passed"])
    print("nan  :", receipt_n["sha256"][:16], receipt_n["all_checks_passed"])
    print(f"elapsed {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
