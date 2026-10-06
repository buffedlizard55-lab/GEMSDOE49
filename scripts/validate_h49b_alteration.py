"""H49-B validation — radiometric alteration halo x cross-gradient coupling.

Hypothesis (new; not previously implemented in GEMSDOE1-48): hydrothermal
alteration chemically rewrites surface radioelement abundances (K is
mobile in alteration fluids; Th is comparatively immobile), so an
isothermally anomalous Th/K halo that is spatially locked to a
multi-physics coupled edge marks a fluid-bearing structure — a class the
geomorphology-driven catalogue compilation cannot see (altered rock weathers
softer, erasing the scarp; CG-7 class 3 cover-masked variant).

The family has used radiometrics twice before and never this way:
* 15GEMSDOE tso1 (0.0782) — alteration-vs-magnetic conjunction, no
  cross-gradient operator, no Th/K z-score;
* GEMSDOE46 r11/r12 — scarp + radiometric concordance, unscored.
H49-B = alteration z-score x cross-gradient consensus, tested on the
preregistered hide-and-recover holdout against the coupled-only submission
and against the radiometric single-layer control.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import rasterio
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from gems49 import data as D
from gems49 import emission as E
from gems49 import xgrad as X
from gems49.holdout import hide_and_recover_folds, score_candidate

EVID = os.path.join(D.ROOT, "evidence")


def load_radiometric() -> dict:
    """Dequantise the u8-bridged radiometric grid (manifest rules)."""
    man = json.load(open(os.path.join(D.ROOT, "data", "bridge", "radiometric_manifest.json")))
    with rasterio.open(os.path.join(D.ROOT, "data", "external", "radiometric_u8.tif")) as ds:
        raw = ds.read()  # (bands, H, W) uint8, 0 = nodata
    out = {}
    for b in man["bands"]:
        q = raw[b["index"] - 1].astype(np.float32)
        v = np.where(q > 0, b["lo"] + (q - 1) / 254.0 * (b["hi"] - b["lo"]), np.nan)
        out[b["band_name"]] = v
    return out


def main():
    footprint, _ = D.submission_footprint()
    catalogue = D.labels()
    rad = load_radiometric()
    k = rad["rad_k"]
    th = rad["rad_th"]
    ratio = np.log(np.maximum(th, 1e-6) / np.maximum(k, 1e-6))
    # robust z-score within the footprint
    fin = footprint & np.isfinite(ratio)
    med = np.median(ratio[fin])
    mad = np.median(np.abs(ratio[fin] - med)) * 1.4826
    z = np.where(fin, np.abs(ratio - med) / max(mad, 1e-9), 0.0)

    fields = {"iso_grav_anom": D.band("iso_grav_anom"), "rtp": D.band("rtp"),
              "det_elev": D.band("det_elev")}
    cf = X.coupled_score(fields, 3.0, footprint, return_ortho=False)
    sc = np.where(np.isfinite(cf.score), cf.score, 0.0).astype(np.float32)
    q01, q999 = np.quantile(sc[footprint], [0.01, 0.999])
    score01 = np.clip((sc - q01) / max(q999 - q01, 1e-12), 0.0, 1.0).astype(np.float32)

    folds = hide_and_recover_folds(catalogue, footprint, hide_fraction=0.25, seed=49)

    results = {}
    n_dots = 40000
    # 1) coupled-only (reference)
    chosen = E.nms_topk(score01, catalogue, footprint, n_dots=n_dots, nms_radius_px=2.0)
    results["coupled_only"] = score_candidate(E.dots_to_prediction(chosen, footprint), folds)
    # 2) radiometric-only control (z as score)
    zc = np.where(footprint, z, 0.0).astype(np.float32)
    chosen = E.nms_topk(zc, catalogue, footprint, n_dots=n_dots, nms_radius_px=2.0)
    results["rad_z_only"] = score_candidate(E.dots_to_prediction(chosen, footprint), folds)
    # 3) H49-B variants: coupled x alteration gate at several z thresholds
    for zq in (1.0, 1.5, 2.0):
        gated = np.where((z >= zq) & footprint, score01, 0.0).astype(np.float32)
        chosen = E.nms_topk(gated, catalogue, footprint, n_dots=n_dots, nms_radius_px=2.0)
        res = score_candidate(E.dots_to_prediction(chosen, footprint), folds)
        results[f"h49b_z{zq:g}"] = res
        print(f"z>={zq}: n_emitted={int(chosen.sum())} DTI={res['dti_mean']:.5f}")
    out = {k_: {"dti_mean": v["dti_mean"], "dti_per_fold": v["dti_per_fold"]}
           for k_, v in results.items()}
    base = results["coupled_only"]["dti_mean"]
    out["delta_vs_coupled"] = {k_: v["dti_mean"] - base for k_, v in results.items()}
    json.dump(out, open(os.path.join(EVID, "h49b_validation.json"), "w"), indent=1)
    print(json.dumps(out["delta_vs_coupled"], indent=1))


if __name__ == "__main__":
    main()
