"""Detector families evaluated against the single-layer gradient baseline.

Every family returns a continuous score in [0, 1] (NaN outside its domain) and is then put
through the *same* metric-aware placement, so the comparison isolates the detector.

Design note on what "agreement as an operator" means here
---------------------------------------------------------
The brief is explicit that agreement must be an operator, not a product of scores.  Three
operators are implemented and one counter-example is kept deliberately as a control:

``balance``    B = 2*min(|G1|,|G2|)/(|G1|+|G2|)
               suppresses an edge that lives in one field only (a buried volcanic edge moves
               the magnetic gradient alone, so B -> 0 there);
``alignment``  A = |G1||G2|cos^2(theta)
               large only when both gradients are large *and* parallel;
``triple``     S = |G_topo| * (|U_g + U_m| / 2) * cos^2(angle(G_topo, U_g + U_m))
               where U are unit gradient vectors.  This asks three questions at once -- is
               there a strong topographic edge, do gravity and magnetics agree with each
               other, and do they agree with the topography -- and it is not a product of
               per-field scores because the middle factor is a *vector sum of directions*.
``product``    control: |G_topo|_norm * max(|G_g|_norm, |G_m|_norm), the naive "product of
               independent scores" that the brief rules out.
"""

from __future__ import annotations

import numpy as np

from .baselines import single_layer_edge
from .crossgrad import (
    horizontal_gradient,
    magnitude,
    pair_coupling,
    robust_unit,
    smooth,
)

NAN32 = np.float32(np.nan)


def _norm(x: np.ndarray, valid: np.ndarray, pct: float = 95.0) -> np.ndarray:
    n, _ = robust_unit(np.nan_to_num(x, nan=0.0), valid & np.isfinite(x), pct)
    return np.where(valid & np.isfinite(x), n, NAN32).astype(np.float32)


def _grad(field: np.ndarray, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    return horizontal_gradient(smooth(field, sigma))


def _coupled(g1, g2, v1, v2, pct=95.0) -> np.ndarray:
    return pair_coupling(g1, g2, v1, v2, pct)["coupled"]


def triple_alignment(g_t, g_g, g_m, valid, pct: float = 95.0) -> tuple[np.ndarray, dict]:
    """Topography-anchored three-field structural agreement.

    Returns ``(score, parts)``.  ``score`` is NaN wherever any of the three gradients is
    undefined.  The middle factor ``|U_g + U_m| / 2`` is 1 when gravity and magnetics point
    the same way and 0 when they are orthogonal, so it *is* the cross-gradient information
    (|U_g + U_m|^2 = 2 + 2 cos(theta)), used as a scalar weight rather than as a product of
    two independent per-field scores.
    """
    gxt, gyt = g_t
    gxg, gyg = g_g
    gxm, gym = g_m
    mt = magnitude(gxt, gyt)
    mg = magnitude(gxg, gyg)
    mm = magnitude(gxm, gym)

    ok = valid & np.isfinite(mt) & np.isfinite(mg) & np.isfinite(mm)
    nt, st = robust_unit(mt, ok, pct)
    ng, sg = robust_unit(mg, ok, pct)
    nm, sm = robust_unit(mm, ok, pct)

    eps = np.float32(1e-12)
    ugx, ugy = gxg / np.maximum(mg, eps), gyg / np.maximum(mg, eps)
    umx, umy = gxm / np.maximum(mm, eps), gym / np.maximum(mm, eps)
    vx, vy = ugx + umx, ugy + umy
    agree_gm = np.hypot(vx, vy) / 2.0            # in [0,1]; 1 == parallel gravity & magnetics

    utx, uty = gxt / np.maximum(mt, eps), gyt / np.maximum(mt, eps)
    mref = np.hypot(vx, vy)
    cos_t = (utx * vx + uty * vy) / np.maximum(mref, eps)
    cos_t = np.clip(cos_t, -1.0, 1.0)

    score = nt * agree_gm * np.square(cos_t)
    score = np.where(ok, score, NAN32).astype(np.float32)
    parts = {
        "topo_mag_norm": np.where(ok, nt, NAN32).astype(np.float32),
        "agree_gm": np.where(ok, agree_gm, NAN32).astype(np.float32),
        "cos2_topo": np.where(ok, np.square(cos_t), NAN32).astype(np.float32),
        "gravity_norm": np.where(ok, ng, NAN32).astype(np.float32),
        "mag_norm": np.where(ok, nm, NAN32).astype(np.float32),
        "scale_topo": st, "scale_grav": sg, "scale_mag": sm,
    }
    return score, parts


def family(
    name: str,
    bands: dict[str, np.ndarray],
    footprint: np.ndarray,
    *,
    sigma: float = 2.0,
    pct: float = 95.0,
) -> np.ndarray:
    """Build one named detector family."""
    g = bands["iso_grav_anom"]
    m = bands["rtp"]
    e = bands["det_elev"]
    s19 = bands["det_elev_slope"]

    if name.startswith("single:"):
        return single_layer_edge(bands[name.split(":", 1)[1]], footprint, sigma)

    # ---- the brief's operator, on the two potential fields --------------------
    if name == "cg_gxm":
        return _coupled(_grad(g, sigma), _grad(m, sigma), footprint, footprint, pct)

    if name == "cg_gxm_s4":
        return _coupled(_grad(g, 4.0), _grad(m, 4.0), footprint, footprint, pct)

    if name == "cg_gxm_hg":
        gp, mp = bands["iso_grav_anom_hg"], bands["tmi_hg"]
        return _coupled(_grad(gp, sigma), _grad(mp, sigma), footprint & np.isfinite(gp),
                        footprint & np.isfinite(mp), pct)

    if name == "cg_gxm_vg":
        gz, mz = bands["iso_grav_anom_vg"], bands["tmi_vg"]
        return _coupled(_grad(gz, sigma), _grad(mz, sigma), footprint & np.isfinite(gz),
                        footprint & np.isfinite(mz), pct)

    # ---- topography-anchored coupling ----------------------------------------
    if name == "tg_g":
        return _coupled(_grad(s19, sigma), _grad(g, sigma), footprint, footprint, pct)

    if name == "tg_m":
        return _coupled(_grad(s19, sigma), _grad(m, sigma), footprint, footprint, pct)

    if name == "tg_gm_max":
        a = np.nan_to_num(_coupled(_grad(s19, sigma), _grad(g, sigma), footprint, footprint, pct), nan=0.0)
        b = np.nan_to_num(_coupled(_grad(s19, sigma), _grad(m, sigma), footprint, footprint, pct), nan=0.0)
        return np.where(footprint, np.maximum(a, b), NAN32)

    if name == "tg_gm_geo":
        a = np.nan_to_num(_coupled(_grad(s19, sigma), _grad(g, sigma), footprint, footprint, pct), nan=0.0)
        b = np.nan_to_num(_coupled(_grad(s19, sigma), _grad(m, sigma), footprint, footprint, pct), nan=0.0)
        return np.where(footprint, np.sqrt(a * b), NAN32)

    if name == "triple":
        score, _ = triple_alignment(_grad(s19, sigma), _grad(g, sigma), _grad(m, sigma), footprint, pct)
        return score

    if name == "triple_s4":
        score, _ = triple_alignment(_grad(s19, 4.0), _grad(g, 4.0), _grad(m, 4.0), footprint, pct)
        return score

    if name == "triple_elev":
        # same operator with the detrended elevation itself as the topographic carrier
        score, _ = triple_alignment(_grad(e, sigma), _grad(g, sigma), _grad(m, sigma), footprint, pct)
        return score

    # ---- controls the brief explicitly rules out -----------------------------
    if name == "product_of_scores":
        t = _norm(magnitude(*_grad(s19, sigma)), footprint, pct)
        gg = _norm(magnitude(*_grad(g, sigma)), footprint, pct)
        mm = _norm(magnitude(*_grad(m, sigma)), footprint, pct)
        return (t * np.maximum(np.nan_to_num(gg, nan=0.0), np.nan_to_num(mm, nan=0.0))).astype(np.float32)

    if name == "cg_align_only":
        out = pair_coupling(_grad(g, sigma), _grad(m, sigma), footprint, footprint, pct)
        return out["align"]

    if name == "cg_confounder_masked":
        out = pair_coupling(_grad(g, sigma), _grad(m, sigma), footprint, footprint, pct)
        align, ortho, bal = out["align"], out["ortho"], out["balance"]
        evidence = align / np.maximum(align + ortho, 1e-12)
        return np.where(out["valid"], bal * evidence, NAN32).astype(np.float32)

    raise KeyError(f"unknown family {name!r}")


CANDIDATE_FAMILIES = (
    "cg_gxm",
    "cg_gxm_s4",
    "cg_gxm_hg",
    "cg_gxm_vg",
    "tg_g",
    "tg_m",
    "tg_gm_max",
    "tg_gm_geo",
    "triple",
    "triple_s4",
    "triple_elev",
    "product_of_scores",
    "cg_align_only",
    "cg_confounder_masked",
)

CORE_BANDS = (
    "iso_grav_anom", "rtp", "det_elev", "det_elev_slope",
    "iso_grav_anom_vg", "iso_grav_anom_hg", "tmi_hg", "tmi_vg",
)
