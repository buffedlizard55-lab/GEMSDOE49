"""The official distance-weighted Tversky index (DOE GEMS Prize / DrivenData #306).

Transcribed from the published equations on
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/

  k(d)   = max(1 - d/R, 0),                R = 300 m
  TP_w   = sum_{g in G} max_{x: d(x,g)<=R} p(x) * k(d(x,g))
  FP_w   = sum_{x: p(x)>0} p(x) * [1 - max_{g in G} k(d(x,g))]
  FN_w   = sum_{g in G} [1 - max_{x: d(x,g)<=R} p(x) * k(d(x,g))]
  DTI    = TP_w / (TP_w + alpha*FP_w + beta*FN_w + eps),   alpha = 0.2, beta = 0.8

Three properties are asserted in the test-suite rather than assumed:

1. the organizers' own worked example  (TP_w 3.00, FP_w 1.89, FN_w 2.00 -> 0.602652),
2. the identity  TP_w + FN_w = |G|,
3. agreement with a literal O(|G| x |X|) double-loop transcription on small rasters.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

from .spec import ALPHA, BETA, R_M

EPS = 1e-9


def kernel_offsets(res_m: float = 100.0, r_m: float = R_M) -> list[tuple[int, int, float]]:
    """Integer pixel offsets (dy, dx) whose Euclidean distance is <= r_m, with their weights.

    Ties are resolved by ``d <= r_m``, which is the published condition
    ``d(x,g) <= R`` (so a pixel at exactly 300 m keeps weight 0.0 and is harmless).
    """
    r_px = r_m / res_m
    n = int(np.floor(r_px))
    out: list[tuple[int, int, float]] = []
    for dy in range(-n, n + 1):
        for dx in range(-n, n + 1):
            d = float(np.hypot(dx, dy)) * res_m
            if d <= r_m + 1e-9:
                out.append((dy, dx, d))
    return out


def triangle_kernel(d_m: np.ndarray | float, r_m: float = R_M) -> np.ndarray | float:
    """k(d) = max(1 - d/R, 0)."""
    return np.maximum(1.0 - np.asarray(d_m, dtype=np.float64) / r_m, 0.0)


def _tp_w(pred: np.ndarray, truth: np.ndarray, offsets, r_m: float) -> float:
    """TP_w = sum over truth pixels of the best weighted prediction within R."""
    idx = np.argwhere(truth)
    if idx.size == 0:
        return 0.0
    h, w = pred.shape
    best = np.zeros(idx.shape[0], dtype=np.float64)
    for dy, dx, d in offsets:
        weight = max(0.0, 1.0 - d / r_m)
        if weight <= 0.0:
            continue
        ys = idx[:, 0] + dy
        xs = idx[:, 1] + dx
        ok = (ys >= 0) & (ys < h) & (xs >= 0) & (xs < w)
        vals = np.zeros(idx.shape[0], dtype=np.float64)
        vals[ok] = pred[ys[ok], xs[ok]] * weight
        np.maximum(best, vals, out=best)
    return float(best.sum())


def _fp_w(pred: np.ndarray, truth: np.ndarray, res_m: float, r_m: float) -> float:
    """FP_w = sum over predicted pixels of p(x) * (1 - max_g k(d(x,g))).

    ``max_g k(d(x,g))`` equals ``max(1 - d_min(x)/R, 0)`` with ``d_min`` the Euclidean
    distance in metres to the nearest truth pixel, so one exact distance transform
    reproduces the published per-pixel maximisation.
    """
    d_min = distance_transform_edt(~truth, sampling=(res_m, res_m))
    best_k = np.maximum(1.0 - d_min / r_m, 0.0)
    return float(np.sum(pred * (1.0 - best_k)))


def dti_components(
    pred: np.ndarray,
    truth: np.ndarray,
    res_m: float = 100.0,
    r_m: float = R_M,
    alpha: float = ALPHA,
    beta: float = BETA,
) -> dict:
    """Return TP_w, FP_w, FN_w and the DTI for a continuous prediction against a truth mask.

    ``pred`` is treated as-is: any pixel whose value is NaN is replaced by 0, and no
    thresholding is applied (the competition scores continuous confidence values).
    """
    pred = np.asarray(pred, dtype=np.float64)
    pred = np.where(np.isfinite(pred), np.clip(pred, 0.0, None), 0.0)
    truth = np.asarray(truth, dtype=bool)
    n_g = int(truth.sum())
    offsets = kernel_offsets(res_m, r_m)
    tp = _tp_w(pred, truth, offsets, r_m)
    fp = _fp_w(pred, truth, res_m, r_m)
    fn = float(n_g) - tp
    denom = tp + alpha * fp + beta * fn + EPS
    return {
        "TPw": tp,
        "FPw": fp,
        "FNw": fn,
        "n_truth": n_g,
        "n_pred_px": int((pred > 0).sum()),
        "DTI": tp / denom if denom > 0 else 0.0,
        "alpha": alpha,
        "beta": beta,
        "R_m": r_m,
    }


def dti_from_components(tp: float, fp: float, n_truth: int, alpha: float = ALPHA, beta: float = BETA) -> float:
    """DTI from (TP_w, FP_w, |G|) using the exact algebraic identity FN_w = |G| - TP_w."""
    fn = n_truth - tp
    return tp / (tp + alpha * fp + beta * fn + EPS)


def marginal_bar(dti: float, alpha: float = ALPHA) -> float:
    """Minimum kernel credit a newly added unit prediction pixel must deliver to raise the score.

    Exact derivation.  Because ``alpha + beta = 1`` and ``FN_w = |G| - TP_w``, the denominator is
    ``D = alpha*(TPw + FPw) + beta*|G|``.  Adding one pixel of unit weight whose credit to its
    nearest truth pixel is ``k`` increases ``TP_w`` by ``k``, its false-positive term by ``1 - k``
    (it is at distance ``R*(1-k)`` from that truth pixel, so ``1 - max_g k(d) = 1 - k``), and both
    ``|G|`` and ``alpha`` are unchanged.  Hence ``Delta(TPw + FPw) = 1`` and ``Delta D = alpha``
    *exactly*, for every ``k``.  Then

        (TPw + k)/(D + alpha) > TPw/D   <=>   k > alpha * TPw/D = alpha * DTI.

    So the bar is ``alpha * DTI`` -- 0.06 at DTI = 0.30, i.e. every pixel within
    ``300*(1-0.06) = 282 m`` of a scored truth pixel pays for itself.

    Note: a sibling repository's analysis states the bar as ``alpha*s/(1 - alpha*s)`` (0.0683 at
    s = 0.3195).  That expression is the correct bar for adding *area* whose credit is delivered
    in bulk; for a single added pixel the exact value is ``alpha*s`` (0.0639 at s = 0.3195).  The
    two imply the same practical distance (281 m vs 280 m) but this module uses the exact one,
    and ``test_marginal_bar_is_exact`` proves it by direct computation on a raster.
    """
    return alpha * dti


def marginal_max_distance(dti: float, res_m: float = 100.0, r_m: float = R_M, alpha: float = ALPHA) -> float:
    """Largest distance (metres) at which a lone added pixel still pays at the given DTI."""
    return r_m * (1.0 - marginal_bar(dti, alpha))


def scoring_domain(pred: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Zero prediction mass on masked pixels.

    DrivenData staff (forum topic 11516) confirmed that pixels corresponding to known
    USGS/INGENIOUS faults are excluded from evaluation; prediction mass there is therefore
    removed *before* the kernel is applied, so it neither earns nor radiates credit.
    """
    out = np.array(pred, dtype=np.float64, copy=True)
    out[mask] = 0.0
    return out
