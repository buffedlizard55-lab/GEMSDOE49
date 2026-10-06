"""Distance-weighted Tversky index (DTI) — exact re-implementation of the
official GEMS Prize metric.

Official specification (verified against the problem description,
https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/):

    k(d)  = max(1 - d/R, 0)                     triangular kernel, R = 300 m
    TP_w  = sum_{g in G} max_{x: d(x,g)<=R} p(x) k(d(x,g))
    FP_w  = sum_{x: p(x)>0} p(x) [1 - max_{g in G} k(d(x,g))]
    FN_w  = sum_{g in G} [1 - max_{x: d(x,g)<=R} p(x) k(d(x,g))]
    DTI   = TP_w / (TP_w + alpha*FP_w + beta*FN_w + eps)

with alpha = 0.2, beta = 0.8, R = 300 m = 3 pixels at 100 m resolution.
Distances are Euclidean between pixel centres.

Implementation notes (pinned against the organiser's reference notebook,
https://github.com/drivendataorg/gems-prize-reference-solution, as extracted
in the family repo GEMSDOE3 `evidence/training-source/metric.py`):

* the kernel offset set is { (dy,dx) : hypot(dy,dx) < R }  (strict `< R`;
  at exactly d = R the weight is 0, so inclusive/exclusive is immaterial);
* FP_w uses 1 - k(d to NEAREST truth pixel) = min(d/R, 1);
* a scoring mask removes catalogue pixels from all three sums (organiser
  staff, forum thread 11516 post 4: the mask is pixel-exact and identical
  to the provided training fault labels; near-known-trace predictions far
  from new-fault truth are fully penalised).

Algebraic identity used for cross-checks:
    TP_w + FN_w = |G_masked|  and  DTI = T / (0.2(T+S-M) + 0.8|G| + eps)
where T = TP_w, S = sum of p over unmasked predicted pixels and
M = sum of p*k(nearest truth) over those pixels.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt

ALPHA = 0.2
BETA = 0.8
R_PX = 3          # 300 m at 100 m resolution (integer pixel radius)
EPS = 1e-7        # same epsilon as the reference implementation


@dataclass
class DTIResult:
    dti: float
    tp_w: float
    fp_w: float
    fn_w: float
    n_truth: int
    n_pred_positive: int
    pred_mass: float

    def as_dict(self) -> dict:
        return {
            "dti": self.dti,
            "tp_w": self.tp_w,
            "fp_w": self.fp_w,
            "fn_w": self.fn_w,
            "n_truth": self.n_truth,
            "n_pred_positive": self.n_pred_positive,
            "pred_mass": self.pred_mass,
        }


def _offsets(radius: int) -> list[tuple[int, int, float]]:
    return [
        (dy, dx, 1.0 - math.hypot(dy, dx) / radius)
        for dy in range(-radius, radius + 1)
        for dx in range(-radius, radius + 1)
        if math.hypot(dy, dx) < radius
    ]


def dti(
    pred: np.ndarray,
    truth: np.ndarray,
    mask: np.ndarray | None = None,
    alpha: float = ALPHA,
    beta: float = BETA,
    radius: int = R_PX,
) -> DTIResult:
    """Compute the distance-weighted Tversky index.

    Parameters
    ----------
    pred : (H, W) float array with finite values in [0, 1].
    truth : (H, W) boolean / 0-1 array — the held-out ground truth.
    mask : (H, W) bool array, True where the pixel is SCORED (unmasked).
        Pixels with mask False are removed from TP/FP/FN entirely.
        Convention matches the competition: catalogue pixels are masked.
    """
    p = np.asarray(pred, dtype=np.float64)
    g = np.asarray(truth).astype(bool)
    if p.shape != g.shape:
        raise ValueError("pred/truth shape mismatch")
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("DTI expects finite predictions in [0,1]")
    m = np.ones(p.shape, bool) if mask is None else np.asarray(mask, bool)
    if m.shape != p.shape:
        raise ValueError("mask shape mismatch")
    # Masked pixels are excluded from evaluation entirely (staff, forum 11516:
    # perfect catalogue reproduction scores zero), so they can neither earn
    # TP credit by proximity nor pay FP cost.
    p = np.where(m, p, 0.0)

    rows, cols = np.nonzero(g)
    select = m[rows, cols]
    rows, cols = rows[select], cols[select]
    n_truth = rows.size

    if n_truth == 0:
        active = m & (p > 0)
        fp = float(p[active].sum())  # no truth => k=0 everywhere => weight 1
        return DTIResult(0.0, 0.0, fp, 0.0, 0, int(active.sum()), float(p[m].sum()))

    # FP kernel uses ONLY the unmasked truth pixels: organiser staff (forum
    # 11516 post 4) confirmed predictions near known faults but far from
    # new-fault truth are FULLY penalised — the residual catalogue does not
    # absorb false-positive mass.
    g_eff = np.zeros_like(g)
    g_eff[rows, cols] = True
    fp_weight = np.minimum(distance_transform_edt(~g_eff) / radius, 1.0)

    h, w = p.shape
    credit = np.zeros(n_truth, dtype=np.float64)
    for dy, dx, k in _offsets(radius):
        rr, cc = rows + dy, cols + dx
        good = (rr >= 0) & (rr < h) & (cc >= 0) & (cc < w)
        credit[good] = np.maximum(credit[good], p[rr[good], cc[good]] * k)

    tp = float(credit.sum())
    fn = float(n_truth) - tp
    active = m & (p > 0)
    fp = float((p[active] * fp_weight[active]).sum())
    dti_val = tp / (tp + alpha * fp + beta * fn + EPS)
    return DTIResult(dti_val, tp, fp, fn, int(n_truth), int(active.sum()), float(p[m].sum()))


def dti_identity_check(res: DTIResult, sum_p_unmasked: float, mp_sum: float, alpha: float = ALPHA, beta: float = BETA) -> float:
    """Return DTI recomputed via T / (alpha(T+S-M) + beta|G|) for auditing."""
    return res.tp_w / (alpha * (res.tp_w + sum_p_unmasked - mp_sum) + beta * res.n_truth + EPS)
