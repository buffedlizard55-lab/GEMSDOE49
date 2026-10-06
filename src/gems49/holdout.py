
"""Hide-and-recover holdout instruments.

Why this file exists
--------------------
The competition's truth is a set of faults that a *reviewing expert added to a map that
already contained the USGS/INGENIOUS catalogue*, and DrivenData staff confirmed on the
forum (topic 11516) that pixels corresponding to the known catalogue are masked out of
evaluation in both prize rounds.  A holdout that scores a model on the catalogue therefore
measures the wrong thing: it rewards re-drawing traces the scorer never looks at.

The instrument built here mirrors the real one as closely as the public data allows:

* the catalogue is cut into connected **segments** and whole segments are removed;
* the retained segments play the role of "the map the expert was given" and are **zeroed
  out of the prediction before scoring**, exactly as the organizer masks them;
* the removed segments play the role of "the faults the map lacked" and are the only truth;
* two variants of the truth set are reported, because the two answer different questions:

  - ``unrestricted`` -- every removed-segment pixel is truth.  A prediction placed on the
    retained catalogue cannot leak (it is zeroed), but a prediction placed *beside* the
    retained catalogue can still earn kernel credit for a removed neighbour.  This variant
    therefore scores the "hidden faults sit next to mapped faults" hypothesis favourably.
  - ``collared`` -- only removed pixels at least ``sep_m`` away from the retained catalogue
    are truth.  This deletes the adjacency channel entirely, so it is the conservative,
    leakage-free instrument and the one used as the promotion gate.

Controls are measured, not assumed: a constant-zero prediction and a prediction equal to
the retained catalogue itself are both scored so that the reader can see the instrument's
floor and its leakage floor.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt, label

from .metric import ALPHA, BETA, R_M, dti_components, kernel_offsets
from dataclasses import dataclass

from .metric import DTIResult, dti  # compatibility layer
from .spec import GRID

EIGHT = np.ones((3, 3), dtype=bool)


def segment_catalogue(catalogue: np.ndarray, min_px: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """8-connected segmentation of the catalogue raster.

    Returns ``(seg_id, sizes)`` where ``seg_id`` is 0 outside the catalogue and 1..n inside,
    and ``sizes[i]`` is the pixel count of segment ``i+1``.  Segments smaller than ``min_px``
    are dropped from the labelling (they are noise specks, not faults).
    """
    seg, n = label(catalogue, structure=EIGHT)
    sizes = np.bincount(seg.ravel(), minlength=n + 1)
    keep = np.zeros(n + 1, dtype=bool)
    keep[1:] = sizes[1:] >= min_px
    remap = np.zeros(n + 1, dtype=np.int32)
    remap[keep] = np.arange(1, int(keep.sum()) + 1)
    seg = remap[seg]
    return seg, np.bincount(seg.ravel(), minlength=int(seg.max()) + 1)[1:]


def fold_assignment(seg: np.ndarray, n_seg: int, n_folds: int, seed: int = 0) -> np.ndarray:
    """Assign whole segments to folds, balancing total segment length across folds."""
    sizes = np.bincount(seg.ravel(), minlength=n_seg + 1)[1:]
    order = np.argsort(-sizes, kind="stable")
    rng = np.random.default_rng(seed)
    # shuffled within equal-size groups so the split is deterministic and balanced
    order = order[rng.permutation(len(order))]
    load = np.zeros(n_folds)
    fold_of_seg = np.zeros(n_seg, dtype=np.int32)
    for s in order:
        f = int(np.argmin(load))
        fold_of_seg[s] = f
        load[f] += sizes[s]
    return fold_of_seg


def build_fold(
    catalogue: np.ndarray,
    footprint: np.ndarray,
    fold_of_seg: np.ndarray,
    seg: np.ndarray,
    k: int,
    *,
    sep_m: float = 600.0,
    res_m: float = GRID["res_m"],
) -> dict:
    """Build one hide-and-recover fold.

    Returns a dict with ``train_mask`` (retained catalogue, zeroed from predictions),
    ``truth_unrestricted``, ``truth_collared``, ``valid`` (prediction domain) and
    ``scored_unrestricted`` / ``scored_collared`` (truth pixels also inside the footprint).
    """
    seg_fold = np.zeros(seg.max() + 1, dtype=np.int32)
    seg_fold[1:] = fold_of_seg
    held = (seg > 0) & (seg_fold[seg] == k)
    train_mask = catalogue & ~held

    dist_train = distance_transform_edt(~train_mask, sampling=(res_m, res_m))
    truth_un = held & footprint
    truth_co = held & footprint & (dist_train >= sep_m)

    return {
        "fold": k,
        "train_mask": train_mask,
        "held_mask": held,
        "truth_unrestricted": truth_un,
        "truth_collared": truth_co,
        "valid": footprint,
        "n_train_px": int(train_mask.sum()),
        "n_held_px": int(held.sum()),
        "n_truth_unrestricted": int(truth_un.sum()),
        "n_truth_collared": int(truth_co.sum()),
        "collared_fraction": float(truth_co.sum() / max(truth_un.sum(), 1)),
    }


def score_fold(
    score_field: np.ndarray,
    fold: dict,
    *,
    truth_key: str = "truth_collared",
    res_m: float = GRID["res_m"],
) -> dict:
    """Score a continuous prediction field on one fold with the official metric.

    Prediction mass on ``train_mask`` is removed before scoring, which is the organizer's
    masking rule; the mask is also removed from the scoreable background so that a zero
    prediction there cannot be charged as a false negative.
    """
    pred = np.array(score_field, dtype=np.float64, copy=True)
    pred[~fold["valid"]] = 0.0
    pred[fold["train_mask"]] = 0.0
    truth = fold[truth_key]
    # pixels in the retained catalogue are not scoreable at all
    return dti_components(pred, truth, res_m=res_m)


class FoldScorer:
    """Cached scorer for one fold and one truth set.

    Everything that depends only on the truth and the mask (the truth pixel indices and the
    per-pixel ``max_g k(d(x,g))`` field) is computed once; scoring a candidate support then
    costs one gather per kernel offset plus one multiply-add over the grid.
    """

    def __init__(self, fold: dict, truth_key: str = "truth_collared", res_m: float = GRID["res_m"]):
        self.fold = fold
        self.truth_key = truth_key
        self.res_m = res_m
        self.pred_mask = ~fold["valid"] | fold["train_mask"]
        self.truth = fold[truth_key] & fold["valid"]
        self.idx = np.argwhere(self.truth).astype(np.int32)
        d_min = distance_transform_edt(~self.truth, sampling=(res_m, res_m))
        self.best_k = np.maximum(1.0 - d_min / R_M, 0.0).astype(np.float32)
        del d_min
        self._truth_mask_size = int(self.truth.sum())
        self.n_truth = int(self.truth.sum())
        self.offsets = [(dy, dx, max(0.0, 1.0 - d / R_M)) for dy, dx, d in kernel_offsets(res_m, R_M)]
        self.h, self.w = fold["valid"].shape
        del self.truth  # only the index list and best_k are needed from here on

    def score(self, support: np.ndarray, alpha: float = ALPHA, beta: float = BETA) -> dict:
        pred = np.where(support > 0, np.asarray(support, dtype=np.float64), 0.0)
        if self.pred_mask.any():
            pred[self.pred_mask] = 0.0
        best = np.zeros(self.idx.shape[0], dtype=np.float64)
        for dy, dx, wgt in self.offsets:
            if wgt <= 0.0:
                continue
            ys = self.idx[:, 0] + dy
            xs = self.idx[:, 1] + dx
            ok = (ys >= 0) & (ys < self.h) & (xs >= 0) & (xs < self.w)
            vals = np.zeros(self.idx.shape[0], dtype=np.float64)
            vals[ok] = pred[ys[ok], xs[ok]] * wgt
            np.maximum(best, vals, out=best)
        tp = float(best.sum())
        fp = float(np.sum(pred * (1.0 - self.best_k.astype(np.float64))))
        fn = float(self.n_truth) - tp
        denom = tp + alpha * fp + beta * fn + 1e-9
        return {"TPw": tp, "FPw": fp, "FNw": fn, "n_truth": self.n_truth,
                "n_pred_px": int((support > 0).sum()), "DTI": tp / denom if denom > 0 else 0.0}


def scorable_domain(fold: dict) -> np.ndarray:
    """Pixels where a prediction is allowed: inside the footprint and off the retained map."""
    return fold["valid"] & ~fold["train_mask"]


def leakage_probe(fold: dict, res_m: float = GRID["res_m"]) -> dict:
    """Score the retained catalogue itself as a prediction.

    Because ``score_fold`` zeroes the retained catalogue, this must return exactly 0.0
    TP_w.  If it does not, the instrument is leaking and every number from it is void.
    """
    field = fold["train_mask"].astype(np.float64)
    return score_fold(field, fold, res_m=res_m)


def null_probe(fold: dict, res_m: float = GRID["res_m"]) -> dict:
    """Score an all-zero prediction: the instrument's floor."""
    return score_fold(np.zeros_like(fold["valid"], dtype=np.float64), fold, res_m=res_m)


# ---------------------------------------------------------------------------
# Compatibility layer: the quadrant hide-and-recover folds from the earlier
# (merged) session, kept so its scripts, receipts and reports still run.
# ---------------------------------------------------------------------------

@dataclass
class Fold:
    fold_id: int
    truth: np.ndarray          # hidden catalogue pixels (bool)
    score_mask: np.ndarray     # True where evaluation is live (residual catalogue masked)
    n_hidden: int


def hide_and_recover_folds(
    catalogue: np.ndarray,
    footprint: np.ndarray,
    hide_fraction: float = 0.25,
    seed: int = 49,
) -> list[Fold]:
    """Build the 4 spatially-blocked hide-and-recover folds.

    Components are 8-connected components of the catalogue.  In each fold we
    hide a deterministic pseudo-random subset of the components whose
    centroids fall in that fold's quadrant, selecting components until at
    least `hide_fraction` of the fold's catalogue pixels are hidden.
    """
    rng = np.random.default_rng(seed)
    lab, n = ndimage.label(catalogue, structure=np.ones((3, 3), dtype=bool))
    comps = np.arange(1, n + 1)
    # component centroids and pixel counts
    sums = ndimage.sum(np.ones_like(lab), lab, index=comps)
    rc = ndimage.center_of_mass(np.ones_like(lab), lab, index=comps)
    rc = np.asarray(rc)
    fold_of_comp = quadrant_fold(catalogue.shape)[rc[:, 0].round().astype(int),
                                                  rc[:, 1].round().astype(int)]
    folds = []
    for f in range(4):
        comp_ids = comps[fold_of_comp == f]
        if len(comp_ids) == 0:
            continue
        sizes = dict(zip(comps, sums))
        total_px = sum(sizes[c] for c in comp_ids)
        order = rng.permutation(comp_ids)
        hidden_px, hidden_set = 0, []
        for c in order:
            hidden_set.append(c)
            hidden_px += sizes[c]
            if hidden_px >= hide_fraction * total_px:
                break
        # all pixels of the hidden components form the recover-truth;
        # components were selected by centroid-in-quadrant, so the truth is
        # (up to boundary-straddling components) inside this fold.
        truth = np.isin(lab, hidden_set)
        score_mask = footprint & ~(catalogue & ~truth)  # residual catalogue masked
        folds.append(Fold(f, truth, score_mask, int(truth.sum())))
    return folds


def score_candidate(pred: np.ndarray, folds: list[Fold]) -> dict:
    """Mean DTI of a candidate prediction raster over the folds."""
    results: list[DTIResult] = [dti(pred, f.truth, mask=f.score_mask) for f in folds]
    per_fold = [r.dti for r in results]
    return {
        "dti_mean": float(np.mean(per_fold)),
        "dti_per_fold": per_fold,
        "tp_w": [r.tp_w for r in results],
        "fp_w": [r.fp_w for r in results],
        "fn_w": [r.fn_w for r in results],
        "n_truth": [r.n_truth for r in results],
        "n_pred_positive": [r.n_pred_positive for r in results],
        "pred_mass": [r.pred_mass for r in results],
    }
