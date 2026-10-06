"""Spatially-blocked hide-and-recover holdout.

Protocol (preregistered in docs/hypotheses.html, following the family's
blocked-holdout convention used by GEMSDOE4..47):

* The only labelled pixels we possess are the public USGS/INGENIOUS
  catalogue.  In live scoring those pixels are MASKED, so they can never
  earn credit directly; the holdout therefore *hides* whole contiguous
  segments of the catalogue, scores candidate emissions against the hidden
  pixels with the official DTI, and masks the residual catalogue.
* Spatial blocking: the footprint is split into 4 quadrants
  (NW/NE/SW/SE).  In each fold, ~25% of the catalogue connected components
  inside the fold's quadrant are hidden (components are kept intact so the
  hidden truth keeps fault-like geometry: "hide-and-recover segments").
* A candidate emission rule must beat the single-layer gradient baseline
  on the mean of the 4 folds before it may be submitted.

Limitation, stated rather than smoothed over: the hidden truth is drawn
from the catalogue, so this instrument measures *placement quality under
the assumption that new faults geometrically resemble catalogue faults*.
It cannot reward a genuinely novel fault family (family registry
IR-32-PROXY-01).  Live-score evidence from the family's scored probes is
combined with this instrument, never replaced by it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from .metric import dti, DTIResult


def quadrant_fold(shape: tuple[int, int]) -> np.ndarray:
    """4-quadrant fold id (0=NW, 1=NE, 2=SW, 3=SE)."""
    h, w = shape
    fold = np.zeros(shape, dtype=np.int8)
    hh, hw = h // 2, w // 2
    fold[:hh, hw:] = 1
    fold[hh:, :hw] = 2
    fold[hh:, hw:] = 3
    return fold


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
