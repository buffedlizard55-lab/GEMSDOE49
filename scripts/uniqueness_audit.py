"""Uniqueness audit — reject relabels before they go near a submission.

Per the session brief: hash the raster and compare rank correlation and
top-k overlap against EVERY prior submission; reject the candidate as a
relabel of an earlier file unless it is clearly different.

Criteria used here (preregistered in docs/hypotheses.html):
  * sha256 must not collide with any prior file (byte-identical);
  * Spearman rank correlation of the candidate's positive-pixel ranking vs
    each prior file's positive-pixel ranking, over the shared positive set,
    must be <= 0.90 for the candidate to be "clearly different" in rank
    structure; AND
  * top-k spatial overlap (Jaccard of the emitted-pixel sets, k matched to
    the smaller emission) must be <= 0.50.
A candidate failing BOTH thresholds against ANY prior file is flagged
REJECT-RELABEL.  Every pairwise number is written to evidence/.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import sys

import numpy as np
import rasterio
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from gems49 import data as D

ROOT = D.ROOT


def load_field(path: str) -> np.ndarray:
    with rasterio.open(path) as ds:
        a = ds.read(1).astype(np.float64)
    a[~np.isfinite(a)] = 0.0
    return a


def topk_set(a: np.ndarray, k: int) -> set[int]:
    """Flat indices of the k largest values (ties broken deterministically)."""
    flat = a.ravel()
    n_pos = int((flat > 0).sum())
    k = min(k, n_pos)
    if k <= 0:
        return set()
    idx = np.argpartition(-flat, k)[:k]
    return set(idx.tolist())


def audit(candidate_path: str, prior_dir: str | None = None,
          rho_max: float = 0.90, jaccard_max: float = 0.50) -> dict:
    prior_dir = prior_dir or os.path.join(ROOT, "data", "prior")
    cand = load_field(candidate_path)
    cand_sha = hashlib.sha256(open(candidate_path, "rb").read()).hexdigest()
    cand_pos = np.argwhere(cand > 0)
    print(f"candidate: {os.path.basename(candidate_path)}  sha256={cand_sha[:16]}  "
          f"positive px={len(cand_pos)}")
    rows = []
    reject = False
    reject_reason = []
    for p in sorted(glob.glob(os.path.join(prior_dir, "*.tif"))):
        name = os.path.basename(p)
        sha = hashlib.sha256(open(p, "rb").read()).hexdigest()
        row = {"file": name, "sha256": sha}
        if sha == cand_sha:
            row.update(byte_identical=True)
            reject = True
            reject_reason.append(f"byte-identical to {name}")
            rows.append(row)
            continue
        row.update(byte_identical=False)
        try:
            prev = load_field(p)
        except Exception as e:
            row.update(error=str(e)[:80])
            rows.append(row)
            continue
        if prev.shape != cand.shape:
            row.update(shape_mismatch=list(prev.shape))
            rows.append(row)
            continue
        # positive-pixel set overlap
        cset = set(np.flatnonzero(cand.ravel() > 0).tolist())
        pset = set(np.flatnonzero(prev.ravel() > 0).tolist())
        k = min(len(cset), len(pset))
        inter = len(cset & pset)
        jacc = inter / max(len(cset | pset), 1)
        jacc_k = inter / max(min(len(cset), len(pset)), 1)  # overlap of smaller set
        row.update(n_candidate=len(cset), n_prior=len(pset),
                   overlap_smaller_set=jacc_k, jaccard_union=jacc)
        # rank correlation over the union of positive pixels
        union = np.array(sorted(cset | pset), dtype=np.int64)
        if len(union) > 10:
            x = cand.ravel()[union]
            y = prev.ravel()[union]
            if np.std(x) > 0 and np.std(y) > 0:
                rho = float(spearmanr(x, y).statistic)
            else:
                rho = 0.0
            row.update(spearman_rho=rho)
            if rho > rho_max and jacc_k > jaccard_max:
                row.update(flag="REJECT-RELABEL")
                reject = True
                reject_reason.append(f"rho={rho:.3f} and top-k overlap={jacc_k:.3f} vs {name}")
            else:
                row.update(flag="distinct")
        rows.append(row)
    out = {
        "candidate": os.path.basename(candidate_path),
        "candidate_sha256": cand_sha,
        "criteria": {"rho_max": rho_max, "jaccard_max": jaccard_max},
        "n_priors_compared": len(rows),
        "verdict": "REJECT-RELABEL" if reject else "UNIQUE",
        "reject_reasons": reject_reason,
        "priors": rows,
    }
    os.makedirs(os.path.join(ROOT, "evidence"), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, "evidence", "uniqueness_audit.json"), "w"), indent=1)
    print(f"verdict: {out['verdict']}  ({len(rows)} priors compared)")
    if reject_reason:
        print("reasons:", reject_reason)
    return out


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else sorted(
        glob.glob(os.path.join(ROOT, "docs", "downloads", "gemsdoe49-*-zeros.tif")))[-1]
    audit(path)
