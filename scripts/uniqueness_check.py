"""Reject a candidate as a relabel: compare it against every prior submission on disk.

The brief requires that a candidate be "clearly different" from every prior submission.  Two
independent measures are used, both computed on the prediction support (the >0 pixels), which is
what the metric actually sees:

* **Spearman rank correlation** between the two rasters over pixels where either is positive
  (ties from binary rasters are unavoidable, so the rank correlation is reported alongside the
  more interpretable overlap measures);
* **top-k Jaccard overlap** for k = 5,000 / 10,000 / 25,000, i.e. the agreement of the two
  files' highest-value pixels;
* **support Jaccard** over all positive pixels.

Run: python scripts/uniqueness_check.py docs/downloads/<candidate>.tif --prior-dir /path/to/tifs
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
TOPK = (5_000, 10_000, 25_000)


def load(path: Path) -> np.ndarray:
    with rasterio.open(path) as s:
        return np.nan_to_num(s.read(1).astype(np.float32), nan=0.0)


def compare(a: np.ndarray, b: np.ndarray) -> dict:
    pa, pb = a > 0, b > 0
    inter = int((pa & pb).sum())
    union = int((pa | pb).sum())
    both = pa | pb
    n = int(both.sum())
    rho = float("nan")
    if n > 3 and a[both].std() > 0 and b[both].std() > 0:
        rho = float(spearmanr(a[both], b[both]).statistic)
    out = {
        "n_pos_a": int(pa.sum()), "n_pos_b": int(pb.sum()),
        "support_intersection": inter, "support_union": union,
        "support_jaccard": round(inter / union, 8) if union else 0.0,
        "spearman_on_union": round(rho, 6) if np.isfinite(rho) else None,
    }
    for k in TOPK:
        ta = np.argpartition(-a.ravel(), min(k, int(pa.sum()) or 1) - 1)[:k] if pa.sum() else np.array([], int)
        tb = np.argpartition(-b.ravel(), min(k, int(pb.sum()) or 1) - 1)[:k] if pb.sum() else np.array([], int)
        sa, sb = np.zeros(a.size, bool), np.zeros(b.size, bool)
        sa[ta] = True
        sb[tb] = True
        u = int((sa | sb).sum())
        out[f"top{k}_jaccard"] = round(int((sa & sb).sum()) / u, 8) if u else 0.0
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate")
    ap.add_argument("--prior-dir", action="append", default=[])
    ap.add_argument("--out", default=None)
    ap.add_argument("--jaccard-threshold", type=float, default=0.10)
    ap.add_argument("--spearman-threshold", type=float, default=0.60)
    args = ap.parse_args()

    cand_path = Path(args.candidate)
    cand = load(cand_path)
    cand_sha = hashlib.sha256(cand_path.read_bytes()).hexdigest()

    priors: list[Path] = []
    for d in args.prior_dir:
        priors.extend(sorted(Path(d).glob("*.tif")))
    priors = [p for p in priors if p.resolve() != cand_path.resolve()]

    rows = []
    for p in priors:
        try:
            b = load(p)
        except Exception:  # noqa: BLE001
            continue
        if b.shape != cand.shape:
            continue
        c = compare(cand, b)
        c["file"] = p.name
        c["sha256_16"] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        rows.append(c)

    if not rows:
        print("no comparable prior rasters found")
        return 2

    # a prior submission that is byte-identical or pixel-identical is an outright relabel
    identical = [r for r in rows if r["support_jaccard"] == 1.0 and r["n_pos_a"] == r["n_pos_b"]]
    worst_top = max(rows, key=lambda r: r["top25000_jaccard"])
    worst_sup = max(rows, key=lambda r: r["support_jaccard"])
    worst_rho = max(rows, key=lambda r: (r["spearman_on_union"] or -1))
    mean_top = float(np.mean([r["top25000_jaccard"] for r in rows]))
    mean_sup = float(np.mean([r["support_jaccard"] for r in rows]))

    verdict = "UNIQUE"
    reasons = []
    if identical:
        verdict = "REJECT"
        reasons.append(f"pixel-identical to {identical[0]['file']}")
    if worst_top["top25000_jaccard"] > args.jaccard_threshold:
        verdict = "REJECT"
        reasons.append(f"top-25000 Jaccard {worst_top['top25000_jaccard']:.4f} exceeds "
                       f"{args.jaccard_threshold} vs {worst_top['file']}")
    if (worst_rho["spearman_on_union"] or 0) > args.spearman_threshold:
        verdict = "REJECT"
        reasons.append(f"Spearman {worst_rho['spearman_on_union']:.4f} exceeds "
                       f"{args.spearman_threshold} vs {worst_rho['file']}")

    top20 = sorted(rows, key=lambda r: -r["top25000_jaccard"])[:20]
    report = {
        "candidate": cand_path.name, "candidate_sha256": cand_sha,
        "n_priors_compared": len(rows),
        "candidate_positives": int((cand > 0).sum()),
        "mean_top25000_jaccard": round(mean_top, 6),
        "mean_support_jaccard": round(mean_sup, 6),
        "max_top25000_jaccard": worst_top["top25000_jaccard"],
        "max_top25000_jaccard_file": worst_top["file"],
        "max_support_jaccard": worst_sup["support_jaccard"],
        "max_support_jaccard_file": worst_sup["file"],
        "max_spearman": worst_rho["spearman_on_union"],
        "max_spearman_file": worst_rho["file"],
        "verdict": verdict, "reasons": reasons,
        "closest_prior_submissions": [
            {k: r[k] for k in ("file", "n_pos_b", "support_jaccard", "top25000_jaccard",
                               "spearman_on_union")} for r in top20],
    }
    print(json.dumps(report, indent=1))
    if args.out:
        Path(args.out).write_text(json.dumps(report, indent=1))
    return 0 if verdict == "UNIQUE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
