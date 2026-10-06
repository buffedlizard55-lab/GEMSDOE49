#!/usr/bin/env python3
"""Check whether rasterized 8-connected label components cross spatial fold cores."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import label as connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe49.holdout import spatial_column_folds  # noqa: E402
from gemsdoe49.raster import load_inputs, sha256_file  # noqa: E402


def summarize_connectivity(labels: np.ndarray, folds: list[dict], connectivity: int) -> dict:
    structure = np.ones((3, 3), dtype=bool) if connectivity == 8 else None
    components, count = connected_components(labels, structure=structure)
    sizes = np.bincount(components.ravel(), minlength=count + 1)[1:]
    hits = np.zeros((count + 1, len(folds)), dtype=bool)
    scored_sizes = np.zeros((count + 1, len(folds)), dtype=np.int64)
    for fold in folds:
        fold_id = int(fold["fold_id"])
        idx = components[fold["domain"] & labels]
        bincounts = np.bincount(idx, minlength=count + 1)
        hits[:, fold_id] = bincounts > 0
        scored_sizes[:, fold_id] = bincounts
    fold_counts = hits[1:].sum(axis=1)
    cross = fold_counts > 1
    return {
        "connectivity": f"{connectivity}-connected" if connectivity == 8 else "4-connected",
        "component_count": int(count),
        "components_gt_one_pixel": int((sizes > 1).sum()),
        "largest_component_pixels": int(sizes.max()) if sizes.size else 0,
        "p99_component_pixels": float(np.percentile(sizes, 99)) if sizes.size else 0.0,
        "components_with_pixels_in_multiple_scored_fold_cores": int(cross.sum()),
        "pixels_from_crossing_components_inside_scored_cores": int(scored_sizes[1:][cross].sum()),
        "fold_membership_counts": {
            str(n): int((fold_counts == n).sum()) for n in range(0, len(folds) + 1)
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / ".cache" / "gems_data")
    parser.add_argument("--out", type=Path, default=ROOT / "research" / "label_fold_component_audit_20261006.json")
    args = parser.parse_args()

    _, labels, footprint, meta = load_inputs(args.data_dir)
    folds = spatial_column_folds(labels.shape, footprint, n_folds=5, boundary_guard_px=3)
    receipt = {
        "audit_id": "GEMSDOE49-label-fold-component-audit",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "label_sha256": meta["labels_sha256"],
        "label_positive_pixels": int(labels.sum()),
        "fold_protocol": "five equal-width east-west cores; three-pixel guard between adjacent cores; three-pixel north/south edge trim",
        "fold_core_bounds": [
            {
                "fold_id": int(f["fold_id"]),
                "row_start": int(f["row_start"]), "row_end": int(f["row_end"]),
                "col_start": int(f["col_start"]), "col_end": int(f["col_end"]),
            }
            for f in folds
        ],
        "interpretation_limit": "Raster connected components are a diagnostic approximation, not geological fault-segment identifiers. Rasterization can split one mapped line or merge crossing/nearby lines.",
        "components": [
            summarize_connectivity(labels, folds, 4),
            summarize_connectivity(labels, folds, 8),
        ],
        "metric_note": "Component pixels are counted only inside scored fold cores; the excluded 300 m boundaries are not assigned to a fold.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt["components"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
