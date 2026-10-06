"""Five-fold spatial hide-and-recover screen on the public known-fault catalogue.

This is a falsification screen, not an estimate of the private test score.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from .metric import ALPHA, BETA, MetricComponents, dti_components, dti_from_components
from .placement import greedy_kernel_suppress


def _allocate_budget(total: int, capacities: np.ndarray) -> np.ndarray:
    if total < 0 or np.any(capacities < 0):
        raise ValueError("budget and capacities must be nonnegative")
    cap_total = int(capacities.sum())
    if total > cap_total:
        raise ValueError(f"Requested total {total} exceeds eligible-cell capacity {cap_total}")
    if total == 0:
        return np.zeros_like(capacities, dtype=np.int64)
    raw = total * capacities.astype(np.float64) / cap_total
    base = np.floor(raw).astype(np.int64)
    left = total - int(base.sum())
    if left:
        order = np.argsort(-(raw - base), kind="mergesort")
        base[order[:left]] += 1
    return base


def spatial_column_folds(
    shape: tuple[int, int],
    footprint: np.ndarray,
    *,
    n_folds: int = 5,
    boundary_guard_px: int = 3,
) -> list[dict]:
    """Create non-overlapping east-west blocks with a 300 m inter-fold guard.

    Fault labels inside a fold are completely held out. We score only within each
    block; no predictions or truth from an adjacent fold can contribute across a
    boundary. The small strip removed at block edges makes the score conservative.
    """
    h, w = shape
    foot = np.asarray(footprint, bool)
    if foot.shape != shape:
        raise ValueError("footprint shape mismatch")
    if n_folds < 2 or boundary_guard_px < 0:
        raise ValueError("n_folds must be at least 2 and boundary_guard_px nonnegative")
    boundaries = np.linspace(0, w, n_folds + 1, dtype=int)
    row_lo = min(boundary_guard_px, h)
    row_hi = max(row_lo, h - boundary_guard_px)
    folds: list[dict] = []
    for fold_id in range(n_folds):
        c0, c1 = int(boundaries[fold_id]), int(boundaries[fold_id + 1])
        if fold_id > 0:
            c0 += boundary_guard_px
        if fold_id < n_folds - 1:
            c1 -= boundary_guard_px
        tile = np.zeros(shape, dtype=bool)
        if c1 > c0 and row_hi > row_lo:
            tile[row_lo:row_hi, c0:c1] = True
        block_foot = tile & foot
        folds.append({
            "fold_id": fold_id,
            "col_start": c0,
            "col_end": c1,
            "row_start": row_lo,
            "row_end": row_hi,
            "domain": tile,
            "footprint_pixels": int(block_foot.sum()),
        })
    if sum(item["footprint_pixels"] for item in folds) == 0:
        raise ValueError("Spatial fold assignment produced no valid footprint cells")
    return folds


def _pool(parts: list[MetricComponents]) -> MetricComponents:
    tp = sum(item.tp_weighted for item in parts)
    fp = sum(item.fp_weighted for item in parts)
    fn = sum(item.fn_weighted for item in parts)
    truth_n = sum(item.truth_pixels for item in parts)
    mass = sum(item.prediction_mass for item in parts)
    return MetricComponents(tp, fp, fn, truth_n, mass, dti_from_components(tp, fp, fn))


def evaluate_spatial_holdout(
    score_fields: Mapping[str, np.ndarray],
    labels: np.ndarray,
    footprint: np.ndarray,
    analysis_support: np.ndarray,
    *,
    total_budget: int = 37_654,
    min_distance_px: float = 3.0,
    n_folds: int = 5,
    boundary_guard_px: int = 3,
) -> dict:
    """Compare each score field at identical mass and placement across held-out blocks."""
    y = np.asarray(labels, bool)
    foot = np.asarray(footprint, bool)
    support = np.asarray(analysis_support, bool)
    if y.ndim != 2 or y.shape != foot.shape or y.shape != support.shape:
        raise ValueError("labels, footprint and analysis_support must share a 2D shape")
    if np.any(y & ~foot):
        raise ValueError("positive labels outside declared footprint")
    if not score_fields:
        raise ValueError("at least one prediction score field is required")
    for name, field in score_fields.items():
        if np.asarray(field).shape != y.shape:
            raise ValueError(f"score field {name!r} has the wrong shape")

    folds = spatial_column_folds(y.shape, foot, n_folds=n_folds, boundary_guard_px=boundary_guard_px)
    capacities = np.asarray([
        int((item["domain"] & foot & support).sum()) for item in folds
    ], dtype=np.int64)
    budgets = _allocate_budget(total_budget, capacities)
    results: dict[str, dict] = {}

    for name, raw in score_fields.items():
        score = np.asarray(raw, dtype=np.float32)
        fold_rows = []
        components = []
        for fold, budget in zip(folds, budgets, strict=True):
            r0, r1 = fold["row_start"], fold["row_end"]
            c0, c1 = fold["col_start"], fold["col_end"]
            sl = np.s_[r0:r1, c0:c1]
            domain = fold["domain"][sl]
            tile_foot = foot[sl]
            tile_support = support[sl]
            tile_score = score[sl]
            valid = domain & tile_foot & tile_support & np.isfinite(tile_score) & (tile_score > 0.0)
            pred, receipt = greedy_kernel_suppress(
                tile_score, valid, int(budget), min_distance_px=min_distance_px
            )
            truth = y[sl] & domain
            comp = dti_components(pred, truth)
            components.append(comp)
            fold_rows.append({
                "fold_id": fold["fold_id"],
                "bounds_pixels": {
                    "row_start": r0, "row_end": r1,
                    "col_start": c0, "col_end": c1,
                },
                "valid_cells": int(capacities[fold["fold_id"]]),
                "prediction_budget": int(budget),
                "truth_pixels": comp.truth_pixels,
                "dti": comp.dti,
                "tp_weighted": comp.tp_weighted,
                "fp_weighted": comp.fp_weighted,
                "fn_weighted": comp.fn_weighted,
                "placement": receipt,
            })
        pooled = _pool(components)
        results[name] = {
            "folds": fold_rows,
            "pooled": pooled.to_dict(),
            "mean_fold_dti": float(np.mean([row["dti"] for row in fold_rows])),
            "fold_dti_std": float(np.std([row["dti"] for row in fold_rows])),
        }

    baseline_names = [name for name in score_fields if name.startswith("single_")]
    if not baseline_names:
        raise ValueError("score_fields must include single_* gradient baselines")
    best_baseline = max(baseline_names, key=lambda name: results[name]["pooled"]["dti"])
    baseline_best_dti = float(results[best_baseline]["pooled"]["dti"])
    candidate = results.get("cross_gradient")
    if candidate is None:
        raise ValueError("score_fields must include a 'cross_gradient' candidate")
    candidate_dti = float(candidate["pooled"]["dti"])
    per_fold_deltas = []
    for i in range(n_folds):
        strongest_baseline_fold = max(results[name]["folds"][i]["dti"] for name in baseline_names)
        per_fold_deltas.append(float(candidate["folds"][i]["dti"] - strongest_baseline_fold))

    controls = [name for name in score_fields if name.startswith("control_")]
    control_wins = {
        name: candidate_dti > float(results[name]["pooled"]["dti"])
        for name in controls
    }
    nonempty_truth = all(row["truth_pixels"] > 0 for row in candidate["folds"])
    positive_fold_wins = sum(delta > 0.0 for delta in per_fold_deltas)
    no_catastrophic_fold = min(per_fold_deltas) >= -0.01
    promotion_pass = bool(
        candidate_dti > baseline_best_dti + 1.0e-6
        and positive_fold_wins >= 4
        and no_catastrophic_fold
        and nonempty_truth
        and all(control_wins.values())
    )
    return {
        "protocol": "five non-overlapping east-west spatial blocks; all catalogue labels in each scored block hidden; 300 m inter-block guard; same mass and 300 m-aware placement for all arms",
        "truth_semantics": "public existing-fault catalogue proxy only; not private expert labels",
        "metric": {"alpha": ALPHA, "beta": BETA, "radius_px": 3.0, "pixel_size_m": 100},
        "mass_budget_total": int(total_budget),
        "min_distance_px": float(min_distance_px),
        "fold_budgets": budgets.tolist(),
        "candidate_id": "cross_gradient",
        "best_single_layer_baseline": best_baseline,
        "pooled_dti_delta_vs_best_single": candidate_dti - baseline_best_dti,
        "per_fold_delta_vs_strongest_single": per_fold_deltas,
        "positive_fold_wins": int(positive_fold_wins),
        "controls_pooled_win": control_wins,
        "nonempty_truth_each_fold": bool(nonempty_truth),
        "no_catastrophic_fold": bool(no_catastrophic_fold),
        "promotion_pass": promotion_pass,
        "results": results,
    }


def apply_analysis_role(result: dict, role: str) -> dict:
    """Prevent sensitivity-only runs on locked folds from authorizing promotion."""
    if role not in {"primary", "sensitivity"}:
        raise ValueError("analysis role must be 'primary' or 'sensitivity'")
    out = dict(result)
    thresholds_met = bool(result.get("promotion_pass", False))
    out["thresholds_met_on_this_run"] = thresholds_met
    out["analysis_role"] = role
    out["promotion_pass"] = bool(role == "primary" and thresholds_met)
    out["promotion_gate_note"] = (
        "Fixed preregistered primary scale; a pass would only authorize uniqueness and format review."
        if role == "primary"
        else "Sensitivity-only result on the same locked labels; numerical thresholds are descriptive and cannot authorize promotion."
    )
    return out
