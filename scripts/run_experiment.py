#!/usr/bin/env python3
"""Run the preregistered H49-XG spatially-blocked public-catalogue screen."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gemsdoe49.features import cross_gradient_alignment, shifted_alignment_control  # noqa: E402
from gemsdoe49.holdout import apply_analysis_role, evaluate_spatial_holdout  # noqa: E402
from gemsdoe49.raster import load_inputs, sha256_file  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / ".cache" / "gems_data")
    parser.add_argument("--out", type=Path, default=ROOT / ".cache" / "analysis" / "h49_xg_holdout.json")
    parser.add_argument("--sigma-px", type=float, default=5.0)
    parser.add_argument(
        "--analysis-role", choices=("primary", "sensitivity"), default="primary",
        help="Only the fixed preregistered 1 km-class primary scale is promotion-eligible; ladder runs are sensitivity-only.",
    )
    parser.add_argument("--budget", type=int, default=37_654)
    parser.add_argument("--random-seed", type=int, default=49_20261006)
    args = parser.parse_args()
    if args.analysis_role == "primary" and not np.isclose(args.sigma_px, 5.0, rtol=0.0, atol=1.0e-12):
        parser.error("the preregistered primary uses --sigma-px 5.0; use --analysis-role sensitivity for other scales")
    if args.analysis_role == "sensitivity" and np.isclose(args.sigma_px, 5.0, rtol=0.0, atol=1.0e-12):
        parser.error("the fixed primary scale must be run with --analysis-role primary")

    fields, labels, footprint, input_meta = load_inputs(args.data_dir)
    common_valid = footprint.copy()
    for field in fields.values():
        common_valid &= np.isfinite(field)
    result = cross_gradient_alignment(
        fields,
        common_valid,
        sigma_px=args.sigma_px,
        pixel_size_m=100.0,
        robust_percentile=99.0,
    )
    shift_score, shift_support = shifted_alignment_control(
        result, shifted_field="rtp", shift_px=(0, 10)
    )
    rng = np.random.default_rng(args.random_seed)
    random_score = rng.random(footprint.shape, dtype=np.float32)
    random_score[~result.valid_support] = 0.0

    score_fields = {
        "cross_gradient": result.score,
        "single_gravity": result.single_layer_baselines["iso_grav_anom"],
        "single_rtp": result.single_layer_baselines["rtp"],
        "single_det_elev": result.single_layer_baselines["det_elev"],
        "control_random": random_score,
        "control_shift_rtp_dx10": shift_score,
    }
    run = evaluate_spatial_holdout(
        score_fields,
        labels,
        footprint,
        result.valid_support,
        total_budget=args.budget,
        min_distance_px=3.0,
        n_folds=5,
        boundary_guard_px=3,
    )
    run = apply_analysis_role(run, args.analysis_role)
    # The orthogonality field is intentionally separate and is not passed into scoring.
    orth_valid = result.valid_support
    orth = result.orthogonality_flag[orth_valid]
    align = result.score[orth_valid]
    implementation_files = [
        "scripts/run_experiment.py",
        "src/gemsdoe49/features.py",
        "src/gemsdoe49/metric.py",
        "src/gemsdoe49/placement.py",
        "src/gemsdoe49/holdout.py",
        "src/gemsdoe49/raster.py",
    ]
    implementation_hashes = {
        path: sha256_file(ROOT / path) for path in implementation_files
    }
    if args.analysis_role == "sensitivity":
        candidate_status = "SENSITIVITY_ONLY_NOT_PROMOTION_ELIGIBLE"
    elif run["promotion_pass"]:
        candidate_status = "PRIMARY_HOLDOUT_PASS_UNIQUENESS_AND_FORMAT_PENDING"
    else:
        candidate_status = "RESEARCH_ONLY_PRIMARY_GATE_FAILED"
    report = {
        "experiment_id": "H49-XG",
        "analysis_role": args.analysis_role,
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "candidate_status": candidate_status,
        "preregistration_file": "research/hypotheses.md",
        "preregistration_sha256": sha256_file(ROOT / "research" / "hypotheses.md"),
        "implementation_sha256": implementation_hashes,
        "input_provenance": "pinned public owner-hosted GitHub mirror; hash proves mirror consistency, not organizer authentication",
        "input_metadata": input_meta,
        "operator": {
            "name": "three-field magnitude-weighted cross-gradient axial alignment",
            "fields": list(fields.keys()),
            "sigma_px": float(args.sigma_px),
            "sigma_m": float(args.sigma_px * 100.0),
            "fwhm_m": float(args.sigma_px * 100.0 * 2.354820045),
            "robust_percentile": 99.0,
            "alignment_p99": float(np.percentile(align, 99)),
            "orthogonality_p99_separate_flag": float(np.percentile(orth, 99)),
            "orthogonality_positive_evidence": False,
            "support_pixels": int(result.valid_support.sum()),
            "shift_null_support_pixels": int(shift_support.sum()),
        },
        "validation": run,
        "interpretation_limits": [
            "This evaluates public existing-fault segments, not the private expert-labelled discoveries.",
            "The candidate is feature-only and does not use labels to construct its score; labels are used only to define the holdout truth.",
            "Five east-west blocks are spatially disjoint with a 300 m guard; they are not statistically independent samples of the Great Basin.",
            "A holdout pass is a screening result, not a leaderboard score or proof of private-set improvement.",
            "The training feature raster is an owner-hosted mirror and has not been authenticated against the login-walled organizer download in this environment.",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "report": str(args.out),
        "candidate_dti": run["results"]["cross_gradient"]["pooled"]["dti"],
        "best_single": run["best_single_layer_baseline"],
        "best_single_dti": run["results"][run["best_single_layer_baseline"]]["pooled"]["dti"],
        "delta": run["pooled_dti_delta_vs_best_single"],
        "positive_fold_wins": run["positive_fold_wins"],
        "controls": run["controls_pooled_win"],
        "thresholds_met_on_this_run": run["thresholds_met_on_this_run"],
        "promotion_pass": run["promotion_pass"],
        "analysis_role": args.analysis_role,
        "truth_by_fold": [row["truth_pixels"] for row in run["results"]["cross_gradient"]["folds"]],
        "prediction_budgets": run["fold_budgets"],
    }, indent=2))
    return 0 if all(row["truth_pixels"] > 0 for row in run["results"]["cross_gradient"]["folds"]) else 2


if __name__ == "__main__":
    raise SystemExit(main())
