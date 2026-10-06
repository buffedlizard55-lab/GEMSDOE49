# GEMSDOE49 — reproducible GEMS fault-mapping research

> **Current decision (2026-10-06): no submission-ready GeoTIFF exists.** The fixed H49-XG primary holdout improved the pooled public-catalogue proxy over its best single-layer baseline, but failed the preregistered fold-stability rule. The 2 km-FWHM sensitivity is promising but cannot authorize promotion because it was evaluated on the same locked labels. No submission slot was used. See [`research/h49_xg_result_20261006.md`](research/h49_xg_result_20261006.md) and the [executive summary site](docs/index.html).

## Mission

**Maximize P(Win)**—the probability of a real, defensible improvement on the organizer's private and final evaluation sets, not a seductive score on public labels. Preserve scientific value: identify credible previously unmapped faults, document uncertainty, and make the work reproducible.

## Core values

- **Own the Outcome.** Close the loop from hypotheses through source verification, implementation, falsification tests, artifact checks, and honest delivery. Do not leave a broken download, an untested command, or an unreported failure.
- **Evidence before narrative.** Separate official facts, mirror observations, model projections, and leaderboard results. Never upgrade a projection or holdout score into a public/private leaderboard claim.
- **Protect the submission.** Do not generate or upload a candidate unless the frozen holdout and control gates pass. Audit a genuinely new raster against all recoverable prior artifacts; filename/hash changes alone do not prove prediction novelty.
- **Be transparent about provenance.** The official DrivenData data tab requires login. Current local rasters are pinned owner-hosted mirrors; checksums prove mirror consistency only, not organizer authentication.
- **Respect the geology.** Treat gravity/magnetic/topographic concordance as a hypothesis, not a fault label. Report native-resolution uncertainty, alternative explanations, and negative results.

## Persistent project prompt — reread at the start of every project session

1. Read this README, [`data/README.md`](data/README.md), [`research/hypotheses.md`](research/hypotheses.md), [`research/source_audit_20261006.md`](research/source_audit_20261006.md), and the latest result/receipt before changing code or selecting a candidate. Check `git status` and preserve existing work.
2. **Maximize P(Win)** while **Owning the Outcome**. Use current official sources and exact evidence; line-check critical claims. Mark unknowns; never hallucinate provenance, board attribution, or score.
3. Keep the preregistered primary scale, metric, fold design, mass budget, baseline, controls, and promotion thresholds fixed. Label later variants as sensitivity/exploratory. Do not select a variant by looking at the same locked holdout.
4. Use only authorized data; report mirror-vs-organizer provenance explicitly. Check source resolution before interpreting 100 m gradients.
5. A submission candidate must beat the best single-layer baseline under the preregistered spatial holdout and satisfy the fold/control criteria. Before a submission slot is used, independently validate GeoTIFF format, hash it, and compare prediction ranks/top-k overlap against every prior recoverable artifact. Never relabel or copy a previous submission.
6. If and only if a candidate is cleared, make the executive summary and one-click GeoTIFF download unmistakable, include exact platform steps, and use a distinct submission name and truthful comment. If no candidate is cleared, say so plainly and do not expose a fake/unsafe submission download.
7. Follow the official rules: disclose generative-AI use in the narrative if applicable; a participant is responsible for the final accuracy and authorship representations. Keep source links, versions, tests, receipts, and unresolved work visible.
8. Review the implementation, data, metric, site, and final claims in multiple passes. Open a PR from this fixed Arena branch when the repo/tooling permits; never change branches.

## Current H49-XG result

- Fixed primary: three-field magnitude-weighted axial gradient alignment from `rtp`, `det_elev`, and `iso_grav_anom`; Gaussian σ=500 m (FWHM 1.177 km); 37,654 prediction mass; five east–west spatial blocks with 300 m guards.
- Pooled public-catalogue DTI: **0.100415**, versus **0.081940** for the best single-layer baseline (`single_det_elev`), a **+0.018475** delta. Four of five paired folds win; random and shifted-RTP controls are beaten.
- Promotion gate: **FAIL**. One paired fold is **−0.011802**, below the preregistered −0.01 stability limit. The 2 km-FWHM sensitivity's numerical screen passed all five folds, but is not promotion-eligible on the same locked labels.
- Spatial-fold caveat: the preregistration asked for whole-segment folds; the implementation used buffered east–west raster blocks because segment IDs were unavailable. An 8-connected diagnostic found 21 raster components spanning two scored blocks (850 pixels). Raster components are not verified geological segments, so no scale is cleared.
- This is a known-fault catalogue proxy, not private-test or leaderboard performance. The public board's 0.2778 entry under `extradr19` is not mapped to the `h33-h33-2-b2` artifact.
- **No candidate TIFF has been generated.** Uniqueness/hash/rank-overlap and final-format receipts remain unrun by design. No DrivenData submission was made.

Detailed evidence: [`research/h49_xg_result_20261006.md`](research/h49_xg_result_20261006.md), [`research/input_spacing_audit_20261006.json`](research/input_spacing_audit_20261006.json), and [`research/source_audit_20261006.md`](research/source_audit_20261006.md).

## Reproduce the research checks

Install the local package and test extras:

```bash
python -m pip install -e '.[test]'
pytest -q
```

Restore the public mirrors if `.cache/gems_data/` is absent, then run the frozen primary and the preregistered sensitivities:

```bash
bash scripts/download_competition_data.sh
python scripts/audit_input_spacing.py
python scripts/run_experiment.py --analysis-role primary
python scripts/run_experiment.py --analysis-role sensitivity --sigma-px 2.5 --out .cache/analysis/h49_xg_sensitivity_500m.json
python scripts/run_experiment.py --analysis-role sensitivity --sigma-px 10 --out .cache/analysis/h49_xg_sensitivity_2000m.json
```

The default run writes its detailed JSON receipt under `.cache/analysis/`, which is intentionally ignored by Git. The primary runner rejects other scales unless they are explicitly labelled sensitivity-only.

## Entry points

- [`docs/index.html`](docs/index.html) — executive summary and evidence dashboard; submission/download gate is visibly closed.
- [`research/hypotheses.md`](research/hypotheses.md) — frozen hypothesis ranking and promotion protocol.
- [`research/h49_xg_result_20261006.md`](research/h49_xg_result_20261006.md) — measured primary and sensitivity results.
- [`research/source_audit_20261006.md`](research/source_audit_20261006.md) — official-source and 0.2778 attribution review.
- [`research/input_spacing_audit_20261006.json`](research/input_spacing_audit_20261006.json) — band metadata and observed lag-correlation audit.
- [`research/receipts/`](research/receipts/) — primary and sensitivity JSON results with input, preregistration and implementation hashes.
- [`data/manifest.json`](data/manifest.json) — pinned mirror refs, byte counts, and hashes.
- [`scripts/run_experiment.py`](scripts/run_experiment.py) and [`src/gemsdoe49/`](src/gemsdoe49/) — metric-aware holdout implementation.

## Official references

- [DrivenData GEMS challenge](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Official problem description, metric, and submission format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [Current public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)
- [September 2026 Official Rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
- [USGS Nevada isostatic-gravity source page](https://pubs.usgs.gov/ds/2006/234/nv_iso.htm) and [complete-Bouguer grid/station details](https://pubs.usgs.gov/ds/2006/234/nv_boug.htm)
