# GEMSDOE49 — auditable GEMS fault-mapping research

> **Read this file first each project session.** It records the standing brief, the current candidate gate, data limitations, and the actions that are safe to take next.

## Executive decision — 2026-10-06

**No GeoTIFF is currently cleared for a DrivenData submission. No weekly slot was used.** Two distinct H49-XG results now exist in this checkout and must not be conflated:

1. **The fixed five-block H49-XG primary produced by this task branch failed its preregistered promotion gate.** Pooled public-catalogue proxy DTI was `0.10041455` versus `0.08193976` for the best single-layer baseline. Four of five blocks won, but one paired delta was `−0.01180243`, beyond the allowed `−0.01`. The 2 km-FWHM result is sensitivity-only. This run generated **no GeoTIFF**. See [`research/h49_xg_result_20261006.md`](research/h49_xg_result_20261006.md).
2. **A separate artifact from merged PR #1 already exists in `docs/downloads/`.** Its NaN-outside twin has been independently checked against the official GeoTIFF geometry/range contract, and its four-fold score reproduces from the committed metric and local input mirrors. It is **not cleared for submission**: its four-fold protocol is not a spatially isolated block holdout, its prediction mask excludes the held-out catalogue components, its selected 300 m gravity smoothing is finer than the 1 km grid spacing described for the relevant official USGS regional product (exact band lineage is unconfirmed), and the full uniqueness corpus is not in the checkout. Audit: [`docs/reports/pr1_candidate_audit_20261006.md`](docs/reports/pr1_candidate_audit_20261006.md).

Do not call either proxy result a public leaderboard score or evidence that a candidate beats `0.2778` on the private test. The 2026-10-06 public board showed `extradr19` at `0.2778` (rank 13, 10 submissions) and `xiaofanhu` at `0.3774` (rank 1); the board does not identify the raster behind the 0.2778 row. The H33-2-B2 project file is marked unscored and its `0.274673...` is a project projection, not an organizer result. The attribution remains unresolved; see [`docs/reports/source_audit_20261006.md`](docs/reports/source_audit_20261006.md).

## Existing PR #1 research artifact — not an upload recommendation

The separate PR #1 candidate was named `GEMSDOE49-XGRAD-ALIGN`. If inspecting its bytes, use the NaN-outside file, not the all-finite zeros twin:

| Item | Independent audit result |
| --- | --- |
| Research artifact | [`gemsdoe49-xgrad-coupled-20261006T203447Z-nan.tif`](docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-nan.tif) |
| SHA-256 | `bef0881aadc8366f1a08e13aac119cfa71492750e93426105b9e5745c7d701bf` |
| Format | Single-band float32, 3,730 × 3,292, EPSG:32611, 100 m; same transform/bounds as the sample grid; `nodata=None`; all 5,167,373 in-footprint cells finite in `[0,1]`; all 7,111,787 outside-footprint cells NaN; 40,000 positive cells; zero positive cells on the catalogue. |
| Recomputed four-fold score | DTI `0.02026648455187942`; fold DTI `0.03650436 / 0.01260255 / 0.01601234 / 0.01594669`. This is a public-catalogue instrument reading, **not** an organizer score. |
| Status | **RESEARCH ONLY — NOT SUBMISSION-CLEARED.** The `-zeros.tif` twin has zero, not null/NaN, outside the footprint and therefore is not the file to use under the official format statement. |

The candidate README's 300 m (`σ=3 px`) setting and 40,000-dot budget were selected from a sweep on the same four-fold proxy. The actual file uses quantile-clipped scores before NMS; the raw-score sweep value `0.0203026` is not the exact final-file score. The independently replayed candidate file matches the reported final-fold receipt. Main's archived 57-file uniqueness receipt reports max *signed* Spearman `ρ = −0.266` and maximum smaller-set overlap `0.290`; it includes the H33-2-B2 raster (reported overlap `0.0210`, `ρ = −0.9796`). The prior raster corpus/manifest is not tracked, three original artifacts were unavailable in that audit, and correlations/overlaps against those absent bytes cannot be independently confirmed here. This is not a complete all-priors clearance.

## Original fixed H49-XG result

- Candidate method: three-field, magnitude-weighted axial horizontal-gradient alignment from `rtp`, `det_elev`, and `iso_grav_anom`; fixed Gaussian `σ=500 m` (FWHM 1.177 km), 37,654 prediction mass, five buffered east–west blocks.
- Pooled public-catalogue DTI: **0.10041455**, versus **0.08193976** for `single_det_elev`, delta **+0.01847479**. Four of five block comparisons win; the failing paired delta is **−0.01180243**.
- Promotion gate: **FAIL** because one block is below the predeclared `−0.01` stability floor. The 2 km-FWHM sensitivity has positive deltas in all five blocks but was measured on the same locked labels and cannot rescue the primary.
- The labels are a public existing-fault proxy, not private expert labels. The data bytes are pinned owner-hosted mirrors, not authenticated against the login-walled competition download. The block audit also found that 8-connected pixels are not reliable geological segment IDs.
- No output raster was generated by this experiment. No leaderboard score or slot use is claimed.

Full receipt and limitations: [`research/h49_xg_result_20261006.md`](research/h49_xg_result_20261006.md), [`research/input_spacing_audit_20261006.json`](research/input_spacing_audit_20261006.json), [`research/receipts/`](research/receipts/).

## Data and source-provenance rules

The official DrivenData data page redirects to login. [`data/manifest.json`](data/manifest.json) pins public, owner-hosted mirrors and exact byte counts/hashes; those hashes prove mirror consistency only, not that the files are byte-identical to organizer downloads. [`data/README.md`](data/README.md) explains restoration. Large rasters and caches are ignored by Git.

The competition page specifies a 100 m, EPSG:32611 feature grid and requires a single-band float32 submission with values in `[0,1]`, same bounds/resolution/CRS, and null/NaN outside the footprint. Staff clarified that the live catalogue mask is pixel-exact, identical to the provided training labels; nearby unmasked predictions remain penalized if far from new-fault truth, and new-fault truth can occur within 300 m of a known trace. Links and checked wording are in the source-audit report.

**Gravity resolution remains a hard caution.** The official USGS data-release README for the regional `GB_iso_grav_anom.tif` says its grid cell size is 1 km. The competition mirror's band 13 has the matching isostatic-gravity name, but its exact provenance/version has not been proved. Do not interpret the 100 m resampled pixels as independent gravity observations or treat a 300 m derivative scale as established. The data-source-to-band mapping must be confirmed or conservatively screened at a 1 km-class scale in a newly frozen spatial validation.

## Persistent project prompt — reread before each session

> Review the repo.
> THE FOLLOWING IS THE HIGHEST URGENCY AND MUST BE FOLLOWED!
> MUST GENERATE A UNIQUE TIF SUBMISSION FOR THE COMPETITION.  DO NOT COPY A PREVIOUS SUBMISSION UNLESS IT'S FOR LEARNING AND EDUCATION.  BUT WE MUST GENERATE A UNIQUE TIF SUBMISSION.
>
> There should be an easy to download submission tif file as described by the prompt.  Read the entire prompt.
>
> Cross-gradient structural coupling of gravity, magnetics and topography. Treat agreement across independent signals as an operator, not a product of scores. Gallardo and Meju (Journal of Geophysical Research, 2004, doi:10.1029/2003JB002716) defined the cross-gradient function, the cross product of two property fields' gradient vectors. It vanishes where both fields change in the same direction, so it measures shared structure without assuming a petrophysical relationship. Fregoso and Gallardo (Geophysics 74(4), L31–L42, 2009) applied it to gravity and magnetics. An ASEG 2012 extended abstract (10.1071/ASEG2012ab273) used the angle and cross product of the horizontal gradients of gravity and reduced-to-pole magnetics directly in the data domain, with no inversion. Compute horizontal-gradient vectors of the isostatic gravity anomaly, the reduced-to-pole magnetic anomaly and a smoothed elevation surface at a common resolution. Score each pixel by magnitude-weighted alignment (both gradients large and parallel). Keep magnitude-weighted orthogonality as a separate confounder flag, not as evidence. A fault juxtaposing rocks of different density and susceptibility puts both edges on one trace, while a buried volcanic edge usually moves the magnetic edge alone, so the coupled score should suppress single-layer edges. Check the native gravity station spacing before trusting 100 m gradients, and expect density steps at basin margins. The candidate must beat the repo's single-layer gradient baseline on hide-and-recover holdout segments before it is allowed near a submission. Normalize to [0,1], apply the repo's metric-aware placement, and write the required float32 GeoTIFF (EPSG:32611, 100 m, NaN only outside the footprint and none inside). Before download, hash the raster and compare rank correlation and top-k overlap against every prior submission, and reject it as a relabel if it isn't clearly different.
>
> We need to study, analyze, and understand the highest score from the GEMSDOE sites where the submission TIF is downloaded from: GEMSDOE32, `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros: 0.2778`. Why and how did this get the highest score and are we able to generate a submission that scores higher than 0.2778? Answer the question using PhD-level experience, knowledge, and judgement. Then use the answer to generate a unique TIF submission into the competition. Must be unique submission unlike any within the GEMSDOE sites. Verify working line by line, no hallucinations.
>
> Before implementing, generate 3–5 candidate geological hypotheses we haven't tried yet, each naming: the specific layer(s) involved, the physical signature being targeted, why it should catch a fault missing from the USGS/INGENIOUS catalogue rather than one already in it, and how it differs from anything already implemented in this repo. Rank them by expected DTI improvement and implementation cost. Validate the top candidate on our spatially-blocked holdout set before touching a weekly submission slot. If a candidate can't be validated without new external data, name the specific free, official source needed and check it's obtainable before proposing the idea as viable.
>
> The goal of this project is to place top of the leaderboard (current top: 0.3195). We need to create a project that can compete and place top of the leaderboard. Create a GitHub page for this repo that has a clean UI, is organized, user friendly, simple and easy to use, with all relevant information in an easy-to-read format and official verified links as sources. The site must make the submission TIF obvious and easy to download (one click), include an executive summary explaining exactly how to make a submission, and fix the portal error "Predicted values must be in range [0, 1]" by guaranteeing every value is in [0,1] with no sentinel nodata. Give the submission a unique name and a short comment for the submission on DrivenData so we can identify it.
>
> Work line by line verifying from official verified trusted sources, provide links for manual review. No manual input — work autonomously. Flag any irregularities for review. No hallucinations. Run this task through multiple passes (implement → review for bugs/missing requirements → re-check against the original request). Create a pull request and merge it onto main. Make suggestions for remaining work and limitations.
>
> Core values: **Maximize P(Win)** — every decision weighed for probability of winning; **Own the Outcome** — own results end to end, treat failure and success as signals.

## Next actions (submission gate stays closed)

1. Build a genuinely spatially isolated, held-out-block protocol for the separate PR #1 operator. Hide every catalogue pixel in each held-out block, allow the predictor to consider those pixels, apply the remaining visible-catalogue mask only, use the same fixed mass/placement for every arm, and compare against all three single-layer baselines. Freeze scale, splits, tie-breaking and gate before looking at the scores.
2. Resolve gravity provenance/native support for competition band 13 or use a conservative scale supported by official metadata. Reconcile this with the fixed H49-XG primary failure; do not select the 2 km sensitivity post hoc.
3. Restore and document the complete prior-raster corpus (or enumerate unavailable artifacts) before claiming uniqueness versus every prior; add a deterministic tie-break to NMS before generating any replacement artifact.
4. Only after validation: build one new NaN-outside GeoTIFF, independently re-read bytes/geometry/range, verify SHA-256 and uniqueness, and update the site with a one-click link, exact upload name/comment, and explicit no-score status. No DrivenData slot is to be touched until then.

## Reproduce tests (no candidate generation)

```bash
python -m pip install -e '.[test]'
python -m pytest tests -q
```

Restoring local inputs, when required:

```bash
bash scripts/download_competition_data.sh
```

The PR #1 `scripts/run_pipeline.py` writes GeoTIFFs; **do not run it as part of tests or before the above validation gate is cleared**.

## Official/manual-review references

- [DrivenData GEMS challenge and official metric/format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [DrivenData public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/)
- [DrivenData staff mask clarification, post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4)
- [Official INGENIOUS GDR 1391](https://gdr.openei.org/submissions/1391)
- [USGS regional geophysical maps release, DOI 10.5066/P9Z6SA1Z](https://doi.org/10.5066/P9Z6SA1Z)
- [USGS Nevada isostatic gravity source](https://pubs.usgs.gov/ds/2006/234/nv_iso.htm) and [Bouguer station/grid notes](https://pubs.usgs.gov/ds/2006/234/nv_boug.htm)
- [Official GEMS rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
