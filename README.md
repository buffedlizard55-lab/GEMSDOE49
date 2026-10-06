# GEMSDOE49 — Cross-gradient structural coupling for the DOE GEMS Prize

> **Read this file first in every session.** It contains the standing brief (below),
> the current state of play, and the ground rules inherited from the GEMSDOE1–48 family.

## What this repo delivers (this session)

A **unique, validator-safe GeoTIFF submission** for the
[DOE GEMS Prize](https://www.drivendata.org/competitions/306/competition-doe-gems/)
(DrivenData #306): magnitude-weighted cross-gradient coupling of isostatic gravity,
reduced-to-pole magnetics and smoothed topography (Gallardo & Meju 2004 operator,
data-domain horizontal-gradient variant), validated on a preregistered
hide-and-recover holdout before being allowed near a submission slot.

| artifact | value |
| --- | --- |
| **Submission file (all-finite, safest)** | [`docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-zeros.tif`](docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-zeros.tif) |
| Submission twin (NaN outside footprint) | [`docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-nan.tif`](docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-nan.tif) |
| Unique submission name | `GEMSDOE49-XGRAD-ALIGN` |
| Holdout DTI (4-fold hide-and-recover) | **0.0203** (fold scores 0.0365 / 0.0126 / 0.0160 / 0.0159) |
| Margin vs best single-layer gradient baseline | **+0.0047 to +0.0058** (+30 % relative, all 3 baselines, all 4 sigmas) |
| Uniqueness audit vs 57 prior family submissions | **UNIQUE** — max Spearman ρ = −0.27, max top-k overlap = 0.29 |
| Format checks | all pass — single band, float32, EPSG:32611, 100 m, 3292×3730, all values in [0,1], no nodata tag |

*The holdout number is an instrument reading against hidden catalogue segments —
it is NOT an organiser score. No organiser score exists for this file.*

Site: open [`docs/index.html`](docs/index.html) (GitHub Pages serves it from `docs/`).

## Repo map

```
data/
  training_features.tif          official 19-band stack (sha256-pinned, reassembled from bridge parts)
  bridge/                        manifest + official rasters (labels, sample submission, parts)
  external/                      lidar scarp field, radiometric u8, topo u8 (sha256-pinned bridges)
  prior/                         57 prior family submission rasters (uniqueness-audit corpus)
src/gems49/
  metric.py                      exact DTI re-implementation (14 tests incl. brute-force oracle)
  data.py                        raster access, band table
  holdout.py                     spatially-blocked hide-and-recover folds
  xgrad.py                       cross-gradient coupling operator
  emission.py                    metric-aware placement (NMS dots), TIF writer, format verifier
scripts/
  run_pipeline.py                end-to-end: diagnostics -> holdout sweep -> emission -> receipts
  uniqueness_audit.py            hash + rank-correlation + top-k-overlap relabel guard
  validate_h49b_alteration.py    next-session hypothesis validation (radiometric alteration gate)
tests/                           pytest suite for metric + emission
evidence/                        every number on this site, as JSON receipts
docs/                            GitHub Pages site (index, executive summary, hypotheses, sources)
```

## Ground rules inherited from the family (verified)

1. **The catalogue is masked in live scoring** (staff, forum thread 11516 post 4:
   pixel-exact mask, identical to the training labels; predictions near known
   faults but far from new-fault truth are fully penalised). Emission on the
   catalogue is pure cost: every emitted file here has **zero** on-catalogue pixels.
2. **The metric is a mass budget.** DTI = T / (0.2(T+S−M) + 0.8|G|): prediction mass
   that does not cover truth costs 0.2 per unit. Dotted, NMS-spaced emission at a
   calibrated budget beats thick surfaces.
3. **Holdout truth is the catalogue, so it cannot reward a genuinely new fault
   family** (family registry IR-32-PROXY-01). Holdout numbers rank *placement
   rules*, not discoveries. Combine with physical plausibility; never substitute.
4. **Portal failure mode:** `"Predicted values must be in range [0, 1]"` is caused
   by values outside [0,1] OR a nodata sentinel (e.g. −3.4e38). All primaries are
   written all-finite with **no nodata tag**.
5. Every claim carries an evidence class; negative results stay in the table.

---

## THE STANDING BRIEF (verbatim from the project owner — read every session)

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
> The goal of this project is to place top of the leaderboard (current top: 0.3195). We need to create a project that can compete and place top of the leaderboard. Create a GitHub page for this repo that has a clean UI, is organized, user friendly, simple and easy to use, with all relevant information in an easy-to-read format and official verified links as sources. The site must make the submission TIF obvious and easy to download (one click), include an executive summary explaining exactly how to make a submission, and fix the portal error "Predicted values must be in range [0, 1]" by guaranteeing every value is in [0,1] with no sentinel nodata. Give the submission a unique name and a short comment for the submit form.
>
> Work line by line verifying from official verified trusted sources, provide links for manual review. No manual input — work autonomously. Flag any irregularities for review. No hallucinations. Run this task through multiple passes (implement → review for bugs/missing requirements → re-check against the original request). Create a pull request and merge it onto main. Make suggestions for remaining work and limitations.
>
> Core values: **Maximize P(Win)** — every decision weighed for probability of winning; **Own the Outcome** — own results end to end, treat failure and success as signals.

(Full session transcript with all site scores, Dropbox mirrors and portal text is
preserved in the git history of this README.)

---

## Known limitations (flagged for review, not smoothed over)

* **No GPU and no DrivenData login in this sandbox.** Training a U-Net and
  downloading `training_features.tif` from the data tab both require an enrolled
  machine; the rasters here were recovered byte-exact from the family's
  sha256-pinned bridge (see `data/bridge/manifest.json`).
* **Gravity station inventory not independently fetched.** The INGENIOUS gravity
  compilation (GDR 1391, DOI 10.15121/1881483) is blocked by sandbox egress; the
  station-spacing check is therefore empirical (gradient-variance retention +
  semivariogram, `evidence/gravity_support_diagnostics.json`). Manual check: open
  <https://gdr.openei.org/submissions/1391> and confirm station spacing.
* **Holdout truth is catalogue geometry** (rule 3 above).
* Direct egress to `*.github.io`, Dropbox, USGS S3 and sciencebase is blocked in
  this sandbox; everything external was fetched via the GitHub Contents API or
  `fetch_page`. Irregularity log: three earliest family submissions
  (`gems-submission-20260925T001403Z-7f00890a`, GEMSDOE4 `20260926T163915Z`,
  5GEMSDOE `20260926T175114Z`) are **not preserved** in any family repo; the
  uniqueness audit substitutes the closest archived files and records this.

## Next session — highest-value actions

1. Run `python3 scripts/validate_h49b_alteration.py` (already written) and adopt
   the alteration gate only if it beats the coupled-only holdout DTI of 0.0203.
2. Extend the budget sweep past 40k dots (holdout DTI was still rising at 40k —
   break-even not yet reached) — but cap emission: α=0.2 keeps marginal cost low,
   yet the family's live evidence says over-emission above ~44–46k px earned
   nothing extra on the leaderboard.
3. Strike-extension of coupled edges past catalogue tips (H49-B v2) — the
   staff-confirmed scoring window is 300 m, catalogue misfit up to 400 m
   (Hermant et al. 2025).
4. En-echelon scarp-pair stepovers (H49-D) using the 1 m lidar scarp field
   already on disk (`data/external/lidar_scarp_features_u8.tif`).
