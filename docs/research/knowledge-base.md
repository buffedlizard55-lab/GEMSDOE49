# Knowledge base — DOE GEMS Prize (DrivenData #306)

Every entry carries an evidence class. **[OFFICIAL]** = read from an organizer/USGS/DOE page or
product, link given. **[MEASURED]** = computed in this checkout from hash-verified bytes, with the
command that produced it. **[BLOCKED]** = could not be obtained here, with the reason.
**[OWNER-REPORT]** = a sibling repository's own unverified ledger.

---

## 1. The task

| # | Fact | Class | Source |
| --- | --- | --- | --- |
| 1.1 | Predict **geological faults** (structures indicative of geothermal resources) as a per-pixel confidence raster over the GeoDAWN region of the northwestern Great Basin, Nevada | OFFICIAL | [page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) |
| 1.2 | The public USGS fault database is *incomplete and may contain inaccurate data*; an expert panel labelled new faults absent from it, and those form the test set | OFFICIAL | ibid. |
| 1.3 | Two prize rounds: an Initial Round against a fixed private set of expert-labelled new faults (top 5 × $10K), then a Final Round re-scored against an **expanded** label set after experts review every team's submission (top 5: $100K/$70K/$40K/$25K/$15K) | OFFICIAL | [competition structure](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#competition-structure) |
| 1.4 | A single submission is selected for scoring across both rounds | OFFICIAL | ibid. |
| 1.5 | Weekly feedback allowance: 3 scored submissions per rolling week | OFFICIAL | [rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf) |
| 1.6 | Public and private test sets are spatial chunks of the GeoDAWN region; the public leaderboard is not the private score | OFFICIAL | [page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) |

## 2. The metric — exact

**[OFFICIAL]** All of section 2 is read from
<https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/>.

```
k(d)  = max(1 - d/R, 0),                                   R = 300 m, 100 m pixels

TP_w  = sum_{g in G} max_{x : d(x,g) <= R}  p(x) * k(d(x,g))
FP_w  = sum_{x : p(x) > 0}          p(x) * [1 - max_{g in G} k(d(x,g))]
FN_w  = sum_{g in G}                 [1 - max_{x : d(x,g) <= R} p(x) * k(d(x,g))]

DTI   = TP_w / (TP_w + alpha*FP_w + beta*FN_w + eps),      alpha = 0.2, beta = 0.8
```

* [OFFICIAL] the organizer's worked example: `TP_w = 3.00, FP_w = 1.89, FN_w = 2.00 -> 0.60`.
  **Reproduced to 0.602652** by `tests/test_metric.py::test_worked_example_from_the_page`.
* [MEASURED] the implementation also agrees with a literal `O(|G| x |X|)` double-loop
  transcription on small rasters, and satisfies `TP_w + FN_w = |G|` exactly.
* [MEASURED] since `alpha + beta = 1` and `FN_w = |G| - TP_w`, the denominator equals
  `0.2*(TP_w + FP_w) + 0.8*|G|`.

### 2.1 Three consequences that drive every design decision here

**(i) Binary support strictly dominates graded support.** [MEASURED, derived]
For a support scaled by `lambda`, `DTI(lambda) = lambda*T / (0.2*lambda*(T+F) + 0.8*|G|)`, and
`dDTI/dlambda = 0.8*|G|*T / (denominator)^2 > 0`. Any graded probability map is dominated by its
own thresholded support. Verified numerically in `tests/test_metric.py`.

**(ii) The exact marginal rule.** [MEASURED, derived] Adding one unit prediction pixel whose
kernel credit is `k` raises `TP_w` by `k`, raises the false-positive term by `1 - k`, so
`TP_w + FP_w` rises by exactly 1 and the denominator by exactly 0.2 — for **every** `k`. Hence the
pixel pays iff

```
k > alpha * DTI
```

At `DTI = 0.28` the bar is 0.056, i.e. every pixel within **283 m** of a scored truth pixel pays.
The metric is a *distance* filter, not a precision filter.

*Note.* A sibling repository states the bar as `alpha*s/(1 - alpha*s)` (0.0683 at s = 0.3195). That
expression is the bar for adding mass whose credit arrives in bulk; for a single added pixel the
exact value is `alpha*s` (0.0639). Both give ~280 m; this repository uses the exact one and proves
it in `tests/test_metric.py::test_marginal_bar_is_exact`.

**(iii) A shadowed pixel always loses.** [MEASURED, derived] A pixel for which no truth pixel takes
it as the best nearby prediction adds false-positive mass and no true-positive mass, so it lowers
the score. This is why the placement rule imposes a minimum separation rather than a threshold.

## 3. The masking rule — the single most decision-relevant fact

* **[OFFICIAL]** DrivenData staff, 2026-09-16: *"pixels corresponding to known USGS/INGENIOUS
  faults are masked / excluded from evaluation, so they do not count towards penalty terms"* — and
  the Final Round re-score masks them too.
  [forum topic 11516](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516)
* **[OFFICIAL]** DrivenData staff, 2026-09-22: *"any fault pixel not already captured by
  USGS/INGENIOUS [counts as new] and can include newly mapped geometry of an existing fault
  system."*
* **[MEASURED]** Independent corroboration: two sibling artifacts whose predictions on the
  catalogue differ enormously (mean 0.95 vs 0.11 on catalogue pixels) reported the identical
  leaderboard score 0.1563. That can only happen if those pixels are excluded pixel-exactly.

* **[OFFICIAL]** DrivenData staff, forum 11516 **post 4**, in more specific wording than the headline
  answer: the live catalogue mask is **pixel-exact and identical to the provided training labels**;
  predictions on unmasked pixels that are far from *new-fault* truth remain **fully penalized**; and
  **new-fault truth can occur within 300 m of a known trace**.
  [forum 11516, post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4)

**Consequences.** (a) Re-drawing the existing catalogue earns nothing. (b) Any local holdout whose
truth *is* the catalogue is an inverted instrument — it rewards exactly what the scorer ignores.
(c) Prediction mass on masked pixels neither earns nor radiates credit, so it must be removed
before scoring, not merely unpenalised. (d) Because new truth can sit within 300 m of a known trace,
a buffer around the catalogue is a **bet**, not a free win: it trades the chance of credit at the
mapped margin for a lower false-positive density everywhere else. This repository measures that bet
rather than assuming it (`scripts/prune_experiment.py`, `scripts/margin_check.py`) and reports the
per-block stability of the result.

## 4. Data inventory — pins re-verified byte-for-byte in this checkout

| file (repo name) | official name | bytes | sha256 | measured content |
| --- | --- | --- | --- | --- |
| `data/training_features.tif` | `gems-geodawn-numerical-features.tif` | 418,912,844 | `4371c82e3b8339b807bdffcf4ef59a225520fe2988d521be208ae33743123bc5` | 19 bands, float32, 3730×3292, EPSG:32611, 100 m, origin (243350, 4508550), nodata sentinel −3.4028235e38 |
| `data/existing_faults.tif` | `existing_faults.tif` (a.k.a. `labels.tif`) | 425,830 | `7ba308ccdc4418b31a178f4f1ef21aaa6e152e4028f2f6f64b01f7eb25ae4093` | int8, values −1 / 0 / 1; 60,988 catalogue pixels; 7,111,787 nodata |
| `data/example_submission.tif` | `example_submission.tif` (a.k.a. `sample_submission.tif`) | 1,599,597 | `2176d08e485aa2cd2860ce8df539db4faf4d76163b38a4dd8c30a40454d35cbc` | float32, NaN outside; **not all-zero** — 1.0 at exactly the 60,988 catalogue pixels |

[MEASURED] Grid: height 3,730, width 3,292, footprint 5,167,373 px, outside 7,111,787 px,
total 12,279,160 px. Command: `python scripts/prepare_data.py`.

**Band inventory**, read from the file's own TIFF tags by
`python scripts/measure_native_resolution.py` and `scripts/prepare_data.py`:

| # | band | category |
| --- | --- | --- |
| 1 | `mag_anom` | magnetic |
| 2 | `rtp` — **reduced-to-pole magnetic anomaly** | magnetic |
| 3 | `tmi_hg` | magnetic |
| 4 | `geod_2ndinv` | geodetic strain |
| 5 | `iso_grav_anom_slope` | gravity |
| 6 | `tc` — "tilt angle **or** total curvature" (ambiguous by its own wording) | magnetic |
| 7 | `geod_shearrate` | geodetic strain |
| 8 | `geod_dilaterate` | geodetic strain |
| 9 | `tmi_vg` | magnetic |
| 10 | `deq_n100a15` | seismic |
| 11 | `iso_grav_anom_vg` | gravity |
| 12 | `det_elev` — **detrended elevation** | topographic |
| 13 | `iso_grav_anom` — **isostatic gravity anomaly** | gravity |
| 14 | `tmi` | magnetic |
| 15 | `depth_to_base_surf` | subsurface |
| 16 | `ieq_n100a15` | seismic |
| 17 | `cond_surf` | subsurface |
| 18 | `iso_grav_anom_hg` | gravity |
| 19 | `det_elev_slope` | topographic |

## 5. Native resolution — the brief's own caution, quantified

[MEASURED] `python scripts/measure_native_resolution.py`. Integral scale = first lag (metres) at
which the along-axis autocorrelation falls below 1/e.

| field | ρ₁ (rows / cols) | integral scale (rows / cols, m) |
| --- | --- | --- |
| `iso_grav_anom` (13) | 0.942 / 0.924 | 5,800 / 8,100 |
| `rtp` (2) | 0.934 / 0.919 | 2,700 / 2,800 |
| `tmi` (14) | 0.935 / 0.919 | 2,900 / 2,800 |
| `det_elev` (12) | 0.940 / 0.923 | 5,400 / 7,400 |
| `det_elev_slope` (19) | **0.814 / 0.810** | **2,700 / 4,000** |
| `cond_surf` (17) | 0.942 / 0.923 | 10,500 / 10,400 |
| `depth_to_base_surf` (15) | 0.942 / 0.925 | 4,800 / 6,600 |

**Reading.** A 100 m horizontal gradient of a field whose integral scale is kilometres resolves
interpolation structure, not a contact edge. `det_elev_slope` is the only supplied field with real
power at the grid scale, and it is exactly the strongest single-layer gradient baseline. This is
the measured answer to the brief's instruction *"Check the native gravity station spacing before
trusting 100 m gradients."*

## 6. The method, from the primary literature

* **[OFFICIAL]** Gallardo & Meju (2004), *Joint two-dimensional DC resistivity and seismic travel
  time inversion with cross-gradients*, JGR 109, B03311,
  [doi:10.1029/2003JB002716](https://doi.org/10.1029/2003JB002716). Defines the cross-gradient
  function `t(x) = grad m1(x) x grad m2(x)`. Key published properties used here: `t` vanishes where
  the two fields change in the same direction; it measures *shared structure*, not shared
  amplitude; and it assumes no petrophysical relationship between the properties.
* **[OFFICIAL]** Fregoso & Gallardo (2009), *Cross-gradients joint 3D inversion with applications
  to gravity and magnetic data*, Geophysics 74(4), L31–L42,
  [doi:10.1190/1.3119263](https://doi.org/10.1190/1.3119263).
* **[OFFICIAL]** ASEG 2012 extended abstract,
  [doi:10.1071/ASEG2012ab273](https://doi.org/10.1071/ASEG2012ab273): the angle and cross product
  of the **horizontal** gradients of gravity and reduced-to-pole magnetics used directly in the
  data domain, with no inversion. This is the construction implemented in
  `src/gems49/crossgrad.py`.

## 7. What was tried, and what it measured

All **[MEASURED]**; receipts in `docs/data/field_screen.json` and `docs/data/final_experiment.json`.

See `docs/research/hypotheses.md` for the ranked hypothesis set with layers, physical signature,
novelty argument and cost. The headline results:

* the pure gravity × magnetics data-domain cross-gradient scores **at the no-information floor** —
  at matched geometry it is statistically indistinguishable from random emission;
* the same occurs for the supplied horizontal-gradient bands (18) and (3) and the
  vertical-gradient pair (11) and (9);
* anchoring the coupling on the topographic edge recovers real signal;
* the cross-gradient used as a **confounder gate** on a multi-scale, coherence-weighted
  topographic carrier is the best-measured family in this repository.

## 8. Irregularities flagged for review

1. **[MEASURED]** The organizer's `example_submission.tif` is **not** "total fault absence",
   although the problem page describes it as such. It carries 1.0 at exactly the 60,988 catalogue
   pixels. Copying the template yields a file that scores zero once the mask is applied.
2. **[MEASURED]** Three coexisting valid-region definitions: the feature stack's nodata sentinel
   (−3.4028235e38, 0 px of it inside the footprint but it defines the outside), the label raster's
   −1, and the template's NaN. This repository uses the template's, because that is what the
   submission format describes.
3. **[OFFICIAL]** The geometry of the organizer's known-fault mask *is* stated in the staff reply's
   post 4: pixel-exact, identical to the provided training labels. The earlier reading on this page
   ("not published") was superseded by that post and is corrected here. What remains unpublished is
   the **new**-fault truth itself, which is the private label set.
4. **[BLOCKED]** No leaderboard feedback. `https://www.drivendata.org/competitions/306/competition-doe-gems/data/`
   redirects unauthenticated requests to `/accounts/login/`, and the terms forbid automated
   monitoring. Every number in this repository is a proxy.
5. **[BLOCKED]** The 1 m DEM is distributed as a link list of ~100 MB tiles. Named source:
   USGS 3DEP staged products, `https://prd-tnm.s3.amazonaws.com/` (obtainable, free, official) —
   but returns HTTP 000 from this environment (measured), and the full region is terabytes.
6. **[MEASURED]** Band 6's own description is ambiguous: *"Tilt angle or total curvature -
   magnetic field derivative for edge detection"*. The band is used only as a baseline, never as
   a required input.
