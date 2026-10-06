# Independent audit of the separate PR #1 candidate — 2026-10-06

## Decision

**Keep the submission gate closed.** The merged PR #1 artifact is a real, non-identical 40,000-pixel GeoTIFF and its NaN-outside twin meets the official raster-format contract. Its four-fold DTI receipt reproduces exactly from the committed method and the locally available, hash-pinned input mirrors. Those checks do **not** establish a submission-ready result: the holdout is component-sampled rather than a true spatially held-out block test, the emission mask is built from the full catalogue including the held-out truth, the chosen gravity smoothing is 300 m while the relevant official regional gravity product is 1 km, and the uniqueness corpus is incomplete and not tracked.

No submission slot was used. This audit did not write, rename, or copy a candidate GeoTIFF into the submission directory. The PR #1 GeoTIFFs remain as pre-existing research artifacts from the already merged PR; they are not a recommendation to upload.

## 1. What was audited

- Remote `origin/main` was fetched at merge commit `2e22802`; its PR #1 candidate was a separate artifact, not the H49-XG raster from this task branch.
- The candidate files are `docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-zeros.tif` and its `-nan.tif` twin. The local prediction inputs were available under `.cache/gems_data/`; the restore receipt identifies public owner-hosted mirrors and explicitly says the hashes authenticate mirror consistency only.
- Candidate execution was reconstructed in `/tmp/gemsdoe49-pr1-review` from `git archive origin/main`. No branch switch was made.
- The candidate file was independently read with Rasterio, its hashes and geometry were checked, the final DTI was recalculated, the score field and NMS output were recomputed without writing a TIFF, and the H33-2-B2 raster was fetched read-only from the pinned `GEMSDOE32` Git commit for a direct pairwise comparison.

## 2. Official competition rules rechecked

### Output format

The official [DrivenData problem and submission-format page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) states that the output must be a single-layer float32 GeoTIFF; predictions are in `[0,1]`; CRS, resolution and bounds match the training grid; and data outside the bounds is **null or NaN**. The 100 m UTM 11N/EPSG:32611 geometry and 300 m triangular DTI kernel are also stated there.

Consequently, the `-nan.tif` twin is the compliant file to inspect. The `-zeros.tif` twin has 7,111,787 zero-valued cells outside the finite sample footprint and does not meet the stated null/NaN outside condition. The PR #1 site's former description of that zeros file as the “safest” primary was incorrect; the reconciled site links the NaN-outside twin and labels it research-only.

### Catalogue mask

The [DrivenData staff clarification, forum thread 11516 post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4), says the mask is pixel-exact and identical to the provided training fault labels; predictions near a known trace but far from new-fault truth remain penalized; and new-fault truth can occur within 300 m of known faults. The PR #1 metric code's scoring mask excludes residual catalogue pixels from TP/FP/FN and computes false-positive proximity from held-out truth only. That scoring arithmetic follows the clarification. The holdout *prediction construction* has a separate issue described below.

### Gravity source scale

The official [INGENIOUS GDR 1391 page](https://gdr.openei.org/submissions/1391) links to the USGS regional geophysical map release [DOI 10.5066/P9Z6SA1Z](https://doi.org/10.5066/P9Z6SA1Z). The attached release README identifies `GB_iso_grav_anom.tif` as a 1 km grid. USGS [Nevada Bouguer notes](https://pubs.usgs.gov/ds/2006/234/nv_boug.htm) also describe a gravity compilation gridded at 1 km. The competition mirror's band 13 is named `iso_grav_anom`, but the organizer download is login-walled and no source crosswalk proves that band 13 is this exact release/version. Therefore 1 km is a material, official scale warning—not a verified lineage claim. The PR #1 primary uses Gaussian `σ=3` pixels = 300 m, finer than the stated cell size of the associated regional product. The empirical variance-retention and row-semivariogram measurements in `gravity_support_diagnostics.json` are measurements of the rasterized mirror, not station-spacing estimates, and do not resolve that mismatch.

## 3. Raster bytes and geometry

The NaN-outside candidate was read from disk and checked against the locally mirrored sample grid:

| Check | Independent result |
| --- | --- |
| SHA-256 / bytes | `bef0881aadc8366f1a08e13aac119cfa71492750e93426105b9e5745c7d701bf` / 419,609 bytes |
| Bands / dtype / dimensions | 1 / float32 / 3,730 × 3,292 |
| CRS / transform / bounds | EPSG:32611; 100 m; `(243350, 4508550)` origin; transform and bounds equal the sample template |
| Nodata metadata | None |
| In-footprint values | 5,167,373 cells; all finite; minimum 0, maximum 1; 40,000 positive cells, all exactly 1 |
| Outside-footprint values | 7,111,787 cells; all NaN, no zeros |
| Catalogue intersection | 0 emitted positive cells on the 60,988 known-label pixels |
| Twin comparison | The NaN and zeros TIFFs have identical values throughout the footprint; their outside values differ as intended |

The zeros twin hash is `154003d9b7220c43a96052364a4380ce8b602a9f7e45c549354b1f73353ef867`; it is not the preferred file under the official outside-bounds requirement.

## 4. Method and holdout review

### What the code does

`src/gems49/xgrad.py` reads the local mirror's band tags (RTP band 2, detrended elevation band 12, isostatic gravity anomaly band 13), smooths each field, computes horizontal gradients, scores sign-invariant pairwise alignment, and requires all three robustly scaled magnitudes via a minimum. The resulting alignment score is a plausible implementation of the requested alignment operator; the local band names/tags were read directly and match the declared bands. The supplied orthogonality diagnostic does **not** match its comment: it computes `cross_gradient² / (|g_i|² |g_j|²) = sin²(theta)`, without a magnitude factor. This diagnostic is excluded from the candidate score, so the defect does not change this raster, but the report/code description should not call it magnitude-weighted until corrected.

`src/gems49/emission.py` applies 2 px NMS, excludes exact catalogue pixels, ranks and emits 40,000 unit-valued dots, then writes both outside-zero and outside-NaN forms. The exact-catalogue exclusion does not add a buffer; that is consistent with the pixel-exact staff clarification. Its `np.argpartition` selection has no stable tie-break although parts of the code/documentation call ordering deterministic. A fresh run under another NumPy version may select different cells where scores tie. Do not claim byte-for-byte reproducibility of the NMS without fixing and testing tie order.

### Holdout construction is not the required spatial-block test

`src/gems49/holdout.py::hide_and_recover_folds` assigns connected-component centroids to four quadrants, then randomly hides about 25% of the components within each quadrant. It does **not** hold out an entire quadrant or an isolated geographic block. `scripts/run_pipeline.py` makes one full-domain 40,000-dot raster and scores that same raster in each fold; unlabelled pixels outside the fold's quadrant still contribute false-positive mass in each fold.

More importantly, the sweep calls `nms_topk(..., catalogue, ...)` with the complete catalogue before scoring. The fold then treats some pixels of that same catalogue as hidden truth and unmasks them in `score_mask`. Thus candidate construction uses the held-out labels to forbid predictions on those exact pixels, even though they are scored as unknown truth. This may suppress direct hits (and is not necessarily optimistic), but it is not a faithful hide-the-label-and-predict protocol and can alter candidate placement and rankings. A valid block holdout must hide the fold labels from the emission mask, use only the remaining visible catalogue as the exact scoring mask, and confine predictions/false-positive accounting to the fold domain with a boundary guard.

The pooled/per-fold arithmetic itself was independently reproduced for the shipped final raster: DTI `0.02026648455187942`; folds `0.036504362325797846`, `0.012602545366713533`, `0.016012338026837405`, `0.015946692488168893`. This equals `evidence/final_holdout.json` exactly. It is a result of the PR #1 protocol, not a spatially independent estimate or leaderboard score.

The code selects `σ=3`, 40k on the same four-fold results. The best-config sweep receipt's raw-score mean is `0.020302617885093063`; after the final q01/q999 clipping and NMS, the actual artifact gives `0.02026648455187942`. The code's “rank-preserving” comment is not strictly true: clipping the top 0.1% to one creates ties. Reconstructing the raw-score mask produced 36,424 of the same 40,000 pixels as the shipped final mask.

The four-fold receipt reports single-layer means `0.0156134 / 0.0144858 / 0.0145133`. A separate replay in this environment (Python 3.11.2, NumPy 2.4.6, SciPy 1.17.1, Rasterio 1.4.4) produced `0.0157704 / 0.0144270 / 0.0145704`; all three remain below the shipped candidate's `0.0202665` on this component-sampled protocol, but exact baseline margins differ. The NMS path uses `argpartition` without a stable tie-break and dependencies are not fully pinned; that is a reproducibility defect to fix before future candidate scoring.

## 5. Uniqueness review

The PR #1 JSON receipt reports 57 comparisons, no recorded byte-hash collision, maximum smaller-set overlap `0.290125`, and maximum *signed* Spearman `−0.26557`. Its largest absolute Spearman magnitude is about `0.98880` and is negative; the signed maximum and absolute maximum are different summaries. The available audit does include the H33-2-B2 file. I fetched that prior from the pinned `GEMSDOE32` commit and independently reproduced:

- H33 SHA-256: `c55bafc470054e8271dcb89347a17e07fefe50de6af6e6ba6c4b169ef7ab6fa9` (37,654 positive cells).
- Candidate/H33 intersection: 792 cells; overlap of the smaller set `0.0210336`; union Jaccard `0.0103042`; Spearman over the union of positive cells `−0.979583`.

That is strong evidence that these two particular rasters are not relabels. It does not close the all-priors gate: the 57 input TIFFs and their manifest are absent from the merged Git tree, and the PR #1 site explicitly says three earliest files are unavailable and substitutes nearby archived files. The uniqueness receipt is therefore a reported audit, not independently reproducible exhaustive evidence. The candidate's unique status against every original submission remains unproven.

## 6. Overall verdict and corrections made

| Gate | Result |
| --- | --- |
| Official NaN-outside format | **PASS for `-nan.tif` only**; zeros twin does not meet the stated outside condition |
| Candidate file hash/geometry/range/catalogue pixels | **PASS**, independently re-read |
| PR #1 final four-fold receipt | **PASS**, exactly reproduced |
| Single-layer comparison on that same component-sampled proxy | Candidate is higher in both the receipt and current replay; exact margin is version/tie-sensitive |
| User-required spatially-blocked holdout | **NOT ESTABLISHED** by PR #1's four-quadrant component protocol |
| Gravity native-scale requirement | **OPEN**; associated official regional product is 1 km, exact competition-band lineage unverified; PR #1 uses 300 m σ |
| Uniqueness against every prior | **PARTIAL**; H33 verified directly, other 56 receipt rows not independently re-read and three original priors are unavailable |
| Organizer score / improvement over 0.2778 | **NONE MEASURED**; no PR #1 score or mapping to the 0.2778 board row |

The site and README have been reconciled to (a) put the NaN-outside twin first as a **research-only** download, (b) state that it is not cleared to upload, (c) remove the unsupported claim that the H33-2-B2 file earned 0.2778, (d) correct the outside-zero format claim, and (e) preserve the separate failed H49-XG primary result. No new submission raster was generated.

## 7. Remaining work before any submission

1. Freeze and hash a new spatial-block holdout that hides full blocks and held-out labels from candidate construction; apply only the visible-catalogue mask and score only the held-out block plus its defined guard. Compare the exact proposed operator and all single-layer baselines at a fixed mass; do not retune on that same holdout.
2. Confirm the competition gravity band's source/version or use an explicitly conservative scale justified by official metadata; then validate the fixed primary on fresh blocks. Do not promote the failed 500 m primary by selecting the already-seen 2 km sensitivity.
3. Recover or enumerate each original prior artifact, commit a lightweight manifest of immutable remote refs/hashes, and make the uniqueness tool's tie-breaking and metrics reproducible.
4. Correct the diagnostic orthogonality weight and the NMS stable tie-breaking; add tests for NaN-only outside semantics and promotion-gate enforcement. The end-to-end PR #1 pipeline currently can continue to file writing after only printing a baseline warning and does not enforce the full holdout/control/uniqueness gate.
5. Only then generate one new unique NaN-outside GeoTIFF and consider using a weekly slot. No board/private performance claim is justified by the current evidence.
