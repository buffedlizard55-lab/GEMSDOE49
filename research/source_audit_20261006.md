# Official-source and prior-result audit — 2026-10-06

## 1. The reported `h33-h33-2-b2` / 0.2778 attribution is unresolved

### What the official live board says

The [DrivenData public leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) was read on 2026-10-06. It displayed **xiaofanhu** at **rank 1** with **0.3774**, and **extradr19** at **rank 13** with best public DW-Tversky **0.2778** and **10 submissions**. The public table names participants and scores, but does not expose submission filenames, candidate slugs, or model recipes.

### What the H33 project says about its file

At the inspected [GEMSDOE32 README revision `b983924b57781edd29b8e249c4923bf33d9902f6`](https://github.com/buffedlizard55-lab/GEMSDOE32/blob/b983924b57781edd29b8e249c4923bf33d9902f6/README.md), the file `gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif` is explicitly marked `UNSCORED`; its README describes **0.2747 as a projected/modelled value**, not an organizer score, and states that no organizer score exists for artifacts in that repository. The repository's audit JSON gives the unrounded projection `0.27467315976932155` and the TIFF SHA-256 `c55bafc470054e8271dcb89347a17e07fefe50de6af6e6ba6c4b169ef7ab6fa9`; streaming the pinned Git blob reproduces that digest. The project identifies `GEMSDOE32-H33-2-B2` as its proposed portal name, but that is not a DrivenData receipt.

### Conclusion

**Do not equate the public 0.2778 board row with the H33-2-B2 file.** The score is a real public-board entry under the same participant name appearing in the owner-maintained project history, but no public record connects that specific score to this specific file. The H33 README's 0.2747 is a projection and cannot be used as evidence for the 0.2778 score. A portal submission receipt or owner-authenticated submission mapping would be needed to resolve attribution; neither is available in this workspace.

The public leaderboard is not the final private score. The [official challenge problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) says the public test labels are only a public portion and that final decisions use a private set of new faults; the current public DTI cannot establish performance on the private test or the later expert-expanded label set.

## 2. Verified official rules and challenge details

The sources below were checked directly on 2026-10-06. The [September 2026 GEMS Prize Official Rules](https://docs.nlr.gov/docs/fy26osti/96647.pdf) are the organizer rules document; the live DrivenData challenge pages are authoritative for current competition-platform details.

| Claim used in this work | Official source and inspected wording |
|---|---|
| Target is faults, including faults absent from the existing public dataset; labels may be incomplete/inaccurate. | [DrivenData problem description](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/): describes publicly available known faults and a manually identified set of new faults as the hidden test. |
| Inputs are a 100 m, EPSG:32611 multiband GeoTIFF. | Same problem description, “Provided features”: `training_features.tif`, projected UTM zone 11N (EPSG:32611), 100 m resolution. |
| Required output is a single-band float32 GeoTIFF, [0,1], same projected CRS/resolution/bounds; outside area null/NaN. | Same problem description, “Submission format”: single layer, float32, confidence/probability in [0,1], same bounds, null/NaN outside. |
| DTI has a 300 m triangular support, α=0.2, β=0.8; false negatives are penalized more. | Same problem description, “Performance metric”; it gives TPw/FPw/FNw formulas and the worked example. The displayed 0.60 is rounded; the formula for (3.00,1.89,2.00) gives 0.6026516673 before rounding. |
| Up to three automated feedback submissions a week, but one final submission across both prize rounds. | [Official rules PDF](https://docs.nlr.gov/docs/fy26osti/96647.pdf), §§3.4–3.5. The final entry is evaluated on the private initial-round labels and later on the expanded expert-reviewed labels. |
| Generative-AI use is allowed but must be disclosed in the narrative. | Same official rules PDF, §3.2. The participant remains responsible for accuracy/authorship representations. |
| Current live public board entry. | [DrivenData leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/), checked 2026-10-06: extradr19 #13, 0.2778, 10 submissions. The board does not identify the specific submitted raster. |
| Official feature-data access is login-gated in this environment. | [DrivenData data tab](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) redirected to the official login page. Restored files therefore remain owner-hosted mirrors, not organizer-authenticated inputs. |

## 3. Gravity spacing: official source and mirror-specific limitation

- [USGS DS 234 Nevada Isostatic Gravity](https://pubs.usgs.gov/ds/2006/234/nv_iso.htm) states that the isostatic residual gravity grid was derived from complete Bouguer gravity by removing the modeled compensating mass for topographic loads. It identifies the source station compilation but does not state the native spacing of the competition mirror's particular band.
- [USGS DS 234 Nevada Bouguer Gravity](https://pubs.usgs.gov/ds/2006/234/nv_boug.htm) states that the complete Bouguer grid used **71,055** station measurements in and adjacent to Nevada and was converted to a **1 km grid** using minimum-curvature gridding.
- The USGS [DS 234 public data directory](https://pubs.usgs.gov/ds/2006/234/data/) exposes `grviso.grd` and `grvcba.grd`, but no evidence here proves the competition's `iso_grav_anom` is byte-for-byte or spatially sourced from that exact DS 234 product.
- Locally observed 100 m autocorrelation and all tested lags are preserved in [`input_spacing_audit_20261006.json`](input_spacing_audit_20261006.json). This measures the smoothness of the resampled 100 m mirror; it is **not** a station-spacing estimate.

## 4. Provenance, source integrity, and irregularities

1. The [official DrivenData data download page](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) redirects to login; this environment has no portal credentials. The local 19-band stack, labels, and template came from [pinned public owner-hosted mirrors](../data/manifest.json). SHA-256 pins prove only that the restored files match those mirrors; they do not authenticate the mirrors against the organizer download.
2. The local feature stack's band tags confirm `rtp` (2), `det_elev` (12), and `iso_grav_anom` (13). Related provided bands include `iso_grav_anom_slope` (5), `iso_grav_anom_vg` (11), `iso_grav_anom_hg` (18), and `det_elev_slope` (19); the experiment intentionally differentiates the named scalar bands and does not treat those derivative bands as independent observations.
3. The official sources do not map `extradr19`'s 0.2778 to `h33-h33-2-b2`. The owner-authored README says the H33 file is unscored. This is a provenance/attribution uncertainty, not proof of misconduct or an irregular leaderboard event.
4. No prediction from this project has been submitted to DrivenData. The runner uses a public existing-fault catalogue proxy, which the official problem page says is incomplete and can contain inaccuracies; it cannot validate private test accuracy.
5. The preregistration requested whole-segment folds, but the implementation uses east–west raster blocks because the local mirror exposes a binary label raster rather than authoritative segment IDs. A follow-up 8-connectivity audit found 21 raster components with pixels in two scored fold cores (850 positive pixels); connectivity is not a geological segment identifier. Reported results are therefore spatial-block results, not verified whole-segment cross-validation. See [`label_fold_component_audit_20261006.json`](label_fold_component_audit_20261006.json).
