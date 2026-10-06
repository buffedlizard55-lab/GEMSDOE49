# GEMSDOE49 — cross-gradient structural coupling for the DOE GEMS Prize

**Competition:** [DOE GEMS Prize, DrivenData #306](https://www.drivendata.org/competitions/306/competition-doe-gems/)
· **Problem description:** [page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
· **Staff mask clarification:** [forum 11516, post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4)
· **Rules PDF:** [docs.nlr.gov/docs/fy26osti/96647.pdf](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
· **Live site:** <https://buffedlizard55-lab.github.io/GEMSDOE49/>
· **Weekly budget:** 3 scored submissions; **one** file is scored in *both* prize rounds.

> ### Core Values — kept central to every decision in this repository
>
> **Maximize P(Win).** *"Maximize the Probability of Winning"* is our decision-making framework.
> In every decision we weigh tradeoffs, assess risk, and choose the path that maximizes the
> probability that we win this prize. We set aside our emotions and make tough decisions in order
> to maximize P(Win). It frees us from constraints and clarifies that we must put the outcome first.
>
> **Own the Outcome.** We own results end to end — not just our individual slice of the work. When
> problems arise and we have the means to act, we act without waiting for permission or assignment.
> We treat failure and success as signals and use them to improve. We stay accountable to the final
> outcome.
>
> Applied here: *Maximize P(Win)* is why the brief's headline construction (the Gallardo–Meju
> cross-gradient of gravity × magnetics) is **reported as a measured rejection** instead of being
> quietly shipped — a submission slot is scarce, and a construction that scores at the random floor
> must not consume one. *Own the Outcome* is why every number below is either **[OFFICIAL]** (read
> from an organizer page, with a link), **[MEASURED]** (computed here from hash-verified bytes, with
> the script that produced it), or **[BLOCKED]** (with the reason), and why a candidate that fails a
> gate is labelled failed rather than renamed.

---

## ⬇ ONE-CLICK SUBMISSION FILE

**[Download the submission](docs/downloads/gems49-gate_ortho_w0.25-40k-20261006T213721Z-nan.tif)**
— `gems49-gate_ortho_w0.25-40k-20261006T213721Z-nan.tif`, 377,092 bytes,
sha256 `c58fe0551e120678fff190239146985d0b2f7255255535fb5fb637a05cb40af3`.

Its twin
[`…-zeros.tif`](docs/downloads/gems49-gate_ortho_w0.25-40k-20261006T213721Z-zeros.tif) contains the
**same predictions** with `0.0` instead of `NaN` outside the footprint — use it if the form rejects
the `-nan` one (see the format section below). A single-file `.zip` of the primary variant sits next
to them. Sizes, hashes and the submission-form note are in
[`docs/data/submission.json`](docs/data/submission.json) and printed on the
[site](https://buffedlizard55-lab.github.io/GEMSDOE49/).

**Submission form fields** (the two boxes to fill):

| Field | Value |
| --- | --- |
| *File to submit* | the `.tif` above (or its `.zip`) |
| *Note* | the `comment` string in [`docs/data/submission.json`](docs/data/submission.json) |

**If the form says `Predicted values must be in range [0, 1]`** — you are hitting one of exactly two
conditions, both reproduced as constructed unit tests in `tests/test_gate.py` and both rejected by
the local gate before upload:

* **mode A** — a `NaN` **inside** the footprint (outside the footprint is fine: the organizer's own
  `example_submission.tif` is `NaN` there). A naive `min >= 0 and max <= 1` check fails on `NaN`.
* **mode B** — a finite value outside `[0, 1]`, e.g. a `-3.4e38` nodata sentinel left in the array.

`scripts/validate_submission.py` checks both plus CRS, resolution, shape, geotransform, dtype, band
count, and byte-level agreement of the `NaN` mask with the official template. The `-zeros` variant
exists so that a form which converts `NaN` to a large negative number still sees a legal file.

---

## THE STANDING BRIEF (read this first, every session)

> Review the repo.
> THE FOLLOWING IS THE HIGHEST URGENCY AND MUST BE FOLLOWED!
> MUST GENERATE A UNIQUE TIF SUBMISSION FOR THE COMPETITION. DO NOT COPY A PREVIOUS SUBMISSION
> UNLESS IT'S FOR LEARNING AND EDUCATION. BUT WE MUST GENERATE A UNIQUE TIF SUBMISSION.
>
> There should be an easy to download submission tif file as described by the prompt. Read the entire
> prompt.
>
> Cross-gradient structural coupling of gravity, magnetics and topography. Treat agreement across
> independent signals as an operator, not a product of scores. Gallardo and Meju (*Journal of
> Geophysical Research*, 2004, [doi:10.1029/2003JB002716](https://doi.org/10.1029/2003JB002716))
> defined the cross-gradient function, the cross product of two property fields' gradient vectors.
> It vanishes where both fields change in the same direction, so it measures shared structure
> without assuming a petrophysical relationship. Fregoso and Gallardo (*Geophysics* 74(4), L31–L42,
> 2009, [doi:10.1190/1.3119263](https://doi.org/10.1190/1.3119263)) applied it to gravity and
> magnetics. An ASEG 2012 extended abstract
> ([doi:10.1071/ASEG2012ab273](https://doi.org/10.1071/ASEG2012ab273)) used the angle and cross
> product of the horizontal gradients of gravity and reduced-to-pole magnetics directly in the data
> domain, with no inversion. Compute horizontal-gradient vectors of the isostatic gravity anomaly,
> the reduced-to-pole magnetic anomaly and a smoothed elevation surface at a common resolution.
> Score each pixel by magnitude-weighted alignment (both gradients large and parallel). Keep
> magnitude-weighted orthogonality as a separate confounder flag, not as evidence. A fault
> juxtaposing rocks of different density and susceptibility puts both edges on one trace, while a
> buried volcanic edge usually moves the magnetic edge alone, so the coupled score should suppress
> single-layer edges. Check the native gravity station spacing before trusting 100 m gradients, and
> expect density steps at basin margins. The candidate must beat the repo's single-layer gradient
> baseline on hide-and-recover holdout segments before it is allowed near a submission. Normalize to
> [0,1], apply the repo's metric-aware placement, and write the required float32 GeoTIFF
> (EPSG:32611, 100 m, NaN only outside the footprint and none inside). Before download, hash the
> raster and compare rank correlation and top-k overlap against every prior submission, and reject
> it as a relabel if it isn't clearly different.
>
> We need to study, analyze, and understand the highest score from the GEMSDOE sites where the
> submission TIF is downloaded from: GEMSDOE32, `h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros:
> 0.2778`. Why and how did this get the highest score and are we able to generate a submission that
> scores higher than 0.2778? Answer the question using PhD-level experience, knowledge, and
> judgement. Then use the answer to generate a unique TIF submission into the competition. Must be
> unique submission unlike any within the GEMSDOE sites. Verify working line by line, no
> hallucinations.
>
> Before implementing, generate 3–5 candidate geological hypotheses we haven't tried yet, each
> naming: the specific layer(s) involved, the physical signature being targeted, why it should catch
> a fault missing from the USGS/INGENIOUS catalogue rather than one already in it, and how it
> differs from anything already implemented in this repo. Rank them by expected DTI improvement and
> implementation cost. Validate the top candidate on our spatially-blocked holdout set before
> touching a weekly submission slot. If a candidate can't be validated without new external data,
> name the specific free, official source needed and check it's obtainable before proposing the idea
> as viable.
>
> The goal of this project is to place top of the leaderboard (current top: 0.3195). We need to
> create a project that can compete and place top of the leaderboard. Create a GitHub page for this
> repo that has a clean UI, is organized, user friendly, simple and easy to use, with all relevant
> information in an easy-to-read format and official verified links as sources. The site must make
> the submission TIF obvious and easy to download (one click), include an executive summary
> explaining exactly how to make a submission, and fix the portal error *"Predicted values must be
> in range [0, 1]"* by guaranteeing every value is in [0,1] with no sentinel nodata. Give the
> submission a unique name and a short comment for the submission on DrivenData so we can identify
> it.
>
> Work line by line verifying from official verified trusted sources, provide links for manual
> review. No manual input — work autonomously. Flag any irregularities for review. No
> hallucinations. Run this task through multiple passes (implement → review for bugs/missing
> requirements → re-check against the original request). Create a pull request and merge it onto
> main. Make suggestions for remaining work and limitations.
>
> Core values: **Maximize P(Win)** — every decision weighed for probability of winning;
> **Own the Outcome** — own results end to end, treat failure and success as signals.

---

## Executive decision — 2026-10-06 (this session, after the merge)

**One GeoTIFF is cleared for a DrivenData slot by this repository's own preregistered gate:
`gems49-gate_ortho_w0.25-40k-20261006T213721Z`.** No slot has been used: this environment has no
DrivenData session (`/data/` redirects to `/accounts/login/`), so nothing here has a leaderboard
score and no claim of one is made.

| Gate | Result |
| --- | --- |
| Spatially blocked holdout, leave-one-fold-out family selection (4 folds, 3,118 catalogue segments removed whole) | candidate `gate_ortho_w0.25` **0.1074** vs best single-layer baseline `single:det_elev_slope` **0.0991** → **PASS** (`docs/data/final_experiment.json`) |
| Same placement rule given to **both** arms at the shipped geometry | collared **0.1131 vs 0.1029**, unrestricted **0.1181 vs 0.1111** (`docs/data/margin_check.json`) |
| Per-block stability floor (no held-out block below −0.01; the floor used by the earlier merged session) | **4/4 blocks won**, worst block delta **+0.00335** collared, **+0.00154** unrestricted → **PASS** |
| Format gate on the committed bytes (13 checks, both encodings) | **11/11 required checks pass** |
| Relabel check against 131 prior submission rasters | max support Jaccard **0.0223**, max top-25k Jaccard **0.0161**, max Spearman **−0.414** → **UNIQUE** |

**Correction carried into this session from the earlier merged work.** The staff clarification is
more specific than the version used earlier in this session: the live catalogue mask is
**pixel-exact** (identical to the provided training labels), *and* **new-fault truth can occur within
300 m of a known trace**, while unmasked predictions far from new-fault truth remain fully penalized
([forum 11516, post 4](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4)).
The shipped 300 m catalogue buffer is therefore a *bet* that pays only if new truth is not
concentrated within 300 m of known traces. It is measured to help on both holdout readings
(`docs/data/prune_experiment.json`) and to be the setting at which the candidate wins all four blocks
(`docs/data/margin_check.json`), but the bet is stated here as a bet, not as a fact.

---

## What this session measured

Everything below is **[MEASURED]** in this checkout from the hash-verified official bytes; receipts
are in `docs/data/*.json` and the site.

### 1. The brief's headline construction does not survive contact with this grid

The pure Gallardo–Meju / ASEG-2012 data-domain cross-gradient of isostatic gravity against
reduced-to-pole magnetics — implemented exactly as described, magnitude-weighted orthogonality, no
inversion — scores **inside the no-information control band** (`cg_gxm` 0.0715 vs `control:random`
0.0718 and `control:uniform_lattice4` 0.0693 at matched geometry; the strict random-support control
in `docs/data/final_experiment.json` reads 0.0549). See `docs/data/field_screen.json` and
`docs/data/final_experiment.json`.

### 2. The reason is measurable, and the brief told us to look for it

The brief says *"Check the native gravity station spacing before trusting 100 m gradients."* Measured
here from the official bytes (`scripts/measure_native_resolution.py`, results in
`docs/data/native_resolution.json`):

| field | lag-1 autocorrelation ρ₁ (rows / cols) | integral scale (m) |
| --- | --- | --- |
| `iso_grav_anom` (13) | 0.942 / 0.924 | 5,800 – 8,100 |
| `rtp` (2) | 0.934 / 0.919 | 2,700 – 2,800 |
| `tmi` (14) | 0.935 / 0.919 | 2,900 – 2,800 |
| `cond_surf` (17) | 0.942 / 0.923 | 10,400 – 10,500 |
| `det_elev` (12) | 0.940 / 0.923 | 5,400 – 7,400 |
| **`det_elev_slope` (19)** | **0.814 / 0.810** | **2,700 – 4,000** |

The measured decorrelation length of the gravity field is ~6–8 km, which independently corroborates
the documentation the earlier merged session surfaced: the official USGS regional
[`GB_iso_grav_anom.tif`](https://pubs.usgs.gov/ds/2006/234/nv_iso.htm) has a **1 km** grid cell, and
competition band 13 carries a matching name, so the 100 m pixels are resampled from a ~1 km product.
A 100 m horizontal gradient of a kilometre-scale field is interpolation structure, not a contact
edge — which is exactly why the pure gravity × magnetics cross-gradient scores at the floor.
`det_elev_slope` is the only supplied field with real power at the grid scale, which is why it is
the strongest single-layer baseline.

### 3. The operator still earns its place — as a *gate*, not as a score

The brief's own reading is the one that works: the coupling's job is to **suppress single-layer
edges**, not to generate the signal. As a multiplicative orthogonality gate on a multi-scale,
coherence-weighted topographic edge it is a **measured no-op**: the gated and ungated variants share
36,602 of 40,000 dots (support Jaccard **0.8434**, `scripts/measure_gate_effect.py`) and differ by
**+0.00003** in holdout score. That is reported as what it is — a confounder filter that does not
damage the carrier — and not as a gain.

### 4. The measurable gain is a placement rule, not a new signal

The artefact the brief asks about — `h33-h33-2-b2` from the sibling project, whose
[audit JSON](docs/research/knowledge-base.md) and local copy were examined in this session — is
described there as *"0.2708 base with every dot at d(catalogue) ≤ 2 px deleted (37,654 dots)"*. The
same rule is reproducible here (`scripts/prune_experiment.py`), and pinning the optimum at the
scoring kernel's own radius (`scripts/margin_check.py`):

| catalogue buffer | candidate collared / unrestricted | baseline collared / unrestricted | margin | blocks won |
| --- | --- | --- | --- | --- |
| none | 0.1074 / 0.1134 | 0.0991 / 0.1077 | +0.0083 / +0.0056 | 3/4, 3/4 |
| 0 px (exact pixels only) | 0.1064 / 0.1127 | 0.0980 / 0.1068 | +0.0083 / +0.0059 | 4/4, 3/4 |
| 200 m | 0.1100 / 0.1166 | 0.1010 / 0.1105 | +0.0090 / +0.0062 | 4/4, 4/4 |
| **300 m (shipped)** | **0.1131 / 0.1181** | **0.1029 / 0.1111** | **+0.0102 / +0.0070** | **4/4, 4/4** |
| 400 m | 0.1156 / 0.1164 | 0.1047 / 0.1086 | +0.0108 / +0.0078 | 4/4, 4/4 |

A placement rule handed to one arm only is not a comparison, so the baseline was given it too before
the margin was believed. Beyond 300 m the *unrestricted* score of both arms falls: the buffer starts
deleting the dots that would have earned credit where a new fault is mapped against a known one.

### 5. What was actually shipped

`gate_ortho_w0.25` — multi-scale topographic edge × orientation coherence × orthogonality gate,
40,000 dots at 300 m separation, 300 m catalogue buffer — selected by leave-one-fold-out on a
spatially blocked, leakage-controlled holdout (leakage probe TP_w = 0, null probe DTI = 0, on-truth
oracle 1.0000, random and uniform-lattice controls at matched geometry).

---

## Artifacts from the earlier merged session (preserved, not deleted)

Two distinct H49-XG results exist in this checkout and must not be conflated with the artifact above:

1. **Fixed five-block H49-XG primary — not cleared by its own preregistered gate.** Pooled
   public-catalogue proxy DTI `0.10041455` vs `0.08193976` for the best single-layer baseline, but
   one paired block delta was `−0.01180243`, beyond the allowed `−0.01`, so the run generated **no
   GeoTIFF**. Receipt: [`research/h49_xg_result_20261006.md`](research/h49_xg_result_20261006.md).
2. **PR #1 candidate `GEMSDOE49-XGRAD-ALIGN`** —
   [`docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-nan.tif`](docs/downloads/gemsdoe49-xgrad-coupled-20261006T203447Z-nan.tif),
   sha256 `bef0881aadc8366f1a08e13aac119cfa71492750e93426105b9e5745c7d701bf`, 40,000 positive cells,
   four-fold proxy DTI `0.02026648455187942`. Status there: **RESEARCH ONLY — NOT
   SUBMISSION-CLEARED** (its four-fold protocol was not a spatially isolated block holdout). Audit:
   [`docs/reports/pr1_candidate_audit_20261006.md`](docs/reports/pr1_candidate_audit_20261006.md).

Both files remain in `docs/downloads/` and are listed on the site. Nothing was overwritten.

---

## Data and source-provenance rules

The official DrivenData data tab requires a registered login and redirects unauthenticated access to
the login page; this environment has no such session. [`data/manifest.json`](data/manifest.json)
pins public, owner-hosted mirrors and exact byte counts/hashes, and
[`data/README.md`](data/README.md) explains restoration. Those hashes prove mirror consistency, not
that the bytes are identical to the organizer downloads; the provenance caveat is repeated in every
receipt. **Gravity resolution remains a hard caution:** do not interpret the 100 m resampled band 13
pixels as independent gravity observations, and do not treat a 300 m derivative scale as established
— the source product is 1 km-class, which is also what this session measured independently.

---

## Repository map

| Path | What it is |
| --- | --- |
| `src/gems49/spec.py` | grid constants and sha256 pins, all measured from the official bytes |
| `src/gems49/metric.py` | the official DTI, verified against the page's worked example and a literal brute force, plus the earlier session's reference implementation |
| `src/gems49/crossgrad.py`, `gating.py`, `xgrad.py` | horizontal gradients, magnitude-weighted alignment, orthogonality flag, the cross-gradient gate |
| `src/gems49/holdout.py`, `emission.py`, `io.py`, `screen_util.py`, `baselines.py`, `methods.py` | folds and leakage probes, metric-aware placement, fail-closed raster IO, screens |
| `src/gemsdoe49/` | the earlier session's package (features, placement, raster, holdout, metric) |
| `scripts/download_competition_data.sh`, `prepare_data.py` | fetch + sha256-verify + audit the official rasters |
| `scripts/final_experiment.py` | the leave-one-fold-out family screen |
| `scripts/prune_experiment.py`, `margin_check.py` | the catalogue-buffer sweep and the both-arms margin |
| `scripts/measure_native_resolution.py` | the "check the native station spacing" measurement |
| `scripts/measure_gate_effect.py` | pixel-level comparison of the gated and ungated supports |
| `scripts/metric_facts.py` | derives the metric algebra the site quotes |
| `scripts/build_submission.py` | normalise → metric-aware placement → write float32 |
| `scripts/validate_submission.py` | the hard format gate (13 checks, both rejection modes as tests) |
| `scripts/verify_committed_submission.py` | re-hashes the committed artifact; runs unattended in CI |
| `scripts/uniqueness_check.py`, `uniqueness_audit.py` | relabel checks against prior submissions |
| `scripts/build_site.py`, `docs/` | the generated GitHub Pages site and its evidence JSON |

## Quickstart

```bash
bash scripts/download_competition_data.sh      # fetch + verify the official rasters (mirrors, hashed)
python3 scripts/prepare_data.py                # audit grid, CRS, dtype, footprint, band tags
python3 -m pytest tests -q                     # 72 tests: metric(s), gate, cross-gradient, holdout, placement
python3 scripts/final_experiment.py            # the LOFO family screen (~20 min, 1 core)
python3 scripts/build_submission.py --family gate_ortho_w0.25 --n-target 40000 --min-sep 3 \
    --exclude-catalogue-px 3 --out-dir docs/downloads     # NEVER pass --purge: other artifacts live here
python3 scripts/validate_submission.py docs/downloads/<file>.tif --json docs/data/submission_gate.json
python3 scripts/uniqueness_check.py docs/downloads/<file>.tif --prior-dir <dir of prior tifs>
python3 scripts/verify_committed_submission.py
python3 scripts/build_site.py                  # regenerate the site from docs/data/*.json
```

## Sources for manual review

| Claim | Source |
| --- | --- |
| Metric, α = 0.2, β = 0.8, R = 300 m, worked example 0.60, submission format | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/> |
| Known-fault mask is pixel-exact; new-fault truth can occur within 300 m of a known trace; unmasked predictions far from new truth stay penalized (staff, post 4) | <https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/4> |
| Competition structure, two prize rounds, expert-labelled new faults | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#competition-structure> |
| Official rules (weekly allowance, one final submission, AI disclosure) | <https://docs.nlr.gov/docs/fy26osti/96647.pdf> |
| Cross-gradient function | Gallardo & Meju (2004), <https://doi.org/10.1029/2003JB002716> |
| Cross-gradients for gravity and magnetics | Fregoso & Gallardo (2009), <https://doi.org/10.1190/1.3119263> |
| Angle/cross product of horizontal gradients in the data domain, no inversion | ASEG 2012, <https://doi.org/10.1071/ASEG2012ab273> |
| Regional isostatic gravity 1 km grid (the gravity-resolution caution) | <https://pubs.usgs.gov/ds/2006/234/nv_iso.htm> |
| GeoDAWN airborne magnetic and radiometric surveys | <https://doi.org/10.5066/P93LGLVQ> |
| INGENIOUS geothermal compilation | <https://gdr.openei.org/submissions/1391> |
| Reference solution (organizer) | <https://github.com/drivendataorg/gems-prize-reference-solution> |
| Official competition data tab (login-gated) | <https://www.drivendata.org/competitions/306/competition-doe-gems/data/> |

## Limitations that are *not* solved

1. **No leaderboard feedback loop.** No DrivenData session exists in this environment and the terms
   forbid automated monitoring. Every local number is a proxy; none is a leaderboard score, and the
   public board's top row changed during this work (0.3195 → 0.3774), which no local instrument can
   track.
2. **The holdout population is not the test population.** Hide-and-recover removes faults from the
   *catalogue* — faults a mapper already found. The competition truth is faults a mapper *missed*,
   plausibly harder and less topographically obvious. The instrument measures placement skill on the
   wrong-population proxy; the sign of the transfer is not guaranteed.
3. **The catalogue-buffer bet.** New-fault truth can sit within 300 m of a known trace (staff, post
   4). The shipped 300 m buffer is measured to help on both holdout readings and to stabilise the
   per-block margin, but a test population concentrated at range fronts hard against mapped traces
   would punish it. This is the largest single interpretive risk in the file.
4. **No 1 m DEM.** The high-resolution elevation is distributed as a link list of ~100 MB tiles
   across terabytes; `prd-tnm.s3.amazonaws.com` returns HTTP 000 from this environment. That is
   probably the largest untapped lever, and it is named and checkable rather than assumed.
5. **The prior-submission corpus is not tracked.** The 131-raster relabel check ran against rasters
   collected into `/tmp/prior/files` during this session; the corpus is not in the repository, so the
   uniqueness verdict is reproducible only by re-fetching it. Within this repository the earlier
   session's 57-file audit reports max signed Spearman −0.266 and max overlap 0.290, giving the same
   qualitative verdict.
6. **No GPU, 4 GB RAM, 2 cores.** Training a competitive U-Net is out of reach here; everything is
   CPU- and memory-lean by construction.

## Next actions

1. Use one of the three weekly slots on the cleared file and record the leaderboard score; every
   local decision so far is unanchored to live feedback.
2. Restore the 1 m / 10 m elevation path (official, free, named) and re-run the topographic carrier
   at that resolution — the only measured lever large enough to move the private score materially.
3. Test the buffer question on live data: submit a 0 px-buffer twin of the same field and compare.
   That is a clean A/B on the one rule whose sign depends on unpublished scoring geometry.
4. Re-run the relabel check with a corpus that is tracked in the repository, so the verdict is
   reproducible by a reviewer without this sandbox.
5. Keep `docs/archive/` and `research/` in sync with the site when pages are regenerated.
