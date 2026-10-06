# GEMSDOE49 — cross-gradient structural coupling for the DOE GEMS Prize

**Competition:** [DOE GEMS Prize, DrivenData #306](https://www.drivendata.org/competitions/306/competition-doe-gems/)
· **Problem description:** [page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
· **Rules PDF:** [docs.nlr.gov/docs/fy26osti/96647.pdf](https://docs.nlr.gov/docs/fy26osti/96647.pdf)
· **Live site:** <https://buffedlizard55-lab.github.io/GEMSDOE49/>
· **Weekly budget:** 3 scored submissions; **one** file is scored in *both* prize rounds.

> ### Our Core Values — kept central to every decision in this repository
>
> **Maximize P(Win).** *"Maximize the Probability of Winning"* is our decision-making framework.
> In every decision we weigh tradeoffs, assess risk, and choose the path that maximizes the
> probability that we win this prize. We set aside our emotions and make tough decisions in order
> to maximize P(Win). It frees us from constraints and clarifies that we must put the outcome
> first.
>
> **Own the Outcome.** We own results end to end — not just our individual slice of the work.
> When problems arise and we have the means to act, we act without waiting for permission or
> assignment. We treat failure and success as signals and use them to improve. We stay
> accountable to the final outcome.
>
> Applied here: *Maximize P(Win)* is why the brief's headline construction (the Gallardo–Meju
> cross-gradient of gravity × magnetics) is **reported as a measured rejection** rather than
> quietly shipped — a submission slot is a scarce resource, and a construction that scores at the
> random floor must not consume one. *Own the Outcome* is why every number below is either
> **[OFFICIAL]** (read from an organizer page, with a link), **[MEASURED]** (computed here from
> hash-verified bytes, with the command), or **[BLOCKED]** (with the reason).

---

## ⬇ ONE-CLICK SUBMISSION FILE

**Download → [`docs/downloads/gems49-gate_ortho_w0.25-40k-20261006T213721Z-nan.tif`](docs/downloads/gems49-gate_ortho_w0.25-40k-20261006T213721Z-nan.tif)**
(377,092 bytes, sha256 `c58fe0551e120678fff190239146985d0b2f7255255535fb5fb637a05cb40af3`).

A twin, [`…-zeros.tif`](docs/downloads/gems49-gate_ortho_w0.25-40k-20261006T213721Z-zeros.tif),
contains the **same predictions** with `0.0` instead of `NaN` outside the footprint; use it if the
form rejects the `-nan` one. A single-file `.zip` of the primary variant is next to them. Sizes,
hashes and the submission-form note are in
[`docs/data/submission.json`](docs/data/submission.json) and printed on the
[site](https://buffedlizard55-lab.github.io/GEMSDOE49/).

**Submission form fields** (the two boxes you have to fill):

| Field | Value |
| --- | --- |
| *File to submit* | the `.tif` (or the `.zip` that contains it) |
| *Note* | see `note_for_submission_form` in [`docs/data/submission.json`](docs/data/submission.json) |

**If the form says `Predicted values must be in range [0, 1]`** — you are hitting one of exactly
two conditions, both of which are reproduced as unit tests in `tests/test_gate.py` and both of
which the local gate rejects before you upload:

* **mode A** — a `NaN` **inside** the footprint (outside the footprint is fine: the organizer's own
  `example_submission.tif` is `NaN` there). A naive `min >= 0 and max <= 1` check fails on `NaN`.
* **mode B** — a finite value outside `[0, 1]`.

`scripts/validate_submission.py` checks both, plus CRS, resolution, shape, geotransform, dtype,
band count, and byte-level agreement of the `NaN` mask with the official template.

---

## THE STANDING BRIEF (read this first, every session)

> Cross-gradient structural coupling of gravity, magnetics and topography. Treat agreement across
> independent signals as an operator, not a product of scores. Gallardo and Meju (*Journal of
> Geophysical Research*, 2004, [doi:10.1029/2003JB002716](https://doi.org/10.1029/2003JB002716))
> defined the cross-gradient function, the cross product of two property fields' gradient vectors.
> It vanishes where both fields change in the same direction, so it measures shared structure
> without assuming a petrophysical relationship. Fregoso and Gallardo (*Geophysics* 74(4),
> L31–L42, 2009, [doi:10.1190/1.3119263](https://doi.org/10.1190/1.3119263)) applied it to
> gravity and magnetics. An ASEG 2012 extended abstract
> ([doi:10.1071/ASEG2012ab273](https://doi.org/10.1071/ASEG2012ab273)) used the angle and cross
> product of the horizontal gradients of gravity and reduced-to-pole magnetics directly in the
> data domain, with no inversion. Compute horizontal-gradient vectors of the isostatic gravity
> anomaly, the reduced-to-pole magnetic anomaly and a smoothed elevation surface at a common
> resolution. Score each pixel by magnitude-weighted alignment (both gradients large and
> parallel). Keep magnitude-weighted orthogonality as a separate confounder flag, not as evidence.
> A fault juxtaposing rocks of different density and susceptibility puts both edges on one trace,
> while a buried volcanic edge usually moves the magnetic edge alone, so the coupled score should
> suppress single-layer edges. Check the native gravity station spacing before trusting 100 m
> gradients, and expect density steps at basin margins. The candidate must beat the repo's
> single-layer gradient baseline on hide-and-recover holdout segments before it is allowed near a
> submission. Normalize to [0,1], apply the repo's metric-aware placement, and write the required
> float32 GeoTIFF (EPSG:32611, 100 m, NaN only outside the footprint and none inside). Before
> download, hash the raster and compare rank correlation and top-k overlap against every prior
> submission, and reject it as a relabel if it isn't clearly different.

**Plus the standing process requirements:** before implementing, generate 3–5 candidate geological
hypotheses not yet tried, each naming the specific layer(s), the physical signature being targeted,
why it should catch a fault missing from the USGS/INGENIOUS catalogue rather than one already in
it, and how it differs from anything already implemented in the sibling repositories; rank them by
expected DTI improvement and implementation cost; validate the top candidate on the spatially
blocked holdout set before touching a weekly submission slot; if a candidate cannot be validated
without new external data, name the specific free official source and check it is obtainable before
proposing the idea as viable. Work line by line, verify from official trusted sources, provide
links for manual review, do the work autonomously, flag irregularities, and hallucinate nothing.

*The ranked hypothesis set for this session is in
[`docs/hypotheses.html`](docs/hypotheses.html) / [`docs/research/hypotheses.md`](docs/research/hypotheses.md);
the knowledge base is in [`docs/research/knowledge-base.md`](docs/research/knowledge-base.md).*

---

## What this repository measured this session

Everything below is **[MEASURED]** in this checkout from the hash-verified official bytes unless
labelled otherwise. Full receipts in `docs/data/*.json`.

### 1. The brief's headline construction does not survive contact with this grid

The pure Gallardo–Meju / ASEG-2012 data-domain cross-gradient of the isostatic gravity anomaly
against the reduced-to-pole magnetic anomaly — implemented exactly as described, magnitude-weighted
orthogonality, no inversion — scores **at the no-information floor** on the hide-and-recover
holdout. See `docs/data/field_screen.json` and `docs/data/final_experiment.json`.

### 2. The reason is measurable, and the brief told us to look for it

The brief says *"Check the native gravity station spacing before trusting 100 m gradients."*
Measured here from the official bytes (`scripts/measure_native_resolution.py`), the integral
scale of each field is:

| field | lag-1 autocorrelation (ρ₁) | integral scale (m) |
| --- | --- | --- |
| `iso_grav_anom` (13) | 0.942 / 0.924 | 5,800 – 8,100 |
| `rtp` (2) | 0.934 / 0.919 | 2,700 – 2,800 |
| `cond_surf` (17) | 0.942 / 0.923 | 10,400 – 10,500 |
| `det_elev_slope` (19) | **0.814 / 0.810** | **2,700 – 4,000** |

A 100 m horizontal gradient of a field whose integral scale is kilometres is interpolation
structure, not a contact edge. `det_elev_slope` is the only supplied field with real power at the
grid scale, which is exactly why it is the strongest single-layer gradient baseline.

### 3. The operator still earns its place — as a *gate*, not as a score

The brief's own reading is the one that works: the coupling's job is to **suppress single-layer
edges**, not to generate the signal. Used as a multiplicative orthogonality gate on a multi-scale,
coherence-weighted topographic edge it is a **measured no-op**: the gated and ungated variants share
36,602 of 40,000 dots (support Jaccard 0.8434, `scripts/measure_gate_effect.py`) and differ by
+0.00003 in holdout score. That is reported as what it is — a documented confounder filter that does
not damage the carrier — and not as a gain. Used as a standalone score it is worthless.

### 4. The placement rule is where the measurable gain was

The highest-scoring artifact in the sibling cohort prunes every dot within 200 m of the catalogue and
claims a live gain from it. Our instrument says the same thing, and pins the optimum at the scoring
kernel's own radius (`scripts/prune_experiment.py`, `scripts/margin_check.py`):

| catalogue buffer | candidate, collared / unrestricted | baseline, collared / unrestricted | margin |
| --- | --- | --- | --- |
| none | 0.1074 / 0.1134 | 0.0991 / 0.1077 | +0.0083 / +0.0056 |
| 0 px (exact pixels, as first shipped) | 0.1064 / 0.1127 | 0.0980 / 0.1068 | +0.0083 / +0.0059 |
| 200 m | 0.1100 / 0.1166 | 0.1010 / 0.1105 | +0.0090 / +0.0062 |
| **300 m (shipped)** | **0.1131 / 0.1181** | **0.1029 / 0.1111** | **+0.0102 / +0.0070** |
| 400 m | 0.1156 / 0.1164 | 0.1047 / 0.1086 | +0.0108 / +0.0078 |

The rule was given to the baseline too before the margin was believed — a placement rule handed to
one arm only is not a comparison. Beyond 300 m the *unrestricted* score of both arms falls: the
buffer starts deleting the dots that would have earned credit where a new fault is mapped right
against a known one.

### 5. What was actually shipped

`gate_ortho_w0.25` — multi-scale topographic edge × orientation coherence × orthogonality gate, 40,000
dots at 300 m separation, 300 m catalogue buffer — under a leave-one-fold-out geometry-selection
protocol on a spatially blocked, leakage-controlled holdout (leakage probe TP_w = 0, null probe
DTI = 0, on-truth oracle 1.0000, random and uniform-lattice controls at matched geometry).

---

## Repository map

| Path | What it is |
| --- | --- |
| `src/gems49/spec.py` | grid constants and sha256 pins, all measured from the official bytes |
| `src/gems49/metric.py` | the official distance-weighted Tversky index, verified against the page's worked example and a literal brute force |
| `src/gems49/crossgrad.py` | horizontal gradients, magnitude-weighted alignment, orthogonality flag, balance operator |
| `src/gems49/gating.py` | the cross-gradient used as a confounder gate (the brief's "operator, not a product of scores") |
| `src/gems49/holdout.py` | hide-and-recover folds, leakage probe, cached scorer |
| `src/gems49/emission.py` | the metric's marginal rule and metric-aware placement |
| `src/gems49/io.py` | fail-closed raster reading and the submission writer |
| `scripts/download_competition_data.sh` | fetch + sha256-verify every official raster |
| `scripts/prepare_data.py` | audit and pin the placement; writes `data/prepared_manifest.json` |
| `scripts/run_validation.py`, `scripts/screen_fields.py`, `scripts/final_experiment.py` | the screens |
| `scripts/measure_native_resolution.py` | the "check the native station spacing" measurement |
| `scripts/metric_facts.py` | derives the metric algebra the site quotes (worked example, marginal bar, binary dominance) |
| `scripts/prune_experiment.py` | the catalogue-buffer sweep |
| `scripts/margin_check.py` | re-measures the margin with the placement rule applied to *both* arms |
| `scripts/measure_gate_effect.py` | pixel-level comparison of the gated and ungated supports |
| `scripts/build_submission.py` | normalise → metric-aware placement → write float32 |
| `scripts/validate_submission.py` | the hard format gate (13 checks, both rejection modes reproduced as tests) |
| `scripts/verify_committed_submission.py` | re-hashes the committed artifact and re-derives the contract; runs unattended in CI |
| `scripts/uniqueness_check.py` | rank correlation and top-k overlap against every prior submission |
| `scripts/build_site.py` | generates the whole GitHub Pages site from the evidence JSON |
| `docs/` | the published site (GitHub Pages) |

## Quickstart

```bash
bash scripts/download_competition_data.sh     # fetch + verify every official raster
python scripts/prepare_data.py                # audit grid, CRS, bands, labels
python -m pytest tests -q                     # metric, gate, holdout, cross-gradient tests
python scripts/final_experiment.py            # the LOFO screen  (~20 min, 1 core)
python scripts/build_submission.py --family gate_ortho_w0.25 --n-target 40000 --min-sep 3 \
    --exclude-catalogue-px 3 --purge        # writes docs/downloads/*.tif + docs/data/submission.json
python scripts/validate_submission.py docs/downloads/<file>.tif --json docs/data/submission_gate.json
python scripts/uniqueness_check.py docs/downloads/<file>.tif --prior-dir <dir of prior tifs>
python scripts/verify_committed_submission.py # hash + contract check of the committed bytes
python scripts/build_site.py                  # regenerate the site from docs/data/*.json
python scripts/prune_experiment.py ; python scripts/margin_check.py   # the placement-rule evidence
```

## Sources for manual review

| Claim | Source |
| --- | --- |
| Metric, α = 0.2, β = 0.8, R = 300 m, worked example 0.60, submission format | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/> |
| Competition structure, two prize rounds, expert-labelled new faults | <https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/#competition-structure> |
| Known USGS/INGENIOUS fault pixels are masked / excluded from evaluation in both rounds (DrivenData staff, 2026-09-16) | <https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516> |
| Official competition data tab (login-gated) | <https://www.drivendata.org/competitions/306/competition-doe-gems/data/> |
| Official rules (weekly allowance, one final submission, blind choice, AI disclosure) | <https://docs.nlr.gov/docs/fy26osti/96647.pdf> |
| Cross-gradient function | Gallardo & Meju (2004), <https://doi.org/10.1029/2003JB002716> |
| Cross-gradients for gravity and magnetics | Fregoso & Gallardo (2009), <https://doi.org/10.1190/1.3119263> |
| Angle/cross product of horizontal gradients in the data domain, no inversion | ASEG 2012 extended abstract, <https://doi.org/10.1071/ASEG2012ab273> |
| GeoDAWN airborne magnetic and radiometric surveys | <https://doi.org/10.5066/P93LGLVQ> |
| INGENIOUS geothermal compilation | <https://gdr.openei.org/submissions/1391> |
| Reference solution (organizer) | <https://github.com/drivendataorg/gems-prize-reference-solution> |

## Limitations that are *not* solved

1. **No leaderboard feedback loop.** There is no DrivenData session in this environment
   (the data tab redirects to `/accounts/login/`, re-verified), and the competition terms forbid
   automated monitoring. Every local number is a proxy; none is a leaderboard score.
2. **The holdout population is not the test population.** Hide-and-recover removes faults from the
   *catalogue*, i.e. faults a mapper already found. The competition truth is faults a mapper
   *missed*, which are plausibly harder and less topographically obvious. The instrument therefore
   measures detection skill on the wrong-population proxy, and it is the best public proxy
   available — but the sign of the transfer is not guaranteed.
3. **No 1 m DEM.** The competition distributes the high-resolution DEM as a link list; the tiles
   are ~100 MB each across a region spanning terabytes, and the USGS 3DEP bucket is not reachable
   from this environment (`prd-tnm.s3.amazonaws.com` → HTTP 000). This is probably the single
   largest untapped lever.
4. **The exact geometry of the organizer's known-fault mask is not published**, only the statement
   that known faults are excluded, and two readings of it are in circulation: catalogue pixels
   dropped from *scoring* entirely, or dropped from *truth* only (so a dot there is a false
   positive). The two readings differ by several hundredths of a point for any submission that emits
   on the catalogue, and they point in opposite directions for buffers wider than 300 m — the
   collared reading keeps rewarding a wider buffer while the unrestricted reading punishes it. This
   repository ships the 300 m buffer because it is the optimum under the *unrestricted* reading,
   which is the faithful analogue of a new fault mapped against a known one, and it is still an
   improvement under the collared reading. **Flagged as the largest remaining interpretive risk.**
5. **No GPU, 4 GB RAM.** Training a full U-Net is out of reach here; everything is CPU and
   memory-lean by construction.
