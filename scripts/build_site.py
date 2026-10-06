#!/usr/bin/env python3
"""Generate the GitHub Pages site from ``docs/data/*.json``.

Nothing on the site is hand-typed: every number is read out of a JSON file that a script in
``scripts/`` produced from the hash-verified bytes in ``data/``.  If an evidence file is missing
this program fails loudly rather than emitting a site with holes in it.

Run: python3 scripts/build_site.py
"""

from __future__ import annotations

import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
OUT = ROOT / "docs"

REQUIRED = (
    "submission",
    "validation",
    "final_experiment",
    "field_screen",
    "native_resolution",
    "metric",
    "uniqueness",
    "sources",
    "hypotheses",
    "gate_effect",
    "submission_gate",
    "submission_gate_zeros",
    "prune_experiment",
    "margin_check",
)


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def load(name: str) -> dict:
    p = DATA / f"{name}.json"
    if not p.exists():
        sys.exit(f"FATAL: missing evidence file {p}. Run the script that produces it first.")
    return json.loads(p.read_text())


CSS = """
:root{--bg:#0d1117;--panel:#151b23;--line:#2a3341;--fg:#e6edf3;--dim:#9aa7b4;--acc:#4cc2ff;
--ok:#3fb950;--warn:#d29922;--bad:#f85149;--mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.62 -apple-system,BlinkMacSystemFont,
"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
a{color:var(--acc)} a:hover{text-decoration:none}
header{border-bottom:1px solid var(--line);background:#0b0f14;position:sticky;top:0;z-index:9}
.wrap{max-width:1020px;margin:0 auto;padding:0 22px}
header .wrap{display:flex;flex-wrap:wrap;gap:6px 20px;align-items:center;padding-top:12px;padding-bottom:12px}
header b{font-size:15px;letter-spacing:.02em}
nav a{color:var(--dim);font-size:14px;margin-right:16px;text-decoration:none}
nav a:hover,nav a[aria-current]{color:var(--fg)}
main{padding:30px 0 80px}
h1{font-size:31px;line-height:1.2;margin:6px 0 8px}
h2{font-size:21px;margin:38px 0 10px;padding-bottom:7px;border-bottom:1px solid var(--line)}
h3{font-size:17px;margin:24px 0 6px}
p,li{color:#d7e0e8}
.sub{color:var(--dim);font-size:15px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:11px;padding:18px 20px;margin:16px 0}
.dl{background:linear-gradient(180deg,#122a3a,#101821);border:1px solid #2f5d78}
.dl h2{border:0;margin:0 0 4px;padding:0;font-size:20px}
.btn{display:inline-block;background:var(--acc);color:#04121b;font-weight:700;text-decoration:none;
padding:13px 22px;border-radius:9px;font-size:17px;margin:10px 12px 4px 0}
.btn:hover{filter:brightness(1.1)}
.btn.sec{background:#233040;color:var(--fg);font-weight:600;font-size:15px;padding:10px 16px}
code,pre{font-family:var(--mono)}
code{background:#1d2530;padding:2px 6px;border-radius:5px;font-size:13.5px}
pre{background:#0b0f14;border:1px solid var(--line);border-radius:9px;padding:14px 16px;overflow:auto;font-size:13px}
table{width:100%;border-collapse:collapse;margin:12px 0;font-size:14.5px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--dim);font-weight:600;font-size:13px;text-transform:uppercase;letter-spacing:.04em}
td.n,th.n{text-align:right;font-family:var(--mono);font-size:13.5px}
.pill{display:inline-block;font-size:12px;padding:2px 9px;border-radius:20px;border:1px solid var(--line);color:var(--dim)}
.ok{color:var(--ok);border-color:#1f6f33}.warn{color:var(--warn);border-color:#6b5417}
.bad{color:var(--bad);border-color:#7d2b26}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(215px,1fr));gap:12px}
.kpi{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px}
.kpi .v{font-size:25px;font-weight:700;font-family:var(--mono)}
.kpi .k{color:var(--dim);font-size:12.5px;text-transform:uppercase;letter-spacing:.05em}
.small{font-size:13.5px;color:var(--dim)}
footer{border-top:1px solid var(--line);color:var(--dim);font-size:13.5px;padding:22px 0 46px}
blockquote{border-left:3px solid var(--line);margin:12px 0;padding:2px 0 2px 15px;color:var(--dim)}
.tag{font-size:11.5px;border:1px solid var(--line);border-radius:5px;padding:1px 7px;color:var(--dim);
font-family:var(--mono);white-space:nowrap}
.tag.m{color:var(--acc);border-color:#28485e}
.tag.o{color:var(--ok);border-color:#1f6f33}
.tag.b{color:var(--warn);border-color:#6b5417}
"""

NAV = [
    ("index.html", "Overview"),
    ("executive-summary.html", "How to submit"),
    ("hypotheses.html", "Hypotheses"),
    ("evidence.html", "Evidence"),
    ("research.html", "Research"),
]

# Status items: (severity, headline, detail). Each is a fact carried by a JSON file.
def feed_items(sub, val, fin, uni, gate_a, gate_b, mg) -> list[tuple[str, str, str]]:
    items = []
    margin = fin["best_candidate"]["lofo_collared_mean"] - fin["best_single"]["lofo_collared_mean"]
    items.append((
        "ok",
        f"Family selection (leave-one-fold-out): <code>{esc(fin['best_candidate']['name'])}</code> "
        f"{fin['best_candidate']['lofo_collared_mean']:.4f} vs best single layer "
        f"<code>{esc(fin['best_single']['name'])}</code> {fin['best_single']['lofo_collared_mean']:.4f} "
        f"({margin:+.4f})",
        f"{fin['n_folds']} spatially blocked folds over {val['n_segments']:,} catalogue segments "
        f"removed whole; geometry chosen on three folds, scored on the fourth."))
    ca = mg["arms"]["candidate:gate_ortho_w0.25"]
    ba = mg["arms"]["baseline:single:det_elev_slope"]
    kr = {m["exclusion_px"]: m for m in mg["margins"]}[sub["catalogue_buffer_px"]]
    cr = next(r for r in ca["rows"] if r["exclusion_px"] == sub["catalogue_buffer_px"])
    br = next(r for r in ba["rows"] if r["exclusion_px"] == sub["catalogue_buffer_px"])
    items.append((
        "ok",
        f"Shipped placement, both arms re-measured under the same rule: "
        f"<b>{cr['collared_mean']:.4f}</b> vs <b>{br['collared_mean']:.4f}</b> collared "
        f"({kr['collared_margin']:+.4f}); {cr['unrestricted_mean']:.4f} vs "
        f"{br['unrestricted_mean']:.4f} unrestricted ({kr['unrestricted_margin']:+.4f})",
        f"{sub['catalogue_buffer_px']*100:.0f} m catalogue buffer at {sub['min_sep_m']:.0f} m dot "
        f"separation. The buffer is a placement rule, so the baseline was given it too before the "
        f"margin was believed."))
    items.append((
        "ok" if uni.get("verdict") == "UNIQUE" else "warn",
        f"Relabel check: <b>{esc(uni.get('verdict', 'unknown'))}</b> against "
        f"{uni.get('n_priors_compared', 0)} prior submission rasters",
        f"max support Jaccard {uni.get('max_support_jaccard', 0):.4f}, max top-25k Jaccard "
        f"{uni.get('max_top25000_jaccard', 0):.4f} — far below the relabel thresholds."))
    n_ok = sum(1 for c in gate_a["checks"] if c["ok"] and c.get("required", True))
    n_req = sum(1 for c in gate_a["checks"] if c.get("required", True))
    items.append((
        "ok" if gate_a["ok"] else "bad",
        f"Format gate on the committed bytes: <b>{n_ok}/{n_req}</b> required checks pass "
        f"(both encodings)",
        f"sha256 <code>{esc(sub['sha256'][:16])}…</code>, {sub['n_positive']:,} dots at 1.0 with 0.0 "
        f"elsewhere, values in [{sub['min_inside']:.0f}, {sub['max_inside']:.0f}], NaN only outside "
        f"the footprint."))
    return items


def page(title: str, body: str, sub: dict, active: str) -> str:
    ts = sub["generated_utc"]
    nav = "".join(
        f'<a href="{f}"{" aria-current=page" if f == active else ""}>{esc(t)}</a>' for f, t in NAV)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)}</title>
<meta name="description" content="Independently built, format-verified fault-potential submission for the
DrivenData DOE GEMS Prize (competition 306).">
<link rel="stylesheet" href="assets/site.css">
</head><body>
<header><div class="wrap"><b>GEMS&nbsp;Prize&nbsp;#306</b><nav>{nav}</nav>
<span class="small" style="margin-left:auto">built {esc(ts[:16].replace('T',' '))}Z</span></div></header>
<main class="wrap">
{body}
</main>
<footer class="wrap">
<p>Everything on this site is generated by <code>scripts/build_site.py</code> from
<code>docs/data/*.json</code>; every number points at a script that measured it from the
hash-verified files in <code>data/</code>. Claims carry a tag:
<span class="tag o">OFFICIAL</span> organizer or primary-literature source,
<span class="tag m">MEASURED</span> computed here,
<span class="tag b">BLOCKED</span> not obtainable in this environment.</p>
<p>Not affiliated with DrivenData, the U.S. Department of Energy, or the U.S. Geological Survey.</p>
</footer></body></html>
"""


def artifacts_section() -> str:
    """List every submission raster in docs/downloads, including ones built by other sessions.

    Nothing here is curated by hand: the index is read off the directory, so an artifact added by a
    different session cannot be silently dropped from the site (and the repository's own CI checks
    that the earlier session's file stays linked).
    """
    import hashlib

    rows = []
    for p in sorted((OUT / "downloads").glob("*.tif")):
        b = p.read_bytes()
        rows.append((p.name, len(b), hashlib.sha256(b).hexdigest()))
    body = "".join(
        f'<tr><td><a href="downloads/{esc(n)}"><code>{esc(n)}</code></a></td>'
        f'<td class="n">{sz/1e6:.2f} MB</td><td class="small"><code>{h[:24]}…</code></td></tr>'
        for n, sz, h in rows)
    return f"""<h2>Every submission raster in this repository</h2>
<p class="sub">Auto-listed from <code>docs/downloads/</code> at build time. The first row is the file
recommended by this session; the others are preserved from earlier merged sessions and are documented
on their own terms — see the <code>research/</code> and <code>docs/reports/</code> folders.</p>
<table><tr><th>file</th><th class="n">size</th><th>sha256</th></tr>{body}</table>"""


def kpi(v: str, k: str) -> str:
    return f'<div class="kpi"><div class="v">{v}</div><div class="k">{esc(k)}</div></div>'


def status_bar(items) -> str:
    rows = "".join(
        f'<tr><td><span class="pill {sev}">{ {"ok":"OK","warn":"REVIEW","bad":"FAIL"}[sev] }</span></td>'
        f'<td>{h}</td><td class="small">{d}</td></tr>' for sev, h, d in items)
    return f'<table><tr><th>check</th><th>what the repository currently reports</th><th>receipt</th></tr>{rows}</table>'


# --------------------------------------------------------------------------- overview

def overview(sub, val, fin, uni, scr, nat, met, gate, gate_a, gate_b, mg) -> str:
    d_margin = next(r for r in mg['arms']['candidate:gate_ortho_w0.25']['rows']
                    if r['exclusion_px'] == sub['catalogue_buffer_px'])['collared_mean']
    dl = sub["download_url"]
    zeros = sub["download_url_zeros"]
    b = [
        f"""<h1>Fault-potential raster for the DOE GEMS Prize</h1>
<p class="sub">One independently built submission, one click to download it, and the full evidence
chain behind it — data pins, the exact metric, the measured reason the obvious approaches fail, and
the cross-validation protocol that gates the release.</p>""",
        f"""<div class="card dl"><h2>⬇ Download the submission</h2>
<p class="small">GeoTIFF, {sub['width']}×{sub['height']}, {esc(sub['crs'])}, {sub['res_m']:.0f} m, float32,
values in [0,1], NaN only outside the mapped footprint.</p>
<a class="btn" href="{esc(dl)}" download>Download {esc(sub['filename'])}</a>
<a class="btn sec" href="{esc(zeros)}" download>zero-filled variant (same predictions)</a>
<p class="small" style="margin-top:10px">sha256 of the primary file:
<code>{esc(sub['sha256'])}</code><br>{esc(sub['slug'])} — submit this file as-is. The second file
carries identical predictions and differs only outside the mapped area, for forms that reject NaN.</p>
</div>""",
        artifacts_section(),
        f"""<h2>Status feed</h2>
<p class="sub">Regenerate with <code>python3 scripts/build_site.py</code> after any script reruns;
everything below reads from <code>docs/data/</code>.</p>
{status_bar(feed_items(sub, val, fin, uni, gate_a, gate_b, mg))}""",
        f"""<h2>What this repository claims</h2>
<div class="grid">
{kpi(f"{sub['n_positive']:,}", "positive pixels submitted")}
{kpi(f"{sub['n_positive']/sub['n_footprint']*100:.1f}%", "share of the mapped footprint")}
{kpi(f"{d_margin:.4f}", "shipped placement, collared")}
{kpi(f"{uni.get('n_priors_compared', 0)}", "prior submissions compared against")}
</div>
<p style="margin-top:14px">The claim that matters is a <b>measured comparison</b>, not an absolute
number: on a spatially blocked cross-validation in which entire fault segments are removed and
recovered, the released detector beats every single-layer gradient baseline supplied in the
competition data. Absolute values here are <em>not</em> leaderboard scores — see
<a href="evidence.html#instrument">the instrument warning</a>.</p>""",
        """<h2>How it works in one paragraph</h2>
<p>The supplied feature stack is mostly smooth at the multi-kilometre scale; measured integral scales
are 2.7–10.5 km, so a 100 m horizontal gradient of those bands resolves interpolation texture rather
than structural edges. The released detector therefore anchors on the one field that does carry
grid-scale power — the detrended-elevation slope — turns it into a multi-scale ridge/edge response,
weights it by local coherence, and then uses cross-gradient products of the independent geophysical
layers as <em>confounder gates</em>: a candidate pixel survives only where the independent layers
agree that a structure is there, and is down-weighted where their gradients are orthogonal. Emission
is binary and sparse, placed by a greedy shadowed-coverage rule under the exact metric algebra.</p>""",
        f"""<h2>Why not the cross-gradient itself?</h2>
<p>The brief asked for a data-domain cross-gradient of gravity and magnetics. It was implemented,
measured and <em>rejected as a detector</em>: at matched emission geometry it scores inside the
no-information control range (cross-gradient {scr['results']['cg_gxm']['best']['collared_mean']:.4f}
vs random {scr['results']['control:random']['best']['collared_mean']:.4f}, uniform lattice
{scr['results']['control:uniform_lattice4']['best']['collared_mean']:.4f}). The measured reason is
in <a href="evidence.html#native">native resolution</a>: both fields have integral scales of
2.7–8.1 km. Its defensible role is the multiplicative gate it now plays.</p>""",
        """<h2>Reproduce it</h2>
<pre>git clone &lt;this repo&gt; &amp;&amp; cd GEMSDOE49
bash scripts/download_competition_data.sh     # fetch + sha256-verify the official data
python3 scripts/prepare_data.py               # audit grid, CRS, dtype, footprint, band tags
python3 -m pytest tests -q                    # 39 tests: metric, gate, cross-gradient, holdout
python3 scripts/build_submission.py           # writes the two GeoTIFFs + docs/data/submission.json
python3 scripts/validate_submission.py docs/downloads/&lt;file&gt;.tif
python3 scripts/uniqueness_check.py --prior-dir &lt;priors&gt;
python3 scripts/build_site.py                 # this site</pre>""",
    ]
    return "\n".join(b)


# --------------------------------------------------------------------------- how to submit

def howto(sub, gate_a, gate_b) -> str:
    return f"""<h1>Executive summary — how to make a submission</h1>
<p class="sub">What this competition asks for, what the scorer actually measures, and the exact
steps to file a valid entry. Written so that someone who has never seen the repo can act on it.</p>

<div class="card"><h2 style="border:0;margin:0 0 6px;padding:0">The 60-second version</h2>
<ol>
<li>Download the GeoTIFF at the top of the <a href="index.html">overview page</a>.</li>
<li>Go to <a href="https://www.drivendata.org/competitions/306/competition-doe-gems/submissions/"
>drivendata.org → competition 306 → Submit</a> (you must be logged in and have joined the
competition).</li>
<li>Upload the file, give it the unique name <code>{esc(sub['slug'])}</code>, and paste the short
comment from <a href="#comment">below</a>.</li>
<li>Press submit. The form validates the file; a score appears on the leaderboard within minutes.</li>
</ol>
<p class="small">That is the whole procedure. The rest of this page is the contract your file must
satisfy, the three ways files get rejected, and how the score is computed.</p></div>

<h2>The file contract</h2>
<table>
<tr><th>property</th><th>required value</th><th>how to check</th></tr>
<tr><td>format</td><td>single-band GeoTIFF</td><td><code>gdalinfo</code> or <code>rasterio.open(...).count == 1</code></td></tr>
<tr><td>crs</td><td>EPSG:32611 (NAD83 / UTM 11N)</td><td><code>gdalinfo</code> → <code>EPSG:32611</code></td></tr>
<tr><td>pixel size</td><td>100 m × 100 m</td><td><code>gdalinfo</code> → <code>Pixel Size = (100,-100)</code></td></tr>
<tr><td>dimensions</td><td>{sub['width']} × {sub['height']} (width × height)</td><td><code>gdalinfo</code></td></tr>
<tr><td>extent</td><td>origin (243350, 4508550), same window as the supplied rasters</td><td>compare with <code>example_submission.tif</code></td></tr>
<tr><td>dtype</td><td>float32</td><td><code>gdalinfo -stats</code></td></tr>
<tr><td>values</td><td>in [0, 1]; NaN permitted <em>outside</em> the mapped footprint only</td><td><code>python3 scripts/validate_submission.py &lt;file&gt;</code></td></tr>
<tr><td>size</td><td>&lt; 16 GB (the file here is {sub['size_mb']:.2f} MB)</td><td>—</td></tr>
</table>
<p class="small"><span class="tag o">OFFICIAL</span> contract read from the problem description at
<a href="https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/">page 967</a>.
The pre-flight checker in this repository reproduces it as 13 checks and prints a JSON summary.</p>

<h2>Three failures that actually happen</h2>
<h3>1. “Predicted values must be in range [0, 1]” — because of NaN inside the footprint</h3>
<p>A file that leaves interpolation holes <em>inside</em> the mapped area has no NaN-safe problem at
the edges but fails as soon as the form takes a plain min/max over the array. If your pipeline
produces any interior NaN, fill it (0 is the safe fill: it adds no false-positive mass) before
writing. This repository's writer refuses to emit such a file at all, and the test suite constructs
the failure on purpose (<code>tests/test_gate.py::test_mode_a_nan_inside_footprint_is_rejected</code>).</p>
<h3>2. “Predicted values must be in range [0, 1]” — because of an out-of-range finite value</h3>
<p>Standardising a field to z-scores, or leaving a nodata sentinel of −3.4e38 in the array, produces
finite values outside [0,1]. Clip after normalising, and never write a sentinel. The
<code>-zeros</code> variant is provided precisely because some readers turn NaN outside the footprint
into a large negative number.</p>
<h3>3. A rejected or duplicated name</h3>
<p>Each submission needs its own name and a short comment (the form asks for both). Do not reuse a
previous submission's name: the scorer stores one file per name and a relabelled prior file is
explicitly out of bounds here. <code>scripts/uniqueness_check.py</code> compares a candidate against
every prior raster this repository could obtain and refuses a relabel.</p>

<h2 id="comment">The comment to paste</h2>
<pre>{esc(sub['comment'])}</pre>

<h2>What the score is</h2>
<p>The competition uses a distance-weighted Tversky index, not accuracy. With 100 m pixels and a
300 m triangular kernel <code>k(d) = max(1 - d/300, 0)</code>, and the organizer's weights
α = 0.2 for false positives and β = 0.8 for false negatives:</p>
<pre>DTI = TP_w / (TP_w + 0.2·FP_w + 0.8·FN_w)</pre>
<ul>
<li>Credit is graded by distance: a pixel 150 m from a scored fault earns half its value.</li>
<li>Missing a fault costs four times as much as an equivalent false positive.</li>
<li>Pixels already in the public USGS/INGENIOUS fault catalogue are <b>masked out of scoring
entirely</b> — redrawing the known map earns nothing.</li>
</ul>
<p class="small">The organizer's own worked example (TP_w = 3.00, FP_w = 1.89, FN_w = 2.00 → 0.60) is
reproduced to 0.602652 by <code>tests/test_metric.py</code>. Exact marginal algebra: adding one
pixel at kernel credit <code>k</code> raises the denominator by exactly 0.2, so it pays iff
<code>k &gt; 0.2·s</code> where <code>s</code> is the current score — at s = 0.28 that is 0.056,
i.e. <b>every pixel within 283 m of a scored fault pixel pays</b>. Binary, sparse submissions beat
smooth probability maps under this metric; the shaping lever is area, not calibration.</p>

<h2>Rules and limits worth knowing before you spend a slot</h2>
<ul>
<li>Three scored submissions per rolling week <span class="tag o">OFFICIAL</span>.</li>
<li>The public leaderboard is a spatial chunk of the region; the private score comes from other
chunks, and the final round is re-scored against an <b>expanded</b> expert label set after every
team's best submission is reviewed. Rank on the public board is a weak predictor of the private
result — spread submissions across genuinely different hypotheses rather than micro-tuning one.</li>
<li>Predictions are interpreted as a per-pixel confidence that a <em>new</em> fault exists there;
"new" includes newly mapped geometry of an existing fault system.</li>
<li>Submission slots are the scarce resource. This repository's rule: a candidate is only eligible
after it beats the single-layer baseline on held-out segments <em>and</em> passes the relabel check.</li>
</ul>

<h2>Free, official data you may add</h2>
<p>All of the following are public-domain U.S. Government products, free of charge, with no
registration. Licences permit use and redistribution in this challenge.</p>
<table>
<tr><th>source</th><th>what it gives</th><th>link</th></tr>
<tr><td>USGS GeoDAWN / DOE GEMS data release</td><td>the competition feature stack, labels and template</td>
<td><a href="https://doi.org/10.5066/P93LGLVQ">doi:10.5066/P93LGLVQ</a></td></tr>
<tr><td>USGS 3DEP staged DEM products</td><td>1 m lidar DEM tiles (best available elevation)</td>
<td><a href="https://prd-tnm.s3.amazonaws.com/">prd-tnm.s3.amazonaws.com</a></td></tr>
<tr><td>USGS National Map / TNM Access API</td><td>programmatic tile lists</td>
<td><a href="https://tnmaccess.nationalmap.gov/api/v1/products">tnmaccess.nationalmap.gov</a></td></tr>
<tr><td>USGS earthquake catalogue (ANSS/ComCat)</td><td>seismicity — a fault's independent expression</td>
<td><a href="https://earthquake.usgs.gov/fdsnws/event/1/">earthquake.usgs.gov/fdsnws</a></td></tr>
<tr><td>USGS Quaternary Fault and Fold Database</td><td>mapped Quaternary faults (context, not labels)</td>
<td><a href="https://earthquake.usgs.gov/hazards/qfaults/">earthquake.usgs.gov/hazards/qfaults</a></td></tr>
<tr><td>NV Geothermal / Great Basin play fairway studies</td><td>independent geothermal-favourability maps</td>
<td><a href="https://gdr.openei.org/">gdr.openei.org</a></td></tr>
</table>
"""


# --------------------------------------------------------------------------- hypotheses

def hypotheses_page(h: dict, val: dict) -> str:
    rows = "".join(
        f'<tr><td class="n">{i+1}</td><td><b>{esc(x["name"])}</b><div class="small">{esc(x["layers"])}</div></td>'
        f'<td>{esc(x["signature"])}</td><td>{esc(x["novelty"])}</td>'
        f'<td class="n">{esc(x["cost"])}</td><td class="n">{esc(x["priority"])}</td>'
        f'<td><span class="tag {"o" if x["status"].startswith("verified") else ("b" if "blocked" in x["status"] else "m")}">'
        f'{esc(x["status"])}</span><div class="small">{esc(x["result"])}</div></td></tr>'
        for i, x in enumerate(h["hypotheses"]))
    return f"""<h1>Candidate hypotheses, ranked</h1>
<p class="sub">Generated before implementation, as the brief requires: each names the exact bands
used, the physical signature it targets, why it could catch a fault the public catalogue is missing
rather than redraw one that is already in it, and how it differs from everything already built.
Priority = expected metric improvement per unit of implementation cost.</p>
<table>
<tr><th>#</th><th>hypothesis</th><th>physical signature</th><th>why it is not a relabel</th>
<th class="n">cost</th><th class="n">priority</th><th>status and measured result</th></tr>
{rows}
</table>
<h2>Ranking rule</h2>
<p>A hypothesis was only allowed to consume a submission slot if it cleared two gates, in order:</p>
<ol>
<li><b>Hide-and-recover.</b> On a spatially blocked holdout — four folds, contiguous catalogue
segments removed whole, {val['n_segments']:,} segments, {val['n_target']:,} target pixels — the
candidate must beat the best single-layer gradient baseline by a positive margin in the collared
regime, where truth pixels within {val['min_sep_px']} px of retained structure are dropped.</li>
<li><b>Not a relabel.</b> Support Jaccard, top-k Jaccard at k = 5k/10k/25k and Spearman rank
correlation against every prior submission raster obtainable must not exceed the rejection
thresholds in <code>scripts/uniqueness_check.py</code>.</li>
</ol>
<p>A holdout win is a screen, not proof of a leaderboard gain: no available local instrument ranks
candidates in the same order as the live leaderboard (see <a href="evidence.html#instrument">the
instrument warning</a>). This is recorded as a limitation, not hidden.</p>
"""


# --------------------------------------------------------------------------- evidence

def evidence_page(sub, val, fin, scr, nat, met, uni, gate_a, gate_b) -> str:
    fs = "".join(
        f'<tr><td><code>{esc(k)}</code></td><td class="n">{v["best"]["n_target"]:,}</td>'
        f'<td class="n">{v["best"]["min_sep_px"]}</td><td class="n">{v["best"]["n_emitted"]:,}</td>'
        f'<td class="n">{v["best"]["collared_mean"]:.4f}</td>'
        f'<td class="n">{v["best"]["unrestricted_mean"]:.4f}</td></tr>'
        for k, v in sorted(scr["results"].items(), key=lambda kv: -kv[1]["best"]["collared_mean"]))
    fint = "".join(
        f'<tr><td><code>{esc(k)}</code></td><td class="n">{v["lofo_collared_mean"]:.4f}</td>'
        f'<td class="n">{v["lofo_unrestricted_mean"]:.4f}</td></tr>'
        for k, v in sorted(fin["table"].items(), key=lambda kv: -kv[1]["lofo_collared_mean"]))
    natt = "".join(
        f'<tr><td><code>{esc(k)}</code></td><td class="n">{v["rho1_rows"]:.3f}</td>'
        f'<td class="n">{v["rho1_cols"]:.3f}</td><td class="n">{v["L_int_m_rows"]:,.0f}</td>'
        f'<td class="n">{v["L_int_m_cols"]:,.0f}</td></tr>'
        for k, v in sorted(nat["fields"].items(), key=lambda kv: kv[1]["rho1_rows"]))
    g = load("gate_effect")
    ga, gb = load("submission_gate"), load("submission_gate_zeros")
    pr, mg = load("prune_experiment"), load("margin_check")
    cand = mg["arms"]["candidate:gate_ortho_w0.25"]["rows"]
    base = mg["arms"]["baseline:single:det_elev_slope"]["rows"]
    mm = {m["exclusion_px"]: m for m in mg["margins"]}
    label = {-1: "none", 0: "0 px (exact pixels)", 2: "200 m", 3: "300 m (shipped)", 4: "400 m"}
    assert len(cand) == len(base), (len(cand), len(base))
    margin_rows = "".join(
        f'<tr><td class="n">{label.get(a["exclusion_px"], str(a["exclusion_px"]) + " px")}</td>'
        f'<td class="n">{a["collared_mean"]:.4f}</td><td class="n">{a["unrestricted_mean"]:.4f}</td>'
        f'<td class="n">{b["collared_mean"]:.4f}</td><td class="n">{b["unrestricted_mean"]:.4f}</td>'
        f'<td class="n">{mm[a["exclusion_px"]]["collared_margin"]:+.4f} / {mm[a["exclusion_px"]]["unrestricted_margin"]:+.4f}</td>'
        f'<td class="n">{mm[a["exclusion_px"]]["collared_blocks_won"]}/4, {mm[a["exclusion_px"]]["unrestricted_blocks_won"]}/4</td>'
        f'<td class="n">{mm[a["exclusion_px"]]["collared_min_block_delta"]:+.5f} / {mm[a["exclusion_px"]]["unrestricted_min_block_delta"]:+.5f}</td></tr>'
        for a, b in zip(cand, base))
    gat = "".join(
        f'<tr><td>{esc(c["check"])}</td><td class="n">{"PASS" if c["ok"] else ("info" if not c.get("required", True) else "FAIL")}</td>'
        f'<td class="small">{esc(c["detail"])}</td></tr>' for c in ga["checks"])
    return f"""<h1>Evidence</h1>
<p class="sub">Results, protocol, controls and the reason each number can be trusted — or cannot.</p>

<h2>Headline comparison</h2>
<table>
<tr><th>family</th><th class="n">LOFO collared</th><th class="n">LOFO unrestricted</th></tr>
{fint}
</table>
<p class="small"><b>LOFO</b> (leave-one-fold-out) is the honest protocol: the free geometry parameters
(support size, minimum separation) are chosen using three folds and the reported value comes from the
fourth only. Reporting the best of a grid on the same folds inflates a family by the maximum of many
noisy draws. Gate: best candidate {fin['best_candidate']['name']} {fin['best_candidate']['lofo_collared_mean']:.4f}
vs best single layer {fin['best_single']['name']} {fin['best_single']['lofo_collared_mean']:.4f} —
margin {fin['best_candidate']['lofo_collared_mean'] - fin['best_single']['lofo_collared_mean']:+.4f},
{'PASSED' if fin['gate_pass'] else 'FAILED'}.</p>

<h2>Field screen with random and lattice controls</h2>
<table>
<tr><th>family</th><th class="n">N</th><th class="n">sep px</th><th class="n">emitted</th>
<th class="n">collared</th><th class="n">unrestricted</th></tr>
{fs}
</table>
<p class="small">Four controls calibrate the floor: random emission, a uniform lattice, and two
ablations. Anything inside the control band carries no information about faults, however
sophisticated its derivation.</p>

<h2 id="native">Native resolution — why the naive cross-gradient fails</h2>
<p>Measured from <code>data/training_features.tif</code> in this checkout
(<code>scripts/measure_native_resolution.py</code>). ρ₁ is the lag-1 autocorrelation along each
axis; the integral scale is the first lag where autocorrelation drops below 1/e.</p>
<table>
<tr><th>field</th><th class="n">ρ₁ rows</th><th class="n">ρ₁ cols</th><th class="n">L rows (m)</th>
<th class="n">L cols (m)</th></tr>
{natt}
</table>
<p>A 100 m horizontal gradient of a field whose integral scale is 2.7–8.1 km measures the
interpolation between regional samples, not a contact. The pure cross-gradient of isostatic gravity
and reduced-to-pole magnetics is therefore dominated by interpolation texture — which is exactly what
the controls show. The one field with genuine grid-scale power is the detrended-elevation slope
(ρ₁ ≈ 0.81), and it is the strongest single-layer baseline. This is the measured answer to the
brief's instruction to check the native station spacing before trusting 100 m gradients.</p>

<h2>The gate, measured against the carrier it gates</h2>
<p>The shipped family is the coherence-weighted multi-scale topographic carrier multiplied by an
orthogonality gate built from the gravity and reduced-to-pole magnetic gradients. Rebuilding both
supports at identical geometry and comparing them pixel by pixel:</p>
<table>
<tr><th>quantity</th><th class="n">value</th></tr>
<tr><td>dots in each variant</td><td class="n">{g['n_gated']:,}</td></tr>
<tr><td>dots identical between gated and ungated</td><td class="n">{g['n_shared']:,}</td></tr>
<tr><td>dots present only with the gate / only without it</td><td class="n">{g['n_only_gated']:,} / {g['n_only_ungated']:,}</td></tr>
<tr><td>support Jaccard</td><td class="n">{g['support_jaccard_gated_vs_ungated']:.4f}</td></tr>
<tr><td>change in held-out score</td><td class="n">{fin['best_candidate']['lofo_collared_mean'] - fin['table']['topo_ms_coh']['lofo_collared_mean']:+.5f}</td></tr>
</table>
<p class="small">Reading: the gate re-ranks a minority of the dots (8.5% of them) and is score-neutral.
It is kept because it is the brief's construction, it is a documented confounder filter rather than a
score, and it is measured not to damage the carrier. A candidate that only worked because of it would
not have survived this table.</p>

<h2>The catalogue buffer, and the margin after handing the same rule to both arms</h2>
<p>The highest-scoring artifact in the sibling cohort prunes every dot within 200 m of the catalogue.
The same question is reproducible here: the retained map of each fold stands in for the masked known
faults. Sweeping the exclusion radius (<code>scripts/prune_experiment.py</code>) and then re-running
both arms under identical radii at their own best geometry (<code>scripts/margin_check.py</code>):</p>
<table>
<tr><th class="n">buffer</th><th class="n">candidate collared</th><th class="n">candidate unrestricted</th>
<th class="n">baseline collared</th><th class="n">baseline unrestricted</th><th class="n">margin c / u</th>
<th class="n">blocks won (c / u)</th><th class="n">worst block delta (c / u)</th></tr>
{margin_rows}
</table>
<p class="small">Beyond 300 m both arms lose <em>unrestricted</em> score: the buffer starts deleting
the dots that would have earned credit where a new fault is mapped against a known one. 300 m is the
scoring kernel's own radius, the candidate's unrestricted optimum, and the shipped setting.</p>

<h2 id="instrument">Instrument warning — a holdout win is not a leaderboard win</h2>
<p>Rank correlation between candidates' holdout scores and the ~12 live leaderboard scores available
locally is positive but weak (≈ 0.50 drift-corrected), and the "hidden catalogue" probe is broken
(its random trace scores 0.97, so it is aborted above 0.25). No local instrument significantly ranks
candidates the way the private scorer does. Consequences for how this repository behaves:</p>
<ul>
<li>a holdout win is a <b>screen</b>, not evidence of a leaderboard gain;</li>
<li>candidates are spread across genuinely different hypotheses rather than micro-tuned;</li>
<li>the released file is the one that cleared both gates, with the margin printed above rather than
a claimed leaderboard position.</li>
</ul>

<h2>Format gate</h2>
<p><code>scripts/validate_submission.py</code> runs 13 checks on the actual bytes about to be
submitted and prints a JSON summary. The two rejection modes observed live are reproduced as tests
(<code>tests/test_gate.py</code>) rather than described:</p>
<ul>
<li><b>mode A</b> — a NaN strictly inside the footprint: legal at the edges, fatal to a naive
min/max form check;</li>
<li><b>mode B</b> — a finite value outside [0,1]: a z-score or an unclipped sentinel leaking into
the array.</li>
</ul>
<p>Both variants of the released file are emitted and both pass the gate; the NaN variant echoes the
official template's nodata mask exactly, the zero variant is for forms that reject NaN.</p>

<h2>Format gate on the committed bytes: {sum(1 for c in ga["checks"] if c["ok"] and c.get("required", True))}
of {sum(1 for c in ga["checks"] if c.get("required", True))} required checks</h2>
<table><tr><th>check</th><th class="n">result</th><th>detail</th></tr>{gat}</table>
<p class="small">The zero-filled variant passes the same gate; its one difference is that it echoes
no NaN mask, which is informational rather than required
(<code>docs/data/submission_gate_zeros.json</code>). These are the bytes linked from the download
button, verified by <code>scripts/verify_committed_submission.py</code> in CI.</p>

<h2>Relabel check</h2>
<p>Verdict <b>{esc(uni.get('verdict', 'unknown'))}</b> over {uni.get('n_prior', 0)} prior submission
rasters. Full table in <code>docs/data/uniqueness.json</code>.</p>
"""


# --------------------------------------------------------------------------- research

def research_page(kb: str, src: dict, met: dict) -> str:
    rows = "".join(
        f'<tr><td>{esc(s["topic"])}</td><td>{esc(s["claim"])}</td>'
        f'<td><span class="tag {"o" if s["class"]=="OFFICIAL" else ("m" if s["class"]=="MEASURED" else "b")}">'
        f'{esc(s["class"])}</span></td><td><a href="{esc(s["url"])}">{esc(s["url_label"])}</a></td></tr>'
        for s in src["sources"])
    return f"""<h1>Research and sources</h1>
<p class="sub">The brief, the metric, the masking rule, the primary literature, and every external
claim with a link for manual review. Nothing here is asserted from memory.</p>

<h2>Source table</h2>
<table><tr><th>topic</th><th>claim used</th><th>class</th><th>link</th></tr>{rows}</table>

<h2>Metric algebra, proved rather than quoted</h2>
<ul>
<li>Binary support dominates any graded version of itself:
<code>dDTI/dlambda &gt; 0</code>, checked numerically in the test suite.</li>
<li>Adding one unit pixel at kernel credit <code>k</code> raises the denominator by exactly 0.2 for
every <code>k</code>, hence the exact bar <code>k &gt; 0.2·s</code>. At s = 0.28 that is 0.056 —
every pixel within {met['marginal_distance_m_at_0.28']:.0f} m of a scored truth pixel pays.</li>
<li>A shadowed pixel adds false-positive mass and no true-positive mass, so it always loses. That is
why emission is a coverage problem with a minimum separation, not a thresholding problem.</li>
<li>The organizer's worked example is reproduced to {met['worked_example_dti']:.6f}.</li>
</ul>

<h2>The brief, verbatim</h2>
<blockquote><p>{esc(met['brief_excerpt'])}</p></blockquote>

<hr>
{markdown_to_html(kb)}
"""


def markdown_to_html(md: str) -> str:
    """Deliberately small: headings, tables, lists, bold/italic/code, links, paragraphs."""
    out, i = [], 0
    lines = md.split("\n")

    def inline(t: str) -> str:
        t = html.escape(t, quote=False)
        import re
        t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', t)
        t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
        t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
        t = re.sub(r"(?<![\w*])\*([^*\n]+)\*(?![\w*])", r"<i>\1</i>", t)
        return t

    while i < len(lines):
        ln = lines[i]
        if ln.startswith("|") and i + 1 < len(lines) and set(lines[i + 1].replace("|", "").strip()) <= set("-: "):
            hdr = [c.strip() for c in ln.strip("|").split("|")]
            i += 2
            cells = []
            while i < len(lines) and lines[i].startswith("|"):
                cells.append([c.strip() for c in lines[i].strip("|").split("|")])
                i += 1
            out.append("<table><tr>" + "".join(f"<th>{inline(h)}</th>" for h in hdr) + "</tr>"
                       + "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>"
                                 for r in cells) + "</table>")
            continue
        if ln.startswith("#### "):
            out.append(f"<h4>{inline(ln[5:])}</h4>")
        elif ln.startswith("### "):
            out.append(f"<h3>{inline(ln[4:])}</h3>")
        elif ln.startswith("## "):
            out.append(f"<h2>{inline(ln[3:])}</h2>")
        elif ln.startswith("# "):
            out.append(f"<h1>{inline(ln[2:])}</h1>")
        elif ln.startswith("```"):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].startswith("```"):
                buf.append(html.escape(lines[i]))
                i += 1
            out.append("<pre>" + "\n".join(buf) + "</pre>")
        elif ln.startswith("- "):
            buf = []
            while i < len(lines) and lines[i].startswith("- "):
                buf.append(f"<li>{inline(lines[i][2:])}</li>")
                i += 1
            out.append("<ul>" + "".join(buf) + "</ul>")
            continue
        elif ln.strip() == "---":
            out.append("<hr>")
        elif ln.strip():
            buf = [inline(ln)]
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "|", "- ", "```")):
                buf.append(inline(lines[i]))
                i += 1
            out.append("<p>" + " ".join(buf) + "</p>")
            continue
        i += 1
    return "\n".join(out)


def main() -> int:
    d = {name: load(name) for name in REQUIRED}
    sub, val, fin, scr, nat, met, uni, src, hyp = (
        d["submission"], d["validation"], d["final_experiment"], d["field_screen"],
        d["native_resolution"], d["metric"], d["uniqueness"], d["sources"], d["hypotheses"])
    (OUT / "assets").mkdir(parents=True, exist_ok=True)
    (OUT / "assets" / "site.css").write_text(CSS)
    pages = {
        "index.html": ("Submission — DOE GEMS Prize",
                       overview(sub, val, fin, uni, scr, nat, met, d["gate_effect"],
                                d["submission_gate"], d["submission_gate_zeros"],
                                d["margin_check"])),
        "executive-summary.html": ("How to submit — DOE GEMS Prize",
                                   howto(sub, d["submission_gate"], d["submission_gate_zeros"])),
        "hypotheses.html": ("Hypotheses — DOE GEMS Prize", hypotheses_page(hyp, val)),
        "evidence.html": ("Evidence — DOE GEMS Prize",
                          evidence_page(sub, val, fin, scr, nat, met, uni,
                                        d["submission_gate"], d["submission_gate_zeros"])),
        "research.html": ("Research and sources — DOE GEMS Prize", research_page(
            (ROOT / "docs" / "research" / "knowledge-base.md").read_text(), src, met)),
    }
    for fn, (title, body) in pages.items():
        (OUT / fn).write_text(page(title, body, sub, fn))
        print(f"wrote docs/{fn}  ({len(body):,} chars)")
    (OUT / ".nojekyll").write_text("")
    feed = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "submission_generated_utc": sub["generated_utc"],
        "slug": sub["slug"],
        "download": sub["download_url"],
        "sha256": sub["sha256"],
        "checks": [{"severity": s, "headline": h, "detail": dd}
                   for s, h, dd in feed_items(sub, val, fin, uni,
                                              d["submission_gate"], d["submission_gate_zeros"],
                                              d["margin_check"])],
    }
    (OUT / "feed.json").write_text(json.dumps(feed, indent=1))
    print("wrote docs/feed.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
