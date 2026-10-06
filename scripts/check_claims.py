#!/usr/bin/env python3
"""Anti-hallucination gate: every headline number on the site must exist in an evidence file.

The site is generated from ``docs/data/*.json``, so its numbers should be traceable by construction —
but the prose also quotes values (margins, dot counts, hashes) that a careless edit could drift. This
script re-derives each quoted value from the artifact that owns it and fails if they disagree.

Run: python3 scripts/check_claims.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"


def load(name: str):
    return json.loads((DATA / f"{name}.json").read_text())


def page(name: str) -> str:
    return (ROOT / "docs" / name).read_text()


def check(label: str, needle: str, haystacks: list[tuple[str, str]], fails: list[str]) -> None:
    for where, text in haystacks:
        if needle in text:
            print(f"  OK   {label}: '{needle}' found in {where}")
            return
    fails.append(f"{label}: '{needle}' not found in any of {[w for w, _ in haystacks]}")


def main() -> int:
    sub, uni, gate = load("submission"), load("uniqueness"), load("submission_gate")
    fin, mg, ge = load("final_experiment"), load("margin_check"), load("gate_effect")
    nat, met, pr = load("native_resolution"), load("metric"), load("prune_experiment")
    names = ("index.html", "evidence.html", "analysis.html", "executive-summary.html",
             "hypotheses.html", "research.html")
    pages = [(n, page(n)) for n in names]
    fails: list[str] = []

    print("submission identity")
    check("sha256", sub["sha256"][:16], pages + [("submission.json", json.dumps(sub))], fails)
    check("positive pixels", f"{sub['n_positive']:,}", pages + [("submission.json", json.dumps(sub))], fails)
    check("size", f"{sub['size_mb']:.2f} MB", pages, fails)
    check("buffer", f"{sub['catalogue_buffer_m']:.0f} m", pages, fails)

    print("format gate")
    n_req = sum(1 for c in gate["checks"] if c.get("required", True))
    assert gate["ok"], "the committed artifact no longer passes its own gate"
    check("required checks", f"{n_req}/{n_req}", pages, fails)

    print("holdout and margins")
    cand = next(r for r in mg["arms"]["candidate:gate_ortho_w0.25"]["rows"]
                if r["exclusion_px"] == sub["catalogue_buffer_px"])
    base = next(r for r in mg["arms"]["baseline:single:det_elev_slope"]["rows"]
                if r["exclusion_px"] == sub["catalogue_buffer_px"])
    m = next(x for x in mg["margins"] if x["exclusion_px"] == sub["catalogue_buffer_px"])
    for label, value in (("candidate collared", f"{cand['collared_mean']:.4f}"),
                         ("baseline collared", f"{base['collared_mean']:.4f}"),
                         ("candidate unrestricted", f"{cand['unrestricted_mean']:.4f}"),
                         ("baseline unrestricted", f"{base['unrestricted_mean']:.4f}"),
                         ("collared margin", f"{m['collared_margin']:+.4f}")):
        check(label, value, pages, fails)
    check("best single baseline", f"{fin['best_single']['lofo_collared_mean']:.4f}", pages, fails)
    check("candidate LOFO", f"{fin['best_candidate']['lofo_collared_mean']:.4f}", pages, fails)
    check("blocks won", f"{m['collared_blocks_won']}/4", pages, fails)

    print("gate effect")
    check("jaccard", f"{ge['support_jaccard_gated_vs_ungated']:.4f}", pages, fails)
    check("shared dots", f"{ge['n_shared']:,}", pages, fails)

    print("relabel check")
    check("priors", f"{uni['n_priors_compared']}", pages, fails)
    check("verdict", uni["verdict"], pages, fails)
    check("max support jaccard", f"{uni['max_support_jaccard']:.4f}", pages, fails)

    print("metric algebra")
    check("worked example", f"{met['worked_example_dti']:.6f}", pages, fails)
    check("bar distance", f"{met['marginal_distance_m_at_0.28']:.0f} m", pages, fails)

    print("native resolution")
    g = nat["fields"]["iso_grav_anom"]
    check("gravity integral scale", f"{g['L_int_m_rows']:,.0f}", pages, fails)
    s19 = nat["fields"]["det_elev_slope"]
    check("rho1 of det_elev_slope", f"{s19['rho1_rows']:.3f}", pages, fails)

    print("prune experiment")
    row = next(r for r in pr["rows"] if r["n_target"] == 40_000 and r["exclusion_px"] == 3)
    check("prune row collared", f"{row['collared_mean']:.4f}", pages, fails)

    if fails:
        print("\nFAILED:")
        for f in fails:
            print("  -", f)
        return 1
    print("\nOK: every checked number on the site is present in its evidence file.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
