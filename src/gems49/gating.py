"""Cross-gradient information used as a CONFOUNDER GATE, not as positive evidence.

The brief is explicit on the role of the potential fields:

    "Score each pixel by magnitude-weighted alignment (both gradients large and parallel).
     Keep magnitude-weighted orthogonality as a separate confounder flag, not as evidence.
     A fault juxtaposing rocks of different density and susceptibility puts both edges on one
     trace, while a buried volcanic edge usually moves the magnetic edge alone, so the coupled
     score should suppress single-layer edges."

Measured on this grid, the *absolute* cross-gradient of gravity and magnetics carries no
detectable fault-location information: on the hide-and-recover instrument it scores 0.0715
against a no-information random control at 0.0718 (see docs/data/field_screen.json).  That is
consistent with the brief's own warning to "check the native gravity station spacing before
trusting 100 m gradients" -- the isostatic gravity field is smooth at the 100 m cell size, so
its gradient is dominated by interpolation structure rather than by contact edges.

What still works is the *relative* geometry: where a strong, linear topographic edge exists,
the cross-gradient angle against the potential-field gradients says whether that edge could be
a contact at all.  So this module uses the cross-gradient product as a multiplicative gate in
[0,1] applied to a sensitive topographic detector, and never as a score on its own.

    gate = 1 - w * (o_g * o_m)          o_f = |sin theta| * |G_f|_norm,  magnitude-weighted
                                        orthogonality of the topographic edge to field f
    gate = w * max(a_g, a_m) + (1-w)    a_f = |cos theta|^2 * |G_f|_norm,
                                        magnitude-weighted alignment

`orthogonality_gate` implements the first (the brief's "confounder flag ... not evidence"
reading); `alignment_gate` implements the second (the "both gradients large and parallel"
reading).  Both are ablations of the same operator and both are reported.
"""

from __future__ import annotations

import numpy as np

from .crossgrad import horizontal_gradient, magnitude, pair_coupling, robust_unit, smooth

NAN32 = np.float32(np.nan)


def _unit_gradients(field: np.ndarray, sigma: float):
    gx, gy = horizontal_gradient(smooth(field, sigma))
    m = magnitude(gx, gy)
    return gx, gy, m


def cross_terms(
    topo: np.ndarray,
    other: np.ndarray,
    footprint: np.ndarray,
    *,
    sigma_topo: float = 2.0,
    sigma_other: float = 2.0,
    pct: float = 95.0,
) -> dict:
    """Magnitude-weighted alignment and orthogonality of a topographic edge vs another field."""
    tx, ty, tm = _unit_gradients(topo, sigma_topo)
    ox, oy, om = _unit_gradients(other, sigma_other)
    ok = footprint & np.isfinite(tm) & np.isfinite(om)
    t_n, _ = robust_unit(tm, ok, pct)
    o_n, _ = robust_unit(om, ok, pct)
    eps = np.float32(1e-12)
    cos = np.clip((tx * ox + ty * oy) / np.maximum(tm * om, eps), -1.0, 1.0)
    sin = (tx * oy - ty * ox) / np.maximum(tm * om, eps)
    align = (t_n * o_n * np.square(cos)).astype(np.float32)          # magnitude-weighted alignment
    ortho = (t_n * o_n * np.square(sin)).astype(np.float32)          # magnitude-weighted orthogonality
    return {
        "topo_n": np.where(ok, t_n, NAN32).astype(np.float32),
        "other_n": np.where(ok, o_n, NAN32).astype(np.float32),
        "cos2": np.where(ok, np.square(cos), NAN32).astype(np.float32),
        "align": np.where(ok, align, NAN32).astype(np.float32),
        "ortho": np.where(ok, ortho, NAN32).astype(np.float32),
        "valid": ok,
    }


def orthogonality_gate(
    topo_edge: np.ndarray,
    gravity: np.ndarray,
    magnetic: np.ndarray,
    footprint: np.ndarray,
    *,
    w: float = 0.5,
    sigma: float = 2.0,
) -> tuple[np.ndarray, dict]:
    """``topo_edge * [1 - w * o_g * o_m]``: suppress edges orthogonal to BOTH potential fields."""
    tg = cross_terms(topo_edge, gravity, footprint, sigma_topo=0.0, sigma_other=sigma)
    tm = cross_terms(topo_edge, magnetic, footprint, sigma_topo=0.0, sigma_other=sigma)
    og = np.nan_to_num(tg["ortho"], nan=0.0)
    om = np.nan_to_num(tm["ortho"], nan=0.0)
    gate = (1.0 - w * og * om).astype(np.float32)
    base = np.nan_to_num(topo_edge, nan=0.0).astype(np.float32)
    score = np.where(footprint, base * np.clip(gate, 0.0, 1.0), NAN32).astype(np.float32)
    return score, {"gate": gate, "ortho_grav": og, "ortho_mag": om,
                   "align_grav": np.nan_to_num(tg["align"], nan=0.0),
                   "align_mag": np.nan_to_num(tm["align"], nan=0.0)}


def alignment_gate(
    topo_edge: np.ndarray,
    gravity: np.ndarray,
    magnetic: np.ndarray,
    footprint: np.ndarray,
    *,
    w: float = 0.5,
    sigma: float = 2.0,
) -> tuple[np.ndarray, dict]:
    """``topo_edge * [(1-w) + w * max(a_g, a_m)]``: keep edges corroborated by a potential field."""
    tg = cross_terms(topo_edge, gravity, footprint, sigma_topo=0.0, sigma_other=sigma)
    tm = cross_terms(topo_edge, magnetic, footprint, sigma_topo=0.0, sigma_other=sigma)
    ag = np.nan_to_num(tg["align"], nan=0.0)
    am = np.nan_to_num(tm["align"], nan=0.0)
    gate = ((1.0 - w) + w * np.maximum(ag, am)).astype(np.float32)
    base = np.nan_to_num(topo_edge, nan=0.0).astype(np.float32)
    score = np.where(footprint, base * np.clip(gate, 0.0, 1.0), NAN32).astype(np.float32)
    return score, {"gate": gate, "align_grav": ag, "align_mag": am}
