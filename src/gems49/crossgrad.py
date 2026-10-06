"""Cross-gradient structural coupling (CGSC) of gravity, magnetics and topography.

Theory
------
Gallardo & Meju (2004), *Joint two-dimensional DC resistivity and seismic travel time
inversion with cross-gradients*, JGR 109, B03311, doi:10.1029/2003JB002716, define the
cross-gradient function between two property fields m1, m2 as the cross product of their
gradients,

    t(x) = grad m1(x)  x  grad m2(x).

The published property that makes it useful is that ``t`` vanishes wherever the two fields
change in the same direction (gradients parallel) AND wherever either field is locally
constant.  It therefore measures *shared structure* rather than shared amplitude, and it
requires no petrophysical relationship between the two properties.  Fregoso & Gallardo
(2009), *Cross-gradients joint 3D inversion with applications to gravity and magnetic data*,
Geophysics 74(4), L31-L42, doi:10.1190/1.3119263, applied it to gravity and magnetics.  An
ASEG 2012 extended abstract (doi:10.1071/ASEG2012ab273) used the angle and cross product of
the *horizontal* gradients of gravity and reduced-to-pole magnetics directly in the data
domain, with no inversion -- which is the construction implemented here.

Implementation in the data domain
---------------------------------
For each field f we form the horizontal gradient vector
``G_f = (df/dx, df/dy)`` on the common 100 m grid after common-scale Gaussian smoothing.
Two scalars are then formed per pixel:

* **alignment energy** (the evidence term)
      A = |G_g| |G_m| cos^2(theta),      cos(theta) = <G_g,G_m>/(|G_g||G_m|)
  which is large only when BOTH gradients are large AND parallel.  cos^2 is used rather
  than cos so that polarity reversals (a density step up against a magnetisation step down,
  which is physically legitimate at a contact) score the same as matching polarity.
* **coupling balance** (the operator that suppresses single-layer edges)
      B = 2 min(|G_g|,|G_m|) / (|G_g| + |G_m|)   in [0,1]
  which is 1 only when the two fields carry comparable edge magnitude.  A buried volcanic
  edge moves the magnetic gradient alone, so B -> 0 there and the coupled score is
  suppressed -- this is the requirement that agreement be an *operator* and not a product
  of independent scores.

* **orthogonality flag** (never evidence; reported separately)
      O = |G_g| |G_m| sin^2(theta) = |t|^2
  the squared cross-gradient magnitude.  It is the confounder indicator: large O with large
  A is impossible, so O marks pixels where a strong edge in one field is unaccompanied in
  the other, or where the two fields genuinely disagree about strike.

The coupled score is  C = B * A_norm, and a confounder-masked variant subtracts the
orthogonality flag from the evidence by weighting with cos^2 computed on the *normalised*
vectors only where both magnitudes exceed a noise floor.

Everything here is label-free: the catalogue is never used to compute a score.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, sobel

from .spec import GRID

RES = GRID["res_m"]


def smooth(field: np.ndarray, sigma_px: float) -> np.ndarray:
    """NaN-aware Gaussian smoothing at a common scale."""
    f = np.asarray(field, dtype=np.float32)
    ok = np.isfinite(f)
    if sigma_px <= 0:
        return np.where(ok, f, np.nan)
    num = gaussian_filter(np.where(ok, f, np.float32(0.0)), sigma_px, mode="nearest")
    den = gaussian_filter(ok.astype(np.float32), sigma_px, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[den < 1e-6] = np.nan
    return out


def horizontal_gradient(field: np.ndarray, res_m: float = RES) -> tuple[np.ndarray, np.ndarray]:
    """Central-difference horizontal gradient (d/dx, d/dy) in units per metre.

    NaN is propagated: a NaN anywhere in the 3x3 stencil taints the pixel, so no gradient
    is ever invented across the footprint boundary.
    """
    f = np.asarray(field, dtype=np.float32)
    ok = np.isfinite(f)
    filled = np.where(ok, f, 0.0)
    gx = sobel(filled, axis=1, mode="nearest") / (8.0 * res_m)
    gy = sobel(filled, axis=0, mode="nearest") / (8.0 * res_m)
    tainted = ~ok
    for _ in range(1):
        t = tainted.copy()
        tainted = (
            t
            | np.roll(t, 1, 0) | np.roll(t, -1, 0)
            | np.roll(t, 1, 1) | np.roll(t, -1, 1)
        )
    gx[tainted] = np.nan
    gy[tainted] = np.nan
    return gx, gy


def magnitude(gx: np.ndarray, gy: np.ndarray) -> np.ndarray:
    return np.hypot(gx, gy)


def robust_unit(x: np.ndarray, valid: np.ndarray, pct: float = 95.0) -> tuple[np.ndarray, float]:
    """Scale a non-negative field so that its ``pct`` percentile over ``valid`` equals 1.

    Potential-field derivatives have heavy tails; a percentile scale keeps a single
    outlier from dominating the magnitude weighting.  Values are CLIPPED at 1, not
    renormalised, so a very large gradient cannot buy an unbounded score.
    """
    x = np.asarray(x, dtype=np.float32)
    sel = valid & np.isfinite(x)
    if not np.any(sel):
        return np.zeros_like(x), 1.0
    scale = float(np.percentile(x[sel], pct))
    if scale <= 0:
        scale = float(np.max(x[sel])) or 1.0
    return np.clip(x / scale, 0.0, 1.0), scale


def cross_gradient(gx1, gy1, gx2, gy2) -> np.ndarray:
    """Scalar cross product t_z = gx1*gy2 - gy1*gx2 (the 2-D cross-gradient)."""
    return gx1 * gy2 - gy1 * gx2


def pair_coupling(
    g1: tuple[np.ndarray, np.ndarray],
    g2: tuple[np.ndarray, np.ndarray],
    v1: np.ndarray,
    v2: np.ndarray,
    pct: float = 95.0,
) -> dict:
    """Cross-gradient coupling between two gradient-vector fields.

    Parameters
    ----------
    g1, g2 : (gx, gy) horizontal gradient pairs (per metre)
    v1, v2 : boolean masks of pixels where each field is defined

    Returns a dict with the normalised magnitudes, the alignment evidence, the balance
    operator, the orthogonality (cross-gradient) confounder flag, and the coupled score.
    """
    gx1, gy1 = g1
    gx2, gy2 = g2
    m1 = magnitude(gx1, gy1)
    m2 = magnitude(gx2, gy2)

    both = v1 & v2 & np.isfinite(m1) & np.isfinite(m2)
    n1, s1 = robust_unit(m1, both, pct)   # both in [0,1], 95th pct == 1
    n2, s2 = robust_unit(m2, both, pct)

    eps = 1e-12
    denom = n1 * n2
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_theta = (gx1 * gx2 + gy1 * gy2) / np.maximum(m1 * m2, eps)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        sin_theta = (gx1 * gy2 - gy1 * gx2) / np.maximum(m1 * m2, eps)
        # inside the domain a vanishing gradient is a real "no coupling" observation, so the
        # angle terms are defined as 0 there rather than NaN; NaN is reserved for undefined data.
        cos_theta = np.where(denom > 0, cos_theta, 0.0)
        sin_theta = np.where(denom > 0, sin_theta, 0.0)
        defined = denom > 0

    align = denom * np.square(cos_theta)          # |G1||G2|cos^2(theta), normalised
    ortho = denom * np.square(sin_theta)          # |G1||G2|sin^2(theta) == |t|^2, normalised
    balance = np.where(defined, 2.0 * np.minimum(n1, n2) / np.maximum(n1 + n2, eps), 0.0)
    coupled = balance * align

    for a in (align, ortho, balance, coupled):
        a[~both] = np.nan
    return {
        "n1": np.where(both, n1, np.nan),
        "n2": np.where(both, n2, np.nan),
        "align": align,          # evidence
        "ortho": ortho,          # confounder flag, NOT evidence
        "balance": balance,      # single-layer-edge suppressor
        "coupled": coupled,      # the operator output
        "scale1": s1,
        "scale2": s2,
        "valid": both,
    }


def cgsc(
    gravity: np.ndarray,
    magnetic: np.ndarray,
    elevation: np.ndarray,
    footprint: np.ndarray,
    *,
    sigma_common: float = 2.0,
    sigma_elevation: float = 3.0,
    pct: float = 95.0,
    pair: str = "gx m",
    third_weight: float = 0.0,
) -> dict:
    """Cross-gradient structural coupling of isostatic gravity, RTP magnetics and topography.

    ``pair`` selects which two fields carry the primary coupling ("gx m" = gravity x
    magnetic, the Fregoso & Gallardo pair).  ``third_weight`` mixes in the pairwise coupling
    of the primary field with the smoothed elevation surface; 0 disables the topographic
    term entirely.
    """
    fg = smooth(gravity, sigma_common)
    fm = smooth(magnetic, sigma_common)
    fe = smooth(elevation, sigma_elevation)

    vg = footprint & np.isfinite(fg)
    vm = footprint & np.isfinite(fm)
    ve = footprint & np.isfinite(fe)

    gg = horizontal_gradient(fg)
    gm = horizontal_gradient(fm)
    ge = horizontal_gradient(fe)

    primary = pair_coupling(gg, gm, vg, vm, pct)
    out = dict(primary)
    out["pair"] = "gx m"
    out["third"] = None
    if third_weight > 0.0:
        second = pair_coupling(gg, ge, vg, ve, pct)
        third = pair_coupling(gm, ge, vm, ve, pct)
        combined = third_weight * np.nan_to_num(third["coupled"], nan=0.0) + (
            (1.0 - third_weight) * np.nan_to_num(primary["coupled"], nan=0.0)
        )
        # only pixels where all three fields are defined keep a topographic contribution
        combined = np.where(np.isfinite(primary["coupled"]), combined, np.nan)
        out["coupled_primary"] = primary["coupled"]
        out["coupled"] = combined
        out["third_balance"] = third["balance"]
        out["second_coupled"] = second["coupled"]
    return out
