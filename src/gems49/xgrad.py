"""Cross-gradient structural coupling of gravity, magnetics and topography.

Method (per the session brief; references):

* Gallardo & Meju, JGR 2004, doi:10.1029/2003JB002716 — the cross-gradient
  function, the cross product of the gradient vectors of two property
  fields, vanishes where both fields change in the same direction: shared
  structure without any petrophysical assumption.
* Fregoso & Gallardo, Geophysics 74(4) L31-L42, 2009 — cross-gradient
  coupling applied jointly to gravity and magnetics.
* ASEG 2012 extended abstract, doi:10.1071/ASEG2012ab273 — angle and cross
  product of the *horizontal* gradients of gravity and reduced-to-pole
  magnetics used directly in the data domain (no inversion).

Implementation here:

1. The three fields — isostatic gravity anomaly (band 13), reduced-to-pole
   magnetic anomaly (band 2) and the detrended elevation surface (band 12)
   — are Gaussian-smoothed to a COMMON scale sigma and differenced with a
   central-difference horizontal gradient, so all three gradient vectors
   live on the same resolution.
2. Agreement across the independent signals is treated as an OPERATOR:
   for each pair (i, j) we compute the 2-D cross-gradient magnitude
       CG_ij = |df_i/dx * df_j/dy - df_i/dy * df_j/dx|
   and the sign-invariant alignment
       c_ij = cos^2(theta_ij) = (grad f_i . grad f_j)^2 / (|grad f_i|^2 |grad f_j|^2)
   The consensus is an order statistic (median of the three pairwise
   alignments = at-least-2-of-3 agreement), NOT a product of three
   normalised score maps.
3. The magnitude weight is a conjunction (min) of the three robustly
   normalised gradient magnitudes: an edge must be expressed in ALL fields
   at once.  A buried volcanic edge that moves only the magnetic field
   therefore gets a near-zero coupled score (single-layer suppression).
4. Magnitude-weighted orthogonality (|grad f_i||grad f_j| sin^2 theta) is
   returned as a SEPARATE confounder flag and is never added to the
   evidence score.

Memory discipline: every intermediate array is float32 (the sandbox has
4 GB RAM; the grid is 3730 x 3292 = 12.28 M px = 49 MB per float32 plane).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter

FIELDS = ("iso_grav_anom", "rtp", "det_elev")
_PAIRS = (
    ("iso_grav_anom", "rtp"),
    ("iso_grav_anom", "det_elev"),
    ("rtp", "det_elev"),
)


def _robust_unit(mag: np.ndarray, footprint: np.ndarray, q_hi: float = 0.99) -> np.ndarray:
    """Normalise a magnitude to [0,1] by its footprint q-quantiles."""
    finite = mag[footprint & np.isfinite(mag)]
    if finite.size == 0:
        return np.zeros_like(mag)
    lo = np.quantile(finite, 0.50)
    hi = np.quantile(finite, q_hi)
    if hi <= lo:
        return np.zeros_like(mag)
    return np.clip((mag - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def _smoothed_grad(field: np.ndarray, sigma_px: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gaussian-smooth (NaN-aware) then central-difference. Returns
    (gx, gy, |grad|), all float32."""
    fin = np.isfinite(field)
    num = gaussian_filter(np.where(fin, field, 0.0).astype(np.float32), sigma_px)
    wgt = gaussian_filter(fin.astype(np.float32), sigma_px)
    sm = np.divide(num, wgt, out=np.full_like(num, np.nan), where=wgt > 1e-3)
    gy, gx = np.gradient(sm)
    mag = np.hypot(gx, gy).astype(np.float32)
    return gx.astype(np.float32), gy.astype(np.float32), mag


@dataclass
class CoupledField:
    score: np.ndarray            # coupled evidence in [0,1] (float32)
    alignment: np.ndarray        # median pairwise cos^2 (0..1)
    magnitude_weight: np.ndarray # min of the three robust units
    ortho_flag: np.ndarray       # magnitude-weighted orthogonality (confounder)
    sigma_px: float


def coupled_score(
    fields: dict[str, np.ndarray],
    sigma_px: float,
    footprint: np.ndarray,
    pairs: tuple = _PAIRS,
    return_ortho: bool = True,
) -> CoupledField:
    """Magnitude-weighted cross-gradient alignment score (see module doc)."""
    grads, mags, units = {}, {}, {}
    for name, f in fields.items():
        gx, gy, mag = _smoothed_grad(f, sigma_px)
        grads[name] = (gx, gy)
        mags[name] = mag
        units[name] = _robust_unit(np.where(footprint, mag, np.nan), footprint)
        del gx, gy, mag

    c2_list = []
    ortho_max = np.zeros(footprint.shape, dtype=np.float32) if return_ortho else None
    for a, b in pairs:
        gax, gay = grads[a]
        gbx, gby = grads[b]
        dot = gax * gbx + gay * gby
        na2 = gax * gax + gay * gay
        nb2 = gbx * gbx + gby * gby
        denom = na2 * nb2
        c2 = np.where(denom > 0, dot * dot / np.maximum(denom, np.float32(1e-30)), 0.0)
        c2_list.append(c2.astype(np.float32))
        if return_ortho:
            cg = np.abs(gax * gby - gay * gbx)          # cross-gradient magnitude
            ortho = np.where(denom > 0, cg * cg / np.maximum(denom, np.float32(1e-30)), 0.0)
            np.maximum(ortho_max, ortho.astype(np.float32), out=ortho_max)
        del dot, na2, nb2, denom, c2
    del grads

    consensus = np.median(np.stack(c2_list, axis=0), axis=0).astype(np.float32)
    del c2_list
    weight = np.min(np.stack([units[n] for n in fields], axis=0), axis=0).astype(np.float32)
    score = (weight * consensus).astype(np.float32)
    score = np.where(footprint, score, np.float32(np.nan)).astype(np.float32)
    return CoupledField(score, consensus, weight, ortho_max, sigma_px)


def single_layer_baseline(name: str, fields: dict[str, np.ndarray], sigma_px: float,
                          footprint: np.ndarray) -> np.ndarray:
    """The repo's single-layer gradient baseline: one field's horizontal
    gradient magnitude alone, robust-normalised to [0,1]."""
    _, _, mag = _smoothed_grad(fields[name], sigma_px)
    u = _robust_unit(np.where(footprint, mag, np.nan), footprint)
    return np.where(footprint, u, np.float32(np.nan)).astype(np.float32)


def gravity_support_diagnostics(grav: np.ndarray, footprint: np.ndarray,
                                seed: int = 49) -> dict:
    """Empirical check of the scale at which the isostatic gravity anomaly
    carries independent information — the brief's "check the native gravity
    station spacing before trusting 100 m gradients".

    Two memory-light measurements:
    (a) the fraction of raw gradient variance retained after smoothing at
        300 m / 500 m / 1 km (a field interpolated from widely spaced
        stations loses little to such smoothing; a 100 m-scale field loses
        most of it);
    (b) a subsampled along-row empirical semivariogram, whose practical
        range is an upper-bound proxy for the native sampling scale.

    The official station inventory itself lives in the INGENIOUS gravity
    compilation (GDR submission 1391, DOI 10.15121/1881483); fetching that
    archive is blocked from this sandbox, so it is listed as a manual
    verification item rather than asserted here.
    """
    fin = footprint & np.isfinite(grav)
    g = np.where(fin, grav, 0.0).astype(np.float32)
    gy0, gx0 = np.gradient(g)
    v0 = float(np.mean((gx0 * gx0 + gy0 * gy0)[fin]))
    out = {"footprint_fraction": float(fin.mean()), "raw_grad_variance": v0}
    del gy0, gx0
    for label, sigma in (("300m", 3.0), ("500m", 5.0), ("1000m", 10.0)):
        s = gaussian_filter(g, sigma)
        gy, gx = np.gradient(s)
        v = float(np.mean((gx * gx + gy * gy)[fin]))
        out[f"grad_variance_retained_{label}"] = v / max(v0, 1e-30)
        del s, gy, gx
    # subsampled semivariogram along rows (lag in pixels)
    rng = np.random.default_rng(seed)
    rows = rng.choice(np.where(fin.any(axis=1))[0], size=min(400, int(fin.any(axis=1).sum())),
                      replace=False)
    max_lag = 60
    semi = np.zeros(max_lag, dtype=np.float64)
    cnt = np.zeros(max_lag, dtype=np.int64)
    vals = grav[rows].astype(np.float64)
    ok = np.isfinite(vals)
    for lag in range(1, max_lag):
        a = vals[:, :-lag]
        b = vals[:, lag:]
        m = ok[:, :-lag] & ok[:, lag:]
        d = (a - b) ** 2
        semi[lag] = d[m].mean() if m.any() else np.nan
        cnt[lag] = m.sum()
    semi[0] = 0.0
    out["semivariogram_lag_px"] = list(range(max_lag))
    out["semivariogram"] = [None if np.isnan(x) else float(x) for x in semi]
    # practical range: first lag where semivariance reaches 90% of the
    # large-lag sill
    valid = ~np.isnan(semi[1:])
    if valid.any():
        sill = np.nanmax(semi[30:]) if np.isfinite(semi[30:]).any() else np.nanmax(semi[1:])
        idx = np.argmax(semi[1:] >= 0.9 * sill) + 1 if np.isfinite(sill) else None
        out["practical_range_px"] = int(idx) if idx else None
    return out
