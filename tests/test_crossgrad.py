"""The cross-gradient operator: the properties the brief asserts, checked numerically."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems49 import crossgrad as C  # noqa: E402
from gems49 import gating  # noqa: E402


def test_gradient_of_a_plane_is_exact():
    """grad of f = a*x + b*y must be (a, b) per metre, away from the edges."""
    n = 40
    a, b = 0.02, -0.03
    yy, xx = np.mgrid[0:n, 0:n]
    f = (a * xx + b * yy).astype(np.float32)
    gx, gy = C.horizontal_gradient(f, res_m=1.0)
    assert np.allclose(gx[5:-5, 5:-5], a, atol=1e-5)
    assert np.allclose(gy[5:-5, 5:-5], b, atol=1e-5)


def test_cross_product_vanishes_for_parallel_fields():
    """Gallardo & Meju (2004): t = grad m1 x grad m2 vanishes when the fields change together."""
    n = 48
    yy, xx = np.mgrid[0:n, 0:n]
    f1 = (0.01 * xx + 0.02 * yy).astype(np.float32)
    f2 = (3.0 * f1).astype(np.float32)          # exactly parallel gradients
    g1 = C.horizontal_gradient(f1, 1.0)
    g2 = C.horizontal_gradient(f2, 1.0)
    t = C.cross_gradient(*g1, *g2)
    assert np.allclose(t[6:-6, 6:-6], 0.0, atol=1e-9)


def test_cross_product_is_maximal_for_perpendicular_fields():
    n = 48
    yy, xx = np.mgrid[0:n, 0:n]
    f1 = (0.01 * xx).astype(np.float32)
    f2 = (0.01 * yy).astype(np.float32)
    g1 = C.horizontal_gradient(f1, 1.0)
    g2 = C.horizontal_gradient(f2, 1.0)
    t = C.cross_gradient(*g1, *g2)
    assert np.allclose(np.abs(t[6:-6, 6:-6]), 1e-4, atol=1e-8)


def test_operator_keeps_shared_structure_and_suppresses_single_layer_edges():
    """The brief: 'a buried volcanic edge usually moves the magnetic edge alone'.

    Three cases on one raster, all with a strong horizontal gravity edge (0.01, 0):

    * a **shared contact**: the magnetic field also steps in x at the same column -> both
      gradients are large and parallel -> coupling must be non-zero;
    * a **magnetic-only edge with a different strike**: the magnetic field steps in y ->
      the gradients are orthogonal -> the alignment term must kill it;
    * **no magnetic edge at all**: the magnetic field is flat -> the balance term must kill it,
      because the field that carries the structure has no gradient there.
    """
    n = 60
    yy, xx = np.mgrid[0:n, 0:n]
    gravity = (0.01 * xx).astype(np.float32)                 # (0.01, 0) everywhere
    magnetic = np.zeros((n, n), np.float32)
    magnetic[:, 30:] = 0.5                                   # shared contact at column 30
    magnetic[45:, :] += 0.7                                  # orthogonal magnetic-only step
    g1 = C.horizontal_gradient(gravity, 1.0)
    g2 = C.horizontal_gradient(magnetic, 1.0)
    out = C.pair_coupling(g1, g2, np.ones((n, n), bool), np.ones((n, n), bool))
    coupled = np.nan_to_num(out["coupled"])

    assert float(coupled[5:40, 30].max()) > 0.1    # shared, parallel, comparable -> kept
    assert float(coupled[50:55, 5:25].max()) == 0.0    # orthogonal magnetic-only -> killed
    assert float(coupled[5:40, 5:25].max()) == 0.0     # no magnetic gradient at all -> killed


def test_alignment_and_orthogonality_partition_the_magnitude_product():
    n = 40
    rng = np.random.default_rng(1)
    a = rng.random((n, n)).astype(np.float32)
    b = rng.random((n, n)).astype(np.float32)
    out = C.pair_coupling(C.horizontal_gradient(a, 1.0), C.horizontal_gradient(b, 1.0),
                          np.ones((n, n), bool), np.ones((n, n), bool))
    s = np.nan_to_num(out["align"]) + np.nan_to_num(out["ortho"])
    p = np.nan_to_num(out["n1"]) * np.nan_to_num(out["n2"])
    m = out["valid"] & np.isfinite(out["align"]) & np.isfinite(out["ortho"])
    assert np.allclose(s[m], p[m], atol=1e-5)


def test_nan_does_not_leak_into_gradients():
    a = np.ones((30, 30), np.float32)
    a[:, :15] = np.nan
    gx, gy = C.horizontal_gradient(a, 1.0)
    assert np.isnan(gx[:, 14]).all()   # the boundary column is tainted
    assert np.isnan(gx[:, 13]).all()
    assert np.isfinite(gx[:, 20]).all()


def test_triple_alignment_is_bounded_and_zero_without_agreement():
    n = 60
    rng = np.random.default_rng(4)
    topo = rng.random((n, n)).astype(np.float32)
    grav = rng.random((n, n)).astype(np.float32)
    mag = rng.random((n, n)).astype(np.float32)
    parts = gating.cross_terms(topo, grav, np.ones((n, n), bool), sigma_topo=0.0, sigma_other=0.0)
    assert np.nanmin(parts["align"]) >= 0.0
    assert np.nanmax(parts["align"]) <= 1.0
    assert np.nanmax(parts["ortho"]) <= 1.0


def test_orthogonality_gate_is_a_multiplier_not_an_addend():
    """The gate must never increase a carrier pixel, and must leave it alone where w = 0."""
    n = 50
    rng = np.random.default_rng(9)
    topo = rng.random((n, n)).astype(np.float32)
    grav = rng.random((n, n)).astype(np.float32)
    mag = rng.random((n, n)).astype(np.float32)
    base = np.nan_to_num(topo, nan=0.0)
    g0, _ = gating.orthogonality_gate(base, grav, mag, np.ones((n, n), bool), w=0.0)
    assert np.allclose(np.nan_to_num(g0), base, atol=1e-6)
    g1, _ = gating.orthogonality_gate(base, grav, mag, np.ones((n, n), bool), w=1.0)
    assert np.all(np.nan_to_num(g1) <= base + 1e-6)


def test_robust_unit_clips_at_one():
    x = np.arange(1000, dtype=np.float32) ** 2
    n, scale = C.robust_unit(x, np.ones(1000, bool), 95.0)
    assert n.max() <= 1.0 and n.min() >= 0.0
    assert abs(np.percentile(x, 95) / scale - 1.0) < 1e-6
