"""Estimated exceedance risk by angular integration and radial quadrature (report Sections 1.1.1 and 1.1.2)."""
from __future__ import annotations

import numpy as np

TWO_PI = 2 * np.pi
_U, _WU = np.polynomial.legendre.leggauss(24)
_U = 0.5 * (_U + 1)
_WU = 0.5 * _WU


def _union_measure(intervals):
    """Total length of a union of arcs given as (start, end) with end - start in [0, 2 pi]."""
    if not intervals:
        return 0.0
    segs = []
    for s, e in intervals:
        L = e - s
        if L >= TWO_PI - 1e-15:
            return TWO_PI
        s = s % TWO_PI
        e = s + L
        if e <= TWO_PI:
            segs.append((s, e))
        else:
            segs.append((s, TWO_PI))
            segs.append((0.0, e - TWO_PI))
    segs.sort()
    total, cs, ce = 0.0, segs[0][0], segs[0][1]
    for s, e in segs[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs
            cs, ce = s, e
    return total + ce - cs


def _circle_arc(z, o, r, A):
    """Angles phi with ||z + A e(phi) - o|| < r (z outside the disc): one arc (centre, half-width) or None."""
    dv = o - z
    d = float(np.hypot(*dv))
    if r <= 0:
        return None
    arg = (A * A + d * d - r * r) / (2 * A * d)
    if arg >= 1:
        return None
    half = np.pi if arg <= -1 else float(np.arccos(arg))
    phi = float(np.arctan2(dv[1], dv[0]))
    return (phi - half, phi + half)


def _halfplane_arc(n, gap, A):
    """Angles with n^T(A e(phi)) < -gap (gap >= 0: current clearance minus the threshold): centred at -n."""
    if gap >= A:
        return None
    half = np.pi if gap <= -A else float(np.arccos(gap / A))
    phi = float(np.arctan2(-n[1], -n[0]))
    return (phi - half, phi + half)


def _rect_arcs_general(z, rect, r, A, rho_body, thr):
    """Exact arcs of {phi : sdist(z + A e(phi), B) - rho_body < thr} by boundary crossings and midpoint tests."""
    from .geometry import sdist_rects
    x0, y0, x1, y1 = rect
    ang = [0.0]
    for X in (x0 - r, x1 + r, x0, x1):
        cth = (X - z[0]) / A
        if -1 < cth < 1:
            t = float(np.arccos(cth))
            ang += [t, -t]
    for Yv in (y0 - r, y1 + r, y0, y1):
        sth = (Yv - z[1]) / A
        if -1 < sth < 1:
            t = float(np.arcsin(sth))
            ang += [t, np.pi - t]
    for cx, cy in ((x0, y0), (x0, y1), (x1, y0), (x1, y1)):
        arc = _circle_arc(z, np.array([cx, cy]), r, A)
        if arc is not None and arc[1] - arc[0] < TWO_PI - 1e-12:
            ang += [arc[0], arc[1]]
    ang = np.unique(np.mod(ang, TWO_PI))
    ends = np.append(ang[1:], ang[0] + TWO_PI)
    mids = 0.5 * (ang + ends)
    pts = z + A * np.column_stack([np.cos(mids), np.sin(mids)])
    inside = sdist_rects(pts, rect[None])[:, 0] - rho_body < thr
    return [(float(s), float(e)) for s, e, k in zip(ang, ends, inside) if k]


def exceedance_given_radius(scene, z, c, A, g=None):
    """Angular risk p_t(c | A) = P(Y > c | ||w|| = A) at the predicted nose point z (report eq. (15))."""
    if g is None:
        g = scene.obstacle_clearances(z)
    psi0 = float(g.min())
    thr = psi0 - c
    nc, nw = len(scene.circles), len(scene.walls)
    arcs = []
    for i in np.flatnonzero(g - A < thr):                          # 1-Lipschitz: other obstacles cannot be reached
        if i < nc:
            x, y, rho = scene.circles[i]
            a = _circle_arc(z, np.array([x, y]), thr + rho + scene.rho_body, A)
            if a is not None:
                arcs.append(a)
        elif i < nc + nw:
            a = _halfplane_arc(scene.walls[i - nc, :2], g[i] - thr, A)
            if a is not None:
                arcs.append(a)
        else:
            arcs += _rect_arcs_general(z, scene.rects[i - nc - nw], thr + scene.rho_body, A, scene.rho_body, thr)
    return _union_measure(arcs) / TWO_PI


def _piece(a0, a1, lo, hi):
    """Nodes and weights on [a0, a1] with A = a0 + (a1 - a0) u^2 (square-root onset at a0), density 1/(hi - lo)."""
    A = a0 + (a1 - a0) * _U ** 2
    return A, _WU * 2 * (a1 - a0) * _U / (hi - lo)


def exceedance_risk(scene, z, cs, beta, c_max):
    """Estimated exceedance risk of each c in cs (report eq. (19)) under the mixture
    (1 - beta) U[0, c_max/2] + beta U[c_max/2, c_max].

    Obstacle i can contribute only for A above its entry magnitude psi_i(z) - psi(z) + c (report eq. (14)), where its
    angular set appears with a square-root onset; each component interval is therefore split at these entry magnitudes
    and every piece uses the u^2 substitution with the 24-point Gauss-Legendre rule (report eqs. (17) and (18);
    error below 4e-5 against adaptive quadrature)."""
    z = np.asarray(z, float)
    g = scene.obstacle_clearances(z)
    psi0 = float(g.min())
    out = np.zeros(len(cs))
    for k, c in enumerate(cs):
        tot = 0.0
        for weight, lo, hi in ((1 - beta, 0.0, c_max / 2), (beta, c_max / 2, c_max)):
            a = max(lo, c)
            if weight <= 0 or hi <= a:
                continue
            br = sorted({a, hi} | {float(gi - psi0 + c) for gi in g if a < gi - psi0 + c < hi})
            for a0, a1 in zip(br[:-1], br[1:]):
                A, w = _piece(a0, a1, lo, hi)
                tot += weight * float(np.dot(w, [exceedance_given_radius(scene, z, c, x, g) for x in A]))
        out[k] = tot
    return out


def monte_carlo_risk(scene, z, cs, beta, c_max, n, rng):
    """Monte Carlo estimate of the exceedance risk, for checking exceedance_risk (not used by the simulation)."""
    comp = rng.random(n) < beta
    A = np.where(comp, rng.uniform(c_max / 2, c_max, n), rng.uniform(0, c_max / 2, n))
    phi = rng.uniform(0, TWO_PI, n)
    w = A[:, None] * np.column_stack([np.cos(phi), np.sin(phi)])
    Y = np.maximum(scene.psi(z) - scene.psi(z + w), 0.0)
    return np.array([np.mean(Y > c) for c in cs])
