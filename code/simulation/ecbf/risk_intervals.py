"""Evaluation-only conditional-risk enclosures.

For a radial bin [a,b], let m=(a+b)/2 and h=(b-a)/2. Since each
clearance and their minimum are 1-Lipschitz, |Y(A,phi)-Y(m,phi)|<=h.
Thus P(Y(m)>c+h) <= P(Y(A)>c) <= P(Y(m)>c-h), pointwise in A.
Weighted bin integration gives an analytical enclosure without a quadrature
error assumption. Largest-uncertainty bisection refines it.

LIMIT: the angular set routine and the affine antiderivative use ordinary
IEEE floating point with explicit padding. These are analytical enclosures
with numerical allowances, NOT outward-rounded machine-certified intervals.
The classification metadata always records this limitation. No controller
module imports this file.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
import heapq
import math
import numpy as np
from .risk import exceedance_given_radius

ANGULAR_PAD = 1e-10
PROBABILITY_PAD = 1e-10

@dataclass(frozen=True)
class RiskInterval:
    lower: float
    upper: float
    evaluations: int
    bins: int
    method: str
    machine_certified: bool = False

    def as_dict(self):
        return asdict(self)


def affine_uniform_risk(c, lo, hi):
    """Integral acos(c/A)/pi under U[lo,hi], in an isolated affine region."""
    if c < 0:
        return 1.0
    if c == 0:
        return 0.5
    if c >= hi:
        return 0.0
    a = max(lo, c)
    def primitive(A):
        if A <= c:
            return 0.0
        r = math.sqrt(max(0.0, (A-c)*(A+c)))
        return (A*math.acos(min(1.0, c/A)) - c*math.asinh(r/c))/math.pi
    return max(0.0, (primitive(hi)-primitive(a))/(hi-lo))


def risk_interval(scene, z, c, beta, c_max, *, threshold=None, max_bins=128,
                  tolerance=1e-5, atom_radius=None, use_affine=True):
    """Enclose P(Y>c). Stop when threshold is separated, width<=tolerance,
    or max_bins is reached. Return an unresolved interval at the work cap.
    beta is the evaluator's true mixture weight, never a policy input.
    """
    z = np.asarray(z, float)
    c, c_max = float(c), float(c_max)
    if z.shape != (2,) or not np.all(np.isfinite(z)) or not math.isfinite(c):
        raise ValueError('risk input must be finite with a two-dimensional nominal point')
    if not math.isfinite(c_max) or c_max <= 0 or max_bins < 2 or tolerance < 0:
        raise ValueError('invalid risk enclosure configuration')
    if atom_radius is None and (not math.isfinite(beta) or not 0 <= beta <= 1):
        raise ValueError('mixture weight must be in [0,1]')
    if atom_radius is not None and not 0 <= atom_radius <= c_max:
        raise ValueError('atom radius violates the declared support')
    if c < 0:
        return RiskInterval(1.0, 1.0, 0, 0, 'nonnegative-loss', True)
    support = c_max if atom_radius is None else atom_radius
    if c >= support:
        return RiskInterval(0.0, 0.0, 0, 0, 'loss-support', True)
    g = scene.feature_clearances(z)
    if not np.all(np.isfinite(g)):
        raise ValueError('nonfinite feature clearance')
    # Explicit allowance in length units. This is not an interval-arithmetic proof.
    length_pad = 256*np.finfo(float).eps*max(1.0, float(np.max(np.abs(z))),
                                          float(np.max(np.abs(g))), c_max)
    evaluations = 0
    def angular(cq, radius):
        nonlocal evaluations
        if cq < 0:
            return 1.0
        if radius <= cq:
            return 0.0
        evaluations += 1
        p = float(exceedance_given_radius(scene, z, cq, radius, g))
        if not math.isfinite(p):
            raise ArithmeticError('nonfinite angular probability')
        return min(1.0, max(0.0, p))
    def padded(lo, hi):
        return max(0.0, lo-PROBABILITY_PAD), min(1.0, hi+PROBABILITY_PAD)
    if atom_radius is not None:
        low = max(0.0, angular(c+length_pad, atom_radius)-ANGULAR_PAD)
        up = min(1.0, angular(c-length_pad, atom_radius)+ANGULAR_PAD)
        low, up = padded(low, up)
        return RiskInterval(low, up, evaluations, 0, 'angular-atom')
    # A whole-ball feature check establishes a single affine clearance.
    if use_affine:
        ok, offset = scene.current_contribution(z, c_max+4*length_pad)
        if ok and math.isinf(offset):
            lo = ((1-beta)*affine_uniform_risk(c+length_pad, 0.0, c_max/2)
                  + beta*affine_uniform_risk(c+length_pad, c_max/2, c_max))
            up = ((1-beta)*affine_uniform_risk(c-length_pad, 0.0, c_max/2)
                  + beta*affine_uniform_risk(c-length_pad, c_max/2, c_max))
            lo, up = padded(lo, up)
            return RiskInterval(lo, up, 0, 0, 'affine-antiderivative')
    heap, counter = [], 0
    def add_bin(a, b, density):
        nonlocal counter
        mid, half = (a+b)/2, (b-a)/2
        mass = density*(b-a)
        low = mass*max(0.0, angular(c+half+length_pad, mid)-ANGULAR_PAD)
        up = mass*min(1.0, angular(c-half-length_pad, mid)+ANGULAR_PAD)
        counter += 1
        heapq.heappush(heap, (-(up-low), counter, a, b, density, low, up))
    for weight, a, b in ((1-beta, 0.0, c_max/2), (beta, c_max/2, c_max)):
        if weight > 0 and b > c:
            add_bin(max(a,c), b, weight/(b-a))
    while True:
        low, up = padded(math.fsum(v[5] for v in heap), math.fsum(v[6] for v in heap))
        separated = threshold is not None and (up < threshold or low > threshold)
        if separated or up-low <= tolerance or len(heap) >= max_bins or not heap:
            return RiskInterval(low, up, evaluations, len(heap), 'radial-lipschitz')
        _, _, a, b, density, _, _ = heapq.heappop(heap)
        mid = (a+b)/2
        if mid in (a,b):
            add_bin(a,b,density)
            return RiskInterval(low, up, evaluations, len(heap), 'radial-lipschitz-resolution-limit')
        add_bin(a, mid, density)
        add_bin(mid, b, density)
