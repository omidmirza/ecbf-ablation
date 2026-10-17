"""Statistical rules, adaptive baselines, and certified-margin selection."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .certification.monitor import RepeatedCertificationMonitor
from .certification.policies import make_policy


# ----------------------------------------------------------------------------------------------- certification
class CertificationMonitor:
    """Certification monitor with the current uncorrected rule ('V').

    The 'VI' label is retained for historical shared-correction compatibility;
    it does not refer to Section VI of the current paper.

    budget_mode='joint_log' requires common periodic overlapping tests. 'joint_log_once' keeps its first allocated
    level but allows no later starts (K=infinity). 'single' uses the full delta for its one initial test.
    With budget_mode='legacy', `restart=False` disables restarts and the historical thresholds are unchanged.
    certs[j] lists (start, end, index) in certification order; a certificate issued by observe() at time t has end t+1
    and is usable from decision time t+1 on.
    """

    def __init__(self, candidates, gains, alpha, delta, rule="V", restart=True, policy=None, budget_mode="legacy",
                 budget_normalizer=2.11):
        """policy: start-schedule spec of ecbf.certification.policies.make_policy, e.g. {"kind": "norestart"} (certify once) or
        {"kind": "fixed", "period": K, "overlap": True} (joint_log).
        Only legacy permits restart-after-certification and other historical start policies."""
        self.grid = np.asarray(candidates, float)
        self.m = len(self.grid)
        self.rule = rule
        if budget_mode != "legacy" and rule != "V":
            raise ValueError("single/joint_log budgets apply only to the uncorrected Section V rule")
        if budget_mode == "legacy":
            if policy is None or not restart:
                policy = {"kind": "restart"} if restart else {"kind": "norestart"}
        elif budget_mode in ("single", "joint_log_once"):
            if policy is None:
                policy = {"kind": "norestart"}
        elif budget_mode == "joint_log":
            if policy is None or not restart:
                raise ValueError("joint_log needs an explicit common periodic policy and restart=True")
        else:
            raise ValueError(f"unknown budget mode {budget_mode!r}")
        self.budget_mode = budget_mode
        self.policy_spec = dict(policy)
        policy = make_policy(policy)
        self.core = RepeatedCertificationMonitor(self.m, gains, alpha, delta, policy,
                                                 {"V": "none", "VI": "shared"}[rule], budget_mode=budget_mode,
                                                 budget_normalizer=budget_normalizer)
        self.certs = [[] for _ in range(self.m)]

    @property
    def t(self):
        return self.core.t

    def observe(self, Y):
        z = Y > self.grid
        new = self.core.observe(z)
        for c in sorted(new, key=lambda c: (c.candidate, c.index)):
            self.certs[c.candidate].append((int(c.start), int(c.end), int(c.index)))
        return new

    def any_certificate(self):
        return any(len(c) for c in self.certs)


def select_certified(certs, D, alpha, alpha0, t):
    """Select a deployable certified margin from the stored certificates.

    The smallest candidate j whose best stored certificate (end <= t) satisfies alpha + d_t^{(r)} <= alpha0 is deployed,
    d_t^{(r)} = max_{s in [start, end)} D[j, s]. The certificate used is the one with
    the smallest current bound, with ties resolved by the largest test index.
    Only observations before the certificate end enter its comparison window.
    Returns (j or m for c_max, r or 0, d of the used certificate or nan, per-candidate best d, nocert, notrans).
    """
    m = len(certs)
    best = np.full(m, np.inf)
    used = (m, 0, np.nan)
    any_cert = False
    for j in range(m):
        cj = [c for c in certs[j] if c[1] <= t]
        if not cj:
            continue
        any_cert = True
        ds = [(float(np.max(D[j, start:end])) if end > start else 1.0, idx) for start, end, idx in cj]
        dmin = min(d for d, _ in ds)
        best[j] = dmin
        if used[0] == m and alpha + dmin <= alpha0 + 1e-15:
            used = (j, max(idx for d, idx in ds if d <= dmin + 1e-15), dmin)
    nocert = not any_cert
    notrans = any_cert and used[0] == m
    return used[0], used[1], used[2], best, nocert, notrans


# ----------------------------------------------------------------------------------------------- historical quantitative monitor
class QuantitativeMonitor:
    """Historical simultaneous confidence-bound comparator retained for compatibility.

    Public experiment entry points do not select this comparator. The internal
    'VII' method label is a historical identifier, not a current section number.

    spacing=None gives the single initial period. New starts at t = K, 2K, ... retain the cumulative counts and store
    L_{sigma}(c), the lower cumulative-count bound at the period start.
    The loss Y_t is observed at time t+1.
    """

    def __init__(self, candidates, delta, xi=(0.5, 1.0, 2.0, 4.0), spacing=None):
        self.grid = np.asarray(candidates, float)
        self.m = len(self.grid)
        self.xi = np.asarray(xi, float)
        self.spacing = spacing
        self.ell = float(np.log(2 * self.m * len(self.xi) / delta))
        self.a = -np.expm1(-self.xi)
        self.b = np.expm1(self.xi)
        self.t = 0
        self.N = np.zeros(self.m, int)
        self.U, self.L = self.prefix(0, self.N)
        self.starts = [0]
        self.start_L = [self.L.copy()]

    def prefix(self, t, N):
        num = self.xi[:, None] * N[None, :]
        U = np.minimum(t, np.min((num + self.ell) / self.a[:, None], axis=0))
        L = np.maximum(0.0, np.max((num - self.ell) / self.b[:, None], axis=0))
        return U, L

    def observe(self, Y):
        self.N += (Y > self.grid)
        self.t += 1
        self.U, self.L = self.prefix(self.t, self.N)
        if self.spacing is not None and self.t % self.spacing == 0:
            self.starts.append(self.t)
            self.start_L.append(self.L.copy())

    def select(self, D, alpha0):
        """Select a candidate from this historical comparator's period-wise risk bounds.

        Returns (j or m, r or 0, B of the used period, per-candidate best B, nocert (no period), notrans, period start).
        """
        t = self.t
        m = self.m
        if t == 0:
            return m, 0, np.nan, np.full(m, np.inf), True, False, -1
        st = np.array([s for s in self.starts if s < t], int)
        Ls = np.array([L for s, L in zip(self.starts, self.start_L) if s < t])
        n = (t - st).astype(float)
        risk = np.minimum(1.0, (self.U[None, :] - Ls) / n[:, None])                 # U_{t,r}(c), shape (R, m)
        suffix = np.cumsum(D[:, ::-1], axis=1)[:, ::-1]                              # sum_{s >= start} D[j, s]
        avg = (suffix[:, st] / n[None, :]).T                                          # mean D over I_{t,r}(c)
        B = risk + avg
        r_best = np.argmin(B, axis=0)                                                 # smallest index minimising B
        Bbest = B[r_best, np.arange(m)]
        ok = np.flatnonzero(Bbest <= alpha0 + 1e-15)
        if ok.size:
            j = int(ok[0])
            return j, int(r_best[j]) + 1, float(Bbest[j]), Bbest, False, False, int(st[r_best[j]])
        return m, 0, np.nan, Bbest, False, True, -1


# ----------------------------------------------------------------------------------------------- baselines
class ACI:
    """Adaptive conformal inference (eq. empirical-baselines): c_t = min{c_max, [c_hat_t]_+}."""

    def __init__(self, alpha0, c_max, lr, window=500, c_init=None):
        self.alpha0, self.c_max, self.lr, self.window = alpha0, c_max, lr, window
        self.a = alpha0
        self.hist = []
        self.c_hat = c_max if c_init is None else c_init

    def margin(self):
        n = len(self.hist)
        if n == 0:
            self.c_hat = self.c_max
        else:
            k = int(np.ceil((n + 1) * (1 - self.a)))
            if k < 1:
                self.c_hat = -np.inf
            elif k > n:
                self.c_hat = np.inf
            else:
                self.c_hat = float(np.sort(self.hist)[k - 1])
        return float(min(self.c_max, max(self.c_hat, 0.0)))

    def observe(self, Y):
        err = float(Y > self.c_hat)
        self.a = self.a + self.lr * (self.alpha0 - err)
        self.hist.append(float(Y))
        if len(self.hist) > self.window:
            self.hist.pop(0)


class QuantileTracking:
    """Quantile tracking (eq. empirical-baselines): c_hat_{t+1} = c_hat_t + eta (1{Y_t > c_hat_t} - alpha0)."""

    def __init__(self, alpha0, c_max, step, c_init=None):
        self.alpha0, self.c_max, self.step = alpha0, c_max, step
        self.c_hat = c_max if c_init is None else c_init

    def margin(self):
        return float(min(self.c_max, max(self.c_hat, 0.0)))

    def observe(self, Y):
        self.c_hat = self.c_hat + self.step * (float(Y > self.c_hat) - self.alpha0)


@dataclass
class Method:
    """A compared method: statistical rule x transfer, or a baseline.

    fraction phi scales the whole deployment bound, d_t^phi(c) = phi d_t(c), i.e. every pairwise entry D_{t,s}(c)
    including the value 1 of a failed condition (stress test; only phi = 1 is valid). policy is the start schedule
    of the certification tests. None retains the historical restart-at-certification
    policy; public repeated-test experiments pass an explicit periodic policy."""
    label: str
    kind: str                 # 'V', 'VI', 'VII', 'unfiltered', 'zero', 'conservative', 'aci', 'qt'
    transfer: str = "none"    # 'uniform' | 'candidate' | 'commondelta' | 'none'
    spacing: int | None = None
    kappa: float = 1.0        # historical geometric multiplier; not the allocation normalizer or budget normalizer
    restart: bool = True
    fraction: float = 1.0
    policy: dict | None = None


def standard_methods(k_list=(None, 400, 200, 100, 50, 25, 10)):
    """Historical method registry; public experiments use the explicit settings registry."""
    out = [Method("Unfiltered", "unfiltered"), Method("Zero margin", "zero"), Method("Conservative", "conservative"),
           Method("V-U", "V", "uniform"), Method("VI-U", "VI", "uniform"),
           Method("V-C", "V", "candidate"), Method("VI-C", "VI", "candidate"),
           Method("V-Cd", "V", "commondelta")]            # candidate-specific sensitivity, common Delta_{t,s}
    out += [Method(f"VII-C(K={'inf' if k is None else k})", "VII", "candidate", spacing=k) for k in k_list]
    out += [Method("ACI", "aci"), Method("QT", "qt")]
    return out


# ----------------------------------------------------------------------------------------------- public method labels
BASELINES = ("Unfiltered", "Zero margin", "Conservative", "ACI", "QT")
BOUNDS = {"V-U": "uniform", "V-D": "candresid", "V-Cd": "commondelta", "V-C": "candidate"}   # V-D = V-Delta, V-Cd = V-L
SV_MAIN = BASELINES + ("V-U", "V-Cd", "V-C")


def make_method(label: str) -> Method:
    """Method from a label. Baselines: BASELINES. Certified methods: 'V-U' | 'V-D' (V-Delta) | 'V-Cd' (V-L) | 'V-C',
    optionally followed by
    '+once' (certify once, no recertification), '+dK' (next test K steps after each certification; '+d0' is the
    default restart at certification), '+pK' (a new test every K steps, sigma_r(c) = (r - 1) K, unfinished tests
    continue in parallel) and '+fPHI' (deployment bound scaled by PHI), e.g. 'V-C+d25', 'V-U+f0.25',
    'V-C+once+f0.5'."""
    kinds = {"Unfiltered": "unfiltered", "Zero margin": "zero", "Conservative": "conservative", "ACI": "aci",
             "QT": "qt"}
    if label in kinds:
        return Method(label, kinds[label])
    parts = label.split("+")
    if parts[0] not in BOUNDS:
        raise ValueError(f"unknown method label {label!r}")
    policy, fraction = None, 1.0
    for p in parts[1:]:
        if p == "once":
            policy = {"kind": "norestart"}
        elif p.startswith("d") and p[1:].isdigit():
            policy = {"kind": "restart", "delay": int(p[1:])}
        elif p.startswith("p") and p[1:].isdigit():            # common periodic overlapping tests: sigma_r = (r - 1) K
            policy = {"kind": "fixed", "period": int(p[1:]), "overlap": True}
        elif p.startswith("f"):
            fraction = float(p[1:])
            if not 0.0 <= fraction <= 1.0:
                raise ValueError(f"fraction must lie in [0, 1]: {label!r}")
        else:
            raise ValueError(f"unknown label part {p!r} in {label!r}")
    return Method(label, "V", BOUNDS[parts[0]], fraction=fraction, policy=policy)
