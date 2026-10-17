"""Sequential tests with periodic starts and allocated error budgets."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .evidence import GainSet
from .budgets import BUDGET_MODES, JOINT_LOG_NORMALIZER, joint_log_delta
from .policies import FixedPeriod, NoRestart


@dataclass(frozen=True)
class RunningMin:
    """Running minimum over a test's monitoring period of `context[name]`.

    per_candidate=False: the context value has shape `shape` and is the same for every test (for example the lower
    component risk bounds of all candidates). per_candidate=True: the context value has shape (m,) + shape and each
    test keeps the entry of its own candidate.
    """
    name: str
    shape: tuple = ()
    per_candidate: bool = False


@dataclass
class Certificate:
    candidate: int          # index j of c^(j+1) in the grid
    index: int              # test index r (order of start)
    start: int              # sigma_r(c)
    end: int                # tau_r(c); monitoring period I_r(c) = [start, end)
    log_correction: np.ndarray   # log A_sigma^lambda per gain
    log_evidence: float     # log V_{tau,r}(c) at certification
    stats: dict = field(default_factory=dict)
    log_w: np.ndarray | None = None   # log E_{tau,r}^lambda(c) per gain at certification
    delta_r: float | None = None      # immutable allocated level of this test
    log_threshold: float | None = None  # log(1 / delta_r), stored when the test starts


class RepeatedCertificationMonitor:
    def __init__(self, n_candidates: int, gains, alpha: float, delta: float, policy, correction=True,
                 trackers=(), delta_weights=None, budget_mode="legacy", budget_normalizer=JOINT_LOG_NORMALIZER):
        if n_candidates < 1 or not 0 < delta < 1:
            raise ValueError("need at least one candidate and 0 < delta < 1")
        self.m = int(n_candidates)
        self.gs = GainSet(gains, alpha)
        self.alpha, self.delta = float(alpha), float(delta)
        self.log_threshold = float(np.log(1 / delta))
        self.policy = policy
        mode = {False: "none", True: "shared", 0: "none", 1: "shared"}.get(correction, correction)
        if mode not in ("none", "shared", "candidate"):
            raise ValueError(f"unknown correction mode {correction!r}")
        if budget_mode not in BUDGET_MODES:
            raise ValueError(f"unknown budget mode {budget_mode!r}")
        if budget_mode != "legacy":
            if mode != "none" or delta_weights is not None:
                raise ValueError("single/joint_log budgets require no historical or candidate correction")
            if budget_mode in ("single", "joint_log_once") and type(policy) is not NoRestart:
                raise ValueError("single/joint_log_once budget requires NoRestart")
            if budget_mode == "joint_log" and not (
                    type(policy) is FixedPeriod and policy.overlap
                    and not policy.restart_on_certification and policy.max_active is None):
                raise ValueError("joint_log requires FixedPeriod(overlap=True, "
                                 "restart_on_certification=False, max_active=None)")
        self.budget_mode = budget_mode
        self.budget_normalizer = float(budget_normalizer) if budget_mode in ("joint_log", "joint_log_once") else None
        if budget_mode in ("joint_log", "joint_log_once"):
            joint_log_delta(self.delta, 1, self.budget_normalizer)  # validates the conservative normalizer
        self.budget_total_upper = (joint_log_delta(self.delta, 1, self.budget_normalizer)
                                   if budget_mode == "joint_log_once" else self.delta * JOINT_LOG_NORMALIZER / self.budget_normalizer
                                   if budget_mode == "joint_log" else self.delta if budget_mode == "single" else np.nan)
        self.mode = mode
        self.correction = mode != "none"
        w = np.full(self.m, 1.0 / self.m) if delta_weights is None else np.asarray(delta_weights, float)
        if mode == "candidate":
            if w.shape != (self.m,) or np.any(w <= 0) or abs(w.sum() - 1) > 1e-12:
                raise ValueError("delta_weights must be positive and sum to one")
            self.delta_c = delta * w
            self.log_thr_c = np.log(1 / self.delta_c)
        else:
            self.delta_c = np.full(self.m, self.delta)
            self.log_thr_c = np.full(self.m, self.log_threshold)
        self.trackers = tuple(trackers)
        self.t = 0
        self.log_A = np.zeros(self.gs.size)
        self.log_A_c = np.zeros((self.m, self.gs.size))
        G = self.gs.size
        self.cand = np.zeros(0, dtype=int)
        self.idx = np.zeros(0, dtype=int)
        self.start = np.zeros(0, dtype=int)
        self.logw = np.zeros((0, G))
        self.loga = np.zeros((0, G))
        self.test_delta = np.zeros(0)
        self.test_log_threshold = np.zeros(0)
        self.stats = {tr.name: np.zeros((0,) + tuple(tr.shape)) for tr in self.trackers}
        self.n_started = np.zeros(self.m, dtype=int)
        self.last_start = np.full(self.m, -1, dtype=int)
        self.certificates = [[] for _ in range(self.m)]
        self.n_discarded = 0
        self.start_tests(range(self.m))           # sigma_1(c) = 0 for every candidate
        self.policy.reset(self)

    # ------------------------------------------------------------------ test bookkeeping
    @property
    def n_active(self) -> int:
        return int(self.cand.size)

    def active_rows(self, candidate: int) -> np.ndarray:
        return np.flatnonzero(self.cand == candidate)

    def start_tests(self, candidates) -> list:
        """Start one new test at the current time for each listed candidate (at most one start per time)."""
        candidates = sorted(set(int(c) for c in candidates))
        if self.budget_mode != "legacy":
            if candidates != list(range(self.m)):
                raise ValueError("single/joint_log tests must start for every candidate together")
            if self.budget_mode in ("single", "joint_log_once") and np.any(self.n_started):
                if self.t == 0:
                    return []
                raise ValueError("single/joint_log_once budget permits only the initial test")
            if self.budget_mode == "joint_log":
                if self.t % self.policy.period:
                    raise ValueError("joint_log tests start only at multiples of K")
                if np.all(self.last_start == self.t):
                    return []
                expected_previous = self.t // self.policy.period
                if np.any(self.n_started != expected_previous):
                    raise ValueError("joint_log start indices must be common and cannot skip or reset")
        started = []
        for j in candidates:
            if not 0 <= j < self.m:
                raise ValueError(f"candidate index {j} out of range")
            if self.last_start[j] == self.t:
                continue                          # finite start times are strictly increasing
            self.n_started[j] += 1
            self.last_start[j] = self.t
            started.append(j)
        if not started:
            return []
        k = len(started)
        G = self.gs.size
        self.cand = np.r_[self.cand, started]
        self.idx = np.r_[self.idx, self.n_started[started]]
        self.start = np.r_[self.start, np.full(k, self.t)]
        if self.budget_mode in ("joint_log", "joint_log_once"):
            levels = np.array([joint_log_delta(self.delta, int(r), self.budget_normalizer) for r in self.n_started[started]])
            thresholds = np.log(1 / levels)
        else:
            thresholds = self.log_thr_c[started].copy()
            levels = self.delta_c[started].copy()
        self.test_delta = np.r_[self.test_delta, levels]
        self.test_log_threshold = np.r_[self.test_log_threshold, thresholds]
        self.logw = np.vstack([self.logw, np.zeros((k, G))])
        if self.mode == "shared":
            a = np.tile(self.log_A, (k, 1))
        elif self.mode == "candidate":
            a = self.log_A_c[started].copy()
        else:
            a = np.zeros((k, G))
        self.loga = np.vstack([self.loga, a])
        for tr in self.trackers:
            self.stats[tr.name] = np.concatenate([self.stats[tr.name], np.full((k,) + tuple(tr.shape), np.inf)])
        return started

    def discard(self, rows) -> None:
        """Discard active tests without certification (allowed: it can only remove certification events)."""
        rows = np.unique(np.asarray(list(rows), dtype=int))
        if rows.size == 0:
            return
        if self.budget_mode in ("joint_log", "joint_log_once"):
            raise ValueError("joint_log/joint_log_once retains every unfinished test")
        keep = np.ones(self.n_active, dtype=bool)
        keep[rows] = False
        self._keep(keep)
        self.n_discarded += int(rows.size)

    def _keep(self, keep) -> None:
        self.cand, self.idx, self.start = self.cand[keep], self.idx[keep], self.start[keep]
        self.logw, self.loga = self.logw[keep], self.loga[keep]
        self.test_delta = self.test_delta[keep]
        self.test_log_threshold = self.test_log_threshold[keep]
        for name in self.stats:
            self.stats[name] = self.stats[name][keep]

    # ------------------------------------------------------------------ evidence
    def log_V(self) -> np.ndarray:
        """log V_{t,r}(c) of every active test (eq. ct-score)."""
        if self.n_active == 0:
            return np.zeros(0)
        return self.gs.log_mixture(self.loga + self.logw)

    def log_E(self) -> np.ndarray:
        """log E_{t,r}(c) of every active test, without the correction."""
        if self.n_active == 0:
            return np.zeros(0)
        return self.gs.log_mixture(self.logw)

    def candidate_log_V(self, reducer=np.max) -> np.ndarray:
        """Per candidate, a summary (default: max) of log V over its active tests; -inf without an active test."""
        out = np.full(self.m, -np.inf)
        lv = self.log_V()
        for j in range(self.m):
            rows = self.cand == j
            if np.any(rows):
                out[j] = float(reducer(lv[rows]))
        return out

    # ------------------------------------------------------------------ one step
    def observe(self, z, context=None, policy_context=None) -> list:
        """Observe exceedances z_j = 1{Y_t > c^(j+1)} of the loss Y_t (at decision time t+1). Returns new certificates."""
        z = np.asarray(z, dtype=bool).ravel()
        if z.size != self.m:
            raise ValueError("need one exceedance indicator per candidate")
        if np.any(np.diff(z.astype(int)) > 0):
            raise ValueError("exceedance indicators must be nonincreasing along the ordered grid")
        context = context or {}
        t = self.t
        if self.n_active:
            self.logw += self.gs.log_factor(z[self.cand])
            for tr in self.trackers:
                value = np.asarray(context[tr.name], dtype=float)
                value = value[self.cand] if tr.per_candidate else np.broadcast_to(value, (self.n_active,) + tuple(tr.shape))
                self.stats[tr.name] = np.minimum(self.stats[tr.name], value)
        if self.mode == "shared":
            self.log_A = self.log_A + self.gs.log_correction_step(z[0])
        elif self.mode == "candidate":
            self.log_A_c = self.log_A_c + np.where(z[:, None], self.gs.log_down[None, :], 0.0)
        new = []
        if self.n_active:
            lv = self.log_V()
            crossed = lv >= self.test_log_threshold
            order = np.flatnonzero(crossed)
            order = order[np.lexsort((self.idx[order], self.cand[order]))]
            for row in order:
                cert = Certificate(candidate=int(self.cand[row]), index=int(self.idx[row]), start=int(self.start[row]),
                                   end=t + 1, log_correction=self.loga[row].copy(), log_evidence=float(lv[row]),
                                   stats={name: self.stats[name][row].copy() for name in self.stats},
                                   log_w=self.logw[row].copy(), delta_r=float(self.test_delta[row]),
                                   log_threshold=float(self.test_log_threshold[row]))
                self.certificates[cert.candidate].append(cert)
                new.append(cert)
            if order.size:
                self._keep(~crossed)
        self.t = t + 1
        starts, discards = self.policy.step(self, new, policy_context or {})
        if len(discards):
            self.discard(discards)
        if len(starts):
            self.start_tests(starts)
        return new

    # ------------------------------------------------------------------ summaries
    def all_certificates(self) -> list:
        return [c for certs in self.certificates for c in certs]

    def summary(self) -> dict:
        certs = self.all_certificates()
        return dict(steps=self.t, correction=self.mode, tests_started=int(self.n_started.sum()), tests_discarded=self.n_discarded,
                    active=self.n_active, certificates=len(certs),
                    certificates_per_candidate=[len(c) for c in self.certificates],
                    log_A=self.log_A.tolist(), budget_mode=self.budget_mode, total_delta=self.delta,
                    budget_normalizer=self.budget_normalizer, budget_total_upper=self.budget_total_upper,
                    policy=self.policy.describe())
