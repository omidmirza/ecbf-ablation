"""Test-start policies used by the certification monitor."""
from __future__ import annotations

import numpy as np


class StartPolicy:
    name = "base"

    def reset(self, monitor) -> None:
        pass

    def step(self, monitor, new_certificates, context):
        return [], []

    def describe(self) -> dict:
        return dict(name=self.name)


class NoRestart(StartPolicy):
    """One test per candidate, started at time 0 (K = infinity)."""
    name = "norestart"


class RestartOnCertification(StartPolicy):
    """sigma_{r+1}(c) = tau_r(c) + delay (budget_mode='legacy' only; not used in the study)."""

    def __init__(self, delay: int = 0):
        if delay < 0:
            raise ValueError("delay must be nonnegative")
        self.delay = int(delay)
        self.name = "restart" if self.delay == 0 else f"delay{self.delay}"

    def reset(self, monitor) -> None:
        self.pending = {}

    def step(self, monitor, new_certificates, context):
        for cert in new_certificates:
            self.pending[cert.candidate] = cert.end + self.delay
        due = [j for j, when in self.pending.items() if when <= monitor.t and monitor.active_rows(j).size == 0]
        for j in due:
            del self.pending[j]
        return due, []

    def describe(self) -> dict:
        return dict(name=self.name, delay=self.delay)


class FixedPeriod(StartPolicy):
    """A new test for every candidate every `period` steps, sigma_r = (r - 1) K. The study uses overlap=True
    (unfinished tests continue in parallel), no restart at certification and no cap on the active tests."""

    def __init__(self, period: int, overlap: bool = False, restart_on_certification: bool = False,
                 max_active: int | None = None):
        if period < 1:
            raise ValueError("period must be at least one step")
        self.period = int(period)
        self.overlap = bool(overlap)
        self.restart_on_certification = bool(restart_on_certification)
        self.max_active = max_active
        tag = "overlap" if overlap else "block"
        self.name = f"fixed{self.period}_{tag}" + ("_restart" if restart_on_certification else "")

    def step(self, monitor, new_certificates, context):
        starts, discards = [], []
        t = monitor.t
        if self.restart_on_certification:
            starts += [c.candidate for c in new_certificates]
        if t % self.period == 0:
            for j in range(monitor.m):
                rows = monitor.active_rows(j)
                if not self.overlap:
                    discards += rows.tolist()
                elif self.max_active is not None and rows.size >= self.max_active:
                    oldest = rows[np.argsort(monitor.start[rows])][: rows.size - self.max_active + 1]
                    discards += oldest.tolist()
                starts.append(j)
        return starts, discards

    def describe(self) -> dict:
        return dict(name=self.name, period=self.period, overlap=self.overlap,
                    restart_on_certification=self.restart_on_certification, max_active=self.max_active)


class EveryStep(StartPolicy):
    """A new test at every step (budget_mode='legacy' only; not used in the study)."""

    def __init__(self, window: int | None = None):
        self.window = window
        self.name = "everystep" + (f"_w{window}" if window else "")

    def step(self, monitor, new_certificates, context):
        discards = []
        if self.window is not None and monitor.n_active:
            discards = np.flatnonzero(monitor.t - monitor.start >= self.window).tolist()
        return list(range(monitor.m)), discards

    def describe(self) -> dict:
        return dict(name=self.name, window=self.window)


class StopLoss(StartPolicy):
    """Restart a test whose evidence falls below `floor` (budget_mode='legacy' only; not used in the study)."""

    def __init__(self, floor: float = 0.05, delay: int = 0):
        if not 0 < floor < 1:
            raise ValueError("floor must lie in (0, 1)")
        self.floor = float(floor)
        self.base = RestartOnCertification(delay)
        self.name = f"stoploss{floor:g}" + (f"_d{delay}" if delay else "")

    def reset(self, monitor) -> None:
        self.base.reset(monitor)

    def step(self, monitor, new_certificates, context):
        starts, _ = self.base.step(monitor, new_certificates, context)
        discards = []
        if monitor.n_active:
            low = np.flatnonzero(monitor.log_E() < np.log(self.floor))
            discards = low.tolist()
            starts = sorted(set(starts) | set(int(monitor.cand[r]) for r in low))
        return starts, discards

    def describe(self) -> dict:
        return dict(name=self.name, floor=self.floor, delay=self.base.delay)


class OnDemand(StartPolicy):
    """Start a test for a certified candidate that is not deployable (budget_mode='legacy' only; not used in the
    study). Needs context['deployable'] (bool per candidate, from the previous decision)."""

    def __init__(self, min_gap: int = 20, max_active: int = 2, delay: int = 0):
        self.min_gap, self.max_active = int(min_gap), int(max_active)
        self.base = RestartOnCertification(delay)
        self.name = f"ondemand{self.min_gap}"

    def reset(self, monitor) -> None:
        self.base.reset(monitor)

    def step(self, monitor, new_certificates, context):
        starts, _ = self.base.step(monitor, new_certificates, context)
        discards = []
        deployable = np.asarray(context.get("deployable", np.ones(monitor.m, bool)), dtype=bool)
        certified = np.asarray([len(c) > 0 for c in monitor.certificates])
        for j in np.flatnonzero(certified & ~deployable):
            if monitor.t - monitor.last_start[j] < self.min_gap or j in starts:
                continue
            rows = monitor.active_rows(int(j))
            if rows.size >= self.max_active:
                discards += rows[np.argsort(monitor.start[rows])][: rows.size - self.max_active + 1].tolist()
            starts.append(int(j))
        return starts, discards

    def describe(self) -> dict:
        return dict(name=self.name, min_gap=self.min_gap, max_active=self.max_active, delay=self.base.delay)


def make_policy(spec: dict) -> StartPolicy:
    """Build a policy from a JSON-serializable spec, e.g. {"kind": "fixed", "period": 50, "overlap": false}."""
    spec = dict(spec)
    kind = spec.pop("kind")
    table = {"norestart": NoRestart, "restart": RestartOnCertification, "fixed": FixedPeriod, "everystep": EveryStep,
             "stoploss": StopLoss, "ondemand": OnDemand}
    if kind not in table:
        raise ValueError(f"unknown policy kind {kind!r}")
    return table[kind](**spec)
