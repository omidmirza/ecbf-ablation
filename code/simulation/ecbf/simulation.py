"""Closed-loop unicycle simulation with conditional-risk evaluation."""
from __future__ import annotations

import json
import time

import numpy as np

from . import core
from .geometry import F, L_geom, p0
from .models import make_model
from .risk import exceedance_given_radius, exceedance_risk

L_MU = 5e-4


def atom_level(mode):
    """'atomA' is the stationary distribution ||w|| = A c_max (near-boundary Theorem-1 setting); None for mixture distributions."""
    return float(mode[4:]) if isinstance(mode, str) and mode.startswith("atom") else None


def beta_at(t, mode="sin", L_mu=L_MU):
    """Mixture weight; published runs use (1 + sin(2 L_mu t))/2 with mode='sin'.

    The low-level API also supports a stationary outer component ('const1') and
    a cosine weight ('cos'). Each mode satisfies d_TV(mu_{t+1}, mu_t) <= L_mu.
    The publication entry points always select 'sin'.
    """
    if atom_level(mode) is not None:
        return np.nan
    if mode == "const1":
        return 1.0
    if mode == "cos":                                        # decreasing risk (recovery), same TV bound
        return 0.5 * (1.0 + np.cos(2.0 * L_mu * t))
    return 0.5 * (1.0 + np.sin(2.0 * L_mu * t))


def sample_w(t, c_max, rng, mode="sin", L_mu=L_MU):
    """Sample the planar disturbance from two radial components with uniform direction.

    Published runs use the sinusoidal mixture weight and uniform magnitudes
    within [0, c_max/2] or [c_max/2, c_max]. An explicit atom mode fixes the radius.
    """
    a = atom_level(mode)
    if a is not None:
        A = a * c_max
        rng.uniform()                                        # keep the stream layout of the mixture distributions
    else:
        comp = int(rng.uniform() < beta_at(t, mode, L_mu))
        A = rng.uniform(comp * c_max / 2, (comp + 1) * c_max / 2)
    ph = rng.uniform(0.0, 2 * np.pi)
    return A * np.array([np.cos(ph), np.sin(ph)])


def risk_row(scene, z, cs, mode, beta, c_max):
    """Estimated p_t(c) under the simulator distribution: mixture (exceedance_risk) or atom ||w|| = a (one radius)."""
    a = atom_level(mode)
    if a is None:
        return exceedance_risk(scene, z, cs, beta, c_max)
    return np.array([exceedance_given_radius(scene, np.asarray(z, float), float(c), a * c_max) if c < a * c_max
                     else 0.0 for c in cs])


class Transfer:
    """History of operating points and the pairwise matrices for one trajectory."""

    def __init__(self, model, kappa=1.0, L_mu=L_MU):
        self.m = model
        self.kappa = kappa
        self.L_mu = L_mu
        c_max = model.p.c_max
        self.grid = model.grid
        self.p0 = p0(self.grid, c_max)
        self.Lg = L_geom(self.grid, model.d_min, c_max)
        offs = list(model.scene.finite_offsets()) + [np.inf]
        self.offsets = np.array(offs)
        self.Ftab = np.array([[float(F(R, model.q_min, c, c_max)) for R in offs] for c in self.grid])  # (m, n_off)
        Kx, Ku = model.pairwise_constants()
        self.Lx = float(np.max(self.Lg)) * Kx
        self.Lu = float(np.max(self.Lg)) * Ku
        self.X, self.U, self.psi_s, self.idx_s, self.ok_s = [], [], [], [], []

    def _off_index(self, R):
        return int(np.argmin(np.abs(self.offsets - R))) if np.isfinite(R) else len(self.offsets) - 1

    def append(self, x, u, z):
        sc = self.m.scene
        psi = float(sc.psi(z))
        self.X.append(np.asarray(x, float).copy())
        self.U.append(np.asarray(u, float).copy())
        self.psi_s.append(psi)
        self.idx_s.append(self._off_index(sc.historical_offset(z)))
        self.ok_s.append(psi >= self.m.q_min - 1e-12)

    def current(self, x, cand_u):
        """Candidate conditions and offsets at the current decision."""
        sc, c_max = self.m.scene, self.m.p.c_max
        Zt = np.array([self.m.z(x, u) if np.all(np.isfinite(u)) else np.full(2, np.nan) for u in cand_u])
        ok_t = np.zeros(len(cand_u), bool)
        idx_t = np.zeros(len(cand_u), int)
        psi_t = np.full(len(cand_u), -np.inf)
        for j, zj in enumerate(Zt):
            if not np.all(np.isfinite(zj)):
                continue
            psi_t[j] = float(sc.psi(zj))
            ok, R = sc.current_contribution(zj, c_max)
            ok_t[j] = ok and psi_t[j] >= self.m.q_min - 1e-12
            idx_t[j] = self._off_index(R) if ok else 0
        return ok_t, idx_t, psi_t

    def matrices(self, t, x, cand_u, kind):
        m = len(self.grid)
        if t == 0:
            return np.zeros((m, 0))
        ok_t, idx_t, psi_t = self.current(x, cand_u)
        ok_s = np.array(self.ok_s)
        idx_s = np.array(self.idx_s)
        eps = np.minimum(1.0, self.kappa * self.L_mu * (t - np.arange(t)))
        kap = np.maximum(self.Ftab[np.arange(m)[:, None], idx_t[:, None]] - self.Ftab[:, idx_s], 0.0)   # (m, t)
        if kind == "candidate":
            gap = np.maximum(psi_t[:, None] - np.array(self.psi_s)[None, :], 0.0)
            K = np.minimum(self.p0[:, None], self.kappa * self.Lg[:, None] * gap + kap)
            D = (1 - eps[None, :]) * K + eps[None, :] * self.p0[:, None]
            D[~ok_t, :] = 1.0
            D[:, ~ok_s] = 1.0
        elif kind == "uniform":
            X, U = np.array(self.X), np.array(self.U)
            dx = np.linalg.norm(x[None, :] - X, axis=1)
            du = np.linalg.norm(np.asarray(cand_u)[:, None, :] - U[None, :, :], axis=2)
            if np.all(ok_t):
                Delta = np.minimum(1.0, eps + np.max(kap, axis=0))
                Delta[~ok_s] = 1.0
            else:
                Delta = np.ones(t)
            D = self.kappa * self.Lx * dx[None, :] + self.kappa * self.Lu * du + Delta[None, :]
        elif kind == "candresid":                    # V-Delta: uniform L_x, L_u, candidate residual Delta_{t,s}(c)
            X, U = np.array(self.X), np.array(self.U)    # (eq. unicycle-candidate-delta)
            dx = np.linalg.norm(x[None, :] - X, axis=1)
            du = np.linalg.norm(np.asarray(cand_u)[:, None, :] - U[None, :, :], axis=2)
            Delta = np.minimum(1.0, eps[None, :] + kap)             # (m, t): own check only
            Delta[~ok_t, :] = 1.0
            Delta[:, ~ok_s] = 1.0
            D = self.kappa * self.Lx * dx[None, :] + self.kappa * self.Lu * du + Delta
        elif kind == "commondelta":                  # candidate-specific sensitivity, common Delta_{t,s}
            gap = np.maximum(psi_t[:, None] - np.array(self.psi_s)[None, :], 0.0)
            K = np.minimum(self.p0[:, None], self.kappa * self.Lg[:, None] * gap)
            if np.all(ok_t):
                Delta = np.minimum(1.0, eps + np.max(kap, axis=0))
                Delta[~ok_s] = 1.0
            else:
                Delta = np.ones(t)
            D = K + Delta[None, :]
        elif kind == "zero":
            D = np.zeros((m, t))
        else:
            raise ValueError(kind)
        return np.clip(D, 0.0, 1.0)


def simulate(bench, method: core.Method, seed, alpha=0.10, alpha0=0.15, delta=0.05, kappa=1.0, record_risk=True,
             restart=True, zero_transfer=False, H=None, model_over=None, eps_num=1e-10, beta_mode="sin",
             switch_at=None, L_mu=L_MU, budget_mode="legacy", budget_normalizer=2.11):
    """switch_at = (t_s, j): from decision t_s on, the disturbance stream is replaced by a fresh stream indexed by j
    (continuation from the saved history at t_s, the fixed-certificate conditional assessment).
    budget_mode='joint_log' requires a V method with common periodic overlapping starts (+pK);
    'joint_log_once' uses the first allocated joint-log level with no later starts (+once, K=infinity);
    'single' instead uses full delta for its initial test; 'legacy' preserves historical thresholds.
    budget_normalizer >= 2.11 is separate from the geometric multiplier kappa. A larger
    normalizer conservatively underspends delta; budget_total_upper records the resulting upper bound.
    Oracle risk values are recorded after the decision and never enter the controller.
    """
    if budget_mode not in ("legacy", "single", "joint_log", "joint_log_once"):
        raise ValueError(f"unknown budget mode {budget_mode!r}")
    if budget_mode != "legacy" and not 0 < alpha < alpha0 < 1:
        raise ValueError("require 0 < alpha < alpha0 < 1")
    if budget_mode != "legacy" and method.kind in ("VI", "VII"):
        raise ValueError("single/joint_log budgets are implemented only for Section V")
    model = make_model(bench, **(model_over or {}))
    p = model.p
    H = H or p.H
    grid = model.grid
    m = len(grid)
    c_max = p.c_max
    rng_x, rng_w = [np.random.Generator(np.random.PCG64(s)) for s in np.random.SeedSequence(seed).spawn(2)]
    x = model.initial(rng_x)
    mon = None
    if method.kind in ("V", "VI"):
        mon = core.CertificationMonitor(grid, (0.2, 0.4, 0.6, 0.8, 1.05), alpha, delta, method.kind,
                                        restart=restart and method.restart, policy=method.policy, budget_mode=budget_mode,
                                        budget_normalizer=budget_normalizer)
    elif method.kind == "VII":
        mon = core.QuantitativeMonitor(grid, delta, spacing=method.spacing)
    aci = core.ACI(alpha0, c_max, 0.02 if bench == "corridor" else 0.005, 500) if method.kind == "aci" else None
    qt = core.QuantileTracking(alpha0, c_max, (0.30 if bench == "corridor" else 0.02) * c_max) if method.kind == "qt" else None
    tr = Transfer(model, kappa * method.kappa, L_mu)
    kind = "zero" if zero_transfer else method.transfer
    rec = {k: [] for k in ("u", "u_nom", "c", "j", "r", "Y", "h", "h_next", "cbf", "unsafe", "clearance", "use", "nocert",
                           "notrans", "feasible_cmax", "p_applied", "p_dep", "a2_excess", "a2_true", "residual",
                           "d_used", "cert_start_used", "cert_end_used")}
    X = [x.copy()]
    viiU, viiL, viiN = [], [], []
    if mon is not None and method.kind == "VII":
        viiU.append(mon.U.copy()); viiL.append(mon.L.copy()); viiN.append(mon.N.copy())
    memory = None
    termination = "horizon"
    failed_filter = {"failed_filter_x": np.zeros((0, len(x))), "failed_filter_u_nom": np.zeros((0, 2)),
                     "failed_filter_margins": np.zeros(0), "failed_filter_inputs": np.zeros((0, 2)),
                     "failed_filter_feasible": np.zeros(0, bool)}
    ctrl_time = 0.0
    for t in range(H):
        t0 = time.perf_counter()
        u_nom, memory = model.nominal(x, memory)
        hx = model.h(x)
        nocert = notrans = False
        r = 0
        d_used = np.nan
        cert_start_used = cert_end_used = -1
        if method.kind == "unfiltered":
            j, u = m, u_nom.copy()
            c = np.nan
            feas_cmax = bool(model.filter(x, u_nom, [c_max])[1][0])
            used_D = None
        else:
            if method.kind in ("V", "VI", "VII"):
                inputs, feas = model.filter(x, u_nom, list(grid) + [c_max])
            else:
                if method.kind == "zero":
                    cz = 0.0
                elif method.kind == "conservative":
                    cz = c_max
                elif method.kind == "aci":
                    cz = aci.margin()
                else:
                    cz = qt.margin()
                inputs, feas = model.filter(x, u_nom, [cz, c_max])
            feas_cmax = bool(feas[-1])
            own_ok = feas_cmax if method.kind in ("V", "VI", "VII", "conservative") else bool(feas[0])
            if not own_ok:                                   # adaptive rules need the c_max fallback (Assumption 1)
                termination = "infeasible"
                failed_filter = dict(failed_filter_x=x[None, :].copy(), failed_filter_u_nom=u_nom[None, :].copy(),
                                     failed_filter_margins=np.asarray(list(grid) + [c_max]
                                         if method.kind in ("V", "VI", "VII") else [cz, c_max]),
                                     failed_filter_inputs=np.asarray(inputs).copy(),
                                     failed_filter_feasible=np.asarray(feas, bool).copy())
                ctrl_time += time.perf_counter() - t0
                break
            used_D = None
            if method.kind in ("V", "VI"):
                D = tr.matrices(t, x, inputs[:m], kind)
                if method.fraction != 1.0:                   # stress test: d_t^phi(c) = phi d_t(c)
                    D = method.fraction * D
                j, r, d, _, nocert, notrans = core.select_certified(mon.certs, D, alpha, alpha0, t)
                used_D = D
                d_used = d
                if j < m:
                    selected_cert = next(cc for cc in mon.certs[j] if cc[2] == r)
                    cert_start_used, cert_end_used = selected_cert[:2]
                    d_used = float(np.max(D[j, cert_start_used:cert_end_used]))
            elif method.kind == "VII":
                D = tr.matrices(t, x, inputs[:m], kind)
                j, r, B, _, nocert, notrans, start = mon.select(D, alpha0)
                used_D = D
            else:
                j = 0 if (method.kind == "zero" or (method.kind in ("aci", "qt") and cz < c_max)) else m
            if method.kind in ("V", "VI", "VII"):
                u = inputs[j].copy() if j < m else inputs[m].copy()
                c = grid[j] if j < m else c_max
            else:
                u = inputs[0].copy()
                c = cz
                j = m if cz >= c_max else 0
        ctrl_time += time.perf_counter() - t0
        z = model.z(x, u)
        beta = beta_at(t, beta_mode, L_mu)
        p_row = risk_row(model.scene, z, grid, beta_mode, beta, c_max) if record_risk else np.full(m, np.nan)
        # deployed-margin risk at the applied input (0 on fallback, Assumption conservative-margin)
        if method.kind in ("V", "VI", "VII"):
            p_dep = float(p_row[j]) if j < m else 0.0
        elif method.kind == "unfiltered":
            p_dep = np.nan
        else:
            p_dep = float(risk_row(model.scene, z, [c], beta_mode, beta, c_max)[0]) if (record_risk and c < c_max) else 0.0
        a2 = a2t = np.nan
        if record_risk and not zero_transfer and method.kind in ("V", "VI", "VII") and j < m and used_D is not None:
            if method.kind == "VII":
                s0, s1 = start, t
            else:
                cert = next(cc for cc in mon.certs[j] if cc[2] == r)
                s0, s1 = cert[0], cert[1]
            past = np.array(rec["p_applied"])[s0:s1, j]
            a2 = float(np.max(p_row[j] - past - used_D[j, s0:s1]))
            a2t = float(p_row[j] - np.min(past) - np.max(used_D[j, s0:s1]))   # Assumption 2 itself (record only)
        if switch_at is not None and t == switch_at[0]:
            rng_w = np.random.Generator(np.random.PCG64(np.random.SeedSequence([int(seed), int(switch_at[1]), 7919])))
        residual = (model.h(model.predict(x, u)) - c - (1 - p.gamma) * hx) if not np.isnan(c) else np.nan
        w = sample_w(t, c_max, rng_w, beta_mode, L_mu)
        xn = model.step(x, u, w)
        Y = max(0.0, float(model.scene.psi(z)) - float(model.scene.psi(z + w)))
        if Y > c_max + 1e-12 or not np.all(np.isfinite(u)):
            termination = "numerical"                        # numerical failure is a recorded outcome
            break
        hn = model.h(xn)
        clear = model.body_clearance(xn)
        for k, v in (("u", u), ("u_nom", u_nom), ("c", c), ("j", j), ("r", r), ("Y", Y), ("h", hx), ("h_next", hn),
                     ("cbf", hn < (1 - p.gamma) * hx - eps_num), ("unsafe", hn < 0), ("clearance", clear),
                     ("use", (not np.isnan(c)) and c < c_max), ("nocert", nocert), ("notrans", notrans),
                     ("feasible_cmax", feas_cmax), ("p_applied", p_row), ("p_dep", p_dep), ("a2_excess", a2), ("a2_true", a2t),
                     ("residual", residual), ("d_used", d_used), ("cert_start_used", cert_start_used),
                     ("cert_end_used", cert_end_used)):
            rec[k].append(v)
        t1 = time.perf_counter()
        tr.append(x, u, z)
        if mon is not None:
            mon.observe(Y)
            if method.kind == "VII":
                viiU.append(mon.U.copy()); viiL.append(mon.L.copy()); viiN.append(mon.N.copy())
        if aci is not None:
            aci.observe(Y)
        if qt is not None:
            qt.observe(Y)
        ctrl_time += time.perf_counter() - t1
        x = xn
        X.append(x.copy())
        if clear <= 0.0:
            termination = "contact"
            break
        if model.at_goal(x):
            termination = "goal"
            break
    N = len(rec["Y"])
    out = {k: np.asarray(v) for k, v in rec.items()}
    out["p_applied"] = np.asarray(rec["p_applied"]).reshape(N, m)
    out["u"] = np.asarray(rec["u"]).reshape(N, 2)
    out["u_nom"] = np.asarray(rec["u_nom"]).reshape(N, 2)
    out["x"] = np.asarray(X)
    out.update(failed_filter)
    certificates = []
    if method.kind in ("V", "VI"):
        # Explicit shared order for assessment rows and their allocated budgets/evidence.
        certificates = sorted(mon.core.all_certificates(), key=lambda cc: (cc.candidate, cc.index))
    out["certs"] = np.array([(cc.candidate, cc.index, cc.start, cc.end) for cc in certificates], int).reshape(-1, 4)
    out["cert_budget_delta"] = np.array([cc.delta_r for cc in certificates], float)
    out["cert_log_threshold"] = np.array([cc.log_threshold for cc in certificates], float)
    out["cert_log_evidence"] = np.array([cc.log_evidence for cc in certificates], float)
    if method.kind == "VII":
        out["vii_U"], out["vii_L"], out["vii_N"] = np.array(viiU), np.array(viiL), np.array(viiN)
    out.update(benchmark=bench, method=method.label, seed=int(seed), T=p.T, H=H, N=N, goal=termination == "goal",
               termination=termination, c_max=c_max, alpha=alpha, alpha0=alpha0, delta=delta, kappa=kappa,
               fraction=float(method.fraction), policy=json.dumps(mon.policy_spec if method.kind in ("V", "VI") else None),
               budget_mode=budget_mode if method.kind in ("V", "VI") else "not_applicable",
               requested_budget_mode=budget_mode, total_delta=float(delta),
               budget_normalizer=mon.core.budget_normalizer if method.kind == "V" and budget_mode in ("joint_log", "joint_log_once") else np.nan,
               budget_total_upper=mon.core.budget_total_upper if method.kind in ("V", "VI") else np.nan,
               model_over=json.dumps(model_over or {}, sort_keys=True),
               beta_mode=str(beta_mode), rho_buf=float(p.rho_buf), L_mu=float(L_mu),
               runtime_ms=1000 * ctrl_time / max(N, 1), grid=list(grid))
    return out
