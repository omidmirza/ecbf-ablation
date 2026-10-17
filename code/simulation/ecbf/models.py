"""Corridor and Maze dynamics, nominal controllers, and safety filters."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import L_geom, corridor_scene, maze_scene

_GLX, _GLW = np.polynomial.legendre.leggauss(12)


def e(th):
    th = np.asarray(th, float)
    return np.stack([np.cos(th), np.sin(th)], axis=-1)


def wrap(a):
    return (np.asarray(a) + np.pi) % (2 * np.pi) - np.pi


def sinc(r):
    r = np.asarray(r, float)
    out = np.ones_like(r)
    k = np.abs(r) > 1e-12
    out[k] = np.sin(r[k]) / r[k]
    return out


# ================================================================================================= corridor
@dataclass
class CorridorParams:
    T: float = 0.1
    H: int = 700
    l: float = 0.15
    rho_body: float = 0.22
    rho_buf: float = 0.05
    gamma: float = 0.25
    v_max: float = 0.6
    omega_max: float = 2.0
    v_cruise: float = 0.5
    k_theta: float = 2.0
    omega_weight: float = 0.0225        # l^2: objective (v - v_nom)^2 + l^2 (omega - omega_nom)^2
    omega_nodes: int = 41
    c_max: float = 0.02
    fractions: tuple = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


class Corridor:
    name = "corridor"
    d = 3

    def __init__(self, p: CorridorParams = None):
        self.p = p or CorridorParams()
        self.scene = corridor_scene(self.p.rho_body)
        self.grid = self.p.c_max * np.asarray(self.p.fractions)
        self.q_min = self.p.l + self.p.rho_buf
        self.d_min = self.q_min + float(self.scene.finite_offsets().min())
        self.omega_grid = np.linspace(-self.p.omega_max, self.p.omega_max, self.p.omega_nodes)
        vb, T, l = self.p.v_max, self.p.T, self.p.l
        self.Kx = float(np.sqrt(1 + (l + vb * T) ** 2))                       # eq. unicycle-first-order-constants
        self.Ku = float(np.sqrt(T ** 2 + (l * T + vb * T ** 2 / 2) ** 2))
        self.goal = np.array(self.scene.meta["goal"], float)
        self.goal_tol = float(self.scene.meta["goal_tol"])

    # ------------------------------------------------------------------------------------------------ dynamics
    def nose(self, x):
        x = np.asarray(x, float)
        return x[..., :2] + self.p.l * e(x[..., 2])

    def predict(self, x, u):
        v, om = u[..., 0], u[..., 1]
        T = self.p.T
        y = x[..., :2] + (T * v * sinc(T * om / 2))[..., None] * e(x[..., 2] + T * om / 2)
        return np.concatenate([y, (x[..., 2] + T * om)[..., None]], axis=-1)

    def z(self, x, u):
        return self.nose(self.predict(x, u))

    def step(self, x, u, w):
        xn = self.predict(x, u)
        xn[:2] = xn[:2] + w
        return xn

    def h(self, x):
        return float(self.scene.psi(self.nose(x))) - self.q_min

    def body_clearance(self, x):
        return float(self.scene.psi(np.asarray(x)[:2]))

    # ------------------------------------------------------------------------------------------------ nominal
    def initial(self, rng):
        x0, x1, y0, y1 = self.scene.meta["start_box"]
        hd = np.deg2rad(self.scene.meta["heading_deg"])
        return np.array([rng.uniform(x0, x1), rng.uniform(y0, y1), rng.uniform(-hd, hd)])

    def nominal(self, x, memory):
        """Goal tracking: heading-error yaw rate, speed reduced when misaligned."""
        dv = self.goal - x[:2]
        err = float(wrap(np.arctan2(dv[1], dv[0]) - x[2]))
        om = float(np.clip(self.p.k_theta * err, -self.p.omega_max, self.p.omega_max))
        v = float(self.p.v_cruise * max(np.cos(err), 0.0))
        return np.array([v, om]), memory

    def at_goal(self, x):
        return float(np.linalg.norm(x[:2] - self.goal)) <= self.goal_tol

    # ------------------------------------------------------------------------------------------------ filter
    def filter(self, x, u_nom, margins):
        """Inputs u(c) for each margin (nan rows when infeasible) and feasibility flags (vectorised over the yaw-rate
        grid; for each omega the admissible speeds are [lo, hi] minus the open root intervals of the cones)."""
        p, sc = self.p, self.scene
        T, l = p.T, p.l
        hx = self.h(x)
        th = x[2]
        oms = np.unique(np.concatenate([self.omega_grid, [np.clip(u_nom[1], -p.omega_max, p.omega_max)]]))
        s = (T * sinc(T * oms / 2))[:, None] * e(th + T * oms / 2)               # dz/dv, (n, 2)
        p0 = x[:2] + l * e(th + T * oms)                                           # z at v = 0, (n, 2)
        a2 = np.sum(s * s, axis=1)                                                 # (n,)
        dp = p0[:, None, :] - sc.circles[None, :, :2]                              # (n, nc, 2)
        b1 = 2 * np.einsum("nd,ncd->nc", s, dp)
        dd = np.sum(dp * dp, axis=2)
        ns = s @ sc.walls[:, :2].T                                                 # (n, nw)
        npz = p0 @ sc.walls[:, :2].T
        out = np.full((len(margins), 2), np.nan)
        feas = np.zeros(len(margins), bool)
        v0 = u_nom[0]
        for k, c in enumerate(margins):
            thr = (1 - p.gamma) * hx + c + self.q_min
            r = thr + sc.circles[:, 2] + sc.rho_body                               # (nc,)
            disc = b1 * b1 - 4 * a2[:, None] * (dd - r[None, :] ** 2)
            has = disc > 0
            sq = np.sqrt(np.where(has, disc, 0.0))
            f_lo = np.where(has, (-b1 - sq) / (2 * a2[:, None]), np.inf)          # forbidden (f_lo, f_hi)
            f_hi = np.where(has, (-b1 + sq) / (2 * a2[:, None]), -np.inf)
            rhs = thr + sc.walls[:, 2] + sc.rho_body - npz                         # need v * ns >= rhs
            lo = np.full(len(oms), -p.v_max)
            hi = np.full(len(oms), p.v_max)
            pos, neg, zer = ns > 1e-14, ns < -1e-14, np.abs(ns) <= 1e-14
            lo = np.maximum(lo, np.max(np.where(pos, rhs / np.where(pos, ns, 1), -np.inf), axis=1))
            hi = np.minimum(hi, np.min(np.where(neg, rhs / np.where(neg, ns, 1), np.inf), axis=1))
            dead = np.any(zer & (rhs > 0), axis=1) | (lo > hi)
            cand = np.concatenate([np.clip(v0, lo, hi)[:, None], lo[:, None], hi[:, None], f_lo, f_hi], axis=1)
            inside = (cand >= lo[:, None] - 1e-15) & (cand <= hi[:, None] + 1e-15)
            blocked = np.any((cand[:, :, None] > f_lo[:, None, :]) & (cand[:, :, None] < f_hi[:, None, :]), axis=2)
            okc = inside & ~blocked & np.isfinite(cand) & ~dead[:, None]
            dist = np.where(okc, np.abs(cand - v0), np.inf)
            bi = np.argmin(dist, axis=1)
            dbest = dist[np.arange(len(oms)), bi]
            J = dbest ** 2 + p.omega_weight * (oms - u_nom[1]) ** 2
            if np.isfinite(dbest).any():
                i = int(np.argmin(np.where(np.isfinite(dbest), J, np.inf)))
                out[k] = (cand[i, bi[i]], oms[i])
                feas[k] = True
        return out, feas

    def pairwise_constants(self):
        return self.Kx, self.Ku


def _closest_admissible(v0, lo, hi, forbidden):
    """Closest point to v0 in [lo, hi] minus the union of open intervals `forbidden`; None if empty."""
    cands = [min(max(v0, lo), hi)]
    for a, b in forbidden:
        for v in (a, b):
            if lo <= v <= hi:
                cands.append(v)
    cands += [lo, hi]
    best = None
    for v in sorted(set(cands), key=lambda v: abs(v - v0)):
        if all(not (a < v < b) for a, b in forbidden):
            best = v
            break
    return best


# ================================================================================================= maze
@dataclass
class MazeParams:
    T: float = 0.2
    H: int = 600
    l: float = 0.15
    rho_body: float = 0.22
    rho_buf: float = 0.05
    gamma: float = 0.25
    a_br: float = 1.2
    a_max: float = 1.5
    nu_max: float = 8.0
    v_max: float = 0.8
    v_rev: float = 0.3
    omega_max: float = 1.5
    v_cruise: float = 0.6
    k_theta: float = 2.5
    k_omega: float = 4.0
    k_v: float = 2.0
    lookahead: float = 0.5
    nu_weight: float = 0.0225           # l^2
    goal_y: float = 8.25                # the scene's y = 8.5 is unreachable for the centre with a nose point
    a_nodes: int = 41
    nu_nodes: int = 25
    c_max: float = 0.01
    fractions: tuple = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


class Maze:
    name = "maze"
    d = 5

    def __init__(self, p: MazeParams = None):
        self.p = p or MazeParams()
        self.scene = maze_scene(self.p.rho_body)
        self.grid = self.p.c_max * np.asarray(self.p.fractions)
        self.q_min = self.p.l + self.p.rho_buf
        self.d_min = self.q_min + float(self.scene.finite_offsets().min())
        T, l, vb = self.p.T, self.p.l, self.p.v_max
        self.Kx = float(np.sqrt(1 + (l + vb * T) ** 2 + T ** 2 + (l * T + vb * T ** 2 / 2) ** 2))  # eq. (unicycle-second-order-constants) of the paper
        self.Ku = float(np.sqrt(T ** 4 / 4 + (l * T ** 2 / 2 + vb * T ** 3 / 6) ** 2))
        self.tau = 0.5 * T * (_GLX + 1)
        self.wtau = 0.5 * T * _GLW
        ref = np.load(__file__.replace("models.py", "maze_reference_v0_4_2.npy"))
        self.path = ref[:, 1:3]
        seg = np.linalg.norm(np.diff(self.path, axis=0), axis=1)
        self.path_s = np.concatenate([[0.0], np.cumsum(seg)])
        self.goal = np.array([self.scene.meta["goal"][0], self.p.goal_y], float)
        self.goal_tol = float(self.scene.meta["goal_tol"])
        self.a_grid = np.linspace(-self.p.a_max, self.p.a_max, self.p.a_nodes)
        self.nu_grid = np.linspace(-self.p.nu_max, self.p.nu_max, self.p.nu_nodes)

    def nose(self, x):
        x = np.asarray(x, float)
        return x[..., :2] + self.p.l * e(x[..., 2])

    def _integrals(self, th, om, nu):
        """I0 = int_0^T e(vartheta), I1 = int_0^T tau e(vartheta), vartheta(tau) = th + om tau + nu tau^2/2."""
        nu = np.atleast_1d(nu)
        ph = th + om * self.tau[None, :] + 0.5 * nu[:, None] * self.tau[None, :] ** 2
        E = e(ph)                                                                  # (n, 12, 2)
        I0 = np.einsum("k,nkd->nd", self.wtau, E)
        I1 = np.einsum("k,nkd->nd", self.wtau * self.tau, E)
        return I0, I1

    def predict(self, x, u):
        return self.predict_state(x, u)

    def _predict_core(self, x, a, nu, I0, I1):
        T = self.p.T
        y, th, v, om = x[:2], x[2], x[3], x[4]
        yn = y + v * I0 + a * I1
        return np.array([yn[0], yn[1], th + om * T + nu * T * T / 2, v + a * T, om + nu * T])

    def z(self, x, u):
        return self.nose(self.predict_state(x, u))

    def predict_state(self, x, u):
        a, nu = float(u[0]), float(u[1])
        I0, I1 = self._integrals(x[2], x[4], nu)
        return self._predict_core(x, a, nu, I0[0], I1[0])

    def step(self, x, u, w):
        xn = self.predict_state(x, u)
        xn[:2] = xn[:2] + w
        return xn

    def h(self, x):
        return float(self.scene.psi(self.nose(x))) - self.q_min - x[3] ** 2 / (2 * self.p.a_br)

    def body_clearance(self, x):
        return float(self.scene.psi(np.asarray(x)[:2]))

    def initial(self, rng):
        sx, sy = self.scene.meta["start"]
        return np.array([sx, sy + rng.uniform(-0.08, 0.08), np.deg2rad(self.scene.meta["heading_deg"]), 0.0, 0.0])

    def nominal(self, x, memory):
        """Path tracking: pure pursuit on the reference path, monotone progress index in `memory`."""
        k0 = 0 if memory is None else memory
        win = slice(k0, min(k0 + 400, len(self.path)))
        k = k0 + int(np.argmin(np.linalg.norm(self.path[win] - x[:2], axis=1)))
        s_t = self.path_s[k] + self.p.lookahead
        kt = min(int(np.searchsorted(self.path_s, s_t)), len(self.path) - 1)
        tgt = self.path[kt]
        dv = tgt - x[:2]
        err = float(wrap(np.arctan2(dv[1], dv[0]) - x[2]))
        om_des = float(np.clip(self.p.k_theta * err, -self.p.omega_max, self.p.omega_max))
        nu = float(np.clip(self.p.k_omega * (om_des - x[4]), -self.p.nu_max, self.p.nu_max))
        dist_goal = float(np.linalg.norm(self.goal - x[:2]))
        v_des = self.p.v_cruise * max(np.cos(err), 0.0) * min(1.0, dist_goal / 0.6 + 0.3)
        a = float(np.clip(self.p.k_v * (v_des - x[3]), -self.p.a_max, self.p.a_max))
        return np.array([a, nu]), k

    def at_goal(self, x):
        return float(np.linalg.norm(x[:2] - self.goal)) <= self.goal_tol

    def filter(self, x, u_nom, margins):
        p, sc = self.p, self.scene
        T, l = p.T, p.l
        hx = self.h(x)
        th, v, om = x[2], x[3], x[4]
        nus = np.unique(np.concatenate([self.nu_grid, [np.clip(u_nom[1], -p.nu_max, p.nu_max)]]))
        as_ = np.unique(np.concatenate([self.a_grid, [np.clip(u_nom[0], -p.a_max, p.a_max)]]))
        I0, I1 = self._integrals(th, om, nus)                                     # (n_nu, 2)
        thT = th + om * T + nus * T * T / 2
        base = x[:2] + v * I0 + l * e(thT)                                         # z at a = 0, (n_nu, 2)
        Z = base[:, None, :] + as_[None, :, None] * I1[:, None, :]                 # (n_nu, n_a, 2)
        vn = v + as_ * T
        omn = om + nus * T
        box = (np.abs(omn) <= p.omega_max + 1e-12)[:, None] & ((vn >= -p.v_rev - 1e-12) & (vn <= p.v_max + 1e-12))[None, :]
        hn = sc.psi(Z) - self.q_min - (vn ** 2)[None, :] / (2 * p.a_br)           # h(f_hat(x, u))
        J = (as_[None, :] - u_nom[0]) ** 2 + p.nu_weight * (nus[:, None] - u_nom[1]) ** 2
        out = np.full((len(margins), 2), np.nan)
        feas = np.zeros(len(margins), bool)
        for k, c in enumerate(margins):
            ok = box & (hn - c >= (1 - p.gamma) * hx)
            if ok.any():
                i, j = np.unravel_index(np.argmin(np.where(ok, J, np.inf)), J.shape)
                out[k] = (as_[j], nus[i])
                feas[k] = True
        return out, feas

    def pairwise_constants(self):
        return self.Kx, self.Ku


def make_model(name, **over):
    if name == "corridor":
        return Corridor(CorridorParams(**over))
    if name == "maze":
        return Maze(MazeParams(**over))
    raise ValueError(name)


def l_geom_all(model):
    return L_geom(model.grid, model.d_min, model.p.c_max)
