"""Obstacle geometry, geometric check and the fractions used by the risk-increase terms ([1, App. D-B])."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


@dataclass
class Scene:
    circles: np.ndarray                     # (nc, 3): x, y, rho
    walls: np.ndarray                       # (nw, 3): n_x, n_y, offset   (feasible side n^T z >= offset)
    rects: np.ndarray                       # (nr, 4): x0, y0, x1, y1
    rho_body: float
    extent: tuple = (0.0, 15.0, 0.0, 2.0)  # x_min, x_max, y_min, y_max (rendering, workspace)
    meta: dict = field(default_factory=dict)

    # ---------------------------------------------------------------------------------------------- clearances
    def obstacle_clearances(self, z):
        """Clearances psi_i(z) of every obstacle i (circles, walls, rectangles), shape z.shape[:-1] + (number of
        obstacles,); psi = min_i psi_i ([1, (49)], [1, (50)], report Section 1.1.1)."""
        z = np.asarray(z, float)
        out = []
        if len(self.circles):
            d = np.linalg.norm(z[..., None, :] - self.circles[:, :2], axis=-1)
            out.append(d - self.circles[:, 2] - self.rho_body)
        if len(self.walls):
            out.append(z @ self.walls[:, :2].T - self.walls[:, 2] - self.rho_body)
        if len(self.rects):
            out.append(sdist_rects(z, self.rects) - self.rho_body)
        return np.concatenate(out, axis=-1)

    def psi(self, z):
        return np.min(self.obstacle_clearances(z), axis=-1)

    def obstacle_kind(self, i):
        """Shape of obstacle i: 'circle', 'wall' or 'rect'."""
        nc, nw = len(self.circles), len(self.walls)
        return "circle" if i < nc else ("wall" if i < nc + nw else "rect")

    # ------------------------------------------------------------------------------------- effective radii and check
    def finite_offsets(self):
        """The set R of finite effective radii ([1, (55)]): cone radius + body radius, and the body radius for
        rectangle corners."""
        R = [float(r + self.rho_body) for r in self.circles[:, 2]]
        if len(self.rects):
            R.append(float(self.rho_body))
        return np.array(sorted(set(R)))

    def current_contribution(self, z, c_max):
        """Geometric check at a predicted nose point z: one obstacle determines psi on the disc of radius c_max about z
        ([1, (58)]) and, for a rectangle, the disc lies in one exterior face or corner region (report eq. (21)).
        Returns (ok, R_t) with R_t the effective radius of that obstacle on the disc (inf: affine)."""
        g = self.obstacle_clearances(z)
        i = int(np.argmin(g))
        others = np.delete(g, i)
        if others.size and np.min(others) < g[i] + 2 * c_max:
            return False, np.nan
        kind = self.obstacle_kind(i)
        if kind == "circle":
            return True, float(self.circles[i, 2] + self.rho_body)
        if kind == "wall":
            return True, np.inf
        x0, y0, x1, y1 = self.rects[i - len(self.circles) - len(self.walls)]
        zx, zy = float(z[0]), float(z[1])
        right, left, up, down = zx - x1 >= c_max, x0 - zx >= c_max, zy - y1 >= c_max, y0 - zy >= c_max
        in_x = zx - x0 >= c_max and x1 - zx >= c_max          # disturbance ball within the face's x-range
        in_y = zy - y0 >= c_max and y1 - zy >= c_max
        if (right or left) and (up or down):
            return True, float(self.rho_body)                  # corner region: circle about the corner
        if (right or left) and in_y or (up or down) and in_x:
            return True, np.inf                                # one exterior face throughout the ball
        return False, np.nan

    def past_offset(self, z):
        """Effective radius R_s of the nearest obstacle at a past nose point z (circle: rho + rho_body, wall: inf,
        rectangle: rho_body about its nearest point), as in the geometric correction ([1, (59)])."""
        g = self.obstacle_clearances(z)
        i = int(np.argmin(g))
        kind = self.obstacle_kind(i)
        if kind == "circle":
            return float(self.circles[i, 2] + self.rho_body)
        if kind == "wall":
            return np.inf
        return float(self.rho_body)


def sdist_rects(z, rects):
    """Signed distance to closed axis-aligned rectangles ([1, (50)]), shape z.shape[:-1] + (nr,)."""
    z = np.asarray(z, float)
    lo, hi = rects[:, :2], rects[:, 2:]
    c, half = (lo + hi) / 2, (hi - lo) / 2
    q = np.abs(z[..., None, :] - c) - half
    outside = np.linalg.norm(np.maximum(q, 0.0), axis=-1)
    inside = np.minimum(np.max(q, axis=-1), 0.0)
    return outside + inside


# -------------------------------------------------------------------------------------------------- scenes
def corridor_scene(rho_body=0.22):
    s = json.loads((HERE / "corridor_scene_v0_3_0.json").read_text())
    circles = np.array([[c["x"], c["y"], s["cone_radius"]] for c in s["cones"]], float)
    W = s["course"]["y_max"]
    walls = np.array([[0.0, 1.0, 0.0], [0.0, -1.0, -W]], float)          # y >= 0 and y <= W
    return Scene(circles, walls, np.zeros((0, 4)), rho_body, (0.0, s["course"]["x_max"], 0.0, W),
                 dict(name="corridor", source="corridor_scene_v0_3_0.json", goal=(s["goal"]["x"], s["goal"]["y"]),
                      goal_tol=s["goal"]["tol"], start_box=(*s["start_box"]["x"], *s["start_box"]["y"]),
                      heading_deg=s["start_box"]["heading_deg"], apron=s["apron"]["width"], raw=s))


def maze_scene(rho_body=0.22):
    s = json.loads((HERE / "maze_scene_v0_4_2.json").read_text())
    X, Y = s["size"]
    walls = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, -X], [0.0, 1.0, 0.0], [0.0, -1.0, -Y]], float)
    return Scene(np.array(s["circles"], float), walls, np.array(s["blocks"], float), rho_body, (0.0, X, 0.0, Y),
                 dict(name="maze", source="maze_scene_v0_4_2.json", goal=tuple(s["goal"]), goal_tol=0.15,
                      start=tuple(s["start"]), heading_deg=s["start_heading_deg"], raw=s))


# ------------------------------------------------------------------------------------------------ exceedance fractions
def p0(c, A):
    """Affine fraction P_inf(c, A) = pi^-1 arccos(c/A) for A > c, else 0 ([1, (54)], report eq. (20)); the fraction
    bound is Pbar(c) = P_inf(c, c_max)."""
    c, A = np.broadcast_arrays(np.asarray(c, float), np.asarray(A, float))
    out = np.zeros(c.shape)
    k = A > c
    out[k] = np.arccos(np.clip(c[k] / A[k], -1, 1)) / np.pi
    return out


def p1(d, c, A):
    """Directional exceedance fraction P_R(g, c, A) of a circular obstacle ([1, (53)]), at the distance d = g + R
    from its centre."""
    d, c, A = np.broadcast_arrays(np.asarray(d, float), np.asarray(c, float), np.asarray(A, float))
    out = np.zeros(d.shape)
    k = A > c
    arg = c[k] / A[k] + (A[k] ** 2 - c[k] ** 2) / (2 * d[k] * A[k])
    out[k] = np.arccos(np.clip(arg, -1, 1)) / np.pi
    return out


def F(R, q, c, A):
    """P_R(q, c, A) for clearance q and effective radius R ([1, (53)]); R = inf gives P_inf ([1, (54)])."""
    if np.isinf(R):
        return p0(c, A)
    return p1(q + R, c, A)


def L_geom(c, d_min, c_max):
    """Clearance sensitivity L^geom(c) ([1, (57)]), with d_min = varrho_min of [1, (55)]."""
    c = np.asarray(c, float)
    return np.sqrt(c_max ** 2 - c ** 2) / (np.pi * d_min * np.sqrt((2 * d_min - c) ** 2 - c_max ** 2))
