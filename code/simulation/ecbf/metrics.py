"""Per-trajectory metrics; decisions are averaged within each trajectory."""
from __future__ import annotations
import numpy as np

EPS_GEO = 0.0      # contact: clearance <= EPS_GEO at a recorded decision time (sampled contact, report Section 1.2.6)

# Transition quantities z_{k,t}: key -> function(rec) giving the length-N array (margins and losses in units of c_max).
TRANSITION = {
    "m_c": lambda r: np.asarray(r["c"], float) / float(r["c_max"]),
    "m_Y": lambda r: np.asarray(r["Y"], float) / float(r["c_max"]),
    "m_int": lambda r: np.sum((np.asarray(r["u"], float) - np.asarray(r["u_nom"], float)) ** 2, axis=1),
    "q_CBF": lambda r: np.asarray(r["cbf"], bool),
    "q_Y": lambda r: np.asarray(r["Y"], float) > np.asarray(r["c"], float),       # I^Y = 1{Y_t > c_t}
    "q_use": lambda r: np.asarray(r["use"], bool),
    "q_unsafe": lambda r: np.asarray(r["unsafe"], bool),
    "q_nocert": lambda r: np.asarray(r["nocert"], bool),
    "q_notrans": lambda r: np.asarray(r["notrans"], bool),
}


def trajectory(rec, eps_geo=EPS_GEO):
    """Per-trajectory values (report eqs. (25), (33), (35), (36)): m_{z,k}(N_k) for every TRANSITION key (nan when
    N_k = 0), T_mis,k = T[R N + (1-R) H],
    goal R_k, sampled contact (clearance <= eps_geo at a recorded position), infeasible termination, tau_use (s, inf
    without use), p_max (max quadrature-estimated p_t(c_t) over used decisions, nan without use), runtime, minimum h and clearance."""
    N, T, H = int(rec["N"]), float(rec["T"]), int(rec["H"])
    R = bool(rec["goal"])
    out = {"seed": int(rec["seed"]), "N": N, "goal": R, "T_mis": T * (N if R else H)}
    for key, fz in TRANSITION.items():
        out[key] = float(np.mean(fz(rec))) if N > 0 else np.nan
    clear = np.asarray(rec["clearance"], float)
    out["contact"] = bool(np.any(clear <= eps_geo))
    out["infeas"] = rec["termination"] == "infeasible"
    used = np.flatnonzero(np.asarray(rec["use"], bool))
    out["tau_use"] = T * used[0] if used.size else np.inf
    pd = np.asarray(rec["p_dep"], float)[used]
    out["p_max"] = float(np.nanmax(pd)) if used.size and np.any(np.isfinite(pd)) else np.nan   # missing != zero
    out["runtime_ms"] = float(rec["runtime_ms"])
    h = np.concatenate([np.asarray(rec["h"], float), np.asarray(rec["h_next"], float)])
    out["h_min"] = float(np.min(h)) if h.size else np.nan
    out["clearance_min"] = float(np.min(clear)) if clear.size else np.nan
    return out
