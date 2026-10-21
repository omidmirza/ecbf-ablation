#!/usr/bin/env python3
"""N200 table fragments of the paper (arXiv and extended versions) and of Section 2 of the ablation report.

Printed fragments: arXiv Table I (main_maze_arxiv_v8), Table IV (main_corridor_arxiv_v8), Table V
(stress_representative_arxiv_percent) and the totals of arXiv App. F-D (assumption_checks_arxiv, a macro file);
report Section 2.1 (table_key) and Tables 3-24 (design_overview, default_{scen}, topfive_full_{scen}, star_maze,
delta_{scen}, period_{scen}, b_{scen}, phi_{scen}, stress_{scen}_cases, checks_{scen}_{L0005,L005},
variation_{scen}). The paper source also includes paper/table_key, table_key_arxiv, design_overview,
default_{scen}, topfive_full_{scen}, topfive_full_{scen}_onepage, star_maze, phi_{scen}, checks_{scen}_{L0005,L005},
variation_{scen} and main_{scen}_extended_v8 (extended version). The other fragments are generated for reference and
checked like the printed ones (KEY_LINES.md lists all of them with their caption notes).

Layout (uniform layout, Layout(uniform=True), baseline_box_u, eq_tex, uniform_stretch): settings row, header, one-line
block row per setting (and per L_mu where it varies), centred rows and values at the full line width with equal
value-column widths, balanced baseline boxes, rows of one line stretched to fill their page, body size in the
appendix tables; multi-page tables break with \\ntablebreak. The catalog fragments (catalog_tables.py) keep the
boxed layout of Layout(uniform=False), except default_ and topfive_full_, which follow the editions.

Setting tuples: every setting tuple prints on one line. delta/period/b/phi, topfive_full, n30_record and the catalog
all-variant tables carry a full-width block row per setting (block_row, block_text); n30_record groups rows by the
changed parameter.

checks_{scen}_arxiv and stress_{scen}_arxiv (not printed; stress_{scen}_cases prints the same cells in the report):
chosen variant at S1-S5, both L_mu; checks: K of the setting and K = inf, six event columns; stress: phi x L_mu,
m_c/c_max, q_CBF^+, q_dep, q_A2, q_pw, dagger for phi < 1. Columns that are 0 or N-A in every row are omitted and
named in the fragment's "% OMITTED" line.

Baseline boxes: every baseline box is centred and spans the full line width (its Method column absorbs the value
columns the box drops) and its first-column Method header restates the left box rule.

Fragment set: table_key and table_key_arxiv (the full keys; besides them only sec6_short, sec6_combined and
sec6_maze_arxiv carry a key box, a two-line key); sec6_short (chosen variant per selected setting, centred baseline
box); topfive_full_{scen} (all four variants, all metrics with q_dep and q_cert, three value lines per row); checks,
detailed events and the delta / K / b / phi parameter tables at the chosen variant only (paper: Default and S3 for
delta, K and b, all six for phi; report: all six); row labels 'Method (alpha0, alpha, delta, K)'; report fragments at
505.89 pt (geometry margin 0.75in), 11 pt. Any other .tex file in the output folders is moved to _superseded/.

The generator extends build_tables_v10_original.py (kept in this folder with input paths and docstrings adapted) with
  * the N200 selected settings S1-S5 and the default setting, all four variants at each;
  * the new-anchor sweeps (b, phi x L_mu, K in {K, inf} x L_mu at S1-S5; N200_numbers/new_anchors);
  * one count per event cell: confirmed plus possible events (the conservative count), "a (u)" with the
    one-sided 95% Clopper-Pearson upper bound u in % rounded up, and \\checkmark when u <= the guarantee level
    (delta at finite K, delta_1 = delta/(2 b ln^2 2) at K = infinity) for the four variants;
  * T_s for the sampling period in the settings lines (T is the decision horizon);
  * a report layout at the report text width (505.89 pt, 11 pt Times): three value lines per method
    and five value columns instead of the paper's two lines and eight columns.
Fragments carry no caption, label or \\clearpage; page blocks are separated by \\ntablebreak. Every printed number
is recorded in tables_values.json.gz; check_tables.py re-derives it.

Inputs (not included): the N200 summary exports in $ECBF_N200_NUMBERS (default <repo>/N200_numbers) and the
candidate review in $ECBF_CANDIDATES (default <repo>/candidate_review.json).
Usage: python3 build_tables.py [--catalog]; then compile_tests.py (layout) and check_tables.py (values).
--catalog also writes the catalog fragments (catalog_tables.py); the edition fragments and
tables_values.json.gz are the same with or without it.
Writes appendices/generated_n200/{paper,report}/*.tex and tools/n200_tables/tables_values.json.gz.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from decimal import ROUND_CEILING, ROUND_HALF_UP, Decimal
from pathlib import Path

from paper_source import unmark   # text read from the paper source

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1]
GEN = SRC / "appendices" / "generated_n200"
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_tables/<this file>)
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports (not included)
FILES = {
    "final": NUM / "final/common/metrics.json",
    "new": NUM / "new_anchors/common/metrics.json",
    "final_temporal": NUM / "final/temporal/temporal_metrics.json",
}
CANDIDATES = Path(os.environ.get("ECBF_CANDIDATES", REPO / "candidate_review.json"))   # candidate review of the unified grid


def rel_input(p):
    """An input path as recorded in the values files and manifests: 'N200_numbers/<path relative to NUM>' for an
    export, the file name for the candidate review, otherwise the path unchanged."""
    p = Path(p)
    try:
        return "N200_numbers/" + p.relative_to(NUM).as_posix()
    except ValueError:
        return p.name if p == CANDIDATES else str(p)


TARGETS = ("paper", "report")
NMC = 200
LN2SQ = math.log(2) ** 2

# ----------------------------------------------------------------------------------------------- notation
SCEN = {"corridor": dict(name="Corridor", cmax="0.02", Ts="0.1"),
        "maze": dict(name="Maze", cmax="0.01", Ts="0.2")}
SCENS = ("corridor", "maze")
VARIANTS = ["V-U", "V-D", "V-Cd", "V-C"]
VNAME = {"V-U": r"\textcolor{vu}{V-U}", "V-D": r"\textcolor{vd}{V-$\Delta$}",
         "V-Cd": r"\textcolor{vl}{V-L}", "V-C": r"\textcolor{vc}{V-C}"}
VSHORT = {"V-U": "V-U", "V-D": "V-$\\Delta$", "V-Cd": "V-L", "V-C": "V-C"}
BASELINES = ["Unfiltered", "Zero margin", "Conservative", "ACI", "QT"]
BNAME = {"Unfiltered": "Unfiltered", "Zero margin": "$c=0$", "Conservative": "$c_{\\max}$", "ACI": "ACI", "QT": "QT"}
DELTAS = [0.025, 0.05, 0.10, 0.15, 0.20]
KS = [10, 50, 100, 200, "infinity"]
ALPHA0S = [0.10, 0.15, 0.20, 0.25, 0.30]
BASE = dict(L_mu=0.0005, bound_scale=1.0, budget_normalizer=2.11)

# The six settings: name, tuple, variant of the selection, data source of the anchor sweeps.
SIX = [
    ("Default", dict(alpha0=0.15, alpha=0.10, delta=0.05, K=50), "V-C", ("final", "current_default")),
    ("S1", dict(alpha0=0.30, alpha=0.25, delta=0.20, K=50), "V-C", ("new", "n200_role1_")),
    ("S2", dict(alpha0=0.30, alpha=0.25, delta=0.15, K=50), "V-C", ("new", "n200_role2_")),
    ("S3", dict(alpha0=0.30, alpha=0.20, delta=0.20, K=50), "V-Cd", ("new", "n200_role3_")),
    ("S4", dict(alpha0=0.30, alpha=0.25, delta=0.025, K=100), "V-D", ("new", "n200_role4_")),
    ("S5", dict(alpha0=0.30, alpha=0.25, delta=0.10, K=10), "V-C", ("new", "n200_role5_")),
]

# metric id -> (symbol, arrow, unit, scale, decimals, boldable), as in build_tables_v10_original.py
MET = {
    "m_c": (r"\widehat m_c/c_{\max}", "down", None, 1, 2, True),
    "m_Y": (r"\widehat m_Y/c_{\max}", "down", None, 1, 2, True),
    "m_int": (r"\widehat m_{\|u-u^{\rm nom}\|^2}", "down", None, 1, 2, True),
    "n_certificates": (r"\widehat m_{\rm cert}", None, None, 1, 2, False),
    "q_CBF": (r"\widehat q_{\rm CBF}^+", "down", "%", 100, 2, True),
    "q_Y": (r"\widehat q_Y^+", "down", "%", 100, 2, True),
    "q_use": (r"\widehat q_{\rm use}^+", "up", "%", 100, 1, True),
    "q_unsafe": (r"\widehat q_{\rm unsafe}^+", "down", "%", 100, 2, True),
    "q_nocert": (r"\widehat q_{\rm nocert}^+", "down", "%", 100, 1, True),
    "q_notrans": (r"\widehat q_{\rm nodep}^+", "down", "%", 100, 1, True),
    "h_min": (r"\widehat m_{h^{\min}}", "up", "m", 1, 2, True),
    "clearance_min": (r"\widehat m_{\psi^{\min}}", "up", "m", 1, 2, True),
    "tau_use": (r"\widehat m_{\tau^{\rm use}}", "down", "s", 1, 2, False),
    "p_max": (r"\widehat m_{p^{\max}}", "down", "%", 100, 1, False),
    "T_mis": (r"\widehat T_{\rm mis}", "down", "s", 1, 2, True),
    # pk: the peak estimated exceedance risk statistics.p_max as a fraction (3 decimals, like m_c/c_max) in the
    # Section VI tables (sec6_*), plain (no arrow, no bold, no colour); mean and SD over the trajectories that use a
    # margin below c_max (available_mean, available_sd), N-A when no trajectory does
    "pk": (r"\widehat m_{p^{\max}}", None, None, 1, 3, False),
}
STAT_SRC = {"pk": "p_max"}            # metric id -> statistics key of the exports
STAT_AVAILABLE = {"pk"}               # these print available_mean / available_sd
# event id -> (symbol, source, guarantee-checked)
EVT = {
    "qcert": (r"\widehat q_{\rm cert}", ("trajectory_events", "qcert"), True),
    "qdep": (r"\widehat q_{\rm dep}", ("trajectory_events", "qdep"), True),
    "rCBF": (r"\widehat q_{\rm CBF}", ("realized_transition_events", "q_CBF"), False),
    "runsafe": (r"\widehat q_{\rm unsafe}", ("realized_transition_events", "q_unsafe"), False),
    "qA2": (r"\widehat q_{\rm A2}", ("risk_increase_checks", "direct"), False),
    "qpw": (r"\widehat q_{\rm pw}", ("risk_increase_checks", "pairwise"), False),
}
# temporal id -> (symbol, quantity, statistic, scale)
TMP = {
    "Vc": (r"10^3\widehat V_{c/c_{\max}}", "normalized_margin", "temporal_variance", 1000),
    "Jc": (r"100\widehat J_{c/c_{\max}}", "normalized_margin", "step_rms", 100),
    "Vp": (r"10^3\widehat V_{\widehat p}", "estimated_applied_risk", "temporal_variance", 1000),
    "Jp": (r"100\widehat J_{\widehat p}", "estimated_applied_risk", "step_rms", 100),
    "nK": (r"|\mathcal K_{c/c_{\max}}|,|\mathcal K_{\widehat p}|", None, None, None),
}
TMP_DEC = 3

PAIRS_P = [["m_c", "m_Y"], ["m_int", "n_certificates"], ["q_CBF", "q_Y"], ["q_use", "q_unsafe"],
           ["q_nocert", "q_notrans"], ["h_min", "clearance_min"], ["tau_use", "p_max"], ["T_mis"]]
TRIPLES_R = [["m_c", "m_Y", "m_int"], ["q_CBF", "q_Y", "q_unsafe"], ["q_use", "tau_use", "p_max"],
             ["q_nocert", "q_notrans", "n_certificates"], ["h_min", "clearance_min", "T_mis"]]
RULE_METHOD = "0.25pt"


def full_cols(target):
    return (PAIRS_P, 2) if target == "paper" else (TRIPLES_R, 3)


# ----------------------------------------------------------------------------------------------- formatting
def half_up(v: float, d: int) -> str:
    q = Decimal(f"{v:.{d + 9}f}").quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    out = f"{q:.{d}f}"
    return out[1:] if out.startswith("-") and q == 0 else out


def upper_pct(u: float) -> str:
    """One-sided 95% upper bound in %, rounded up to one decimal."""
    q = (Decimal(repr(u)) * 100).quantize(Decimal("0.1"), rounding=ROUND_CEILING)
    return f"{q:.1f}"


def cp_upper(x: int, n: int, conf: float = 0.95) -> float:
    """Exact one-sided Clopper-Pearson upper bound (same rule as the stored upper95 fields)."""
    if x >= n:
        return 1.0

    def cdf(p):
        if p <= 0:
            return 1.0
        if p >= 1:
            return 0.0
        lp, lq = math.log(p), math.log1p(-p)
        return sum(math.exp(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * lp + (n - k) * lq)
                   for k in range(x + 1))

    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(mid) > 1 - conf:
            lo = mid
        else:
            hi = mid
    return hi


def num(x) -> str:
    """Setting value: at least two decimals, no trailing zeros beyond."""
    if x == "infinity":
        return r"\infty"
    s = f"{float(x):.4f}".rstrip("0")
    if len(s.split(".")[1]) < 2:
        s = f"{float(x):.2f}"
    return s


def kfmt(K):
    return r"\infty" if K == "infinity" else str(int(float(K)))


def key(x):
    return x if isinstance(x, str) else round(float(x), 6)


def guarantee_level(cfg):
    """delta at finite K (arXiv Theorem 3); delta_1 = delta/(2 b ln^2 2) at K = infinity (arXiv Section V-B)."""
    if cfg["K"] == "infinity":
        return cfg["delta"] / (2 * cfg["budget_normalizer"] * LN2SQ)
    return cfg["delta"]


def delta1_pct(delta, b=2.11):
    return half_up(100 * delta / (2 * b * LN2SQ), 2)


# The budget normalizer b equals kappa (arXiv Section V-B, kappa ~ 2.11; the simulations use 2.11)
# by default. Displays print the symbol b only where b is compared; the default value prints as \kappa.
DEFAULT_B = 2.11


def is_default_b(x):
    try:
        return x is not None and abs(float(x) - DEFAULT_B) < 1e-9
    except (TypeError, ValueError):
        return False


def b_value(x):
    """Printed value of b (math mode, no $): \\kappa for the default, else the shortest plain form."""
    return r"\kappa" if is_default_b(x) else f"{float(x):g}"


# ----------------------------------------------------------------------------------------------- data
class Data:
    def __init__(self):
        self.idx, self.sha, self.anchor_ids = {}, {}, {}
        for src in ("final", "new"):
            path = FILES[src]
            print(f"loading {path} ...", flush=True)
            raw = path.read_bytes()
            self.sha[src] = hashlib.sha256(raw).hexdigest()
            d = json.loads(raw)
            del raw
            for v in d["views"]:
                df = v["definition"]
                fam = df["family"]
                assert len(v["cells"]) == 1
                c = v["cells"][0]
                for s in c.get("statistics", {}).values():
                    s.pop("observed_seeds", None)
                seeds = c.pop("seeds", None)
                assert c["n"] == NMC and seeds == list(range(970001, 970001 + NMC)), v["view_id"]
                k = (src, fam) + tuple(sorted((a, key(b)) for a, b in c["configuration"].items() if b is not None))
                assert k not in self.idx, k
                self.idx[k] = dict(src=src, view_id=v["view_id"], family=fam, anchor_id=df.get("anchor_id"),
                                   n100_view_id=df.get("n100_view_id"), cell=c)
                if df.get("anchor_id"):
                    self.anchor_ids.setdefault(src, set()).add(df["anchor_id"])
            del d
        self.tmp = {}
        path = FILES["final_temporal"]
        print(f"loading {path} ...", flush=True)
        raw = path.read_bytes()
        self.sha["final_temporal"] = hashlib.sha256(raw).hexdigest()
        d = json.loads(raw)
        del raw
        for v in d["views"]:
            c = v["cells"][0]
            k = ("final_temporal", v["definition"]["family"]) + tuple(
                sorted((a, key(b)) for a, b in c["configuration"].items() if b is not None))
            st = c["statistics"]
            keep = {q: {m: st[q][m] for m in ("temporal_variance", "step_rms")} | {"assessed_n": st[q]["assessed_n"],
                                                                                   "total_n": st[q]["total_n"]}
                    for q in ("normalized_margin", "estimated_applied_risk")}
            self.tmp[k] = dict(src="final_temporal", view_id=v["view_id"], family=v["definition"]["family"],
                               cell=dict(configuration=c["configuration"], configuration_id=c["configuration_id"],
                                         n=c["n"], temporal=keep))
        del d
        print(f"  {len(self.idx)} common cells, {len(self.tmp)} temporal cells", flush=True)

    def get(self, src, fam, **cfg):
        k = (src, fam) + tuple(sorted((a, key(b)) for a, b in cfg.items() if b is not None))
        return (self.tmp if src == "final_temporal" else self.idx).get(k)

    def anchor(self, src, prefix):
        ids = [a for a in self.anchor_ids[src] if a.startswith(prefix)]
        assert len(ids) == 1, (src, prefix, ids)
        return ids[0]


# ----------------------------------------------------------------------------------------------- record
class Recorder:
    def __init__(self):
        self.entries = []

    def add(self, **kw):
        self.entries.append(kw)


REC = Recorder()


def prov(found):
    c = found["cell"]
    return dict(src=found["src"], view_id=found["view_id"], n100_view_id=found.get("n100_view_id"),
                family=found["family"], configuration_id=c.get("configuration_id"), configuration=c["configuration"])


# ----------------------------------------------------------------------------------------------- cells
def split_spec(spec):
    """'metric[#decimals][@sub]' -> (metric, sub); sub selects the cell (scenario, K or infinity)."""
    m, _, sub = spec.partition("@")
    return m.partition("#")[0], sub


def spec_decimals(spec):
    m = spec.partition("@")[0]
    return int(m.partition("#")[2]) if "#" in m else None


def stat_cell(found, metric, dec_override=None):
    if found is None:
        return "--", dict(kind="not_simulated")
    src = STAT_SRC.get(metric, metric)
    s = found["cell"]["statistics"][src]
    sym, arrow, unit, scale, dec, boldable = MET[metric]
    dec = dec_override or dec
    info = dict(prov(found), source=f"statistics.{src}", scale=scale, decimals=dec)
    mk, sk = ("available_mean", "available_sd") if metric in STAT_AVAILABLE else ("mean", "sd")
    if not s["applicable"] or s["observed_n"] == 0 or s[mk] is None:
        return r"\na", dict(info, kind="NA")
    m, sd = half_up(s[mk] * scale, dec), half_up(s[sk] * scale, dec)
    if metric in ("h_min", "clearance_min") and s["mean"] < 0 and Decimal(m) == 0:   # keep the sign
        m = "-" + m
    info.update(kind="meanSD", mean_str=m, sd_str=sd, conditional=s["conditional"],
                boldable=bool(boldable and not s["conditional"]), arrow=arrow)
    return None, info


def event_cell(found, ev, certified):
    """One conservative count: possible events (they include the confirmed ones)."""
    if found is None:
        return "--", dict(kind="not_simulated")
    grp, sub = EVT[ev][1]
    e = found["cell"][grp][sub]
    info = dict(prov(found), source=f"{grp}.{sub}")
    if grp == "trajectory_events":
        n, cnt, u = e["applicable_n"], e["possible"], e["upper95_possible"]
    elif grp == "realized_transition_events":
        n, cnt, u = e["assessed_n"], e["count"], None   # CBF / safe-set violation counts: printed without a bound
    else:
        n, cnt = e["assessed_trajectory"], e["any_excess"]
        u = cp_upper(cnt, NMC) if n else None   # q_A2, q_pw divide by all N_MC trajectories (report Section 1.4)
    if not n:
        return r"\na", dict(info, kind="NA", assessed=n)
    up = upper_pct(u) if u is not None else None
    s = (f"{cnt}" if n == NMC else f"{cnt}/{n}") + (f" ({up})" if up is not None else "")
    check = False
    if certified and EVT[ev][2]:
        level = guarantee_level(found["cell"]["configuration"])
        check = Decimal(up) <= Decimal(repr(100 * level))
        info["level"] = level
    if check:
        s += r"\,\checkmark"
    return s, dict(info, kind="event", count=cnt, n=n, upper95=u, upper_str=up, check=check, certified=certified)


# Submission ranges ("the submission": the earlier submitted version of the paper): the baseline cells that its
# Table I prints as a range print the same range "confirmed--possible (bound of possible)" in every table that prints
# that baseline cell (sec6_combined; default_maze and checks_maze_L0005 of the paper and the report; the catalog's
# checks_all_maze_L0005 and topfive_events; the baselines at its L_mu); every other event cell keeps the one
# conservative count. check_tables.py re-derives it.
SUBMISSION_RANGES = {("maze", "QT", 0.15, "qdep")}
SUBMISSION_L_MU = 0.0005


def apply_submission_ranges(rows, scen, alpha0=None):
    """scen: the scenario of single-scenario rows; for multi-scenario rows the spec's sub names the scenario.
    alpha0: the rows' alpha0 (default: each row's recorded alpha0)."""
    for r in rows:
        if abs(r["record"].get("L_mu", SUBMISSION_L_MU) - SUBMISSION_L_MU) > 1e-12:
            continue
        a0 = r["record"]["alpha0"] if alpha0 is None else alpha0
        for spec, cell in r.get("cells", {}).items():
            met, sub = split_spec(spec)
            sc = sub if sub in SCENS else scen
            if (sc, r["method"], a0, met) not in SUBMISSION_RANGES:
                continue
            info = cell[1]
            e = r["found"][sub if sub in r["found"] else ""]["cell"]["trajectory_events"][met]
            assert e["applicable_n"] == NMC and e["possible"] == info["count"]
            if e["confirmed"] != e["possible"]:
                cell[0] = f"{e['confirmed']}--{e['possible']} ({info['upper_str']})"
                info.update(count_rule="range", confirmed=e["confirmed"])


def temporal_cell(found, tid):
    if found is None:
        return "--", dict(kind="not_simulated")
    t = found["cell"]["temporal"]
    info = dict(prov(found), source=f"temporal.{tid}")
    if tid == "nK":
        a, b = t["normalized_margin"]["assessed_n"], t["estimated_applied_risk"]["assessed_n"]
        return f"{a}, {b}", dict(info, kind="count", counts=[a, b])
    sym, q, stat, scale = TMP[tid]
    s = t[q][stat]
    if not s or s.get("n", 0) < 2 or s.get("mean") is None:
        return r"\na", dict(info, kind="NA")
    m, sd = half_up(s["mean"] * scale, TMP_DEC), half_up(s["sd"] * scale, TMP_DEC)
    return f"${m}\\vpm {sd}$", dict(info, kind="tmeanSD", mean_str=m, sd_str=sd, scale=scale, decimals=TMP_DEC)


def fill_rows(rows, specs):
    """rows: dicts with 'found' (sub -> found), 'is_variant', 'group'. Adds rows[i]['cells'][spec]."""
    groups = {}
    for r in rows:
        r.setdefault("cells", {})
        for spec in specs:
            m, sub = split_spec(spec)
            f = r["found"].get(sub)
            if m in MET:
                tex, info = stat_cell(f, m, spec_decimals(spec))
            elif m in EVT:
                tex, info = event_cell(f, m, r["is_variant"])
            else:
                tex, info = temporal_cell(f, m)
            r["cells"][spec] = [tex, info]
        if r["is_variant"] and r.get("group") is not None:
            groups.setdefault(r["group"], []).append(r)
    for g, rs in groups.items():
        for spec in specs:
            m, _ = split_spec(spec)
            if m not in MET:
                continue
            infos = [r["cells"][spec][1] for r in rs]
            if len(infos) != 4 or not all(i.get("kind") == "meanSD" and i["boldable"] for i in infos):
                continue
            vals = [Decimal(i["mean_str"]) for i in infos]
            best = max(vals) if infos[0]["arrow"] == "up" else min(vals)
            for i, v in zip(infos, vals):
                i["bold"] = v == best
    for r in rows:
        for spec, (tex, info) in r["cells"].items():
            if tex is None:
                mm = rf"\mathbf{{{info['mean_str']}}}" if info.get("bold") else info["mean_str"]
                r["cells"][spec][0] = f"${mm}\\vpm {info['sd_str']}$"
                info.setdefault("bold", False)


# ----------------------------------------------------------------------------------------------- headers
def head(spec):
    m, _ = split_spec(spec)
    if m in MET:
        sym, arrow, unit, *_ = MET[m]
    elif m in EVT:
        sym, arrow, unit = EVT[m][0], "down", None
        if m in ("qA2", "qpw"):
            arrow = None
    else:
        sym, arrow, unit = TMP[m][0], ("down" if m != "nK" else None), None
    a = {"down": r"\!\downarrow", "up": r"\!\uparrow", None: ""}[arrow]
    u = {None: "", "%": r"\,(\%)", "m": r"\,(m)", "s": r"\,(s)"}[unit]
    return f"${sym}{a}${u}"


KEY_MSD = (r"$\pm$: one sample SD over the trajectories that enter the mean; $^+$: decision fraction; "
           r"$\downarrow$/$\uparrow$: lower/higher is better; bold: lowest ($\downarrow$) or highest ($\uparrow$) mean of the four variants at one setting")
KEY_NA = r"; N-A: not applicable"
KEY_DASH = r"; --: not simulated at this $\alpha_0$"
KEY_EVT = (r"Counts out of $N_{\rm MC}=200$ trajectories, confirmed and possible events together (possible: not excluded "
           r"by the numerical risk bounds); parentheses: $95\%$ upper confidence bound (\%); "
           r"\checkmark: bound $\le\delta$ at finite $K$ and $\le\delta_1$ in \eqref{eq:empirical-first-test-budget} at $K=\infty$")
KEY_EVT_NOCHECK = (r"Counts out of $N_{\rm MC}=200$ trajectories, confirmed and possible events together (possible: not "
                   r"excluded by the numerical risk bounds); parentheses: $95\%$ upper confidence bound (\%)")
KEY_ASSESSED = r"; $a/m$: $a$ trajectories with an excess above $\varepsilon_{\rm chk}$ among the $m$ with a checked decision (\appsecref{app:empirical-assumption-checks}; $a$ alone if $m=N_{\rm MC}$); the parenthesized bound uses all $N_{\rm MC}$ trajectories"
KEY_QNOPLUS = r"; $\widehat q$ without $^+$: trajectories with the event"


# ----------------------------------------------------------------------------------------------- LaTeX tables
class Widths:
    def __init__(self, gid):
        self.gid = gid
        self.cols = {}

    def add(self, col, s):
        if s and "multirow" not in s:
            self.cols.setdefault(col, set()).add(s)

    def name(self, col):
        return f"{self.gid}{col}"

    def tex(self):
        out = []
        for col, ss in self.cols.items():
            out.append(rf"\vnnew{{{self.name(col)}}}")
            for s in sorted(ss):
                out.append(rf"\vnfit{{{self.name(col)}}}{{{s}}}")
        return "\n".join(out)


def stack(lines):
    return r"\begin{tabular}{@{}c@{}}" + r"\\".join(lines) + r"\end{tabular}"


class Layout:
    """Column layout shared by the main table(s) and the baseline boxes of one fragment.

    level_cols: list of (col_id, header); value_cols: list of spec lists (1..nlines entries);
    groups: optional list of (label, first value column, span) printed above the metric headers."""

    def __init__(self, gid, level_cols, value_cols, nlines, rule_weights, method_header="", groups=None,
                 align="c", max_sep=None, inner_rules=True, stretch=None, uniform=False):
        self.gid, self.level_cols, self.value_cols, self.nlines = gid, level_cols, value_cols, nlines
        self.method_header, self.rule_weights, self.groups = method_header, rule_weights, groups or []
        # align "r" right-aligns the value columns (headers centred; a cell without a check mark in a column with
        # check marks ends in \vnnock, the check mark's width, so that the bounds and the +- signs line up); max_sep
        # caps the column separation (the table keeps its natural width, centred, instead of spreading to
        # \linewidth). No fragment uses either: every fragment keeps align "c" and max_sep None. inner_rules=False
        # drops the rules between the rows of one innermost level block (rules only between level blocks).
        self.align, self.max_sep, self.inner_rules = align, max_sep, inner_rules
        self.stretch = stretch   # \arraystretch of this fragment (local to \ntable's group); None: 0.93
        # uniform layout of every edition table (see the uniform-layout section below): equal value-column
        # widths at the full line width, centred cells, \vnnock padding in check-mark columns, balanced baseline boxes.
        self.uniform = uniform
        # size: appendix tables at the body size in every context (a table float would set them at 8 pt in the
        # extended version); Section VI tables keep the size of their table* float (normalsize=False)
        self.normalsize = uniform
        self.heads = []          # TeX of independent baseline boxes (widths, equalization, separation)
        self.pad_all = None      # check-mark columns of the whole fragment (\vnnock padding on every page)
        self.W = Widths(gid)
        self.ncol = len(level_cols) + 1 + len(value_cols)
        self.hl = max(len(v) for v in value_cols)
        for cid, h in level_cols:
            self.W.add(cid, h)
        self.W.add("meth", method_header)
        if not level_cols and not uniform:  # the baseline box shares the Method column (uniform: only when aligned)
            for b in BNAME.values():
                self.W.add("meth", b)
        for j, vc in enumerate(value_cols):
            for m in vc:
                self.W.add(f"v{j}", head(m))
        self.specs = [m for vc in value_cols for m in vc]

    def w(self, col):
        return rf"\vnw{{{self.W.name(col)}}}"

    def sep_tex(self):
        widths = " - ".join([self.w(c) for c, _ in self.level_cols] + [self.w("meth")] +
                            [self.w(f"v{j}") for j in range(len(self.value_cols))])
        n_inner = len(self.value_cols)
        out = (rf"\vnnewsep{{{self.gid}}}{{\dimexpr(\linewidth - {widths} - {n_inner}\arrayrulewidth"
               rf" - 2\vnbox)/{2 * self.ncol}\relax}}")
        if self.max_sep:
            out += rf"\ifdim\vnsep{{{self.gid}}}>{self.max_sep}\vnnewsep{{{self.gid}}}{{{self.max_sep}}}\fi"
        return out

    def inner_tex(self):
        """Width of a full-width header row's p column (settings line): the line width minus the box rules, or, when
        the separation is capped, the table's natural inner width."""
        if not self.max_sep:
            return r"\linewidth-2\vnbox-2\tabcolsep"
        widths = " + ".join([self.w(c) for c, _ in self.level_cols] + [self.w("meth")] +
                            [self.w(f"v{j}") for j in range(len(self.value_cols))])
        return rf"{widths} + {2 * self.ncol - 2}\tabcolsep + {len(self.value_cols)}\arrayrulewidth"


def rule_between(lay, prev, cur, nlev, method_col):
    for k in range(nlev):
        if prev["levels"][k] != cur["levels"][k]:
            w = lay.rule_weights[k]
            if k == 0:
                return rf"\Xhline{{{w}}}"
            return rf"\Xcline{{{k + 1}-{lay.ncol}}}{{{w}}}"
    if not lay.inner_rules:
        return None
    return rf"\Xcline{{{method_col}-{lay.ncol}}}{{{RULE_METHOD}}}"


def level_spans(rows, nlev, nlines):
    spans = []
    for k in range(nlev):
        sp = [0] * len(rows)
        i = 0
        while i < len(rows):
            j = i
            while (j + 1 < len(rows) and rows[j + 1]["levels"][:k + 1] == rows[i]["levels"][:k + 1]
                   and rows[j + 1].get("block_head") == rows[i].get("block_head")):   # no span across blocks
                j += 1
            sp[i] = (j - i + 1) * nlines
            i = j + 1
        spans.append(sp)
    return spans


UNIFORM_PAD = "4pt"   # fixed side padding of full-width rows (settings and block rows), whatever the column separation


def prow(ncol, width, text, rule=r"\hline", pad=None):
    """Full-width wrapped header row (settings line or key) of a boxed table of inner width `width`.
    pad: fixed side padding instead of \\tabcolsep (the text starts at the box rule whatever the separation)."""
    if pad:
        spec = (rf"!{{\vrule width\vnbox}}@{{\hskip {pad}}}p{{\dimexpr {width} + 2\tabcolsep - 2\dimexpr {pad}\relax\relax}}"
                rf"@{{\hskip {pad}}}!{{\vrule width\vnbox}}")
    else:
        spec = rf"!{{\vrule width\vnbox}}p{{\dimexpr {width}\relax}}!{{\vrule width\vnbox}}"
    return rf"\multicolumn{{{ncol}}}{{{spec}}}{{\raggedright {text}}}\\{rule}"


BLOCK_RULE = "1.0pt"   # rule above a block row (the block's own lines keep their usual separation)


def block_row(ncol, text, pad=None):
    """Full-width one-line block row naming the setting of the rows below it, e.g.
    'S1: V-C (alpha0, alpha, delta, K) = (0.30, 0.25, 0.20, 50)'; an l column, so the tuple never wraps.
    pad: fixed left padding instead of \\tabcolsep."""
    spec = (rf"!{{\vrule width\vnbox}}@{{\hskip {pad}}}l!{{\vrule width\vnbox}}" if pad
            else r"!{\vrule width\vnbox}l!{\vrule width\vnbox}")
    return rf"\multicolumn{{{ncol}}}{{{spec}}}{{{text}}}\\\Xhline{{0.4pt}}"


def dn(name):
    """Displayed setting name: the defined term 'default setting' replaces the label 'Default'."""
    return "default setting" if name == "Default" else name


FIELDS4 = ("alpha0", "alpha", "delta", "K")
FSYM = dict(alpha0=r"\alpha_0", alpha=r"\alpha", delta=r"\delta", K="K")
NO_K = ("alpha0", "alpha", "delta")   # block-row fields when K varies within the block


def fields_sym(fields=FIELDS4):
    """The symbol tuple '(alpha_0, alpha, delta, K)' in math mode, restricted to `fields`."""
    return "$(" + ",".join(FSYM[f] for f in fields) + ")$"


def block_text(name, t, v=None, extra="", fields=FIELDS4):
    """Block-row text: bold setting name, the variant (if one per block), the tuple of `fields` on one line
    (a swept field is never printed in the block row)."""
    vv = f"{VNAME[v]} " if v else ""
    return rf"\textbf{{{dn(name)}}}: {vv}{fields_sym(fields)} = {tuple_plain(t, fields)}{extra}"


def value_header_lines(lay, cols, colw, last_sep):
    """Header lines for value columns. cols: list of (j, specs) shown; colw(j) gives the column width macro."""
    hl = max(len(k) for _, k in cols)
    lines = []
    for jline in range(hl):
        cells = []
        for n_, (vj, k) in enumerate(cols):
            if len(k) == hl or (len(k) > 1 and jline < len(k)):
                c = head(k[jline])
            elif len(k) == 1:
                c = (rf"\multirow{{{hl}}}{{*}}{{{head(k[0])}}}" if hl > 1 else head(k[0])) if jline == 0 else ""
            else:
                c = ""
            last = n_ == len(cols) - 1
            cells.append(c if (last or not last_sep) else rf"\multicolumn{{1}}{{w{{c}}{{{colw(vj)}}}|}}{{{c}}}")
        lines.append(cells)
    return lines, hl


def value_line_cells(lay, r, cols, jline, L, fragment, target, table_id, wcol, pad_check=()):
    cells = []
    for vj, k in cols:
        if jline < len(k) and (len(k) > 1 or jline == 0):
            m = k[jline]
            tex, info = r["cells"][m]
            REC.add(target=target, fragment=fragment, table=table_id, row=r["record"], method=r["method"], spec=m,
                    printed=tex, **info)
            if m in pad_check and r"\checkmark" not in tex:
                tex += r"\vnnock"
            lay.W.add(wcol(vj), tex)
            cells.append(rf"\multirow{{{L}}}{{*}}{{{tex}}}" if (len(k) == 1 and L > 1) else tex)
        else:
            cells.append("")
    return cells


def group_header_line(lay, cols, nleft, last_sep):
    """Group labels (e.g. Corridor / Maze) over value columns; cols: kept (j, specs) in order."""
    if not lay.groups:
        return None, None
    out, i, starts = [""] * nleft, 0, []
    pos = {vj: n for n, (vj, _) in enumerate(cols)}
    cells = []
    for label, j0, span in lay.groups:
        kept = [j for j in range(j0, j0 + span) if j in pos]
        if not kept:
            continue
        starts.append((pos[kept[0]] + nleft + 1, pos[kept[-1]] + nleft + 1))
        cells.append((len(kept), label, pos[kept[-1]] == len(cols) - 1))
    line = out[:]
    for span, label, is_last in cells:
        sep = r"!{\vrule width\vnbox}" if is_last else "|"   # the last group restates the right box rule
        line.append(rf"\multicolumn{{{span}}}{{c{sep}}}{{{label}}}")
    clines = "".join(rf"\Xcline{{{a}-{b}}}{{{RULE_METHOD}}}" for a, b in starts)
    return " & ".join(line) + r"\\" + clines, len(starts)


def main_table(lay, rows, settings, keyrow, fragment, target, table_id):
    gid, nlev, L = lay.gid, len(lay.level_cols), lay.nlines
    cs = [rf"w{{c}}{{{lay.w(c)}}}" for c, _ in lay.level_cols] + [rf"w{{c}}{{{lay.w('meth')}}}"]
    vs = [rf"w{{{lay.align}}}{{{lay.w(f'v{j}')}}}" for j in range(len(lay.value_cols))]
    spec = r"!{\vrule width\vnbox}" + "".join(cs) + "|" + "|".join(vs) + r"!{\vrule width\vnbox}"
    inner = lay.inner_tex()
    out = [rf"{{\setlength{{\tabcolsep}}{{\vnsep{{{gid}}}}}\noindent\begin{{tabular}}{{{spec}}}\Xhline{{\vnbox}}"]
    pad = UNIFORM_PAD if lay.uniform else None
    if keyrow:
        out.append(prow(lay.ncol, inner, keyrow, pad=pad))
    if settings:
        out.append(prow(lay.ncol, inner, settings, pad=pad))
    cols = list(enumerate(lay.value_cols))
    nleft = nlev + 1
    gl, ng = group_header_line(lay, cols, nleft, True)
    hlines, hl = value_header_lines(lay, cols, lambda j: lay.w(f"v{j}"), False)
    if lay.align != "c":    # value headers stay centred over right-aligned value columns
        nv = len(cols)
        right = ["|"] * (nv - 1) + [r"!{\vrule width\vnbox}"]
        hlines = [[(rf"\multicolumn{{1}}{{c{right[n]}}}{{{c}}}" if c else c) for n, c in enumerate(cells)]
                  for cells in hlines]
    pad_check = ({m for m in lay.specs if any(r"\checkmark" in r["cells"][m][0] for r in rows)}
                 if (lay.align != "c" or lay.uniform) else set())   # uniform: centred cells of one width
    if lay.pad_all is not None:
        pad_check = lay.pad_all
    tot = hl + (1 if gl else 0)
    if gl:
        left = [(rf"\multirow{{{tot}}}{{*}}{{{h}}}" if h else "") for _, h in lay.level_cols]
        mh = lay.method_header
        left.append(rf"\multirow{{{tot}}}{{*}}{{{mh}}}" if mh else "")
        gparts = gl.split(" & ")
        out.append(" & ".join(left + gparts[nleft:]))
    for j, cells in enumerate(hlines):
        left = []
        if not gl:
            for _, h in lay.level_cols:
                left.append((rf"\multirow{{{hl}}}{{*}}{{{h}}}" if hl > 1 else h) if j == 0 else "")
            mh = lay.method_header
            left.append((rf"\multirow{{{hl}}}{{*}}{{{mh}}}" if hl > 1 else mh) if j == 0 else "")
        else:
            left = [""] * nleft
        out.append(" & ".join(left + cells) + r"\vnsd\\")
    out.append(r"\Xhline{0.6pt}")
    spans = level_spans(rows, nlev, L)
    for i, r in enumerate(rows):
        bh = r.get("block_head")
        if bh is not None and (i == 0 or rows[i - 1].get("block_head") != bh):
            # a full-width block row names the setting, its tuple on one line (never split)
            if i:
                out.append(rf"\Xhline{{{BLOCK_RULE}}}")
            out.append(block_row(lay.ncol, bh, pad=pad))
        elif i:
            rb = rule_between(lay, rows[i - 1], r, nlev, nlev + 1)
            if rb is not None:
                out.append(rb)
        for j in range(L):
            cells = []
            for k in range(nlev):
                if j == 0 and spans[k][i]:
                    cells.append(rf"\multirow{{{spans[k][i]}}}{{*}}{{{r['level_disp'][k]}}}" if spans[k][i] > 1
                                 else r["level_disp"][k])
                    lay.W.add(lay.level_cols[k][0], r["level_disp"][k])
                else:
                    cells.append("")
            if r.get("name_lines"):   # variant on line 1, (alpha0, alpha, delta, K) on line 2
                nl = r["name_lines"]
                cells.append(nl[j] if j < len(nl) else "")
                if j < len(nl):
                    lay.W.add("meth", nl[j])
            else:
                cells.append((rf"\multirow{{{L}}}{{*}}{{{r['name']}}}" if L > 1 else r["name"]) if j == 0 else "")
                if j == 0:
                    lay.W.add("meth", r["name"])
            cells += value_line_cells(lay, r, cols, j, L, fragment, target, table_id, lambda vj: f"v{vj}",
                                      pad_check)
            out.append(" & ".join(cells) + r"\\")
    out.append(r"\Xhline{\vnbox}\end{tabular}}\par")
    return "\n".join(out)


def baseline_table(lay, brows, settings, blevels, keyrow, fragment, target, table_id):
    """Boxed baseline table aligned with the main table (as in build_tables_v10_original.py). Columns that are N-A
    or not simulated for every baseline are dropped; a column keeps only its available members."""
    gid, L = lay.gid, lay.nlines
    keep = []
    for j, vc in enumerate(lay.value_cols):
        k = [m for m in vc if any(r["cells"][m][1].get("kind") in ("meanSD", "event", "event_pct_sweep", "tmeanSD", "count")
                                  for r in brows)]
        if k:
            keep.append((j, k))
    nb = len(blevels)
    left_main = " + ".join([lay.w(c) for c, _ in lay.level_cols] + [lay.w("meth")])
    nleft = len(lay.level_cols) + 1
    bl_widths = " - ".join([""] + [lay.w("B" + c) for c, _ in blevels]) if blevels else ""
    # the box spans the full \linewidth like the main table below it. The Method column takes the
    # left part of the main table plus every value column the box drops (width + 2 \tabcolsep + one rule each), so
    # the kept value columns stay aligned with the main table up to the first dropped one.
    dropped = [j for j in range(len(lay.value_cols)) if j not in {k for k, _ in keep}]
    fill = "".join(rf" + {lay.w(f'v{j}')} + 2\tabcolsep + \arrayrulewidth" for j in dropped)
    meth_w = rf"\dimexpr {left_main} + {2 * (nleft - 1 - nb)}\tabcolsep {bl_widths}{fill}\relax"
    lay.W.add("bname", "Method")
    cs = [rf"w{{c}}{{{lay.w('B' + c)}}}" for c, _ in blevels] + [rf"w{{c}}{{{meth_w}}}"]
    for cid, h in blevels:
        lay.W.add("B" + cid, h)
    for r in brows:
        lay.W.add("bname", r["name"])
    sep = r"!{\hskip\arrayrulewidth}"
    vs = [rf"w{{c}}{{{lay.w(f'v{j}')}}}" for j, _ in keep]
    spec = r"!{\vrule width\vnbox}" + "".join(cs) + sep + sep.join(vs) + r"!{\vrule width\vnbox}"
    ncol = nb + 1 + len(keep)
    # inner width of the box: all columns and separators minus the two outer paddings
    wsum = " + ".join([lay.w("B" + c) for c, _ in blevels] + [meth_w] + [lay.w(f"v{j}") for j, _ in keep])
    inner = rf"{wsum} + {2 * (ncol - 1)}\tabcolsep + {len(keep)}\arrayrulewidth - 0.01pt"
    # full width (\linewidth), centred like the main table
    out = [rf"{{\setlength{{\tabcolsep}}{{\vnsep{{{gid}}}}}\vncheckleft{{{meth_w}}}{{{lay.w('bname')}}}{{{table_id}}}"
           rf"\noindent\makebox[\linewidth][c]{{\begin{{tabular}}{{{spec}}}\Xhline{{\vnbox}}"]
    if keyrow:
        out.append(prow(ncol, inner, keyrow))
    out.append(prow(ncol, inner, settings))
    gl, _ = group_header_line(lay, keep, nb + 1, False)
    hlines, hl = value_header_lines(lay, keep, lambda j: lay.w(f"v{j}"), True)
    tot = hl + (1 if gl else 0)
    # a \multicolumn in the first column replaces its left box rule too; restate it (complete Method box)
    lb = "" if nb else r"!{\vrule width\vnbox}"
    meth_cell = lambda c: rf"\multicolumn{{1}}{{{lb}w{{c}}{{{meth_w}}}|}}{{{c}}}"
    if gl:
        gparts = gl.split(" & ")
        left = [(rf"\multirow{{{tot}}}{{*}}{{{h}}}") for _, h in blevels]
        left.append(meth_cell(rf"\multirow{{{tot}}}{{*}}{{Method}}"))
        out.append(" & ".join(left + gparts[nb + 1:]))
    for j, cells in enumerate(hlines):
        if gl:
            left = [""] * nb + [meth_cell("")]
        else:
            left = [(rf"\multirow{{{hl}}}{{*}}{{{h}}}" if hl > 1 else h) if j == 0 else "" for _, h in blevels]
            c = (r"\multirow{%d}{*}{Method}" % hl if hl > 1 else "Method") if j == 0 else ""
            left.append(meth_cell(c))
        out.append(" & ".join(left + cells) + r"\vnsd\\")
    out.append(r"\Xhline{0.6pt}")
    spans = level_spans(brows, nb, L) if nb else []

    def rule_after(i):
        if i == len(brows) - 1:
            return r"\Xhline{\vnbox}"
        if nb and brows[i + 1]["levels"] != brows[i]["levels"]:
            return r"\Xhline{0.6pt}"
        if L > 1:
            return rf"\Xcline{{{nb + 1}-{ncol}}}{{{RULE_METHOD}}}"
        return None

    for i, r in enumerate(brows):
        if i and rule_after(i - 1):
            out.append(rule_after(i - 1))
        for j in range(L):
            cells = []
            for k in range(nb):
                cells.append((rf"\multirow{{{spans[k][i]}}}{{*}}{{{r['level_disp'][k]}}}" if spans[k][i] > 1
                              else r["level_disp"][k]) if (j == 0 and spans[k][i]) else "")
                if j == 0 and spans[k][i]:
                    lay.W.add("B" + blevels[k][0], r["level_disp"][k])
            cells.append((rf"\multirow{{{L}}}{{*}}{{{r['name']}}}" if L > 1 else r["name"]) if j == 0 else "")
            cells += value_line_cells(lay, r, keep, j, L, fragment, target, table_id, lambda vj: f"v{vj}")
            pad = j == L - 1 and "_" in r["name"] and rule_after(i)
            out.append(" & ".join(cells) + (r"\vnsd\\" if pad else r"\\"))
    out.append(r"\Xhline{\vnbox}\end{tabular}}}\par")
    return "\n".join(out)


# ----------------------------------------------------------------------------------------------- uniform layout
# Every edition table follows the layout of topfive_full (report Tables 6 and 7): a settings row, one header, a
# full-width one-line block row per setting ("S3: V-L (alpha0, alpha, delta, K) = (...)", plus L_mu where it varies),
# rows with centred names (the method, or the swept value) and centred values, boxed, at the full line width. Balanced:
# the value columns share one width (their largest natural width) whenever the table then keeps a column separation of
# at least UNIFORM_MINSEP_PT, and the free width is spread evenly over all columns. No level columns, no rotated or
# stacked labels (exception: the Section VI short table keeps its one row per setting, the setting name in its own
# column). A baseline box keeps the grid of its table when it drops no value column (aligned); otherwise it gets its own
# balanced grid at the full line width (independent). Rows of one line get a larger \arraystretch so that a table fills
# its page. Long tables continue on the next page with \ntablebreak ("Table N (continued)."); every page block fits the
# page (compile_tests.py checks the page fill).
UNIFORM_MINSEP_PT = 3.0


def eq_tex(cols, fixed, nrules, ncol, tag):
    """TeX that sets the width macros `cols` to their largest width when the table then keeps a column
    separation of at least UNIFORM_MINSEP_PT at \\linewidth (`fixed`: the other width macros of the row, `nrules`: rules
    between the columns, `ncol`: columns). Otherwise the natural widths stay."""
    m = f"{tag}eqv"
    out = ["% equal value-column widths when they fit", rf"\vnnew{{{m}}}"]
    out += [rf"\ifdim\vnw{{{c}}}>\vnw{{{m}}}\global\vnw{{{m}}}=\vnw{{{c}}}\fi" for c in cols]
    fx = "".join(rf"-\vnw{{{c}}}" for c in fixed)
    out.append(rf"\ifdim\dimexpr\linewidth{fx}-{len(cols)}\vnw{{{m}}}-{nrules}\arrayrulewidth-2\vnbox"
               rf"-{2 * ncol * UNIFORM_MINSEP_PT:g}pt\relax>0pt" +
               "".join(rf"\global\vnw{{{c}}}=\vnw{{{m}}}" for c in cols) + r"\fi")
    return "\n".join(out) + "\n"


def baseline_box_u(lay, brows, settings, fragment, target, table_id):
    """Baseline box of the uniform layout. Value columns that are N-A or not simulated for every baseline are dropped. Without a
    dropped column the box keeps the grid of the main table (aligned: same first column, widths and separation);
    with one, the box gets its own balanced grid at \\linewidth (independent: Method column at its natural width,
    the kept value columns at one common width, its own separation). Block rows (rows with 'block_head', e.g.
    L_mu = 0.0005) as in the main table. Header rules and the outer box as in baseline_table."""
    gid, L = lay.gid, lay.nlines
    keep = []
    for j, vc in enumerate(lay.value_cols):
        k = [m for m in vc if any(r["cells"][m][1].get("kind") in ("meanSD", "event", "event_pct_sweep", "tmeanSD", "count")
                                  for r in brows)]
        if k:
            keep.append((j, k))
    aligned = len(keep) == len(lay.value_cols) and not lay.level_cols
    if aligned:
        for b in ["Method"] + [r["name"] for r in brows]:
            lay.W.add("meth", b)
        mw, colw, wcol, sepn = lay.w("meth"), (lambda j: lay.w(f"v{j}")), (lambda vj: f"v{vj}"), gid
    else:
        bt = "B" + "".join(str(j) for j, _ in keep)   # boxes with the same kept columns share one grid
        tagb = f"{gid}{bt}"
        lay.W.add(f"{bt}m", "Method")
        for r in brows:
            lay.W.add(f"{bt}m", r["name"])
        for j, k in keep:
            for m in k:
                lay.W.add(f"{bt}v{j}", head(m))
        mw, colw, wcol, sepn = lay.w(f"{bt}m"), (lambda j: lay.w(f"{bt}v{j}")), (lambda vj: f"{bt}v{vj}"), tagb
        ncol_b = 1 + len(keep)
        vnames = [lay.W.name(f"{bt}v{j}") for j, _ in keep]
        widths = " - ".join([mw] + [colw(j) for j, _ in keep])
        hd = (f"% independent baseline box {tagb}\n"
                         + eq_tex(vnames, [lay.W.name(f"{bt}m")], len(keep), ncol_b, tagb)
                         + rf"\vnnewsep{{{tagb}}}{{\dimexpr(\linewidth - {widths} - {len(keep)}\arrayrulewidth"
                           rf" - 2\vnbox)/{2 * ncol_b}\relax}}" + "\n")
        if hd not in lay.heads:
            lay.heads.append(hd)
    sep = r"!{\hskip\arrayrulewidth}"
    vs = [rf"w{{c}}{{{colw(j)}}}" for j, _ in keep]
    spec = r"!{\vrule width\vnbox}" + rf"w{{c}}{{{mw}}}" + sep + sep.join(vs) + r"!{\vrule width\vnbox}"
    ncol = 1 + len(keep)
    wsum = " + ".join([mw] + [colw(j) for j, _ in keep])
    inner = rf"{wsum} + {2 * (ncol - 1)}\tabcolsep + {len(keep)}\arrayrulewidth - 0.01pt"
    out = [rf"{{\setlength{{\tabcolsep}}{{\vnsep{{{sepn}}}}}"
           rf"\noindent\makebox[\linewidth][c]{{\begin{{tabular}}{{{spec}}}\Xhline{{\vnbox}}"]
    out.append(prow(ncol, inner, settings, pad=UNIFORM_PAD))
    gl, _ = group_header_line(lay, keep, 1, False)
    hlines, hl = value_header_lines(lay, keep, colw, True)
    tot = hl + (1 if gl else 0)
    meth_cell = lambda c: rf"\multicolumn{{1}}{{!{{\vrule width\vnbox}}w{{c}}{{{mw}}}|}}{{{c}}}"
    if gl:
        gparts = gl.split(" & ")
        out.append(" & ".join([meth_cell(rf"\multirow{{{tot}}}{{*}}{{Method}}")] + gparts[1:]))
    for j, cells in enumerate(hlines):
        if gl:
            left = [meth_cell("")]
        else:
            left = [meth_cell((r"\multirow{%d}{*}{Method}" % hl if hl > 1 else "Method") if j == 0 else "")]
        out.append(" & ".join(left + cells) + r"\vnsd\\")
    out.append(r"\Xhline{0.6pt}")

    def rule_after(i):
        if i == len(brows) - 1:
            return r"\Xhline{\vnbox}"
        if L > 1 and brows[i + 1].get("block_head") == brows[i].get("block_head"):
            return rf"\Xcline{{1-{ncol}}}{{{RULE_METHOD}}}"
        return None

    for i, r in enumerate(brows):
        bh = r.get("block_head")
        if bh is not None and (i == 0 or brows[i - 1].get("block_head") != bh):
            if i:
                out.append(rf"\Xhline{{{BLOCK_RULE}}}")
            out.append(block_row(ncol, bh, pad=UNIFORM_PAD))
        elif i and rule_after(i - 1):
            out.append(rule_after(i - 1))
        for j in range(L):
            cells = [(rf"\multirow{{{L}}}{{*}}{{{r['name']}}}" if L > 1 else r["name"]) if j == 0 else ""]
            cells += value_line_cells(lay, r, keep, j, L, fragment, target, table_id, wcol)
            pad = j == L - 1 and "_" in r["name"] and rule_after(i)
            out.append(" & ".join(cells) + (r"\vnsd\\" if pad else r"\\"))
    out.append(r"\Xhline{\vnbox}\end{tabular}}}\par")
    return "\n".join(out)


def uniform_stretch(target, lines, fixed_pt, fill=0.9, lo=0.93, hi=1.35):
    """\\arraystretch of a one-line-row table whose largest page block has `lines` table rows (each
    \\arraystretch x \\baselineskip tall) and `fixed_pt` of other material (caption, rules, gaps): the stretch that
    fills `fill` of the text height, within [lo, hi], rounded down to 0.05."""
    H, bs = (696.0, 12.0) if target == "paper" else (686.565, 13.192)
    s = (fill * H - fixed_pt) / (lines * bs)
    return max(lo, min(hi, math.floor(s * 20) / 20))


def uniform_pages(groups_per_page, box_keys=None, same_sheet=()):
    """Explicit page plan. groups_per_page: list of (key, rows) per page block; box_keys: keys whose first page
    block carries the baseline box (default: every new key); same_sheet: indices of page blocks that continue the
    sheet of the previous block (joined by \\vspace instead of \\ntablebreak)."""
    pages, seen = [], set()
    for p, (k, rows) in enumerate(groups_per_page):
        box = (k not in seen) if box_keys is None else (k in box_keys and k not in seen)
        seen.add(k)
        pages.append(dict(key=k, rows=list(rows), box=box, first=p == 0, same_sheet=p in same_sheet))
    return pages


# ----------------------------------------------------------------------------------------------- row builders
def settings_line(scen, tag=None, with_scen_consts=True, **kv):
    """Bold scenario name (or a given title), optional bold tag ("Baselines"), then the settings in math."""
    order = ["alpha0", "alpha", "delta", "K", "b", "phi", "L_mu", "cmax", "Ts", "NMC"]
    sym = dict(alpha0=r"\alpha_0", alpha=r"\alpha", delta=r"\delta", K="K", L_mu=r"L_\mu", cmax=r"c_{\max}",
               Ts="T_s", NMC=r"N_{\rm MC}", b="b", phi=r"\phi")
    kv = dict(kv)
    if is_default_b(kv.get("b")):   # b = kappa by default (arXiv App. F-B); b is printed only where compared
        kv.pop("b")
    if scen in SCEN and with_scen_consts:
        kv.setdefault("cmax", SCEN[scen]["cmax"] + r"\,\mathrm{m}")
        kv.setdefault("Ts", SCEN[scen]["Ts"] + r"\,\mathrm{s}")
    kv.setdefault("NMC", "200")
    parts = []
    for k in order:
        if k in kv and kv[k] is not None:
            v = kv[k]
            parts.append(f"{sym[k]}={v if isinstance(v, str) else (kfmt(v) if k == 'K' else num(v))}")
    title = SCEN[scen]["name"] if scen in SCEN else scen
    name = rf"\textbf{{{title}}}" + (rf"\quad\textbf{{{tag}}}" if tag else "")
    return rf"{name}\quad ${',\\ '.join(parts)}$"


def baseline_found(D, scen, method, alpha0, L_mu, temporal=False):
    src = "final_temporal" if temporal else "final"
    f = D.get(src, "matched_baselines_200", benchmark=scen, method=method, alpha0=alpha0, L_mu=L_mu, rho_buf=0.0)
    if f is None and method == "ACI" and abs(L_mu - 0.0005) < 1e-9:
        f = D.get(src, "matched_ACI_200", benchmark=scen, method=method, alpha0=alpha0, L_mu=L_mu, rho_buf=0.0)
    return f


def baseline_rows(D, scens, alpha0, L_mu, level=None, temporal=False):
    """scens: scenario or {sub: scenario} for multi-scenario rows."""
    rows = []
    for b in BASELINES:
        if isinstance(scens, dict):
            found = {sub: baseline_found(D, s, b, alpha0, L_mu, temporal) for sub, s in scens.items()}
            sc = "+".join(scens.values())
        else:
            found = {"": baseline_found(D, scens, b, alpha0, L_mu, temporal)}
            sc = scens
        rows.append(dict(found=found, is_variant=False, group=None, name=BNAME[b], method=b,
                         levels=level[0] if level else (), level_disp=level[1] if level else (),
                         record=dict(scenario=sc, alpha0=alpha0, L_mu=L_mu, baseline=True)))
    return rows


def vrow(found, v, levels, disp, record, group, is_variant=True):
    return dict(found=found, is_variant=is_variant, group=group, name=VNAME[v], method=v, levels=levels,
                level_disp=disp, record=record)


def cell_or_fail(D, src, fam, **cfg):
    f = D.get(src, fam, **cfg)
    if f is None:
        raise KeyError((src, fam, cfg))
    return f


METHOD_HEAD = r"Method $(\alpha_0,\alpha,\delta,K)$"


def tuple_plain(t, fields=FIELDS4):
    """The setting as '(alpha0, alpha, delta, K)', e.g. (0.30, 0.25, 0.20, 50), or the given fields only
    (a swept field is left out)."""
    vals = [(r"$\infty$" if t[f] == "infinity" else kfmt(t[f])) if f == "K" else num(t[f]) for f in fields]
    return "(" + ", ".join(vals) + ")"


def set_label(r, t, L):
    """Row label 'Method (alpha0, alpha, delta, K)': one line, or variant over tuple when a row has L >= 2 lines."""
    tp = tuple_plain(t)
    r["record"]["label"] = tp
    if L == 1:
        r["name"] = f"{VNAME[r['method']]} {tp}"
    else:
        r["name_lines"] = [VNAME[r["method"]], tp] + [""] * (L - 2)
    return r


def paired(si):
    """The variant named with the setting (Default V-C, S1 V-C, S2 V-C, S3 V-L, S4 V-Delta, S5 V-C)."""
    return SIX[si][2]


# ----------------------------------------------------------------------------------------------- fragment shell
def fragment_head(lay):
    nock = (r"\normalsize" + "\n") if lay.normalsize else ""   # body size before any width is measured
    nock += (r"\providecommand{\vnnock}{\phantom{\,\checkmark}}" + "\n"
             if lay.align != "c" or lay.uniform else "")   # right-aligned layout; also the centred uniform layout
    if lay.stretch:
        nock += rf"\renewcommand{{\arraystretch}}{{{lay.stretch}}}" + "\n"
    eq = ""
    if lay.uniform:   # equal value-column widths (balanced), when the table still fits
        fixed = [lay.W.name(c) for c, _ in lay.level_cols] + [lay.W.name("meth")]
        eq = eq_tex([lay.W.name(f"v{j}") for j in range(len(lay.value_cols))], fixed, len(lay.value_cols),
                    lay.ncol, lay.gid)
    return (f"% widths and column separation for layout {lay.gid}\n{nock}{lay.W.tex()}\n{eq}{lay.sep_tex()}\n"
            + "".join(lay.heads))


FIRST_CUT = 0   # catalog mode only: value lines kept free on a fragment's first page for headings and caption


def paginate(groups, budget, has_baselines, box_each_page, same_sheet=True):
    """Greedy packing of (key, nlines, rows) groups into pages. A page holds groups of one baseline key.
    budget = (first page, later page with a baseline box, later page without one), in value lines.
    The paper repeats the baseline box on every page; the report shows it once per key."""
    # No half-empty pages: a new baseline key starts on the SAME sheet, with its own baseline box,
    # when its first group and the box fit in the space left; otherwise on a new sheet.
    BOX_LINES, SLACK = 10, 9
    pages, cur = [], None
    for bkey, n, rows in groups:
        if (same_sheet and cur is not None and bkey != cur["key"]
                and cur["sheet_used"] + (BOX_LINES if has_baselines else 2) + n + SLACK <= cur["sheet_cap"]):
            used0 = cur["sheet_used"] + (BOX_LINES if has_baselines else 2)
            cur = dict(key=bkey, rows=[], used=0, box=has_baselines, first=False, cap=cur["sheet_cap"],
                       same_sheet=True, sheet_used=used0, sheet_cap=cur["sheet_cap"])
            pages.append(cur)
        elif cur is None or bkey != cur["key"] or cur["sheet_used"] + n > cur["sheet_cap"]:
            first = cur is None
            box = has_baselines and (box_each_page or first or bkey != cur["key"])
            cap = budget[0] - FIRST_CUT if first else (budget[1] if box else budget[2])
            cur = dict(key=bkey, rows=[], used=0, box=box, first=first, cap=cap, same_sheet=False,
                       sheet_used=0, sheet_cap=cap)
            pages.append(cur)
        cur["rows"] += rows
        cur["used"] += n
        cur["sheet_used"] += n
        assert cur["sheet_used"] <= cur["sheet_cap"] or len(pages[-1]["rows"]) == len(rows), "group larger than a page"
    return pages


def paginate_rows(groups, first_cap, cap_box, cap_nobox, L, min_rows=1):
    """Report delta_corridor and period_corridor: pack the rows of the setting blocks into pages, a block splitting
    at a row boundary (its block row is repeated on the next page), a new baseline key starting on the same sheet with
    its box when at least min_rows rows fit behind the box. Caps are in value lines (L per row)."""
    BOXC, HDRC = cap_nobox - cap_box, 6   # HDRC: header lines of the table that starts behind the box
    pages, cur = [], None
    for bkey, n, rows in groups:
        rows = list(rows)
        while rows:
            free = 0
            if cur is not None and bkey == cur["key"]:
                free = cur["cap"] - cur["used"]
            elif cur is not None and (cur["cap"] - cur["used"] - BOXC - HDRC) >= min_rows * L:
                free = cur["cap"] - cur["used"] - BOXC - HDRC
                cur = dict(key=bkey, rows=[], used=0, box=True, first=False, cap=free, same_sheet=True,
                           sheet_used=0, sheet_cap=free)
                pages.append(cur)
            if free // L < 1:
                first = cur is None
                box = first or bkey != cur["key"]
                cap = first_cap if first else (cap_box if box else cap_nobox)
                cur = dict(key=bkey, rows=[], used=0, box=box, first=first, cap=cap, same_sheet=False,
                           sheet_used=0, sheet_cap=cap)
                pages.append(cur)
                free = cap
            k = min(len(rows), free // L)
            cur["rows"] += rows[:k]
            cur["used"] += k * L
            cur["sheet_used"] += k * L
            rows = rows[k:]
    return pages


def balance_pages(groups, cap):
    """Pages per baseline key with the fewest pages that fit cap, the groups spread as evenly as possible."""
    keys = []
    for g in groups:
        if not keys or keys[-1][0] != g[0]:
            keys.append((g[0], []))
        keys[-1][1].append(g)
    pages = []
    for key, gs in keys:
        size = max(n for _, n, _ in gs)
        per = max(1, cap // size)
        m = math.ceil(len(gs) / per)
        q, r = divmod(len(gs), m)
        i = 0
        for j in range(m):
            take = q + (1 if j < r else 0)
            chunk = gs[i:i + take]
            i += take
            pages.append(dict(key=key, rows=[row for _, _, rows in chunk for row in rows],
                              used=sum(n for _, n, _ in chunk), box=(j == 0), first=not pages, cap=cap,
                              same_sheet=False, sheet_used=sum(n for _, n, _ in chunk), sheet_cap=cap))
    return pages


def key_box(text):
    """Full-width boxed key above the first box of a fragment (notation and number formats)."""
    return (r"%KEYBOX" "\n" r"{\setlength{\tabcolsep}{4pt}\noindent\begin{tabular}{!{\vrule width\vnbox}"
            r"p{\dimexpr\linewidth-2\vnbox-2\tabcolsep\relax}!{\vrule width\vnbox}}\Xhline{\vnbox}"
            rf"\raggedright {text}\tabularnewline\Xhline{{\vnbox}}\end{{tabular}}}}\par\vspace{{2pt}}")


def emit(lay, pages, frag, target, key, main_sets, base_fn=None, blevels=(), box_alone=False, box_alone_first=False,
         split_first_box=False):
    """Page blocks: [baseline box] + main table; key row in the first box of the fragment.
    box_alone: the baseline box takes a page of its own (report L_mu table, ten two-L_mu baselines)."""
    blocks, joins = [], []
    if lay.uniform:   # \vnnock padding decided over the whole fragment; balanced baseline boxes
        lay.pad_all = {m for m in lay.specs if any(r"\checkmark" in r["cells"][m][0] for pg in pages for r in pg["rows"])}
        baseline_table = lambda lay_, brows, bsets, blevels, keyrow, frag_, target_, tid: baseline_box_u(
            lay_, brows, bsets, frag_, target_, tid)
    else:
        baseline_table = globals()["baseline_table"]
    gap = r"\nopagebreak\vspace{5pt}" if lay.uniform else r"\vspace{5pt}"   # uniform: box and rows on one page
    for p, pg in enumerate(pages):
        joins.append(r"\vspace{6pt}" if (pg.get("same_sheet") and not box_alone)
                     or (box_alone_first and p > 0 and pg["box"]) else r"\ntablebreak")
        parts = [key_box(key)] if (p == 0 and key) else []
        keyrow = None
        if pg["box"] and base_fn is not None:
            brows, bsets = base_fn(pg["key"])
            if split_first_box and p == 0:   # first half of the box under the caption, second half
                h = len(brows) // 2          # above the first rows (report phi, Corridor)
                blocks.append("\n".join(parts + [baseline_table(lay, brows[:h], bsets, list(blevels), keyrow, frag,
                                                                target, f"{lay.gid}-base-p{p + 1}a")]))
                joins.append(r"\ntablebreak")
                parts = [baseline_table(lay, brows[h:], bsets, list(blevels), keyrow, frag, target,
                                        f"{lay.gid}-base-p{p + 1}b"), r"\vspace{5pt}",
                         main_table(lay, pg["rows"], main_sets(pg), keyrow, frag, target, f"{lay.gid}-p{p + 1}")]
                blocks.append("\n".join(parts))
                continue
            parts += [baseline_table(lay, brows, bsets, list(blevels), keyrow, frag, target,
                                     f"{lay.gid}-base-p{p + 1}"), gap]
            keyrow = None
            if box_alone or box_alone_first:   # report delta: each box closes the page above its rows
                blocks.append("\n".join(parts[:-1]))
                joins.append(r"\ntablebreak")
                parts = []
        if pg["rows"]:   # a first page may hold the baseline box only (report period_corridor)
            parts.append(main_table(lay, pg["rows"], main_sets(pg), keyrow, frag, target, f"{lay.gid}-p{p + 1}"))
        else:
            parts = parts[:-1] if parts and parts[-1] == gap else parts
        blocks.append("\n".join(parts))
    joins = joins[1:len(blocks)]
    if lay.uniform:   # a break may open a new group (a float of its own in the extended version): restate size and stretch
        setup = ((r"\normalsize" + "\n") if lay.normalsize else "") + r"\providecommand{\vnnock}{\phantom{\,\checkmark}}"
        if lay.stretch:
            setup += "\n" + rf"\renewcommand{{\arraystretch}}{{{lay.stretch}}}"
        joins = [(j + "\n" + setup) if j == r"\ntablebreak" else j for j in joins]
    out = blocks[0]
    for j, blk in zip(joins, blocks[1:]):
        out += "\n" + j + "\n" + blk
    return fragment_head(lay) + out + "\n"


def key_msd(rows_all, brows_all=()):
    k = KEY_MSD
    kinds = {i.get("kind") for r in list(rows_all) + list(brows_all) for _, i in r["cells"].values()}
    if "NA" in kinds:
        k += KEY_NA
    if "not_simulated" in kinds:
        k += KEY_DASH
    return k + "."


def key_events(rows_all, check=True, na=True):
    k = KEY_EVT if check else KEY_EVT_NOCHECK
    k += KEY_QNOPLUS
    if any(i.get("kind") == "event" and i["n"] != NMC for r in rows_all for _, i in r["cells"].values()):
        k += KEY_ASSESSED
    kinds = {i.get("kind") for r in rows_all for _, i in r["cells"].values()}
    if "NA" in kinds and na:
        k += KEY_NA
    return k + "."


# value-line budgets (first page, later page with box, later page without box)
# report delta_corridor and period_corridor: value lines on the first page, on a page with a box, on one without
REPORT_SPLIT = (15, 24, 36)
B_MULTI = {"paper": (40, 40, 40), "report": (24, 24, 36)}


def boxes(target):
    return target == "paper"


# ----------------------------------------------------------------------------------------------- fragments
def tag3(a0):
    return f"{int(round(a0 * 100)):03d}"


def six_cfg(t, **over):
    return dict(dict(t, **BASE), **over)


def anchor_cell(D, si, scen, v, fam, **over):
    """Cell of a sweep family at one of the six settings (final old-anchor export for the default, new_anchors
    for S1-S5)."""
    name, t, _, (src, aid) = SIX[si]
    cfg = six_cfg(t, **over)
    f = cell_or_fail(D, src, fam, benchmark=scen, method=v, rho_buf=0.0, **cfg)
    want = aid if src == "final" else D.anchor(src, aid)
    assert f["anchor_id"] == want, (f["anchor_id"], want)
    return f, cfg


def sweep_fragment(D, scen, alpha0, target, gid):
    cols, L = full_cols(target)
    lay = Layout(gid, [("a", r"$\alpha$"), ("d", r"$\delta$"), ("k", r"$K$")], cols, L, ["1.6pt", "1.0pt", "0.6pt"])
    alphas = [x / 100 for x in range(5, int(round(alpha0 * 100)), 5)]
    groups = []
    for a in alphas:
        for d in DELTAS:
            for K in KS:
                lv = (num(a), num(d), kfmt(K))
                disp = (rf"\rotatebox[origin=c]{{90}}{{$\alpha={lv[0]}$}}", rf"\rotatebox[origin=c]{{90}}{{$\delta={lv[1]}$}}",
                        f"${lv[2]}$")
                cfg = dict(alpha0=alpha0, alpha=a, delta=d, K=K, **BASE)
                rows = [vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)},
                             v, lv, disp, dict(scenario=scen, table="sweep", **cfg), (scen, a, d, K)) for v in VARIANTS]
                groups.append((0, len(rows) * L, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    brows = baseline_rows(D, scen, alpha0, 0.0005)
    fill_rows(brows, lay.specs)
    frag = f"sweep_{scen}_a{tag3(alpha0)}"
    pages = paginate(groups, B_MULTI[target], True, boxes(target))
    return emit(lay, pages, frag, target, key_msd(allrows, brows),
                lambda pg: settings_line(scen, alpha0=alpha0, b="2.11", phi="1", L_mu="0.0005"),
                lambda k: (brows, settings_line(scen, tag="Baselines", alpha0=alpha0, L_mu="0.0005")))


def param_fragment(D, scen, kind, target, gid):
    """kind: 'delta' / 'period' (main grid), 'b', 'phi' (phi x L_mu; anchor sweeps); chosen variant only.
    Paper: Default and S3 (phi: all six settings); report: all six settings. A full-width block row per setting
    names the setting, its chosen variant and its tuple (alpha0, alpha, delta, K) on one line; the rows below give the
    swept value (delta, K, b, or phi) in the first column, so no tuple is split over two lines.
    Uniform layout; phi: one block row per (setting, L_mu), e.g.
    'S3: V-L (alpha0, alpha, delta, K) = (0.30, 0.20, 0.20, 50), L_mu = 0.0005', rows phi (no rotated L_mu column);
    the baseline box gives both L_mu as block rows; paper phi: one page per setting, the box on the first page of each
    alpha0."""
    cols, L = full_cols(target)
    paper = target == "paper"
    sel = [0, 3] if paper and kind != "phi" else list(range(6))   # phi: all six settings in both editions
    sym = {"delta": r"$\delta$", "period": r"$K$", "b": r"$b$", "phi": r"$\phi$"}[kind]
    lay = Layout(gid, [], cols, L, [], method_header=sym, uniform=True)
    if kind == "delta" and 1 in sel and 2 in sel:   # S1 and S2 share (alpha0, alpha, K) and V-C: one block
        merged = [[m] for m in sel if m not in (1, 2)]
        merged.insert(1, [1, 2])
    else:
        merged = [[m] for m in sel]
    bf = tuple(f for f in FIELDS4 if f != {"delta": "delta", "period": "K"}.get(kind))   # no swept field
    groups = []
    for members in merged:
        si = members[0]
        name, t, _, _ = SIX[si]
        v = paired(si)
        assert all(paired(m) == v for m in members)
        if len(members) > 1:   # S1 and S2 coincide once delta is dropped: one name pair, one tuple
            assert len({tuple_plain(SIX[m][1], bf) for m in members}) == 1
            head_txt = (" and ".join(rf"\textbf{{{dn(SIX[m][0])}}}" for m in members)
                        + rf": {VNAME[v]} {fields_sym(bf)} = {tuple_plain(t, bf)}")
        else:
            head_txt = block_text(name, t, v, fields=bf)
        if kind == "phi":
            vals = [(lm, p) for lm in (0.0005, 0.005) for p in (0.0, 0.25, 0.5, 0.75, 1.0)]
        else:
            vals = {"delta": DELTAS, "period": KS, "b": [2.11, 2.5, 3.0, 4.0]}[kind]
        rows = []
        for val in vals:
            bh = head_txt
            if kind == "delta":
                cfg = six_cfg(t, delta=val)
                f = cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
                vname = f"${num(val)}$"
            elif kind == "period":
                cfg = six_cfg(t, K=val)
                f = cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
                vname = f"${kfmt(val)}$"
            elif kind == "b":
                f, cfg = anchor_cell(D, si, scen, v, "budget_normalizer_200", budget_normalizer=val)
                vname = f"${b_value(val)}$"   # the default b prints as kappa
            else:
                lm, p = val
                f, cfg = anchor_cell(D, si, scen, v, "stress_drift_200", bound_scale=p, L_mu=lm)
                bh = head_txt + (r", $L_\mu=0.0005$" if lm == 0.0005 else r", $L_\mu=0.005$")
                vname = f"${p:g}$"
            r = vrow({"": f}, v, (), (), dict(scenario=scen, table=kind, setting="+".join(SIX[m][0] for m in members),
                                              **cfg), (scen, si, str(val)))
            tt = dict(alpha0=cfg["alpha0"], alpha=cfg["alpha"], delta=cfg["delta"], K=cfg["K"])
            r["record"]["label"] = tuple_plain(tt)
            r["record"]["block"] = [tuple_plain(SIX[m][1], bf) for m in members]   # tuples of the block row
            r["record"]["block_fields"] = bf   # the swept field is not in the block row
            r["record"]["swept"] = kind
            r["name"] = vname
            r["block_head"] = bh
            rows.append(r)
        if kind == "phi":   # one group per (setting, L_mu) block: its block row and five phi rows
            groups += [((t["alpha0"],), 5 * L + 1, rows[i:i + 5]) for i in (0, 5)]
        else:
            groups.append(((t["alpha0"],), len(rows) * L, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    frag = f"{kind}_{scen}"
    budget = {("phi", "paper"): (32, 32, 32), ("phi", "report"): (48, 48, 48),
              ("delta", "report"): (36, 36, 36) if scen == "corridor" else B_MULTI[target]}.get((kind, target), B_MULTI[target])   # report delta Corridor: first box alone
    if paper and kind != "phi":   # paper parameter tables (not in an edition): one page per alpha0, no overflow
        budget = (24, 24, 24)
    pages = paginate(groups, budget, True, boxes(target), same_sheet=kind != "phi" and not paper)
    split_rows = kind in ("delta", "period") and not paper and scen == "corridor"
    if split_rows:   # report delta_corridor and period_corridor split at rows so that every page but the last is full
        pages = paginate_rows(groups, REPORT_SPLIT[0] if kind == "delta" else 0, REPORT_SPLIT[1], REPORT_SPLIT[2], L)
    if kind == "phi" and paper:   # one page per setting (both L_mu blocks); the box on the first page of its alpha0
        blocks = [((SIX[m[0]][1]["alpha0"],), [r for g in groups if g[2][0]["record"]["setting"] == SIX[m[0]][0]
                                               for r in g[2]]) for m in merged]
        packed = []   # the first setting of each alpha0 shares its page with the box; the others go two per page
        for k, rows_ in blocks:
            if packed and packed[-1][0] == k and packed[-1][2] == 1 and not packed[-1][3]:
                packed[-1][1].extend(rows_)
                packed[-1][2] = 2
            else:
                packed.append([k, list(rows_), 1, not any(q[0] == k for q in packed)])
        pages = uniform_pages([(k, rows_) for k, rows_, _, _ in packed])
    if kind == "phi" and not paper:   # spread each key's blocks evenly (3,3,2,2, not 3,3,3,1)
        pages = balance_pages(groups, budget[0])
    bcache = {}

    def base_fn(k):
        a0 = k[0]
        if k not in bcache:
            if kind == "phi":   # both L_mu in one box, each a block row
                brows = []
                for lm, disp in ((0.0005, "0.0005"), (0.005, "0.005")):
                    rs = baseline_rows(D, scen, a0, lm)
                    for r in rs:
                        r["block_head"] = rf"$L_\mu={disp}$"
                    brows += rs
                bsets = settings_line(scen, tag="Baselines", alpha0=a0)
            else:
                brows = baseline_rows(D, scen, a0, 0.0005)
                bsets = settings_line(scen, tag="Baselines", alpha0=a0, L_mu="0.0005")
            fill_rows(brows, lay.specs)
            bcache[k] = (brows, bsets)
        return bcache[k]

    fixed = {"delta": dict(b="2.11", phi="1", L_mu="0.0005"), "period": dict(b="2.11", phi="1", L_mu="0.0005"),
             "b": dict(phi="1", L_mu="0.0005"), "phi": dict(b="2.11")}[kind]

    def sets(pg):   # the settings are named in the block rows
        return settings_line(scen, **fixed)

    # report phi: the ten two-L_mu baselines (30 value lines) take a page of their own
    return emit(lay, pages, frag, target, None, sets, base_fn, (), box_alone=(kind == "phi" and not paper),
                box_alone_first=(kind == "delta" and not paper and scen == "corridor" and not split_rows),
                split_first_box=(kind == "phi" and paper and scen == "corridor"))
    # box_alone_first would overflow for report phi Corridor (2 blocks + half box). Paper phi Corridor: the first page
    # block is the L_mu = 0.0005 half of the baseline box (it starts under the text of its subsection); the
    # L_mu = 0.005 half and the default setting follow after \ntablebreak.


# the columns of topfive_full (report Tables 6 and 7), three lines per method
TOP_COLS = TRIPLES_R + [["qdep", "qcert"]]


def default_fragment(D, scen, target, gid):
    """The default setting in the layout of topfive_full (report Tables 4 and 5): settings row, block row
    'default setting: (alpha0, alpha, delta, K) = (0.15, 0.10, 0.05, 50)', the four variants with three value lines
    each (TOP_COLS, with q_dep and q_cert), the baselines at alpha0 = 0.15 in a box above."""
    lay = Layout(gid, [], TOP_COLS, 3, [], method_header="Method", uniform=True)
    t = SIX[0][1]
    cfg = six_cfg(t)
    rows = []
    for v in VARIANTS:
        r = vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)}, v, (),
                 (), dict(scenario=scen, table="default", setting="Default", **cfg), (scen, "default"))
        r["record"]["label"] = tuple_plain(t)
        r["record"]["block"] = [tuple_plain(t)]
        r["block_head"] = block_text("Default", t)
        rows.append(r)
    fill_rows(rows, lay.specs)
    brows = baseline_rows(D, scen, 0.15, 0.0005)
    fill_rows(brows, lay.specs)
    apply_submission_ranges(brows, scen, 0.15)   # the submission's range cell (SUBMISSION_RANGES), all editions
    pages = uniform_pages([((0.15,), rows)])
    return emit(lay, pages, f"default_{scen}", target, None,
                lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"),
                lambda k: (brows, settings_line(scen, tag="Baselines", alpha0=0.15, L_mu="0.0005")))


def topfive_fragment(D, scen, target, gid, with_default=True):
    """Catalog: S1-S5 (+ the default), all four variants, all metrics."""
    cols, L = full_cols(target)
    lay = Layout(gid, [], cols, L, [])
    groups = []
    for si in range(0 if with_default else 1, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        rows = [vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)}, v,
                     (), (), dict(scenario=scen, table="topfive", setting=name, **cfg), (scen, si))
                for v in VARIANTS]
        for r in rows:   # block row per setting, tuple on one line
            r["block_head"] = block_text(name, t)
        groups.append(((t["alpha0"],), len(rows) * L + 2, rows))   # + block row
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    pages = paginate(groups, B_MULTI[target], True, boxes(target))
    bc = {}

    def base_fn(k):
        if k not in bc:
            br = baseline_rows(D, scen, k[0], 0.0005)
            fill_rows(br, lay.specs)
            bc[k] = (br, settings_line(scen, tag="Baselines", alpha0=k[0], L_mu="0.0005"))
        return bc[k]

    allb = [r for pg in pages if pg["box"] for r in base_fn(pg["key"])[0]]
    return emit(lay, pages, f"topfive_{scen}", target, key_msd(allrows, allb),
                lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"), base_fn)


def event_percent_rows(rows):
    """Display trajectory events in percent, keeping the original confidence limits and flags."""
    for row in rows:
        for spec, (tex, info) in row["cells"].items():
            if info["kind"] == "event":
                pct = half_up(100 * info["count"] / NMC, 1)
                tex = f"{pct} ({info['upper_str']})" + (r"\,\checkmark" if info["check"] else "")
                row["cells"][spec] = [tex, dict(info, kind="event_pct_sweep", pct_str=pct)]


def event_percent_headers(tex):
    for spec in EVT:
        tex = tex.replace(head(spec), head(spec) + r"\,(\%)")
    return tex


def topfive_full_fragment(D, scen, target, gid, percent=False):
    """report Tables 6 and 7 (and the paper fragments of the same name): S1-S5, all four variants, all metrics of
    build_tables_v10_original.py plus (q_dep, q_cert); baselines at alpha0 = 0.30 once; block row
    'S1: (alpha0, alpha, delta, K) = (...)'. The reference of the uniform layout."""
    lay = Layout(gid, [], TOP_COLS, 3, [], method_header="Method", uniform=True)
    L = 3
    groups = []
    for si in range(1, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        rows = []
        for v in VARIANTS:   # block row 'S1: (alpha0, alpha, delta, K) = (...)', one line per setting
            r = vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)},
                     v, (), (), dict(scenario=scen, table="topfive_full", setting=name, **cfg), (scen, si))
            r["record"]["label"] = tuple_plain(t)
            r["block_head"] = block_text(name, t)
            rows.append(r)
        groups.append(((t["alpha0"],), len(rows) * L, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    budget = {"paper": (24, 36, 36), "report": (24, 36, 36)}[target]
    if percent and scen == "corridor":
        budget = (12, 48, 48)  # a short first page (baselines and S1), then S2-S5
        lay.stretch = "0.93"  # slightly tighter rows
    pages = paginate(groups, budget, True, False)
    brows = baseline_rows(D, scen, 0.30, 0.0005)
    fill_rows(brows, lay.specs)
    frag = f"topfive_full_{scen}" + ("_arxiv_percent" if percent else "")
    if percent:
        event_percent_rows(allrows + brows)
    tex = emit(lay, pages, frag, target, None,
               lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"),
               lambda k: (brows, settings_line(scen, tag="Baselines", alpha0=0.30, L_mu="0.0005")))
    return event_percent_headers(tex) if percent else tex


# topfive_full_{scen}_onepage: one page per scenario, baselines and S1-S5 with all four variants, all metrics of
# build_tables_v10_original.py plus (q_dep, q_cert), two metrics per cell (PAIRS_P), one page block, in the size of
# its float (\footnotesize in arxiv_restored_figures.tex of the paper source).
FULL_ONEPAGE_COLS = PAIRS_P + [["qdep", "qcert"]]


def topfive_full_onepage_fragment(D, scen):
    lay = Layout(f"fop{scen[0]}", [], FULL_ONEPAGE_COLS, 2, [], method_header="Method", uniform=True)
    lay.normalsize = False   # the size of its float (\footnotesize in arxiv_restored_figures.tex of the paper source)
    lay.stretch = "0.93"
    groups = []
    for si in range(1, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        rows = []
        for v in VARIANTS:
            r = vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)},
                     v, (), (), dict(scenario=scen, table="topfive_full_onepage", setting=name, **cfg), (scen, si))
            r["record"]["label"] = tuple_plain(t)
            r["block_head"] = block_text(name, t)
            rows.append(r)
        groups.append(((t["alpha0"],), len(rows) * 2, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    pages = paginate(groups, (1000, 1000, 1000), True, False)   # one page block: baselines and S1-S5
    assert len(pages) == 1, len(pages)
    brows = baseline_rows(D, scen, 0.30, 0.0005)
    fill_rows(brows, lay.specs)
    frag = f"topfive_full_{scen}_onepage"
    event_percent_rows(allrows + brows)
    tex = emit(lay, pages, frag, "paper", None,
               lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"),
               lambda k: (brows, settings_line(scen, tag="Baselines", alpha0=0.30, L_mu="0.0005")))
    return event_percent_headers(tex)


# secondary_{scen}_arxiv (not printed): the secondary and temporal metrics for S1-S5, in fragments of their own.
SECONDARY_ARXIV_COLS = [
    ["m_Y", "m_int", "n_certificates"],
    ["q_Y", "q_unsafe", "q_use"],
    ["q_nocert", "q_notrans", "tau_use"],
    ["h_min", "clearance_min", "Vc@temporal"],
    ["Jc@temporal", "Vp@temporal", "Jp@temporal"],
]


def secondary_arxiv_fragment(D, scen):
    """Eleven secondary performance metrics and four temporal metrics, all variants and baselines.

    Decision fractions already use percent mean +- SD. The certificate metric is a mean number of certificates
    per trajectory, not a frequency. Temporal populations are all N_MC for applicable methods (asserted).
    """
    target, frag, gid = "paper", f"secondary_{scen}_arxiv", f"sa{scen[0]}"
    lay = Layout(gid, [], SECONDARY_ARXIV_COLS, 3, [], method_header="Method", uniform=True)
    groups = []
    for si in range(1, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        rows = []
        for v in VARIANTS:
            found = {sub: cell_or_fail(D, src, "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
                     for sub, src in (("", "final"), ("temporal", "final_temporal"))}
            r = vrow(found, v, (), (), dict(scenario=scen, table=frag, setting=name, **cfg), (scen, si))
            r["record"]["label"] = tuple_plain(t)
            r["block_head"] = block_text(name, t)
            rows.append(r)
        groups.append(((t["alpha0"],), len(rows) * 3, rows))
    allrows = [r for _, _, rows in groups for r in rows]
    brows = baseline_rows(D, scen, 0.30, 0.0005)
    bt = {r["method"]: r for r in baseline_rows(D, scen, 0.30, 0.0005, temporal=True)}
    for r in brows:
        r["found"]["temporal"] = bt[r["method"]]["found"][""]
    for r in allrows + brows:
        for quantity in ("normalized_margin", "estimated_applied_risk"):
            n = r["found"]["temporal"]["cell"]["temporal"][quantity]["assessed_n"]
            assert n in (0, NMC), (scen, r["method"], quantity, n)
    fill_rows(allrows, lay.specs)
    fill_rows(brows, lay.specs)
    pages = paginate(groups, (24, 36, 36), True, False)
    return emit(lay, pages, frag, target, None,
                lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"),
                lambda k: (brows, settings_line(scen, tag="Baselines", alpha0=0.30, L_mu="0.0005")))


def table_key_arxiv_fragment():
    return key_box(
        r"\textbf{Table key.} $\pm$: one sample SD over applicable trajectories. "
        r"$\widehat q^+$: decision fraction; $\widehat q$: trajectory frequency, both in percent. "
        r"Parentheses: $95\%$ upper confidence bound from possible events. "
        r"$\checkmark$: this bound meets its guarantee level. "
        r"$\downarrow$/$\uparrow$: lower/higher is better. "
        r"Bold: best mean of the four variants at one setting, including ties, except certificate number, "
        r"first-use time and peak risk. N-A or ---: not applicable. "
        r"$\widehat m_{\rm cert}$ is the mean number of certificates per trajectory.") + "\n"


# Section VI keys (at most two lines at the text width; the full key is the table_key fragment).
SEC6_KEY = (r"$\pm$: SD; $^+$: decision fraction; $\downarrow$/$\uparrow$: lower/higher is better; counts: trajectories "
            r"of $N_{\rm MC}=200$ with a possible event (\appsecref{app:empirical-risk-classification}), $95\%$ upper "
            r"confidence bound (\%) in parentheses; \checkmark: bound meets its level; N-A: not applicable.")
# sec6_combined: the key names the count exceptions (S* confirmed counts, QT range).
SEC6_KEY_COMBINED = SEC6_KEY.replace(r"with a possible event (\appsecref{app:empirical-risk-classification})",
                                     r"with an event (\appsecref{app:sim-display})")
# "Key to all tables" (table_key): one branch per edition (\versionchoice of the paper source), then the common part.
KEY_CHECK = r"bound at most $\delta$, or $\delta_1$ for $K=\infty$"
KEY_REPORT = (r"$\pm$: SD; $^+$: decision fraction; $\downarrow$/$\uparrow$: lower/higher is better; counts: trajectories of "
              r"$N_{\rm MC}$ with a possible event, or the range from confirmed to possible events "
              r"(\appsecref{app:empirical-risk-classification}), $95\%$ upper confidence bound (\%) in parentheses; "
              rf"\checkmark: {KEY_CHECK}; N-A: not applicable")
KEY_ARXIV = (r"$\pm$: SD; $^+$: decision fraction; $\downarrow$/$\uparrow$: lower/higher is better; counts: trajectories of "
             r"$N_{\rm MC}$ with a possible event (\appsecref{app:empirical-risk-classification}), $95\%$ upper confidence "
             rf"bound (\%) in parentheses; \checkmark: {KEY_CHECK}; N-A: not applicable")
KEY_EXTENDED = (r"$\pm$: SD; $^+$: decision fraction; $\downarrow$/$\uparrow$: lower/higher is better; counts: trajectories "
                r"of $N_{\rm MC}$ with a possible event, or the range from confirmed to possible events, except "
                r"$\widehat q_{\rm cert}$ of S$^\ast$, which gives confirmed events (\appsecref{app:sim-display}); "
                rf"$95\%$ upper confidence bound (\%) in parentheses; \checkmark: {KEY_CHECK}; N-A: not applicable")
KEY_COMMON = (r"; bold: best mean of the four variants at one setting (ties all bold; never SDs, $\widehat m_{\rm cert}$, "
              r"$\widehat m_{\tau^{\rm use}}$, $\widehat m_{p^{\max}}$); $a/m$: $a$ trajectories with a positive bound "
              r"excess among the $m$ with a checked decision ($a$ alone if $m=N_{\rm MC}$), bound over all $N_{\rm MC}$.")


def table_key_fragment(target="paper"):
    """The "Key to all tables". The report target names the S* exception, the bound "from possible events" and
    eps_chk (defined only in the report) in its report branch."""
    rep, common = KEY_REPORT, KEY_COMMON
    if target == "report":
        rep = rep.replace(r"(\appsecref{app:empirical-risk-classification}), $95\%$ upper confidence bound (\%) in "
                          r"parentheses",
                          r"(\appsecref{app:empirical-risk-classification}), except $\widehat q_{\rm cert}$ of S$^\ast$, "
                          r"which gives confirmed events (Table~\ref{tab:report-star-maze}), $95\%$ upper confidence "
                          r"bound (\%) from possible events in parentheses")
        common = common.replace("a positive bound excess", r"an excess above $\varepsilon_{\rm chk}$")
    key = (r"\textbf{Key to all tables.} \versionchoice{report}{" + rep + r"}{\versionchoice{arxiv}{" + KEY_ARXIV
           + "}{" + KEY_EXTENDED + "}}" + common)
    return key_box(key) + "\n"


STAR = HERE / "starred_setting.json"


def starred_row(lay, frag):
    """Row S* (sec6_combined) with the values of starred_setting.json, copied from an N_MC = 200 table outside the
    N200 exports; check_tables.py checks these cells against the JSON and its source lines."""
    S = json.loads(STAR.read_text())
    cfg = dict(S["configuration"], method=S["method"], benchmark=S["scenario"])
    r = dict(found={}, is_variant=True, group=None, method=S["method"], levels=("S*",), level_disp=(S["label"],),
             record=dict(scenario=S["scenario"], table=frag, setting="S*", starred=True,
                         **{k: v for k, v in S["configuration"].items()}))
    t = {k: S["configuration"][k] for k in ("alpha0", "alpha", "delta", "K")}
    set_label(r, t, 1)
    cells = {}
    for c in S["cells"]:
        if c["kind"] == "meanSD":
            tex = f"${c['mean']}\\vpm {c['sd']}$"
        else:
            tex = f"{c['count']} ({c['upper_pct']})"
            if Decimal(c["upper_pct"]) <= Decimal(repr(100 * S["configuration"]["delta"])):
                tex += r"\,\checkmark"
        cells[c["spec"]] = [tex, dict(kind="starred", src="starred_setting.json", configuration=cfg,
                                      source_file=S["source_file"], source_line=c["line"], json_spec=c["spec"])]
    assert set(cells) == set(lay.specs), (sorted(cells), lay.specs)
    r["cells"] = cells
    return S["scenario"], r


def sec6_short_fragment(D, gid):
    """Section VI table (not printed; the reference of sec6_combined and sec6_maze_arxiv): per scenario the baselines
    at alpha0 = 0.30 in a box, then one row per selected setting S1-S5 with its chosen variant.
    Uniform layout: centred cells of equal width at the full line width, a balanced baseline box and the fixed
    padding of the uniform layout; one row per setting (each setting has one variant here), its name in the first
    column and 'Method (alpha0, alpha, delta, K)' on one line, the scenario line above the box."""
    target, frag = "paper", "sec6_short"
    lay = Layout(gid, [("n", "")], SEC6_COLS, 1, ["0.6pt"], method_header=METHOD_HEAD, uniform=True)
    lay.normalsize = False   # the size of its table* float
    per = []
    for scen in SCENS:
        rows = []
        for si in range(1, 6):
            name, t, _, _ = SIX[si]
            cfg = six_cfg(t)
            v = paired(si)
            rows.append(set_label(vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v,
                                                         rho_buf=0.0, **cfg)}, v, (name,), (dn(name),),
                                       dict(scenario=scen, table=frag, setting=name, **cfg), (scen, si)), t, 1))
        fill_rows(rows, lay.specs)
        brows = baseline_rows(D, scen, 0.30, 0.0005)
        fill_rows(brows, lay.specs)
        per.append((scen, rows, brows))
    # \nsecsixkey prints the key box; a document that sets this table after another one with the same key can
    # define it to gobble its argument.
    parts = [r"\providecommand{\nsecsixkey}[1]{#1}\nsecsixkey{" + key_box(SEC6_KEY) + "}"]
    for n, (scen, rows, brows) in enumerate(per):
        parts.append(rf"\vspace{{{2 if n == 0 else 6}pt}}\noindent "
                     + settings_line(scen, b="2.11", phi="1", L_mu="0.0005") + r"\par\vspace{2pt}")
        parts.append(baseline_box_u(lay, brows, r"\textbf{Baselines}\quad $\alpha_0=0.30$", frag, target,
                                    f"{gid}-{scen}-base"))
        parts += [r"\vspace{3pt}", main_table(lay, rows, None, None, frag, target, f"{gid}-{scen}")]
    return fragment_head(lay) + "\n".join(parts) + "\n"


BASE_A0 = (0.15, 0.30)   # alpha0 of the ACI and QT rows in sec6_combined


def sec6_combined_fragment(D, gid):
    """Section VI table with all settings (not printed): per scenario the baselines once (Unfiltered, c = 0 and c_max
    do not depend on alpha0; ACI and QT at alpha0 = 0.15 and 0.30, labelled), then the default setting with V-C (the
    Table I setting of the submission, the earlier submitted version of the paper), S1-S5 with their chosen variants
    and, in the Maze, S* (starred_setting.json). Layout of sec6_short (uniform cells, balanced baseline box); the key
    box once (SEC6_KEY_COMBINED, no bold: one row per setting). The all-method default tables are default_{scen}."""
    target, frag = "paper", "sec6_combined"
    lay = Layout(gid, [("n", "")], SEC6_COLS, 1, ["0.6pt"], method_header=METHOD_HEAD, uniform=True)
    lay.normalsize = False   # the size of its table* float
    S = json.loads(STAR.read_text())
    per = []
    for scen in SCENS:
        rows = []
        for si in range(0, 6):
            name, t, _, _ = SIX[si]
            cfg = six_cfg(t)
            v = paired(si)
            rows.append(set_label(vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v,
                                                         rho_buf=0.0, **cfg)}, v, (name,), (dn(name),),
                                       dict(scenario=scen, table=frag, setting=name, **cfg), (scen, si)), t, 1))
        fill_rows(rows, lay.specs)
        if S["scenario"] == scen:
            rows.append(starred_row(lay, frag)[1])
        bset = {a0: baseline_rows(D, scen, a0, 0.0005) for a0 in BASE_A0}
        for a0 in BASE_A0:
            fill_rows(bset[a0], lay.specs)
            apply_submission_ranges(bset[a0], scen, a0)
        brows = []
        for b in BASELINES:
            same = b not in ("ACI", "QT")
            for a0 in (BASE_A0[:1] if same else BASE_A0):
                r = next(x for x in bset[a0] if x["method"] == b)
                if same:   # these baselines do not depend on alpha0: the printed values must agree
                    o = next(x for x in bset[BASE_A0[1]] if x["method"] == b)
                    assert all(r["cells"][m][0] == o["cells"][m][0] for m in lay.specs), (scen, b)
                else:
                    r["name"] = rf"{BNAME[b]}, $\alpha_0={num(a0)}$"
                brows.append(r)
        per.append((scen, rows, brows))
    parts = [key_box(SEC6_KEY_COMBINED)]
    for n, (scen, rows, brows) in enumerate(per):
        parts.append(rf"\vspace{{{2 if n == 0 else 6}pt}}\noindent "
                     + settings_line(scen, b="2.11", phi="1", L_mu="0.0005") + r"\par\vspace{2pt}")
        parts.append(baseline_box_u(lay, brows, r"\textbf{Baselines}", frag, target, f"{gid}-{scen}-base"))
        parts += [r"\vspace{3pt}", main_table(lay, rows, None, None, frag, target, f"{gid}-{scen}")]
    return fragment_head(lay) + "\n".join(parts) + "\n"


# sec6_maze_arxiv: Maze-only Section VI table (not printed), S1-S5 with the four variants each, six key metrics.
SEC6_MAZE_COLS = [["m_c#3"], ["q_CBF"], ["pk"], ["qdep"], ["qcert"], ["T_mis"]]
SEC6_KEY_MAZE = SEC6_KEY.replace(r"better; counts",
                                 r"better; bold: best mean of the four variants at one setting; counts")


def sec6_maze_arxiv_fragment(D, gid):
    """Maze-only Section VI table (not printed). The key box (SEC6_KEY_MAZE), the baselines at alpha0 = 0.30 in a
    box, then one main table with the Maze settings row and a block row per selected setting S1-S5
    ('S1: (alpha0, alpha, delta, K) = (...)', as topfive_full), the four variants as rows, one line each. Columns SEC6_MAZE_COLS (the
    Section VI precision: m_c/c_max and m_p^max three decimals); bold: best mean of the four variants at one setting.
    The cells are those of the same exports as topfive_full_maze and sec6_short; check_tables.py compares them with
    the Maze cells of sec6_short (chosen variant rows and baselines). sec6_short stays generated (reference)."""
    target, frag, scen = "paper", "sec6_maze_arxiv", "maze"
    lay = Layout(gid, [], SEC6_MAZE_COLS, 1, [], method_header="Method", uniform=True)
    lay.normalsize = False   # the size of its table* float
    rows = []
    for si in range(1, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        for v in VARIANTS:
            r = vrow({"": cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)},
                     v, (), (), dict(scenario=scen, table=frag, setting=name, **cfg), (scen, si))
            r["record"]["label"] = tuple_plain(t)
            r["block_head"] = block_text(name, t)
            rows.append(r)
    fill_rows(rows, lay.specs)
    brows = baseline_rows(D, scen, 0.30, 0.0005)
    fill_rows(brows, lay.specs)
    pages = [dict(key=(0.30,), rows=rows, box=True, first=True, same_sheet=False)]
    return emit(lay, pages, frag, target, SEC6_KEY_MAZE,
                lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"),
                lambda k: (brows, settings_line(scen, tag="Baselines", alpha0=0.30, L_mu="0.0005")))


# The compact one-column tables: main_maze_arxiv_v8 (arXiv Table I), main_corridor_arxiv_v8 (arXiv Table IV, App. F-C)
# and main_{scen}_extended_v8 (extended version). booktabs rules only, at the 245.7 pt \columnwidth, no scaling: the
# baselines (arXiv: alpha0 = 0.30; extended: c = 0, c_max and Unfiltered once, ACI and QT at alpha0 = 0.15 and 0.30),
# then a bold setting row between \cmidrule lines and one line per method: S1-S5 with the four variants (extended:
# also the default setting and, in the Maze, S* with V-C). Columns MAIN_SPECS: m_c/c_max, T_mis, q_CBF^+ and
# m_{p^max} as mean +- SD with the decimals of MAIN_DEC (m_{p^max}: statistics.p_max available mean and SD, a fraction,
# over the trajectories that use a margin below c_max); q_dep and q_cert as the percentage of the N_MC trajectories
# with a possible event (arXiv) or as counts (extended), with the one-sided 95% upper bound in %; '---' where not
# applicable. Bold: lowest printed mean of the four variants of one setting in the MAIN_BOLD columns (ties all bold).
# Cell: mean | {\scriptsize$\vpm SD$} and count | \,(bound) in r@{}l halves, '---' over both halves; header arrows
# '{\downarrow}' set tight, units in the caption; an event percentage or bound equal to 100 prints '100' (exact for
# N_MC = 200), every other one keeps one decimal; the SD of m_c/c_max, T_mis and m_{p^max} drops its leading zero
# (every such SD is below 1; asserted). check_tables.py re-derives every cell, bold mark and formatting rule.
NA_DASH = "---"
MAIN_SPECS = ("m_c", "T_mis", "q_CBF", "pk", "qdep", "qcert")
MAIN_DEC = {"m_c": 2, "T_mis": 1, "q_CBF": 1, "pk": 2}
MAIN_BOLD = ("m_c", "T_mis", "q_CBF", "pk")
MAIN_FONT, MAIN_SEP = r"\footnotesize", "1.5pt"
MAIN_NOLEAD = ("m_c", "T_mis", "pk")
MAIN_HEAD = (r"$\widehat m_c/c_{\max}{\downarrow}$", r"$\widehat T_{\rm mis}{\downarrow}$",
             r"$\widehat q^{+}_{\rm CBF}{\downarrow}$", r"$\widehat m_{p^{\max}}{\downarrow}$",
             r"$\widehat q_{\rm dep}{\downarrow}$", r"$\widehat q_{\rm cert}{\downarrow}$")


def pct_event(found, ev, certified):
    """'p (u)' with p = 100 * possible / N_MC (one decimal, exact for N_MC = 200) and u the stored bound (ceil 0.1 %),
    '100' for a full count or bound."""
    e = found["cell"]["trajectory_events"][ev]
    info = dict(prov(found), source=f"trajectory_events.{ev}")
    n, cnt, u = e["applicable_n"], e["possible"], e["upper95_possible"]
    if not n or (ev == "qcert" and not certified):
        return NA_DASH, dict(info, kind="NA_dash", assessed=n)
    assert n == NMC, (ev, n)
    p = half_up(100 * cnt / NMC, 1)
    up = upper_pct(u)
    level = guarantee_level(found["cell"]["configuration"]) if certified else None
    if certified:
        assert Decimal(up) <= Decimal(repr(100 * level)), ("bound of a variant above its level", ev, up, level)
    info = dict(info, kind="event_pct", count=cnt, n=n, upper95=u, upper_str=up, pct_str=p, level=level,
                certified=certified)
    full = lambda s: "100" if s == "100.0" else s
    p, u = full(info["pct_str"]), full(info["upper_str"])
    return f"{p} ({u})", dict(info, pct_str=p, upper_str=u)


def plain_stat(found, metric):
    """mean +- SD of one statistic with the decimals of MAIN_DEC; '---' when not applicable."""
    src = STAT_SRC.get(metric, metric)
    s = found["cell"]["statistics"][src]
    scale = MET[metric][3]
    dec = MAIN_DEC[metric]
    mk, sk = ("available_mean", "available_sd") if metric in STAT_AVAILABLE else ("mean", "sd")
    info = dict(prov(found), source=f"statistics.{src}", scale=scale, decimals=dec, stat_keys=[mk, sk])
    if not s["applicable"] or s["observed_n"] == 0 or s[mk] is None:
        return NA_DASH, dict(info, kind="NA_dash")
    m, sd = half_up(s[mk] * scale, dec), half_up(s[sk] * scale, dec)
    return f"${m}\\vpm {sd}$", dict(info, kind="meanSD_plain", mean_str=m, sd_str=sd)


def main_split(spec, tex, info):
    """The printed cell (two halves of an r@{}l column) from the one-cell string `tex` and its record."""
    if tex == NA_DASH:
        return r"\multicolumn{2}{c}{---}"
    if info.get("kind") == "meanSD_plain" or info.get("compact_mean"):
        m = f"\\mathbf{{{info['mean_str']}}}" if info.get("bold") else info["mean_str"]
        sd = info["sd_str"]
        if spec in MAIN_NOLEAD:
            assert sd.startswith("0."), (spec, sd)
            sd = sd[1:]
        return f"${m}$&{{\\scriptsize$\\vpm {sd}$}}"
    p, u = re.fullmatch(r"(\S+) \((\S+)\)", tex).groups()   # '100' already for a full count or bound
    return f"{p}&\\,({u})"


# table ids "mmh-*" (main_maze_arxiv_v8), "mch-*" (main_corridor_arxiv_v8), "mex-*" and "cex-*" (extended)
MAIN_TID = {"maze": "mmh", "corridor": "mch"}


def compact_count_event(found, ev, certified):
    tex, info = event_cell(found, ev, certified)
    if info["kind"] == "NA":
        return NA_DASH, dict(info, kind="NA_dash")
    e = found["cell"]["trajectory_events"][ev]
    cfg = found["cell"]["configuration"]
    count = str(info["count"])
    if (cfg["benchmark"], cfg["method"], cfg["alpha0"], ev) in SUBMISSION_RANGES:
        count = f"{e['confirmed']}--{e['possible']}"
    up = "100" if info["upper_str"] == "100.0" else info["upper_str"]
    return f"{count} ({up})", dict(info, kind="event_count_compact", count_str=count, upper_str=up)


# The q_dep and q_cert columns, the events bounded by arXiv Theorem 3 and its proof (every setting of these tables
# has finite K), are marked with a black box and the label "Theorem 3" above the header in all four compact tables.
# \guaranteeboxrow (label row) and \guaranteeboxend (a \noalign before \bottomrule) are defined in release_helper.tex
# of the paper source (TikZ marks); the fallback line at the top of the fragment makes both empty elsewhere (compile
# tests, catalog). No cell or data row changes. Each table has its own TikZ mark id.
GUARANTEE_BOX_ID = {("maze", False): "T1", ("corridor", False): "T4", ("corridor", True): "TXX", ("maze", True): "TIext"}
GUARANTEE_BOX_ROW = r"\guaranteeboxrow{%s}"
GUARANTEE_BOX_END = r"\guaranteeboxend{%s}"
# In main_maze_extended_v8 the "default setting" row (a defined term, not shortened) reaches into the boxed columns;
# the box's left line is interrupted at that row (\guaranteeboxgap, \guaranteeboxendgap in release_helper.tex).
GUARANTEE_BOX_GAP = {("maze", True)}


def main_arxiv_v8_fragment(D, scen, extended=False):
    target, frag = "paper", (f"main_{scen}_extended_v8" if extended else f"main_{scen}_arxiv_v8")
    tp = ("mex" if scen == "maze" else "cex") if extended else MAIN_TID[scen]
    out = [rf"{{{MAIN_FONT}\setlength{{\tabcolsep}}{{{'1.2pt' if extended else MAIN_SEP}}}\renewcommand{{\arraystretch}}{{1.0}}%",
           r"\begin{tabular}{@{}l" + "r@{}l" * 6 + r"@{}}\toprule",
           "Method & " + " & ".join(rf"\multicolumn{{2}}{{c}}{{{h}}}" for h in MAIN_HEAD) + r"\\\midrule"]

    def cells_of(found, certified):
        res = []
        for spec in MAIN_SPECS:
            if spec in EVT:
                res.append((spec,) + (compact_count_event(found, spec, certified) if extended else pct_event(found, spec, certified)))
            else:
                res.append((spec,) + plain_stat(found, spec))
        return res

    def emit(name, cells, method, record, tid):
        printed = []
        for spec, tex, info in cells:
            tex = main_split(spec, tex, info)
            REC.add(target=target, fragment=frag, table=tid, row=record, method=method, spec=spec, printed=tex, **info)
            printed.append(tex)
        return " & ".join([name] + printed) + r"\\"

    for a0 in ((0.15, 0.30) if extended else (0.30,)):
        if extended:
            out.append(rf"\multicolumn{{13}}{{@{{}}l}}{{\textbf{{Baselines}}: $\alpha_0={num(a0)}$}}\\\cmidrule{{1-13}}")
        for b in BASELINES:
            if extended and a0 == 0.30 and b not in ("ACI", "QT"):
                continue
            f = baseline_found(D, scen, b, a0, 0.0005)
            assert f is not None, b
            out.append(emit(BNAME[b], cells_of(f, False), b,
                            dict(scenario=scen, alpha0=a0, L_mu=0.0005, baseline=True), f"{tp}-base-{a0}" if extended else f"{tp}-base"))
    for k, si in enumerate(range(0 if extended else 1, 6)):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        out.append(r"\midrule" if k == 0 else r"\cmidrule{1-13}")
        gap = r"\guaranteeboxgap{%s}" % GUARANTEE_BOX_ID[(scen, extended)] if (scen, extended) in GUARANTEE_BOX_GAP and si == 0 else ""
        out.append(rf"\multicolumn{{13}}{{@{{}}l}}{{{gap}\textbf{{{dn(name)}}}: $(\alpha_0,\alpha,\delta,K)=$ {tuple_plain(t)}}}\\"
                   r"\cmidrule{1-13}")
        rows = []
        for v in (("V-C",) if extended and si == 0 else VARIANTS):
            f = cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
            rows.append((v, cells_of(f, True)))
        for j, spec in enumerate(MAIN_SPECS):
            if spec not in MAIN_BOLD or len(rows) != 4:
                continue
            infos = [cs[j][2] for _, cs in rows]
            if not all(i.get("kind") == "meanSD_plain" for i in infos):
                continue
            best = min(Decimal(i["mean_str"]) for i in infos)
            for i in infos:
                i["bold"] = Decimal(i["mean_str"]) == best
        for v, cs in rows:
            out.append(emit(VSHORT[v], cs, v, dict(scenario=scen, table=frag, setting=name, **cfg), f"{tp}-{name}"))
    if extended and scen == "maze":
        S = json.loads(STAR.read_text())
        out.append(r"\cmidrule{1-13}\multicolumn{13}{@{}l}{\textbf{S$^\ast$}: $(\alpha_0,\alpha,\delta,K)=$ (0.25, 0.225, 0.15, 50)}\\\cmidrule{1-13}")
        cs = []
        for spec in MAIN_SPECS:
            rawspec = "m_c#3" if spec == "m_c" else spec
            c = next(x for x in S["cells"] if x["spec"] == rawspec)
            info = dict(kind="starred_compact", json_spec=rawspec, source_file=S["source_file"], source_line=c["line"])
            if c["kind"] == "meanSD":
                d = MAIN_DEC[spec]
                m, sd = half_up(float(c["mean"]), d), half_up(float(c["sd"]), d)
                info.update(mean_str=m, sd_str=sd, decimals=d, compact_mean=True)
                cs.append((spec, f"${m}\\vpm {sd}$", info))
            else:
                cs.append((spec, f"{c['count']} ({c['upper_pct']})", info))
        out.append(emit("V-C", cs, "V-C", dict(scenario=scen, setting="S*", label="(0.25, 0.225, 0.15, 50)",
                                             starred=True, **S["configuration"]), "mex-star"))
    if (scen, extended) in GUARANTEE_BOX_ID:   # all four compact tables (Theorem 3)
        gid = GUARANTEE_BOX_ID[(scen, extended)]
        if extended:   # the extended tables keep 1.2 pt column gaps; 2 pt more on each side of the box line clears the text
            out[1] = (r"\begin{tabular}{@{}l" + "r@{}l" * 4 + r"@{\hspace{\dimexpr\tabcolsep+2pt\relax}}>{\hspace{2pt}}r@{}l"
                      + "r@{}l" + r"@{}}\toprule")
        out[1] = out[1].removesuffix(r"\toprule") + "\n" + GUARANTEE_BOX_ROW % gid + r"\toprule"
        out.insert(0, r"\providecommand{\guaranteeboxrow}[1]{}\providecommand{\guaranteeboxend}[1]{}"
                      r"\providecommand{\guaranteeboxgap}[1]{}\providecommand{\guaranteeboxendgap}[1]{}%")
        end = r"\guaranteeboxendgap{%s}" if (scen, extended) in GUARANTEE_BOX_GAP else GUARANTEE_BOX_END
        out.append(end % gid + r"\bottomrule\end{tabular}}")
    else:
        out.append(r"\bottomrule\end{tabular}}")
    return "\n".join(out) + "\n"


# The Section VI tables (sec6_*) print m_c/c_max with three decimals and q_use^+ with two; the peak estimated
# exceedance risk pk follows q_unsafe^+.
SEC6_COLS = [["m_c#3"], ["q_use#2"], ["q_CBF"], ["q_unsafe"], ["pk"], ["qdep"], ["qcert"], ["T_mis"]]


# star_maze (report Table 8): the S* setting with all four variants. The setting is not on the N200 main grid; its
# cells are copied from the N_MC = 200 table of the paper source named in starred_setting.json ("source_file") at the
# lines STAR_VARIANT_LINES and the columns of STAR_MAP. Counts: the qcert count is that of the confirmed events (lower end of the source range)
# with the source bound (as for row S* of sec6_combined), every other count the upper end with its bound.
# check_tables.py re-reads every cell from the source (check_starred_variant).
# The full-metric layout of topfive_full (TOP_COLS, three lines per variant). Each source mean and SD
# is converted to the unit of MET (p_max: fraction -> %, h_min: mm -> m) and printed with the decimals of MET when
# every value the source digits allow rounds (half up) to the same string, else with the source digits; bold as in
# fill_rows (lowest/highest printed mean of the four variants, ties all bold). clearance_min is not in the source: "--".
STAR_VARIANT_LINES = {"V-U": (13, 27, 41, 55), "V-D": (14, 28, 42, 56), "V-Cd": (15, 29, 43, 57),
                      "V-C": (16, 30, 44, 58)}
STAR_SOURCE_NAME = {"V-U": "V-U", "V-D": r"V-$\Delta$", "V-Cd": "V-L", "V-C": "V-C"}
# spec -> (source block, source column, factor from the source unit to the MET unit)
STAR_MAP = {"m_c": (0, 2, "1"), "q_use": (0, 3, "1"), "p_max": (0, 4, "100"), "tau_use": (0, 5, "1"),
            "T_mis": (0, 1, "1"), "q_CBF": (1, 2, "1"), "q_unsafe": (1, 3, "1"), "h_min": (1, 5, "0.001"),
            "qdep": (2, 1, None), "qcert": (2, 3, None), "n_certificates": (2, 5, "1"),
            "m_Y": (3, 1, "1"), "q_Y": (3, 2, "1"), "m_int": (3, 3, "1"), "q_nocert": (3, 4, "1"),
            "q_notrans": (3, 5, "1")}
STAR_MISSING = {"clearance_min"}


def star_number(txt, factor, dec):
    """One source number ('0.610', '{<}0.01') in the MET unit: (source digits, rounded string or None, sort key).
    The source prints the value x half up with step u, so the value lies in [x - u/2, x + u/2); the rounded string
    (`dec` decimals, half up) exists when that whole interval rounds to one string."""
    lt = txt.startswith("{<}")
    raw = Decimal(txt[3:] if lt else txt)
    f = Decimal(factor)
    u = Decimal(1).scaleb(raw.as_tuple().exponent) * f
    val = raw * f
    src_dec = max(0, -u.normalize().as_tuple().exponent)
    src = f"{val:.{src_dec}f}"
    if lt:   # an upper bound only: kept as printed; it sorts just below its bound
        return "{<}" + src, ("{<}" + src if src_dec <= dec else None), val - u / 10
    if src_dec <= dec:
        return src, src, val
    q = Decimal(1).scaleb(-dec)
    lo, hi = val - u / 2, val + u / 2
    a = lo.quantize(q, rounding=ROUND_HALF_UP)
    b = hi.quantize(q, rounding=ROUND_HALF_UP)
    if hi == b - q / 2:   # hi itself is excluded: a value just below hi rounds down
        b -= q
    if a != b:
        return src, None, val
    out = f"{a:.{dec}f}"
    return src, (out[1:] if out.startswith("-") and a == 0 else out), val


def star_source_cell(text, spec, delta):
    """TeX of one S* cell from its source text (marks removed by unmark), and its record (digits and bold are set
    by star_fragment)."""
    if spec in ("qdep", "qcert"):
        m = re.fullmatch(r"(?:\\textbf\{)?(\d+)(?:--(\d+))?/200(?:\})? \((\d+\.\d)\\%\)", text)
        assert m, (spec, text)
        cnt = int(m.group(1)) if spec == "qcert" else int(m.group(2) or m.group(1))
        tex = f"{cnt} ({m.group(3)})"
        if Decimal(m.group(3)) <= Decimal(repr(100 * delta)):
            tex += r"\,\checkmark"
        return tex, dict(kind="starred_variant", count=cnt, upper_str=m.group(3), source_text=text,
                         count_rule="confirmed" if spec == "qcert" else "upper")
    if text == "---":   # not defined for this variant (no use of a smaller margin)
        return r"\na", dict(kind="starred_variant", na=True, source_text=text)
    src = re.sub(r"\\mathbf\{([^}]*)\}", r"\1", text)   # bold of the source table
    m = re.fullmatch(r"\$(.+?)\\pm\s*(\S+)\$", src)
    assert m, (spec, text)
    _, arrow, _, _, dec, boldable = MET[spec]
    factor = STAR_MAP[spec][2]
    return None, dict(kind="starred_variant", num=(star_number(m.group(1), factor, dec),
                                                   star_number(m.group(2), factor, dec)),
                      boldable=boldable, arrow=arrow, source_text=text, factor=factor, decimals=dec)


def star_fragment(gid, target="paper"):
    """star_maze (report Table 8): S* in the Maze with the four variants, in the layout of topfive_full."""
    S = json.loads(STAR.read_text())
    lines = (SRC / S["source_file"]).read_text().splitlines()
    cfg = S["configuration"]
    lay = Layout(gid, [], TOP_COLS, 3, [], method_header="Method", uniform=True)
    t = {k: cfg[k] for k in ("alpha0", "alpha", "delta", "K")}
    rows = []
    for v in VARIANTS:
        r = vrow({}, v, (), (), dict(scenario=S["scenario"], table="star", setting="S*", starred=True, **cfg),
                 None)
        r["record"]["label"] = tuple_plain(t)
        r["block_head"] = rf"\textbf{{S$^\ast$}}: $(\alpha_0,\alpha,\delta,K)$ = {tuple_plain(t)}"
        r["cells"] = {}
        for spec in lay.specs:
            if spec in STAR_MISSING:
                r["cells"][spec] = ["--", dict(kind="starred_missing", src="starred_variant",
                                               source_file=S["source_file"],
                                               configuration=dict(cfg, method=v, benchmark=S["scenario"]))]
                continue
            blk, col, _ = STAR_MAP[spec]
            ln = STAR_VARIANT_LINES[v][blk]
            parts = [x.strip() for x in lines[ln - 1].rstrip().rstrip("\\").split(" & ")]
            assert parts[0] == STAR_SOURCE_NAME[v], (v, ln, parts[0])
            tex, info = star_source_cell(unmark(parts[col]), spec, cfg["delta"])
            info.update(src="starred_variant", source_file=S["source_file"], source_line=ln, source_column=col,
                        configuration=dict(cfg, method=v, benchmark=S["scenario"]))
            r["cells"][spec] = [tex, info]
        rows.append(r)
    for spec in lay.specs:
        infos = [r["cells"][spec][1] for r in rows if r["cells"][spec][1].get("num")]
        # digits: the MET decimals when every mean and SD of the column rounds unambiguously, else the source digits
        rounded = all(x[1] is not None for i in infos for x in i["num"])
        for i in infos:
            (ms, mr, mk), (ss, sr, _) = i.pop("num")
            i.update(mean_str=mr if rounded else ms, sd_str=sr if rounded else ss,
                     sort_key=str(Decimal(mr) if rounded and not mr.startswith("{<}") else mk),
                     digits="rounded" if rounded else "source")
    for spec in lay.specs:   # bold as in fill_rows: lowest/highest printed mean of the four variants
        infos = [r["cells"][spec][1] for r in rows]
        if not all(i.get("mean_str") is not None and i["boldable"] for i in infos):
            continue
        vals = [Decimal(i["sort_key"]) for i in infos]
        best = max(vals) if infos[0]["arrow"] == "up" else min(vals)
        for i, x in zip(infos, vals):
            i["bold"] = x == best
    for r in rows:
        for spec, (tex, info) in r["cells"].items():
            if tex is None:
                mm = rf"\mathbf{{{info['mean_str']}}}" if info.get("bold") else info["mean_str"]
                r["cells"][spec][0] = f"${mm}\\vpm {info['sd_str']}$"
                info.setdefault("bold", False)
    sets = settings_line(S["scenario"], b=str(cfg["budget_normalizer"]), phi="1", L_mu=str(cfg["L_mu"]))
    sets += r"\quad --: not available; a column keeps its recorded digits where rounding to those of Table~\ref{tab:sim-topfive-maze} is ambiguous"
    pages = uniform_pages([((cfg["alpha0"],), rows)])
    return emit(lay, pages, "star_maze", target, None, lambda pg: sets)


def two_lines(a, b):
    return stack([a, b])


def topfive_events_fragment(D, target, gid, with_default=True):
    cols = [["qcert@corridor"], ["qdep@corridor"], ["qcert@maze"], ["qdep@maze"]]
    lay = Layout(gid, [("n", ""), ("a0", r"$\alpha_0$"), ("a", r"$\alpha$"), ("d", r"$\delta$"), ("k", r"$K$")],
                 cols, 1, ["1.0pt"] * 5, method_header="Method",
                 groups=[(r"\textbf{Corridor}", 0, 2), (r"\textbf{Maze}", 2, 2)])
    frag = "topfive_events" if with_default else "topfive_events_nodef"
    groups = []
    for si in range(0 if with_default else 1, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        lv = (name,) + tuple(num(t[k]) for k in ("alpha0", "alpha", "delta")) + (kfmt(t["K"]),)
        disp = (dn(name),) + tuple(f"${x}$" for x in lv[1:])
        rows = []
        for v in VARIANTS:
            found = {s: cell_or_fail(D, "final", "main_grid_200", benchmark=s, method=v, rho_buf=0.0, **cfg)
                     for s in SCENS}
            rows.append(vrow(found, v, lv, disp, dict(scenario="corridor+maze", table=frag, setting=name, **cfg), None))
        groups.append((0, 4, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    budget = {"paper": (40, 40, 40), "report": (24, 40, 40)}[target]
    pages = paginate(groups, budget, True, boxes(target))
    a0s = [0.15, 0.30] if with_default else [0.30]
    brows = []
    for a0 in a0s:
        brows += baseline_rows(D, {"corridor": "corridor", "maze": "maze"}, a0, 0.0005,
                               level=((num(a0),), (f"${num(a0)}$",)) if len(a0s) > 1 else None)
    fill_rows(brows, lay.specs)
    apply_submission_ranges(brows, None)   # the submission's range cell (SUBMISSION_RANGES)
    bsets = settings_line("Corridor and Maze", tag="Baselines", alpha0=None if len(a0s) > 1 else 0.30, L_mu="0.0005")
    return emit(lay, pages, frag, target, key_events(allrows + brows),
                lambda pg: settings_line("Corridor and Maze", b="2.11", phi="1", L_mu="0.0005"),
                lambda k: (brows, bsets), [("a0", r"$\alpha_0$")] if len(a0s) > 1 else ())


def kgroup_fragment(D, scen, target, gid, frag, table, sel, lms, k_specs, inf_specs, drop, k_head, inf_head,
                    block_extra=None, stretch=None, sets_extra=""):
    """The chosen variant of each setting in `sel` with the setting's K and K = infinity side by side (column
    groups `k_head`, `inf_head`); a block row per setting 'S1: V-C (alpha0, alpha, delta, K) = (...)' and one row
    per L_mu in `lms` (first column L_mu). drop (arXiv): event specs that are 0 or N-A in every cell of both
    groups are omitted and listed in the '% OMITTED' line. Used by checks_{scen}_arxiv, checks_{scen}_L0005_nodef and
    detailed_events_{scen}."""
    rows = []
    for si in sel:
        name, t, _, _ = SIX[si]
        v = paired(si)
        bh = block_text(name, t, v, fields=NO_K) + (block_extra(t) if block_extra else "")   # K is a column group
        for lm in lms:
            fK, cfg = anchor_cell(D, si, scen, v, "detailed_200", L_mu=lm)
            fI, _ = anchor_cell(D, si, scen, v, "detailed_200", L_mu=lm, K="infinity")
            r = vrow({"K": fK, "inf": fI}, v, (), (), dict(scenario=scen, table=table, setting=name, **cfg), None)
            r["record"]["label"] = tuple_plain(t)
            r["record"]["block"] = [tuple_plain(t, NO_K)]
            r["record"]["block_fields"] = NO_K
            r["record"]["swept"] = "period"
            r["name"] = "$0.0005$" if lm == 0.0005 else "$0.005$"
            r["block_head"] = bh
            rows.append(r)
    specs = [f"{m}@K" for m in k_specs] + [f"{m}@inf" for m in inf_specs]
    fill_rows(rows, specs)
    gone = []
    if drop:
        base = [m for m in ARX_EVTS if m in k_specs or m in inf_specs]
        gone = [m for m in base if all(is_zero(r["cells"][s][1]) for r in rows
                                       for s in (f"{m}@K", f"{m}@inf") if s in r["cells"])]
        k_specs = [m for m in k_specs if m not in gone]
        inf_specs = [m for m in inf_specs if m not in gone]
        for r in rows:
            for m in gone:
                for s in (f"{m}@K", f"{m}@inf"):
                    r["cells"].pop(s, None)
    cols = [[f"{m}@K"] for m in k_specs] + [[f"{m}@inf"] for m in inf_specs]
    groups = [(k_head, 0, len(k_specs)), (inf_head, len(k_specs), len(inf_specs))]
    lay = Layout(gid, [], cols, 1, [], method_header=r"$L_\mu$", groups=groups, stretch=stretch, uniform=True)
    pages = uniform_pages([(0, rows)], box_keys=())
    sets = lambda pg: settings_line(scen, b="2.11", phi="1") + sets_extra
    return (omitted_line(gone) if drop else "") + emit(lay, pages, frag, target, None, sets)


def detailed_events_fragment(D, scen, target, gid):
    """Chosen variant at the six settings, L_mu in {0.0005, 0.005}; K of the setting and K = infinity.
    Block row per setting with delta_1, rows L_mu, column groups K of the setting / K = infinity."""
    st = uniform_stretch(target, 1 + 2 + 6 * 3, 60)
    return kgroup_fragment(D, scen, target, gid, f"detailed_events_{scen}", "detailed_events", range(6),
                           (0.0005, 0.005), ["qcert", "qdep"], ["qcert", "qdep"], False,
                           two_lines(r"$K$ of the setting", r"(level $\delta$)"),
                           two_lines(r"$K=\infty$", r"(level $\delta_1$)"),
                           block_extra=lambda t: rf", $\delta_1={delta1_pct(t['delta'])}\%$", stretch=st,
                           sets_extra=r"\quad $\delta_1=\delta/(2\kappa\ln^2 2)$")   # b = kappa


def checks_fragment(D, scen, L_mu, target, gid, nodef=False):
    """Chosen variant only: Default and S1-S5 at K of the setting and K = infinity (nodef: S1-S5, arXiv).
    Block row per setting 'S1: V-C (alpha0, alpha, delta, K) = (...)', rows K (the setting's K
    and infinity), six event columns; the baselines of each alpha0 in a balanced box above its settings; one page,
    rows stretched to fill it. nodef (not used by an edition): kgroup_fragment at L_mu = 0.0005."""
    lmd = "0.0005" if L_mu == 0.0005 else "0.005"
    frag = f"checks_{scen}_{'L0005' if L_mu == 0.0005 else 'L005'}" + ("_nodef" if nodef else "")
    if nodef:
        return kgroup_fragment(D, scen, "paper", gid, frag, "checks_nodef", range(1, 6), (L_mu,),
                               ["qcert", "qdep", "rCBF", "runsafe", "qA2", "qpw"], ["qcert", "qdep"], False,
                               r"$K$ of the setting", r"$K=\infty$",
                               stretch=uniform_stretch("paper", 1 + 2 + 5 * 2, 60))
    cols = [["qcert"], ["qdep"], ["rCBF"], ["runsafe"], ["qA2"], ["qpw"]]
    # one page: box (1 + 1 + 5) + default block (1 + 1 + 1 + 2), box (7) + S1-S5 (1 + 1 + 5 x 3) = 36 table rows
    fixed = 70 if target == "paper" else 55
    lay = Layout(gid, [], cols, 1, [], method_header=r"$K$", stretch=uniform_stretch(target, 36, fixed), uniform=True)
    keyed = []
    for si in range(6):
        name, t, _, _ = SIX[si]
        v = paired(si)
        rows = []
        for K in (t["K"], "infinity"):
            f, cfg = anchor_cell(D, si, scen, v, "detailed_200", L_mu=L_mu, K=K)
            r = vrow({"": f}, v, (), (), dict(scenario=scen, table="checks", setting=name, **cfg), None)
            r["record"]["label"] = tuple_plain(dict(t, K=K))
            r["record"]["block"] = [tuple_plain(t, NO_K)]
            r["record"]["block_fields"] = NO_K
            r["record"]["swept"] = "period"   # the K = infinity row differs from its block tuple in K only
            r["name"] = f"${kfmt(K)}$"
            r["block_head"] = block_text(name, t, v, fields=NO_K)   # K is the row variable
            rows.append(r)
        keyed.append(((t["alpha0"],), rows))
    allrows = [r for _, rs in keyed for r in rs]
    fill_rows(allrows, lay.specs)
    plan = {}
    for k, rs in keyed:
        plan.setdefault(k, []).extend(rs)
    pages = uniform_pages(list(plan.items()), same_sheet=set(range(1, len(plan))))
    sets = lambda pg: settings_line(scen, b="2.11", phi="1", L_mu=lmd)
    bc = {}

    def base_fn(k):
        if k not in bc:
            br = baseline_rows(D, scen, k[0], L_mu)
            fill_rows(br, lay.specs)
            apply_submission_ranges(br, scen, k[0])   # the submission's range cell (SUBMISSION_RANGES), all editions
            bc[k] = (br, settings_line(scen, tag="Baselines", alpha0=k[0], L_mu=lmd))
        return bc[k]

    return emit(lay, pages, frag, target, None, sets, base_fn)


# ----------------------------------------------------------------------------------------------- arXiv
ARX_EVTS = ["qcert", "qdep", "rCBF", "runsafe", "qA2", "qpw"]
ARX_STRESS = ["m_c", "q_CBF", "qdep", "qA2", "qpw"]
OMITTED_TAG = "% OMITTED (0 or N-A in every row):"


def is_zero(info):
    """A printed cell that carries nothing: N-A, not simulated, a zero count, or 0.00 +- 0.00."""
    k = info.get("kind")
    if k in ("NA", "not_simulated"):
        return True
    if k == "event":
        return info["count"] == 0
    if k == "meanSD":
        return Decimal(info["mean_str"]) == 0 and Decimal(info["sd_str"]) == 0
    return False


def drop_zero(rows, specs):
    """arXiv (no all-zero columns): (kept specs, omitted specs)."""
    gone = [m for m in specs if all(is_zero(r["cells"][m][1]) for r in rows)]
    return [m for m in specs if m not in gone], gone


def omitted_line(gone):
    return f"{OMITTED_TAG} {', '.join(gone) if gone else 'none'}\n"


def checks_arxiv_fragment(D, scen, gid):
    """checks_{scen}_arxiv (not printed): chosen variant at S1-S5, L_mu in {0.0005, 0.005}, K of the setting and
    K = infinity; q_cert, q_dep, q_CBF, q_unsafe, q_A2, q_pw with the one-count rule; columns that are 0 or N-A in
    every row are omitted (listed in the '% OMITTED' line).
    Uniform layout: block row per setting, rows L_mu, the setting's K and K = infinity as two column groups; full
    width, centred, rows stretched so that the pair of tables fills its page."""
    # the Corridor and Maze tables share one page: 2 x (settings + 2 header lines + 5 block rows + 10 rows)
    return kgroup_fragment(D, scen, "paper", gid, f"checks_{scen}_arxiv", "checks_arxiv", range(1, 6),
                           (0.0005, 0.005), ARX_EVTS, ARX_EVTS, True, r"$K$ of the setting", r"$K=\infty$",
                           stretch=uniform_stretch("paper", 36, 75))


def stress_arxiv_fragment(D, scen, gid, target="paper", percent=False):
    """Stress table: chosen variant at S1-S5, phi in {0, 0.25, 0.5, 0.75, 1} x L_mu in {0.0005, 0.005};
    m_c/c_max, q_CBF^+, q_dep, q_A2, q_pw; phi < 1 marked with a dagger (scaled bound, Assumption 2 may fail);
    all-zero columns omitted. Uniform layout: block row per (setting, L_mu) 'S1: V-C (alpha0, alpha, delta, K) =
    (...), L_mu = 0.0005', rows phi; two pages (five blocks each, 'Table N (continued).'), rows stretched to fill
    them. Target "paper": stress_{scen}_arxiv (not printed); target "report": the same cells as stress_{scen}_cases
    (report Tables 17 and 18, Section 2.7)."""
    frag = f"stress_{scen}_arxiv" if target == "paper" else f"stress_{scen}_cases"
    if percent:
        frag += "_percent"
    rows, blocks = [], []
    for si in range(1, 6):
        name, t, _, _ = SIX[si]
        v = paired(si)
        for lm in (0.0005, 0.005):
            bh = block_text(name, t, v) + (r", $L_\mu=0.0005$" if lm == 0.0005 else r", $L_\mu=0.005$")
            blk = []
            for p in (0.0, 0.25, 0.5, 0.75, 1.0):
                f, cfg = anchor_cell(D, si, scen, v, "stress_drift_200", bound_scale=p, L_mu=lm)
                r = vrow({"": f}, v, (), (), dict(scenario=scen, table="stress_arxiv", setting=name, **cfg), None)
                r["name"] = f"${p:g}" + (r"^\dagger$" if p < 1 else "$")
                r["record"]["label"] = tuple_plain(t)
                r["record"]["block"] = [tuple_plain(t)]
                r["block_head"] = bh
                blk.append(r)
            rows += blk
            blocks.append(blk)
    fill_rows(rows, ARX_STRESS)
    kept, gone = drop_zero(rows, ARX_STRESS)
    if percent:
        for row in rows:
            for spec in kept:
                tex, info = row["cells"][spec]
                if info["kind"] == "event":
                    pct = half_up(100 * info["count"] / NMC, 1)
                    tex = f"{pct} ({info['upper_str']})" + (r"\,\checkmark" if info["check"] else "")
                    row["cells"][spec] = [tex, dict(info, kind="event_pct_sweep", pct_str=pct)]
    # a page: settings + header + 5 x (block row + 5 rows) = 32 table rows
    lay = Layout(gid, [], [[m] for m in kept], 1, [], method_header=r"$\phi$",
                 stretch=uniform_stretch(target, 32, 55, hi=1.4), uniform=True)
    pages = uniform_pages([(0, [r for b in blocks[:5] for r in b]), (0, [r for b in blocks[5:] for r in b])], box_keys=())
    sets = lambda pg: (settings_line(scen, b="2.11") +
                       r"\quad $^\dagger$: $d_t^{\mathrm V}(c;\sigma,\tau)$ scaled by $\phi<1$")   # one line (the caption names Assumption 2)
    return omitted_line(gone) + emit(lay, pages, frag, target, None, sets)


def stress_representative_arxiv_fragment(D):
    """S3 with V-L in both scenarios: all five bound factors at both recorded drift bounds.

    The four displayed metrics retain the stress-table precision and the percentage event/UCB convention.
    No rows are selected by their outcomes. Full setting comparisons remain in the report fragments.
    """
    frag, target, gid = "stress_representative_arxiv_percent", "paper", "strp"
    specs = ["m_c", "q_CBF", "qdep", "qcert"]
    si = 3
    name, setting, variant, _ = SIX[si]
    assert name == "S3" and variant == "V-Cd"
    rows = []
    for scen in SCENS:
        for lm in (0.0005, 0.005):
            for phi in (1.0, 0.75, 0.5, 0.25, 0.0):
                found, cfg = anchor_cell(D, si, scen, variant, "stress_drift_200", bound_scale=phi, L_mu=lm)
                row = vrow({"": found}, variant, (), (),
                           dict(scenario=scen, table=frag, setting=name, **cfg), None)
                row["name"] = f"${phi:g}$"
                row["record"]["label"] = tuple_plain(setting)
                row["record"]["block"] = [tuple_plain(setting)]
                row["block_head"] = rf"\textbf{{{SCEN[scen]['name']}}}: $L_\mu={lm:g}$"
                rows.append(row)
    fill_rows(rows, specs)
    event_percent_rows(rows)
    for row in rows:
        for spec in ("qdep", "qcert"):
            row["cells"][spec][0] = row["cells"][spec][0].replace(r"\,\checkmark", "")
    lay = Layout(gid, [], [[spec] for spec in specs], 1, [], method_header=r"$\phi$",
                 stretch=1.0, inner_rules=False, uniform=True)
    pages = uniform_pages([(0, rows)], box_keys=())
    settings = block_text(name, setting, variant) + r"\quad $N_{\rm MC}=200$"
    return event_percent_headers(emit(lay, pages, frag, target, None, lambda pg: settings))


def variation_fragment(D, scen, target, gid, arxiv=False):
    """Block row per setting 'S1: (alpha0, alpha, delta, K) = (...)', the four variants as
    rows; the baselines at both alpha0 in an aligned box with a block row per alpha0; one page, rows stretched."""
    cols = [["Vc"], ["Jc"], ["Vp"], ["Jp"], ["nK"]]
    # one page: box (1 + 1 + 2 + 10) + table (1 + 1 + 6 x 5) = 46 table rows
    fixed = 60 if target == "paper" else 45
    lay = Layout(gid, [], cols, 1, [], method_header="Method", stretch=uniform_stretch(target, 46, fixed), uniform=True)
    rows = []
    for si in range(1 if arxiv else 0, 6):
        name, t, _, _ = SIX[si]
        cfg = six_cfg(t)
        for v in VARIANTS:
            r = vrow({"": cell_or_fail(D, "final_temporal", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0,
                                       **cfg)}, v, (), (), dict(scenario=scen, table="variation", setting=name, **cfg),
                     None)
            r["record"]["label"] = tuple_plain(t)
            r["record"]["block"] = [tuple_plain(t)]
            r["block_head"] = block_text(name, t)
            rows.append(r)
    fill_rows(rows, lay.specs)
    brows = []
    for a0 in ((0.30,) if arxiv else (0.15, 0.30)):
        rs = baseline_rows(D, scen, a0, 0.0005, temporal=True)
        for r in rs:
            r["block_head"] = rf"$\alpha_0={num(a0)}$"
        brows += rs
    fill_rows(brows, lay.specs)
    bsets = settings_line(scen, tag="Baselines", L_mu="0.0005")
    pages = uniform_pages([(0, rows)])
    return emit(lay, pages, f"variation_{scen}" + ("_arxiv" if arxiv else ""), target, None,
                lambda pg: settings_line(scen, b="2.11", phi="1", L_mu="0.0005"), lambda k: (brows, bsets))


def old_candidates(D):
    cr = json.loads(CANDIDATES.read_text())
    out = []
    for i, c in enumerate(cr["candidates"], 1):
        t = dict(alpha0=c["alpha0"], alpha=c["alpha"], delta=c["delta"], K=c["K"])
        aid = f"candidate_{i}"
        assert aid in D.anchor_ids["final"]
        out.append((aid, t, c["variant"]))
    return out


def n30_fragment(D, scen, target, gid):
    """A block row per setting (variant and tuple on one line); the rows are grouped by the
    changed parameter (base; L_mu; phi; b) in a 'Change' column, value in the next, thin rules between groups."""
    if target == "paper":
        cols, L = [["m_c"], ["q_CBF"], ["q_unsafe"], ["T_mis"], ["qcert"], ["qdep"], ["qA2"], ["qpw"]], 1
    else:
        cols, L = [["m_c", "T_mis"], ["q_CBF", "q_unsafe"], ["qcert", "qdep"], ["qA2", "qpw"]], 2
    lay = Layout(gid, [("c", "Change")], cols, L, ["0.6pt"], method_header="Value", inner_rules=L > 1)
    groups = []
    for aid, t, var in old_candidates(D):
        tp = tuple_plain(t)
        head_txt = rf"{VNAME[var]} $(\alpha_0,\alpha,\delta,K)$ = {tp}"
        changes = [("base", "--", "stress_drift_200", {}),
                   (r"$L_\mu$", "$0.005$", "stress_drift_200", dict(L_mu=0.005))]
        changes += [(r"$\phi$", f"${p:g}$", "stress_drift_200", dict(bound_scale=p)) for p in (0.0, 0.25, 0.5, 0.75)]
        changes += [("$b$", f"${b:g}$", "budget_normalizer_200", dict(budget_normalizer=b)) for b in (2.5, 3.0, 4.0)]
        rows = []
        for par, val, fam, over in changes:
            cfg = dict(dict(t, **BASE), **over)
            f = cell_or_fail(D, "final", fam, benchmark=scen, method=var, rho_buf=0.0, **cfg)
            assert f["anchor_id"] == aid, (f["anchor_id"], aid)
            lab = "base" if par == "base" else f"{par[:-1]}={val[1:]}"
            r = vrow({"": f}, var, (par,), (par,), dict(scenario=scen, table="n30", anchor=aid, change=lab, label=tp,
                                                       **cfg), None)
            r["name"] = val
            r["block_head"] = head_txt
            rows.append(r)
        groups.append((0, len(rows) * L, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_rows(allrows, lay.specs)
    budget = {"paper": (45, 45, 45), "report": (18, 36, 36)}[target]
    pages = paginate(groups, budget, False, False)
    return emit(lay, pages, f"n30_record_{scen}", target, None,
                lambda pg: settings_line(scen) + r"\quad base: $b=\kappa$, $\phi=1$, $L_\mu=0.0005$")


# ----------------------------------------------------------------------------------------------- design tables
# No Count column; columns Comparison | Varied values | Fixed values; one row per
# comparison; every varied parameter on its own aligned line; fixed parameters one per line; names centred vertically.
# The printed value sets are derived from the exports (design_sets) and re-derived by check_tables.py.
def design_sets(D):
    """Parameter sets of each comparison, taken from the exports (final: main grid, baselines; new: S1-S5 studies)."""

    def cfgs(src, fams):
        return [f["cell"]["configuration"] for k, f in D.idx.items() if k[0] == src and k[1] in fams]

    def vals(cs, k):
        return sorted({key(c[k]) for c in cs if c.get(k) is not None}, key=lambda x: (isinstance(x, str), x))

    main = cfgs("final", {"main_grid_200"})
    pairs = {}
    for c in main:
        pairs.setdefault(key(c["alpha0"]), set()).add(key(c["alpha"]))
    base = cfgs("final", {"matched_baselines_200", "matched_ACI_200"})
    S = dict(main=dict(pairs={a: sorted(v) for a, v in sorted(pairs.items())}, delta=vals(main, "delta"),
                       K=vals(main, "K"), methods=sorted({c["method"] for c in main})),
             base=dict(methods=sorted({c["method"] for c in base})))
    S["base"]["alpha0"] = {m: vals([c for c in base if c["method"] == m], "alpha0") for m in S["base"]["methods"]}
    S["base"]["L_mu"] = vals(base, "L_mu")
    for fam in ("budget_normalizer_200", "stress_drift_200", "detailed_200"):
        cs = cfgs("new", {fam})
        S[fam] = dict(b=vals(cs, "budget_normalizer"), phi=vals(cs, "bound_scale"), L_mu=vals(cs, "L_mu"),
                      K=vals(cs, "K"), methods=sorted({c["method"] for c in cs}))
    return S


def _set(vs, fmt=None):
    fmt = fmt or num
    return r"\{" + ",".join(fmt(v) for v in vs) + r"\}"


def _plain(x):
    """b and phi: shortest plain form (2.11, 2.5, 3, 0.25)."""
    return f"{float(x):g}"


def lines_cell(lines):
    """One aligned block: each entry a line; an entry (left, right) goes into the aligned two-column part.
    Struts (top 2.6ex, bottom 1.1ex) sit inside the first and last line (a \\multicolumn must open its cell)."""
    out = []
    n = len(lines)
    for i, ln in enumerate(lines):
        top = r"\rule{0pt}{2.6ex}" if i == 0 else ""
        bot = r"\rule[-1.1ex]{0pt}{0pt}" if i == n - 1 else ""
        if isinstance(ln, tuple):
            out.append(rf"{top}{ln[0]} & {ln[1]}{bot}")
        else:
            out.append(rf"\multicolumn{{2}}{{@{{}}l@{{}}}}{{{top}{ln}{bot}}}")
    return r"\begin{tabular}[c]{@{}r@{\ }l@{}}" + r"\\".join(out) + r"\end{tabular}"


def design_rows(D, arxiv, target="paper"):
    """Rows (name, varied lines, fixed lines, record) of the sweep-design table."""
    S = design_sets(D)
    m = S["main"]
    assert m["methods"] == sorted(VARIANTS), m["methods"]
    four = r"methods: V-U, V-$\Delta$, V-L, V-C"
    pair_lines = [(rf"$\alpha_0={num(a0)}$:", rf"$\alpha={num(al[0])}$" if len(al) == 1 else rf"$\alpha\in{_set(al)}$")
                  for a0, al in m["pairs"].items()]
    main_fixed = [r"$\phi=1$", r"$L_\mu=0.0005$"]   # b = kappa by default (arXiv App. F-B), not listed
    rows = [(r"Main grid", pair_lines + [rf"$\delta\in{_set(m['delta'])}$", rf"$K\in{_set(m['K'], kfmt)}$", four],
             main_fixed, dict(study="Main grid", sets=m))]
    b = S["base"]
    assert set(b["methods"]) == set(BASELINES), b["methods"]
    shared = sorted(set.intersection(*(set(v) for v in b["alpha0"].values())))
    extra = sorted(set(b["alpha0"]["ACI"]) - set(shared))
    assert all(set(v) == set(shared) for k_, v in b["alpha0"].items() if k_ != "ACI"), b["alpha0"]
    rows.append((r"Baselines",
                 [r"methods: Unfiltered, $c=0$, $c_{\max}$, ACI, QT", rf"$\alpha_0\in{_set(shared)}$",
                  *((rf"ACI also at $\alpha_0\in{_set(extra)}$", r"with $L_\mu=0.0005$ (in the repository)") if target == "report"
                    else (rf"ACI also at $\alpha_0\in{_set(extra)}$ with $L_\mu=0.0005$",)),
                  rf"$L_\mu\in{_set(b['L_mu'])}$"],
                 [r"no certification parameters"], dict(study="Baselines", sets=b)))
    if not arxiv:
        t = SIX[0][1]
        rows.append((r"Default setting",
                     [rf"$(\alpha_0,\alpha,\delta,K)=({num(t['alpha0'])},{num(t['alpha'])},{num(t['delta'])},"
                      rf"{kfmt(t['K'])})$", four],
                     main_fixed + [r"(part of the main grid)"], dict(study="Default setting")))
    top = [(rf"\Tname{{{s[0][1:]}}}:", rf"$({num(s[1]['alpha0'])},{num(s[1]['alpha'])},{num(s[1]['delta'])},{kfmt(s[1]['K'])})$, "
                         rf"{VSHORT[s[2]]}") for s in SIX[1:]]
    rows.append((r"Selected settings", [r"settings $(\alpha_0,\alpha,\delta,K)$:"] + top + [four],
                 main_fixed + [r"(part of the main grid)"], dict(study="Selected settings")))
    at = r"settings: \Tname{1}--\Tname{5}" if arxiv else r"settings: default setting, \Tname{1}--\Tname{5}"
    nb, nd, nk = S["budget_normalizer_200"], S["stress_drift_200"], S["detailed_200"]
    for s_ in (nb, nd, nk):
        assert s_["methods"] == sorted(VARIANTS), s_["methods"]
    rows += [
        (r"Budget normalizer $b$", [rf"$b\in{_set(nb['b'], b_value)}$", at, four],
         [r"$K$ of the setting", r"$\phi=1$", r"$L_\mu=0.0005$"], dict(study="Budget normalizer", sets=nb)),
        (r"Bound factor $\phi$ and distribution-drift bound $L_\mu$",
         [rf"$\phi\in{_set(nd['phi'], _plain)}$", rf"$L_\mu\in{_set(nd['L_mu'])}$", r"($\phi$ and $L_\mu$ crossed)",
          at, four],
         [r"$K$ of the setting"], dict(study="Bound factor", sets=nd)),
        (r"Test period $K$ (detailed checks)",
         [r"$K$ of the setting and $K=\infty$", rf"$L_\mu\in{_set(nk['L_mu'])}$", r"($K$ and $L_\mu$ crossed)", at, four],
         [r"$\phi=1$"], dict(study="Test period", sets=nk)),
    ]
    return rows


def design_fragment(D, target, gid, arxiv=False):
    rows = design_rows(D, arxiv, target)
    frag = "arxiv_design_overview" if arxiv else "design_overview"
    head_ = [r"Comparison", r"Varied values", r"Fixed values"]
    # p-columns as fractions of the inner width; every body cell is a [c]-box, so the cells of a row share their
    # vertical centre (comparison names centred vertically against the varied-value block)
    fr = [0.22, 0.52, 0.26]
    n = len(fr)
    inner = rf"(\linewidth-2\vnbox-{2 * n}\tabcolsep-{n - 1}\arrayrulewidth)"
    spec = r"!{\vrule width\vnbox}" + "|".join(
        rf">{{\raggedright\arraybackslash}}p{{\dimexpr {f}\dimexpr{inner}\relax\relax}}"
        for f in fr) + r"!{\vrule width\vnbox}"
    out = [rf"{{\setlength{{\tabcolsep}}{{4pt}}\noindent\begin{{tabular}}{{{spec}}}\Xhline{{\vnbox}}"]
    out.append(" & ".join(head_[:-1]) + rf" & {head_[-1]}\rule[-1.1ex]{{0pt}}{{3.7ex}}" + r"\\\Xhline{0.6pt}")
    for i, (a, varied, fixed, rec) in enumerate(rows):
        if i:
            out.append(rf"\Xhline{{{RULE_METHOD}}}")
        vc, fc = lines_cell(varied), lines_cell(fixed)
        out.append(rf"\parbox[c]{{\linewidth}}{{\raggedright {a}}} & {vc} & {fc}" + r"\\")
        REC.add(target=target, fragment=frag, table=f"{gid}", row=dict(study=rec["study"]), method=None,
                spec="design", printed=vc, kind="design_values", fixed=fc,
                sets=json.loads(json.dumps(rec.get("sets"), default=str)) if rec.get("sets") else None)
    out.append(r"\Xhline{\vnbox}\end{tabular}}\par")
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------------------------------- build
# The totals of the assumption checks and the main-grid events that arXiv App. F-D states, as a macro file (no
# table). Scope of the totals: the four variants of every
# configuration with phi = 1 in the paper and the report, i.e. main_grid_200 (all), detailed_200,
# budget_normalizer_200 and stress_drift_200 at the default anchor (definition.anchor_id == "current_default") of the
# final export, and every view of new_anchors (new_top5_role1..5); deduplicated by configuration_id (equal ids carry
# equal check fields, asserted). The old-anchor sweeps candidate_1..5 of the final export are not reported and are
# excluded. The text states that every failure count is zero: asserted here. check_tables.py (check_assumption_totals)
# re-derives every number from the two JSON files with its own code.
AC_FAMILIES = ("detailed_200", "budget_normalizer_200", "stress_drift_200")
AC_FLAGS = ("support_excess_observed", "fallback_infeasibility_observed", "filter_residual_failure_observed",
            "terminal_fallback_feasibility_unknown", "missing_tail")
AC_MACROS = (("acCfg", "cfg"), ("acTraj", "traj"), ("acExec", "exec"), ("acChecked", "checked"),
             ("acNoChk", "nochk"), ("acUnresDep", "unres_dep"), ("acUnresCert", "unres_cert"), ("acCert", "cert"))
GC_MACROS = (("gcDep", "dep"), ("gcMissFin", "miss_fin"), ("gcMissInf", "miss_inf"), ("gcMissInfLow", "miss_inf_low"),
             ("gcMissInfCorr", "miss_inf_corr"))


def thousands(n):
    """An integer with the LaTeX thousands separator {,} (3240 -> 3{,}240)."""
    return f"{int(n):,}".replace(",", "{,}")


def ac_scope(D):
    """configuration_id -> cell of every configuration of a variant in the assumption-check scope (phi = 1 and
    phi < 1)."""
    keep = {}
    for k, f in D.idx.items():
        src, fam = k[0], k[1]
        c = f["cell"]
        if c["configuration"].get("method") not in VARIANTS:
            continue
        if not (src == "new" or (src == "final" and (fam == "main_grid_200" or
                                                            (fam in AC_FAMILIES and f["anchor_id"] == "current_default")))):
            continue
        cid = c["configuration_id"]
        if cid in keep:
            for fld in ("risk_increase_checks", "numerical_ambiguity", "assumption_observed_flags", "n",
                        "assumption_assessed_population"):
                assert keep[cid][fld] == c[fld], (cid, fld)
            continue
        keep[cid] = c
    return keep


def ac_totals(cells):
    t = dict(cfg=0, traj=0, exec=0, checked=0, checked_pw=0, exc_a2=0, exc_pw=0, exc_traj=0, nochk=0, nochk_vc=0,
             flags=0, num_err=0, unres_dep=0, unres_cert=0, cert=0, goal=0, contact=0, infeas=0)
    for c in cells:
        r = c["risk_increase_checks"]
        na = c["numerical_ambiguity"]
        t["cfg"] += 1
        t["traj"] += c["n"]
        t["exec"] += c["assumption_assessed_population"]["executed_decisions"]
        t["checked"] += r["direct"]["assessed_decisions"]
        t["checked_pw"] += r["pairwise"]["assessed_decisions"]
        t["exc_a2"] += r["direct"]["excess_decisions"]
        t["exc_pw"] += r["pairwise"]["excess_decisions"]
        t["exc_traj"] += r["direct"]["any_excess"] + r["pairwise"]["any_excess"]
        nochk = r["direct"]["assessed_decisions"] == 0
        t["nochk"] += nochk
        t["nochk_vc"] += nochk and c["configuration"]["method"] == "V-C"
        t["flags"] += sum(int(c["assumption_observed_flags"][x]) for x in AC_FLAGS)
        t["num_err"] += na["risk_numerical_errors"]
        t["unres_dep"] += na["deployment_intervals_unresolved"]
        t["unres_cert"] += na["certificate_events_unresolved"]
        stat = c["statistics"]["n_certificates"]
        assert stat["observed_n"] == c["n"] and not stat["missing_n"], c["configuration_id"]
        certificates = stat["mean"] * stat["observed_n"]
        assert abs(certificates - round(certificates)) < 1e-6, certificates
        t["cert"] += round(certificates)
        for o in ("goal", "contact", "infeas"):
            t[o] += c["outcomes"][o]["count"]
    return t


def main_grid_events(D):
    """The 3000 main-grid cells (375 settings x 4 variants x 2 scenarios, L_mu = 0.0005): cells with a
    deployment event (possible count > 0), and cells where the count of q_cert or q_dep misses its level (printed
    one-sided 95% bound above 100 x level, the rule of event_cell), split by finite K and K = infinity."""
    g = dict(cells=0, dep=0, dep_a0_max=0.0, miss_fin=0, miss_inf=0, miss_inf_low=0, miss_inf_corr=0)
    for k, f in D.idx.items():
        if k[:2] != ("final", "main_grid_200"):
            continue
        c = f["cell"]
        cf = c["configuration"]
        if cf["method"] not in VARIANTS:
            continue
        assert cf["L_mu"] == 0.0005 and cf["bound_scale"] == 1.0 and cf["budget_normalizer"] == 2.11, cf
        g["cells"] += 1
        te = c["trajectory_events"]
        if te["qdep"]["possible"] > 0:
            g["dep"] += 1
            g["dep_a0_max"] = max(g["dep_a0_max"], cf["alpha0"])
        level = guarantee_level(cf)
        miss = []
        for ev in ("qcert", "qdep"):
            e = te[ev]
            if e["applicable_n"] and Decimal(upper_pct(e["upper95_possible"])) > Decimal(repr(100 * level)):
                miss.append((ev, e["possible"], e["applicable_n"]))
        if not miss:
            continue
        if cf["K"] != "infinity":
            g["miss_fin"] += 1
            assert [(m[0], m[1], m[2]) for m in miss] == [("qcert", 2, NMC)] and cf["benchmark"] == "corridor" \
                and (cf["alpha"], cf["delta"], cf["K"]) == (0.15, 0.025, 10), (cf, miss)
        else:
            g["miss_inf"] += 1
            if cf["delta"] == 0.025:
                g["miss_inf_low"] += 1
            else:
                assert cf["benchmark"] == "corridor" and cf["alpha"] == 0.05 and cf["delta"] in (0.05, 0.10), (cf, miss)
                g["miss_inf_corr"] += 1
    n_inf_low = sum(1 for k, f in D.idx.items() if k[:2] == ("final", "main_grid_200")
                    and f["cell"]["configuration"]["method"] in VARIANTS
                    and f["cell"]["configuration"]["K"] == "infinity" and f["cell"]["configuration"]["delta"] == 0.025)
    assert g["miss_inf_low"] == n_inf_low, (g, n_inf_low)   # every K = infinity cell with delta = 0.025 misses
    assert g["cells"] == 3000 and g["dep_a0_max"] <= 0.25, g
    return g


def assumption_checks_fragment(D):
    """appendices/generated_n200/paper/assumption_checks_arxiv.tex: the assumption-check totals and the main-grid
    event counts of arXiv App. F-D as \\providecommand/\\renewcommand pairs (no table, prints nothing)."""
    frag = "assumption_checks_arxiv"
    scope = ac_scope(D)
    all_t = ac_totals(scope.values())
    phi1 = [c for c in scope.values() if c["configuration"]["bound_scale"] == 1.0]
    t = ac_totals(phi1)
    fails = {x: t[x] for x in ("exc_a2", "exc_pw", "exc_traj", "nochk_vc", "flags", "num_err", "contact", "infeas")}
    assert all(v == 0 for v in fails.values()), fails   # the text states these zeros
    assert t["goal"] == t["traj"] and t["checked"] == t["checked_pw"], t
    g = main_grid_events(D)
    print(f"assumption-check totals (phi = 1): {t}; all phi: {all_t['cfg']} configurations; main-grid events: {g}")
    out = ["% Generated by tools/n200_tables/build_tables.py; do not edit by hand. Assumption-check totals (arXiv "
           "App. F-D).",
           "% Data: N200_numbers/{final,new_anchors}; values traced in tables_values.json.gz."]
    for mac, fld in AC_MACROS:
        printed = thousands(t[fld])
        REC.add(target="paper", fragment=frag, kind="a6_total", name=mac, field=fld, value=t[fld], printed=printed,
                scope="phi=1, variants, main_grid_200 + default-anchor detailed/budget/stress + new_anchors")
        out.append(rf"\providecommand{{\{mac}}}{{}}\renewcommand{{\{mac}}}{{{printed}}}")
    out.append("% Main grid (final main_grid_200, 3000 grid cells, L_mu = 0.0005): deployment events and missed levels.")
    for mac, fld in GC_MACROS:
        printed = thousands(g[fld])
        REC.add(target="paper", fragment=frag, kind="a6_total", name=mac, field=fld, value=g[fld], printed=printed,
                scope="main_grid_200, variants, L_mu=0.0005")
        out.append(rf"\providecommand{{\{mac}}}{{}}\renewcommand{{\{mac}}}{{{printed}}}")
    return "\n".join(out) + "\n"


RAW_FRAGMENTS = {"assumption_checks_arxiv"}   # written without HEADER (the fragment carries its own header)


def build_all(D):
    """All fragments of both targets (main() moves any other .tex file in the output folders to _superseded/)."""
    frags = {t: {} for t in TARGETS}
    for target in TARGETS:
        T = target[0]
        F = frags[target]
        F["table_key"] = table_key_fragment(target)
        F["design_overview"] = design_fragment(D, target, f"dz{T}")
        for scen in SCENS:
            s = scen[0]
            F[f"default_{scen}"] = default_fragment(D, scen, target, f"df{s}{T}")
            F[f"topfive_full_{scen}"] = topfive_full_fragment(D, scen, target, f"tg{s}{T}")
            for kind, tag in (("delta", "dl"), ("period", "pr"), ("b", "bn"), ("phi", "ph")):
                F[f"{kind}_{scen}"] = param_fragment(D, scen, kind, target, f"{tag}{s}{T}")
            F[f"detailed_events_{scen}"] = detailed_events_fragment(D, scen, target, f"de{s}{T}")
            for lm, lt in ((0.0005, "L0005"), (0.005, "L005")):
                F[f"checks_{scen}_{lt}"] = checks_fragment(D, scen, lm, target, f"ck{s}{lt[1:]}{T}")
            F[f"variation_{scen}"] = variation_fragment(D, scen, target, f"va{s}{T}")
            F[f"n30_record_{scen}"] = n30_fragment(D, scen, target, f"nt{s}{T}")   # not printed
    # paper only: Section VI and arXiv
    P = frags["paper"]
    P["sec6_short"] = sec6_short_fragment(D, "sxs")   # reference of sec6_combined (not printed)
    P["sec6_maze_arxiv"] = sec6_maze_arxiv_fragment(D, "sxm")   # Maze-only Section VI table (not printed)
    P["main_maze_arxiv_v8"] = main_arxiv_v8_fragment(D, "maze")   # arXiv Table I: one column at 245.7 pt (units in the caption)
    P["main_corridor_arxiv_v8"] = main_arxiv_v8_fragment(D, "corridor")   # arXiv Table IV (App. F-C)
    P["table_key_arxiv"] = table_key_arxiv_fragment()
    for scen in SCENS:
        P[f"secondary_{scen}_arxiv"] = secondary_arxiv_fragment(D, scen)
        P[f"topfive_full_{scen}_arxiv_percent"] = topfive_full_fragment(D, scen, "paper", f"tfp{scen[0]}", percent=True)
        P[f"topfive_full_{scen}_onepage"] = topfive_full_onepage_fragment(D, scen)
        P[f"variation_{scen}_arxiv"] = variation_fragment(D, scen, "paper", f"var{scen[0]}", arxiv=True)
    # Section VI table with the default setting, S1-S5 and S* (not printed)
    P["sec6_combined"] = sec6_combined_fragment(D, "sxc")
    P["main_maze_extended_v8"] = main_arxiv_v8_fragment(D, "maze", extended=True)
    P["main_corridor_extended_v8"] = main_arxiv_v8_fragment(D, "corridor", extended=True)
    frags["report"]["star_maze"] = star_fragment("sstr", target="report")
    P["star_maze"] = star_fragment("sst")   # extended version
    P["arxiv_design_overview"] = design_fragment(D, "paper", "dza", arxiv=True)
    for scen in SCENS:
        P[f"checks_{scen}_L0005_nodef"] = checks_fragment(D, scen, 0.0005, "paper", f"ckn{scen[0]}", nodef=True)
        P[f"checks_{scen}_arxiv"] = checks_arxiv_fragment(D, scen, f"ckx{scen[0]}")      # not printed
        P[f"stress_{scen}_arxiv"] = stress_arxiv_fragment(D, scen, f"stx{scen[0]}")      # not printed
        P[f"stress_{scen}_arxiv_percent"] = stress_arxiv_fragment(D, scen, f"stp{scen[0]}", percent=True)
        # the same cells for report Tables 17 and 18 (report target only)
        frags["report"][f"stress_{scen}_cases"] = stress_arxiv_fragment(D, scen, f"sk{scen[0]}r", target="report")
    P["stress_representative_arxiv_percent"] = stress_representative_arxiv_fragment(D)
    P["assumption_checks_arxiv"] = assumption_checks_fragment(D)   # assumption-check totals of arXiv App. F-D (macros)
    return frags


def collisions(D):
    """(environment, method name as printed) of every method with contact in at least one trajectory of any cell
    (outcomes.contact.count of the exports); the colour switch marks their names."""
    names = dict(BNAME) | VSHORT
    out = set()
    for f in D.idx.values():
        c = f["cell"]
        cnt = ((c.get("outcomes") or {}).get("contact") or {}).get("count")
        if cnt:
            cf = c["configuration"]
            meth = cf["method"]
            assert meth in names, f"collision of a method without a table name: {meth}"
            out.add((SCEN[cf["benchmark"]]["name"], re.sub(r"\s+", "", names[meth])))
    return out


def write_gz(path, text):
    """A deterministic gzip file (no file name, mtime 0, level 9): the same text gives the same bytes."""
    import gzip
    import io
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, compresslevel=9, mtime=0) as fh:
        fh.write(text.encode())
    Path(path).write_bytes(buf.getvalue())


HEADER = ("% Generated by tools/n200_tables/build_tables.py ({target}, {width}); do not edit by hand.\n"
          "% Data: N200_numbers/{{final,new_anchors}}; values traced in tables_values.json.gz.\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", type=Path, help="optional pickle of the loaded exports, written on the first run and "
                                               "read on later runs (it must match the exports; faster rebuilds)")
    ap.add_argument("--catalog", action="store_true",
                    help="after the edition fragments (unchanged), also write the catalog fragments "
                         "to appendices/generated_n200/catalog/{paper,report}/ (catalog_tables.py)")
    a = ap.parse_args()
    import pickle
    if a.cache and a.cache.exists():
        D = pickle.loads(a.cache.read_bytes())
        for k, p in FILES.items():  # the cache is valid only for unchanged exports
            assert hashlib.sha256(p.read_bytes()).hexdigest() == D.sha[k], f"{p} changed; delete the cache"
    else:
        D = Data()
        if a.cache:
            a.cache.write_bytes(pickle.dumps(D, protocol=5))
    chk = 0
    for k, f in list(D.idx.items())[:300]:
        e = f["cell"]["realized_transition_events"]["q_CBF"]
        if e["assessed_n"]:
            assert abs(cp_upper(e["count"], e["assessed_n"]) - e["upper95"]) < 1e-9
            chk += 1
    print(f"Clopper-Pearson check against stored upper95: {chk} cells OK")
    frags = build_all(D)
    widths = {"paper": "505.89pt \\textwidth, 10pt", "report": "505.89pt \\textwidth, 11pt Times"}
    import shutil
    import table_colour   # colour switch: coloured and plain variant per table (\ifcolortables, default plain)
    table_colour.COLLISIONS = collisions(D)
    print("colour switch: collisions " + ", ".join(f"{e}/{m}" for e, m in sorted(table_colour.COLLISIONS)))
    colour_stats, verdicts = {}, {}
    for target in TARGETS:   # first pass: guarantee verdicts, shared by all fragments
        for name in table_colour.COLOUR_TABLES:
            if name in frags[target]:
                table_colour.collect(frags[target][name], verdicts, name)
    for target in TARGETS:
        out = GEN / target
        out.mkdir(parents=True, exist_ok=True)
        # .tex files of fragments not in this build: moved (never deleted) into _superseded/
        for old in sorted(out.glob("*.tex")):
            if old.stem not in frags[target]:
                sup = out / "_superseded"
                sup.mkdir(exist_ok=True)
                dst = sup / old.name
                if dst.exists():
                    dst = sup / f"{old.stem}_{hashlib.sha256(old.read_bytes()).hexdigest()[:8]}.tex"
                shutil.move(str(old), str(dst))
                print(f"  superseded: {target}/{old.name}")
        for name, tex in frags[target].items():
            if name in table_colour.COLOUR_TABLES:
                tex = table_colour.wrap(name, tex, colour_stats.setdefault(f"{target}/{name}", {}), target,
                                        verdicts, records=REC.entries)
            head_ = "" if name in RAW_FRAGMENTS else HEADER.format(target=target, width=widths[target])
            (out / f"{name}.tex").write_text(head_ + tex)
    values = dict(generator="tools/n200_tables/build_tables.py", inputs={k: rel_input(v) for k, v in FILES.items()},
                  sha256=D.sha, candidates=rel_input(CANDIDATES),
                  rounding=dict(rule="half-up after rounding the stored double to d + 9 decimals",
                                event_upper="ceil to 0.1 %",
                                temporal_decimals=TMP_DEC),
                  count_rule="possible events (confirmed plus numerically undecided); n shown when assessed != 200",
                  check_rule="variants, qcert/qdep only: printed (rounded-up) bound <= 100*level, level = "
                             "delta (finite K) or delta/(2 b ln^2 2) (K = infinity)",
                  fragments={t: sorted(frags[t]) for t in TARGETS}, cells=REC.entries)
    write_gz(HERE / "tables_values.json.gz", json.dumps(values, default=str, separators=(",", ":")))
    print("colour switch: " + ", ".join(f"{k} {v.get('painted', 0)}+{v.get('collision', 0)}x"
                                        for k, v in sorted(colour_stats.items())))
    print(f"{len(REC.entries)} printed cells recorded; fragments: paper {len(frags['paper'])}, report {len(frags['report'])}")
    if a.catalog:   # separate fragments, values file and manifest; the edition outputs above are final
        import catalog_tables
        REC.entries = []
        global FIRST_CUT
        FIRST_CUT = catalog_tables.FIRST_CUT
        cat = catalog_tables.build(sys.modules[__name__], D)
        catalog_tables.write(sys.modules[__name__], D, cat, REC.entries)


if __name__ == "__main__":
    main()
