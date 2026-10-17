#!/usr/bin/env python3
"""v10 table generator for the N200 result tables.

v10 = v9 with these formatting rules: baseline headers read "<Scenario>  Baselines  <settings>";
method names, the alpha0 label and all values centred; thin rules between methods in the paired baseline box;
a small extra depth (\vnsd) only on rows whose subscripts sit on a rule. Data: final/common/metrics.json
(N200 view IDs; each printed cell also records definition.n100_view_id).

Reads the N200 common metrics export and writes LaTeX fragments in the boxed layout:
  (1) variant sweep tables per (scenario, alpha0), nested alpha -> delta -> K -> [V-U, V-Delta, V-L, V-C],
      seven paired metric columns plus T_mis, two value lines per method;
  (2) a baseline box on every table page (Unfiltered, c=0, c_max, ACI, QT), same column pairs, all-N-A columns dropped;
  (3) L_mu, stress-factor (phi) and budget-normalizer (b) tables at the default and the five N30 anchors;
  (4) the unified checks table from detailed_200;
  (5) the default plus the five candidate settings selected on the first 30 trajectories (arXiv summary).
Every printed cell is listed in tables_values.json with its view_id, configuration and source statistic;
check_tables.py re-reads metrics.json and asserts the printed strings.

Rounding (half-up on the stored double, as papertables.fmt_half_up): normalized quantities m_c/c_max, m_Y/c_max,
control intervention and certificate counts 2 decimals; percentages (fractions x 100) 2 decimals; seconds 2 decimals;
metres 3 decimals (c_max is 0.01-0.02 m, so 2 decimals would hide margin-scale differences); event upper bounds in %
rounded UP to 1 decimal (paper Table I, e.g. 0/200 (1.5)).

Usage:
  python3 build_tables.py                       # sample build (alpha0 = 0.15, alpha in {0.05,0.10}, delta in {0.05,0.10})
  python3 build_tables.py --sweep-full 0.30     # also emit the complete sweep table(s) for alpha0 = 0.30
  python3 build_tables.py --no-compile          # fragments and JSON only
Writes only inside its own folder (or --out). Inputs: ECBF_N200_NUMBERS (default <repository>/N200_numbers) and
ECBF_CANDIDATES (default <repository>/candidate_review.json), both not included.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from decimal import ROUND_HALF_UP, ROUND_CEILING, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_tables/<this file>)
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports (not included)
DEF_METRICS = NUM / "final/common/metrics.json"
DEF_CANDIDATES = Path(os.environ.get("ECBF_CANDIDATES", REPO / "candidate_review.json"))   # candidate review


def recorded_path(p):
    """Path written into the values file: an export as "N200_numbers/<path relative to NUM>", the candidate review
    as "candidate_review.json", any other path as given."""
    p = Path(p).resolve()
    if p == DEF_CANDIDATES.resolve():
        return "candidate_review.json"
    try:
        return "N200_numbers/" + p.relative_to(NUM.resolve()).as_posix()
    except ValueError:
        return str(p)

# ----------------------------------------------------------------------------------------------- notation
SCEN = {
    "corridor": dict(name="Corridor", cmax="0.02", T="0.1"),
    "maze": dict(name="Maze", cmax="0.01", T="0.2"),
}
VARIANTS = ["V-U", "V-D", "V-Cd", "V-C"]
VNAME = {"V-U": r"\textcolor{vu}{V-U}", "V-D": r"\textcolor{vd}{V-$\Delta$}",
         "V-Cd": r"\textcolor{vl}{V-L}", "V-C": r"\textcolor{vc}{V-C}"}
BASELINES = ["Unfiltered", "Zero margin", "Conservative", "ACI", "QT"]
BNAME = {"Unfiltered": "Unfiltered", "Zero margin": "$c=0$", "Conservative": "$c_{\\max}$", "ACI": "ACI", "QT": "QT"}
DEFAULT = dict(alpha0=0.15, alpha=0.10, delta=0.05, K=50)
ANCHORS = ["current_default", "candidate_1", "candidate_2", "candidate_3", "candidate_4", "candidate_5"]
DELTAS = [0.025, 0.05, 0.10, 0.15, 0.20]
KS = [10, 50, 100, 200, "infinity"]
ALPHA0S = [0.10, 0.15, 0.20, 0.25, 0.30]

# metric id -> (symbol, arrow, unit, scale, decimals, boldable)
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
}
# event id -> (symbol, source)
EVT = {
    "qcert": (r"\widehat q_{\rm cert}", ("trajectory_events", "qcert")),
    "qdep": (r"\widehat q_{\rm dep}", ("trajectory_events", "qdep")),
    "rCBF": (r"\widehat q_{\rm CBF}", ("realized_transition_events", "q_CBF")),
    "runsafe": (r"\widehat q_{\rm unsafe}", ("realized_transition_events", "q_unsafe")),
    "qA2": (r"\widehat q_{\rm A2}", ("risk_increase_checks", "direct")),
    "qpw": (r"\widehat q_{\rm pw}", ("risk_increase_checks", "pairwise")),
}
PAIRS = [["m_c", "m_Y"], ["m_int", "n_certificates"], ["q_CBF", "q_Y"], ["q_use", "q_unsafe"],
         ["q_nocert", "q_notrans"], ["h_min", "clearance_min"], ["tau_use", "p_max"], ["T_mis"]]
RULE_METHOD = "0.25pt"


# ----------------------------------------------------------------------------------------------- formatting
def half_up(v: float, d: int) -> str:
    q = Decimal(f"{v:.{d + 9}f}").quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    out = f"{q:.{d}f}"
    return "0" + out[2:] if out.startswith("-0") and q == 0 else out


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


# ----------------------------------------------------------------------------------------------- data
class Data:
    def __init__(self, path: Path):
        print(f"loading {path} ...", flush=True)
        raw = path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        d = json.loads(raw)
        del raw
        self.idx = {}
        for v in d["views"]:
            fam = v["definition"]["family"]
            for c in v["cells"]:
                for s in c.get("statistics", {}).values():
                    s.pop("observed_seeds", None)
                c.pop("seeds", None)
                cfg = c["configuration"]
                k = (fam,) + tuple(sorted((a, key(b)) for a, b in cfg.items() if b is not None))
                assert k not in self.idx, k
                self.idx[k] = (v["view_id"], v["definition"], c)
        self.anchor = {}
        for (fam, *kv), (vid, df, c) in self.idx.items():
            if df.get("anchor_id"):
                cfg = c["configuration"]
                t = (key(cfg["alpha0"]), key(cfg["alpha"]), key(cfg["delta"]))
                self.anchor.setdefault(df["anchor_id"], set()).add((fam, t, key(cfg["K"])))
        print(f"  {len(self.idx)} cells", flush=True)

    def get(self, fam, **cfg):
        k = (fam,) + tuple(sorted((a, key(b)) for a, b in cfg.items() if b is not None))
        return self.idx.get(k)


# ----------------------------------------------------------------------------------------------- cells
class Recorder:
    def __init__(self):
        self.entries = []

    def add(self, **kw):
        self.entries.append(kw)


REC = Recorder()


def stat_cell(found, metric, ctx):
    """Return (latex, info) for a mean +- SD cell; info carries provenance."""
    if found is None:
        return "--", dict(kind="not_simulated")
    vid, df, c = found
    s = c["statistics"][metric]
    sym, arrow, unit, scale, dec, boldable = MET[metric]
    info = dict(view_id=vid, n100_view_id=df.get("n100_view_id"), family=df["family"],
                configuration_id=c.get("configuration_id"), configuration=c["configuration"],
                source=f"statistics.{metric}", scale=scale, decimals=dec)
    if not s["applicable"] or s["observed_n"] == 0 or s["mean"] is None:
        return r"\na", dict(info, kind="NA", applicable=s["applicable"], observed_n=s["observed_n"])
    m, sd = half_up(s["mean"] * scale, dec), half_up(s["sd"] * scale, dec)
    info.update(kind="meanSD", mean_str=m, sd_str=sd, conditional=s["conditional"],
                boldable=bool(boldable and not s["conditional"]), arrow=arrow)
    return None, info


def event_cell(found, ev):
    if found is None:
        return "--", dict(kind="not_simulated")
    vid, df, c = found
    grp, sub = EVT[ev][1]
    e = c[grp][sub]
    info = dict(view_id=vid, n100_view_id=df.get("n100_view_id"), family=df["family"],
                configuration_id=c.get("configuration_id"), configuration=c["configuration"], source=f"{grp}.{sub}")
    if grp == "trajectory_events":
        n, lo, hi, u = e["applicable_n"], e["confirmed"], e["possible"], e["upper95_possible"]
    elif grp == "realized_transition_events":
        n, lo, hi, u = e["assessed_n"], e["count"], e["count"], e["upper95"]
    else:
        n, lo, hi = e["assessed_trajectory"], e["any_excess"], e["any_excess"]
        u = cp_upper(hi, n) if n else None
    if not n:
        return r"\na", dict(info, kind="NA", assessed=n)
    cnt = f"{lo}" if lo == hi else f"{lo}--{hi}"
    s = f"{cnt}/{n} ({upper_pct(u)})"
    return s, dict(info, kind="event", count_lo=lo, count_hi=hi, n=n, upper95=u)


def fill_stat_rows(rows):
    """rows: list of dicts with 'found', 'is_variant', 'group' (bold group key). Adds rows[i]['cells'][metric]."""
    groups = {}
    for r in rows:
        r.setdefault("cells", {})
        for pair in PAIRS:
            for m in pair:
                tex, info = stat_cell(r["found"], m, r)
                r["cells"][m] = [tex, info]
        if r["is_variant"]:
            groups.setdefault(r["group"], []).append(r)
    for g, rs in groups.items():
        for m in MET:
            infos = [r["cells"][m][1] for r in rs]
            if not all(i.get("kind") == "meanSD" and i["boldable"] for i in infos) or len(infos) < 2:
                continue
            vals = [Decimal(i["mean_str"]) for i in infos]
            best = max(vals) if infos[0]["arrow"] == "up" else min(vals)
            for i, v in zip(infos, vals):
                i["bold"] = v == best
    for r in rows:
        for m, (tex, info) in r["cells"].items():
            if tex is None:
                mm = rf"\mathbf{{{info['mean_str']}}}" if info.get("bold") else info["mean_str"]
                r["cells"][m][0] = f"${mm}\\vpm {info['sd_str']}$"
                info.setdefault("bold", False)


def fill_event_rows(rows, evs):
    for r in rows:
        r.setdefault("cells", {})
        for ev in evs:
            tex, info = event_cell(r["found"], ev)
            r["cells"][ev] = [tex, info]


# ----------------------------------------------------------------------------------------------- headers
def head(m):
    if m in MET:
        sym, arrow, unit, *_ = MET[m]
    else:
        sym, arrow, unit = EVT[m][0], "down", None
    a = {"down": r"\!\downarrow", "up": r"\!\uparrow", None: ""}[arrow]
    u = {None: "", "%": r"\,(\%)", "m": r"\,(m)", "s": r"\,(s)"}[unit]
    return f"${sym}{a}${u}"


# ----------------------------------------------------------------------------------------------- LaTeX tables
class Widths:
    """Collects every string per column so LaTeX can measure the widest one (fixed, aligned column widths)."""

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
    """Column layout shared by a main table and its baseline boxes in one fragment.

    level_cols: list of (col_id, header); value_cols: list of metric-id lists (1 or 2 entries)."""

    def __init__(self, gid, level_cols, value_cols, nlines, rule_weights, method_header="Method"):
        self.gid, self.level_cols, self.value_cols, self.nlines = gid, level_cols, value_cols, nlines
        self.method_header = method_header
        self.rule_weights = rule_weights
        self.W = Widths(gid)
        self.ncol = len(level_cols) + 1 + len(value_cols)
        # header strings
        for cid, h in level_cols:
            self.W.add(cid, h)
        self.W.add("meth", method_header)
        for j, vc in enumerate(value_cols):
            for m in vc:
                self.W.add(f"v{j}", head(m))

    def w(self, col):
        return rf"\vnw{{{self.W.name(col)}}}"

    def sep_tex(self):
        """\\tabcolsep that makes the main table exactly \\linewidth wide."""
        widths = " - ".join([self.w(c) for c, _ in self.level_cols] + [self.w("meth")] +
                            [self.w(f"v{j}") for j in range(len(self.value_cols))])
        n_inner = len(self.value_cols)  # rule after Method and between value columns
        return (rf"\vnnewsep{{{self.gid}}}{{\dimexpr(\linewidth - {widths} - {n_inner}\arrayrulewidth"
                rf" - 2\vnbox)/{2 * self.ncol}\relax}}")


def rule_between(lay, prev, cur, nlev, method_col):
    """Horizontal rule between two method rows according to the level hierarchy."""
    for k in range(nlev):
        if prev["levels"][k] != cur["levels"][k]:
            w = lay.rule_weights[k]
            if k == 0:
                return rf"\Xhline{{{w}}}"
            return rf"\Xcline{{{k + 1}-{lay.ncol}}}{{{w}}}"
    return rf"\Xcline{{{method_col}-{lay.ncol}}}{{{RULE_METHOD}}}"


def level_spans(rows, nlev, nlines):
    spans = []
    for k in range(nlev):
        sp = [0] * len(rows)
        i = 0
        while i < len(rows):
            j = i
            while j + 1 < len(rows) and rows[j + 1]["levels"][:k + 1] == rows[i]["levels"][:k + 1]:
                j += 1
            sp[i] = (j - i + 1) * nlines
            i = j + 1
        spans.append(sp)
    return spans


def main_table(lay: Layout, rows, settings, keyrow, fragment, table_id):
    """Boxed main table. rows: dicts with levels (tuple of latex), level_disp (latex per level), name, cells."""
    gid, nlev, L = lay.gid, len(lay.level_cols), lay.nlines
    cs = [rf"w{{c}}{{{lay.w(c)}}}" for c, _ in lay.level_cols] + [rf"w{{c}}{{{lay.w('meth')}}}"]
    vs = [rf"w{{c}}{{{lay.w(f'v{j}')}}}" for j in range(len(lay.value_cols))]
    spec = r"!{\vrule width\vnbox}" + "".join(cs) + "|" + "|".join(vs) + r"!{\vrule width\vnbox}"
    out = [rf"{{\setlength{{\tabcolsep}}{{\vnsep{{{gid}}}}}\noindent\begin{{tabular}}{{{spec}}}\Xhline{{\vnbox}}"]
    full = rf"!{{\vrule width\vnbox}}l!{{\vrule width\vnbox}}"
    out.append(rf"\multicolumn{{{lay.ncol}}}{{{full}}}{{{settings}}}\\\hline")
    if keyrow:
        out.append(rf"\multicolumn{{{lay.ncol}}}{{{full}}}{{{keyrow}}}\\\hline")
    # header
    hl = 2 if any(len(v) == 2 for v in lay.value_cols) else 1
    for j in range(hl):
        cells = []
        for cid, h in lay.level_cols:
            cells.append((rf"\multirow{{{hl}}}{{*}}{{{h}}}" if hl > 1 else h) if j == 0 else "")
        mh = lay.method_header
        cells.append((rf"\multirow{{{hl}}}{{*}}{{{mh}}}" if hl > 1 else mh) if j == 0 else "")
        for vc in lay.value_cols:
            if len(vc) == 2:
                cells.append(head(vc[j]))
            else:
                cells.append((rf"\multirow{{{hl}}}{{*}}{{{head(vc[0])}}}" if hl > 1 else head(vc[0])) if j == 0 else "")
        out.append(" & ".join(cells) + r"\vnsd\\")  # header subscripts (c_max, h^min, unsafe) need depth
    out.append(r"\Xhline{0.6pt}")
    spans = level_spans(rows, nlev, L)
    for i, r in enumerate(rows):
        if i:
            out.append(rule_between(lay, rows[i - 1], r, nlev, nlev + 1))
        for j in range(L):
            cells = []
            for k in range(nlev):
                if j == 0 and spans[k][i]:
                    cells.append(rf"\multirow{{{spans[k][i]}}}{{*}}{{{r['level_disp'][k]}}}" if spans[k][i] > 1
                                 else r["level_disp"][k])
                    lay.W.add(lay.level_cols[k][0], r["level_disp"][k])
                else:
                    cells.append("")
            cells.append((rf"\multirow{{{L}}}{{*}}{{{r['name']}}}" if L > 1 else r["name"]) if j == 0 else "")
            if j == 0:
                lay.W.add("meth", r["name"])
            for vj, vc in enumerate(lay.value_cols):
                if len(vc) == 2 or j == 0:
                    m = vc[j] if len(vc) == 2 else vc[0]
                    tex, info = r["cells"][m]
                    lay.W.add(f"v{vj}", tex)
                    REC.add(fragment=fragment, table=table_id, row=r["record"], method=r["method"], metric=m,
                            printed=tex, **info)
                    cells.append(rf"\multirow{{{L}}}{{*}}{{{tex}}}" if (len(vc) == 1 and L > 1) else tex)
                else:
                    cells.append("")
            out.append(" & ".join(cells) + r"\\")
    out.append(r"\Xhline{\vnbox}\end{tabular}}\par")
    return "\n".join(out)


def baseline_table(lay: Layout, brows, settings, blevels, fragment, table_id):
    """Boxed baseline table aligned with the main table: left block = level columns + Method of the main table.

    blevels: list of (col_id, header) shown in the baseline box (e.g. alpha0 or L_mu); brows: rows with
    levels/level_disp of the same length. Value columns that are N-A for every baseline are dropped;
    a pair with one all-N-A member becomes a single centred value."""
    gid, L = lay.gid, lay.nlines
    # which metrics survive
    keep = []
    for vc in lay.value_cols:
        k = [m for m in vc if any(r["cells"][m][1].get("kind") in ("meanSD", "event") for r in brows)]
        if k:
            keep.append((lay.value_cols.index(vc), k))
    nb = len(blevels)
    # left block width: all main left columns + their separators, minus the baseline level columns
    left_main = " + ".join([lay.w(c) for c, _ in lay.level_cols] + [lay.w("meth")])
    nleft = len(lay.level_cols) + 1
    bl_widths = " - ".join([""] + [lay.w("B" + c) for c, _ in blevels]) if blevels else ""
    meth_w = rf"\dimexpr {left_main} + {2 * (nleft - 1 - nb)}\tabcolsep {bl_widths}\relax"
    cs = [rf"w{{c}}{{{lay.w('B' + c)}}}" for c, _ in blevels] + [rf"w{{c}}{{{meth_w}}}"]
    for cid, h in blevels:
        lay.W.add("B" + cid, h)
    for r in brows:
        lay.W.add("bname", r["name"])
    sep = r"!{\hskip\arrayrulewidth}"
    vs = [rf"w{{c}}{{{lay.w(f'v{j}')}}}" for j, _ in keep]
    spec = r"!{\vrule width\vnbox}" + "".join(cs) + sep + sep.join(vs) + r"!{\vrule width\vnbox}"
    ncol = nb + 1 + len(keep)
    out = [rf"{{\setlength{{\tabcolsep}}{{\vnsep{{{gid}}}}}\vncheckleft{{{meth_w}}}{{{lay.w('bname')}}}{{{table_id}}}"
           rf"\noindent\begin{{tabular}}{{{spec}}}\Xhline{{\vnbox}}"]
    full = rf"!{{\vrule width\vnbox}}l!{{\vrule width\vnbox}}"
    out.append(rf"\multicolumn{{{ncol}}}{{{full}}}{{{settings}}}\\\hline")
    hl = 2 if any(len(k) == 2 for _, k in keep) else 1
    for j in range(hl):
        cells = []
        for idx, (cid, h) in enumerate(blevels):
            c = (rf"\multirow{{{hl}}}{{*}}{{{h}}}" if hl > 1 else h) if j == 0 else ""
            cells.append(c)
        c = (r"\multirow{%d}{*}{Method}" % hl if hl > 1 else "Method") if j == 0 else ""
        cells.append(rf"\multicolumn{{1}}{{w{{c}}{{{meth_w}}}|}}{{{c}}}")
        for n_, (vj, k) in enumerate(keep):
            if len(k) == 2:
                c = head(k[j])
            else:
                c = (rf"\multirow{{{hl}}}{{*}}{{{head(k[0])}}}" if hl > 1 else head(k[0])) if j == 0 else ""
            last = n_ == len(keep) - 1
            cells.append(c if last else rf"\multicolumn{{1}}{{w{{c}}{{{lay.w(f'v{vj}')}}}|}}{{{c}}}")
        out.append(" & ".join(cells) + r"\vnsd\\")
    out.append(r"\Xhline{0.6pt}")
    spans = level_spans(brows, nb, L) if nb else []
    def rule_after(i):
        """Rule printed below baseline row i (None: no rule)."""
        if i == len(brows) - 1:
            return r"\Xhline{\vnbox}"
        if nb and brows[i + 1]["levels"] != brows[i]["levels"]:
            return r"\Xhline{0.6pt}"
        if L > 1:  # paired box: thin rule between methods, as in the variant table
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
            for vj, k in keep:
                if len(k) == 2 or j == 0:
                    m = k[j] if len(k) == 2 else k[0]
                    tex, info = r["cells"][m]
                    lay.W.add(f"v{vj}", tex)
                    REC.add(fragment=fragment, table=table_id, row=r["record"], method=r["method"], metric=m,
                            printed=tex, **info)
                    cells.append(rf"\multirow{{{L}}}{{*}}{{{tex}}}" if (len(k) == 1 and L > 1) else tex)
                else:
                    cells.append("")
            pad = j == L - 1 and "_" in r["name"] and rule_after(i)  # e.g. the c_max label above a rule
            out.append(" & ".join(cells) + (r"\vnsd\\" if pad else r"\\"))
    out.append(r"\Xhline{\vnbox}\end{tabular}}\par")
    dropped = [m for vc in lay.value_cols for m in vc if all(m not in k for _, k in keep)]
    return "\n".join(out), dropped


# ----------------------------------------------------------------------------------------------- row builders
def settings_line(scen, tag=None, **kv):
    """Bold scenario name, optional bold tag ("Baselines"), then the settings in math."""
    s = SCEN[scen]
    parts = []
    order = ["alpha0", "alpha", "delta", "K", "L_mu", "cmax", "T", "NMC", "b", "phi"]
    sym = dict(alpha0=r"\alpha_0", alpha=r"\alpha", delta=r"\delta", K="K", L_mu=r"L_\mu", cmax=r"c_{\max}", T="T",
               NMC=r"N_{\rm MC}", b="b", phi=r"\phi")
    kv = dict(kv)
    kv.setdefault("cmax", s["cmax"] + r"\,\mathrm{m}")
    kv.setdefault("T", s["T"] + r"\,\mathrm{s}")
    kv.setdefault("NMC", "200")
    for k in order:
        if k in kv and kv[k] is not None:
            v = kv[k]
            parts.append(f"{sym[k]}={v if isinstance(v, str) else num(v)}")
    joined = ",\\ ".join(parts)
    name = rf"\textbf{{{s['name']}}}" + (rf"\quad\textbf{{{tag}}}" if tag else "")
    return rf"{name}\quad ${joined}$"


def baseline_rows(D, scen, alpha0, L_mu, record_extra=None, level=None):
    rows = []
    for b in BASELINES:
        f = None
        if abs(alpha0 - 0.15) < 1e-9 or abs(alpha0 - 0.30) < 1e-9:
            f = D.get("matched_baselines_200", benchmark=scen, method=b, alpha0=alpha0, L_mu=L_mu, rho_buf=0.0)
        if f is None and b == "ACI" and abs(L_mu - 0.0005) < 1e-9:
            f = D.get("matched_ACI_200", benchmark=scen, method=b, alpha0=alpha0, L_mu=L_mu, rho_buf=0.0)
        rows.append(dict(found=f, is_variant=False, group=None, name=BNAME[b], method=b,
                         levels=level[0] if level else (), level_disp=level[1] if level else (),
                         record=dict(scenario=scen, alpha0=alpha0, L_mu=L_mu, **(record_extra or {}))))
    return rows


def variant_rows(D, fam, scen, levels, level_disp, record, **cfg):
    rows = []
    for v in VARIANTS:
        f = D.get(fam, benchmark=scen, method=v, rho_buf=0.0, **cfg)
        if f is None:
            raise KeyError((fam, scen, v, cfg))
        rows.append(dict(found=f, is_variant=True, group=(scen,) + tuple(sorted(cfg.items(), key=str)),
                         name=VNAME[v], method=v, levels=levels, level_disp=level_disp, record=dict(record, **cfg)))
    return rows


def anchor_tuple(D, anchor, fam="stress_drift_200"):
    ts = {(t, K) for f, t, K in D.anchor[anchor] if f == fam and K != "infinity"}
    assert len(ts) == 1, (anchor, ts)
    (a0, a, d), K = ts.pop()
    return dict(alpha0=a0, alpha=a, delta=d, K=K)


def tuple_rot(t):
    """Rotated, self-labelled setting tuple for the anchor tables (saves about 30pt of width)."""
    return (r"\rotatebox[origin=c]{90}{$(\alpha_0,\alpha,\delta,K)=(" +
            ",".join([num(t["alpha0"]), num(t["alpha"]), num(t["delta"]), kfmt(t["K"])]) + ")$}")


def tuple_disp(t, with_K=True):
    lines = [rf"$\alpha_0={num(t['alpha0'])}$", rf"$\alpha={num(t['alpha'])}$", rf"$\delta={num(t['delta'])}$"]
    if with_K:
        lines.append(rf"$K={kfmt(t['K'])}$")
    return stack(lines)


# ----------------------------------------------------------------------------------------------- fragments
CAP_MSD = r"Mean$\,\pm\,$SD over $N_{\rm MC}$ trajectories; bold: best V-variant mean per setting."


def caption(text):
    return rf"\vncaption{{{text}}}"


def page(parts):
    return "\n".join(parts) + "\n\\clearpage\n"


def fragment_head(lay):
    return f"% widths and column separation for layout {lay.gid}\n{lay.W.tex()}\n{lay.sep_tex()}\n"


def pack(groups, budget):
    """Greedy packing of (key, nlines, rows) groups into pages; a page holds groups with one baseline key."""
    pages, cur, used, ck = [], [], 0, None
    for bkey, n, rows in groups:
        if cur and (used + n > budget or bkey != ck):
            pages.append((ck, cur))
            cur, used = [], 0
        cur.extend(rows)
        used += n
        ck = bkey
    if cur:
        pages.append((ck, cur))
    return pages


def sweep_fragment(D, scen, alpha0, alphas, deltas, Ks, fname, gid, budget, sample_note):
    lay = Layout(gid, [("a", r"$\alpha$"), ("d", r"$\delta$"), ("k", r"$K$")], PAIRS, 2,
                 ["1.6pt", "1.0pt", "0.6pt"], method_header="")
    groups = []
    for a in alphas:
        for d in deltas:
            rows = []
            for K in Ks:
                lv = (num(a), num(d), kfmt(K))
                disp = (rf"\rotatebox[origin=c]{{90}}{{$\alpha={lv[0]}$}}", rf"\rotatebox[origin=c]{{90}}{{$\delta={lv[1]}$}}",
                        f"${lv[2]}$")
                rows += variant_rows(D, "main_grid_200", scen, lv, disp,
                                     dict(scenario=scen, table="sweep"), alpha0=alpha0, alpha=a, delta=d, K=K,
                                     L_mu=0.0005, bound_scale=1.0, budget_normalizer=2.11)
            groups.append((0, len(rows) * 2, rows))
    allrows = [r for g in groups for r in g[2]]
    fill_stat_rows(allrows)
    brows = baseline_rows(D, scen, alpha0, 0.0005)
    fill_stat_rows(brows)
    valid_alpha = [x / 100 for x in range(5, int(round(alpha0 * 100)), 5)]
    keyrow = (rf"\fbox{{$\alpha\in\{{{','.join(num(x) for x in valid_alpha)}\}};\ "
              rf"\delta\in\{{{','.join(num(x) for x in DELTAS)}\}};\ K\in\{{10,50,100,200,\infty\}}$}}")
    sets = settings_line(scen, alpha0=alpha0, L_mu="0.0005", b="2.11", phi="1")
    bsets = settings_line(scen, tag="Baselines", alpha0=alpha0, L_mu="0.0005")
    body = []
    pages = pack(groups, budget)
    for p, (_, rows) in enumerate(pages):
        btex, dropped = baseline_table(lay, brows, bsets, [], fname, f"{gid}-base-p{p + 1}")
        mtex = main_table(lay, rows, sets, keyrow, fname, f"{gid}-p{p + 1}")
        cap = CAP_MSD
        if any(r["found"] is None for r in brows):
            cap += r" --: baseline not simulated at this $\alpha_0$."
        if sample_note:
            cap += " " + sample_note
        body.append(page([btex, r"\vspace{5pt}", mtex, caption(cap)]))
    return fragment_head(lay) + "\n".join(body), lay


def anchor_param_fragment(D, scen, kind, fname, gid, budget, L_mu=None):
    """kind: 'lmu' (L_mu sweep), 'phi' (stress factor at L_mu), 'b' (budget normalizer)."""
    if kind == "lmu":
        pcol, fam, vals = ("p", r"$L_\mu$"), "stress_drift_200", [0.0005, 0.005]
    elif kind == "phi":
        pcol, fam, vals = ("p", r"$\phi$"), "stress_drift_200", [0.0, 0.25, 0.5, 0.75, 1.0]
    else:
        pcol, fam, vals = ("p", r"$b$"), "budget_normalizer_200", [2.11, 2.5, 3.0, 4.0]
    lay = Layout(gid, [("t", ""), pcol], PAIRS, 2, ["1.6pt", "0.6pt"], method_header="")
    groups = []
    for an in ANCHORS:
        t = anchor_tuple(D, an, fam)
        rows = []
        for val in vals:
            cfg = dict(t, L_mu=0.0005, bound_scale=1.0, budget_normalizer=2.11)
            if kind == "lmu":
                cfg["L_mu"] = val
                disp = "0.0005" if val == 0.0005 else "0.005"
            elif kind == "phi":
                cfg["L_mu"], cfg["bound_scale"] = L_mu, val
                disp = f"{val:g}"
            else:
                cfg["budget_normalizer"] = val
                disp = f"{val:g}"
            rows += variant_rows(D, fam, scen, (an, disp), (tuple_rot(t), f"${disp}$"),
                                 dict(scenario=scen, table=kind, anchor=an), **cfg)
        groups.append(((t["alpha0"],), len(rows) * 2, rows))
    fill_stat_rows([r for g in groups for r in g[2]])
    body = []
    for p, (bkey, rows) in enumerate(pack(groups, budget)):
        a0 = bkey[0]
        if kind == "lmu":
            brows = []
            for lm, disp in ((0.0005, "0.0005"), (0.005, "0.005")):
                brows += baseline_rows(D, scen, a0, lm,
                                       level=((disp,), (rf"\rotatebox[origin=c]{{90}}{{$L_\mu={disp}$}}",)))
            blev = [pcol]
            bsets = settings_line(scen, tag="Baselines", alpha0=a0)
            sets = settings_line(scen, b="2.11", phi="1")
        elif kind == "phi":
            brows = baseline_rows(D, scen, a0, L_mu)
            blev = []
            lmd = "0.0005" if L_mu == 0.0005 else "0.005"
            bsets = settings_line(scen, tag="Baselines", alpha0=a0, L_mu=lmd)
            sets = settings_line(scen, L_mu=lmd, b="2.11")
        else:
            brows = baseline_rows(D, scen, a0, 0.0005)
            blev = []
            bsets = settings_line(scen, tag="Baselines", alpha0=a0, L_mu="0.0005")
            sets = settings_line(scen, L_mu="0.0005", phi="1")
        fill_stat_rows(brows)
        btex, _ = baseline_table(lay, brows, bsets, blev, fname, f"{gid}-base-p{p + 1}")
        mtex = main_table(lay, rows, sets, None, fname, f"{gid}-p{p + 1}")
        cap = CAP_MSD
        if kind == "phi":
            cap += r" $\phi<1$ scales the derived bound below its proven value."
        body.append(page([btex, r"\vspace{5pt}", mtex, caption(cap)]))
    return fragment_head(lay) + "\n".join(body), lay


def checks_fragment(D, scen, L_mu, fname, gid, budget):
    evs = ["qcert", "qdep", "rCBF", "runsafe", "qA2", "qpw"]
    lay = Layout(gid, [("t", ""), ("k", r"$K$")], [[e] for e in evs], 1, ["1.6pt", "0.6pt"])
    groups = []
    for an in ANCHORS:
        t = anchor_tuple(D, an, "detailed_200")
        rows = []
        for K in (t["K"], "infinity"):
            cfg = dict(alpha0=t["alpha0"], alpha=t["alpha"], delta=t["delta"], K=K, L_mu=L_mu, bound_scale=1.0,
                       budget_normalizer=2.11)
            rows += variant_rows(D, "detailed_200", scen, (an, kfmt(K)), (tuple_disp(t, False), f"${kfmt(K)}$"),
                                 dict(scenario=scen, table="checks", anchor=an), **cfg)
        groups.append((0, len(rows), rows))
    fill_event_rows([r for g in groups for r in g[2]], evs)
    lmd = "0.0005" if L_mu == 0.0005 else "0.005"
    body = []
    for p, (_, rows) in enumerate(pack(groups, budget)):
        brows = []
        for a0 in (0.15, 0.30):
            brows += baseline_rows(D, scen, a0, L_mu, level=((num(a0),), (f"${num(a0)}$",)))
        fill_event_rows(brows, evs)
        btex, _ = baseline_table(lay, brows, settings_line(scen, tag="Baselines", L_mu=lmd),
                                 [("t", r"$\alpha_0$")], fname, f"{gid}-base-p{p + 1}")
        mtex = main_table(lay, rows, settings_line(scen, L_mu=lmd, b="2.11", phi="1"), None, fname,
                          f"{gid}-p{p + 1}")
        cap = (r"Trajectory counts/$N$ (one-sided 95\% upper bound, \%): $\widehat q_{\rm cert}$, $\widehat q_{\rm dep}$ "
               r"confirmed--possible events; $\widehat q_{\rm CBF}$, $\widehat q_{\rm unsafe}$ realized transitions; "
               r"$\widehat q_{\rm A2}$, $\widehat q_{\rm pw}$ any $e_{{\rm A2},t}$, $e_{{\rm pw},t}$ excess among "
               r"assessed trajectories.")
        body.append(page([btex, r"\vspace{5pt}", mtex, caption(cap)]))
    return fragment_head(lay) + "\n".join(body), lay


def topfive_fragment(D, cand_path, fname, gid):
    cr = json.loads(Path(cand_path).read_text())
    cands = cr["candidates"]
    cols = [["T_mis"], ["m_c"], ["q_CBF"], ["q_unsafe"], ["qcert"], ["qdep"]]
    lay = Layout(gid, [("a0", r"$\alpha_0$"), ("a", r"$\alpha$"), ("d", r"$\delta$"), ("k", r"$K$")], cols, 1,
                 ["1.0pt", "1.0pt", "1.0pt", "1.0pt"])
    body = []
    parts = []
    nomination = []
    for scen in ("corridor", "maze"):
        t0 = anchor_tuple(D, "current_default", "detailed_200")
        assert (t0["alpha0"], t0["alpha"], t0["delta"], t0["K"]) == (0.15, 0.1, 0.05, 50)
        lv = tuple(num(t0[k]) for k in ("alpha0", "alpha", "delta")) + (kfmt(t0["K"]),)
        rows = variant_rows(D, "detailed_200", scen, tuple(f"0:{x}" for x in lv), tuple(f"${x}$" for x in lv),
                            dict(scenario=scen, table="topfive", anchor="current_default"),
                            L_mu=0.0005, bound_scale=1.0, budget_normalizer=2.11, **t0)
        for i, c in enumerate(cands, 1):
            t = anchor_tuple(D, f"candidate_{i}", "detailed_200")
            assert (t["alpha0"], t["alpha"], t["delta"], key(t["K"])) == (c["alpha0"], c["alpha"], c["delta"], key(c["K"])), (i, t, c)
            f = D.get("detailed_200", benchmark=scen, method=c["variant"], rho_buf=0.0, L_mu=0.0005, bound_scale=1.0,
                      budget_normalizer=2.11, **t)
            lv = tuple(num(t[k]) for k in ("alpha0", "alpha", "delta")) + (kfmt(t["K"]),)
            rows.append(dict(found=f, is_variant=False, group=None, name=VNAME[c["variant"]], method=c["variant"],
                             levels=tuple(f"{i}:{x}" for x in lv), level_disp=tuple(f"${x}$" for x in lv),
                             record=dict(scenario=scen, table="topfive", anchor=f"candidate_{i}", **t)))
            if scen == "corridor":
                nomination.append(dict(candidate=i, role=c["role"], variant=c["variant"], **t))
        fill_stat_rows(rows)
        fill_event_rows(rows, ["qcert", "qdep"])
        brows = []
        for a0 in (0.15, 0.30):
            brows += baseline_rows(D, scen, a0, 0.0005, level=((num(a0),), (f"${num(a0)}$",)))
        fill_stat_rows(brows)
        fill_event_rows(brows, ["qcert", "qdep"])
        btex, _ = baseline_table(lay, brows, settings_line(scen, tag="Baselines", L_mu="0.0005"),
                                 [("a0", r"$\alpha_0$")], fname, f"{gid}-{scen}-base")
        mtex = main_table(lay, rows, settings_line(scen, L_mu="0.0005", b="2.11", phi="1"), None, fname,
                          f"{gid}-{scen}")
        parts += [btex, r"\vspace{3pt}", mtex, r"\vspace{8pt}"]
    cap = (r"First block: default setting; next five rows: candidates 1--5 of the existing N30 nomination, each with "
           r"its nominating variant (not re-selected at $N_{\rm MC}=200$). " + CAP_MSD +
           r" Events: counts/$N_{\rm MC}$ (one-sided 95\% upper bound, \%).")
    body.append(page(parts[:-1] + [caption(cap)]))
    return fragment_head(lay) + "\n".join(body), lay, nomination


# ----------------------------------------------------------------------------------------------- preamble/sample
PREAMBLE = r"""% v10 table preamble (generated by build_tables.py). Requires array (w columns), multirow, makecell, xcolor.
\usepackage{array,multirow,makecell,xcolor,graphicx}
\definecolor{vu}{HTML}{56B4E9}\definecolor{vd}{HTML}{0072B2}\definecolor{vl}{HTML}{000000}\definecolor{vc}{HTML}{E69F00}
\newcommand{\na}{\textup{N-A}}
\newcommand{\vpm}{{\pm}}
\newlength{\vnbox}\setlength{\vnbox}{0.8pt}
\setlength{\arrayrulewidth}{0.4pt}
\renewcommand{\arraystretch}{0.93}
\makeatletter
\newlength{\vn@tmp}
\newcommand{\vnw}[1]{\csname vnw@#1\endcsname}
\newcommand{\vnnew}[1]{\expandafter\newlength\csname vnw@#1\endcsname\global\csname vnw@#1\endcsname=0pt}
\newcommand{\vnfit}[2]{\settowidth{\vn@tmp}{#2}\ifdim\vn@tmp>\csname vnw@#1\endcsname\global\csname vnw@#1\endcsname=\vn@tmp\fi}
\newcommand{\vnsep}[1]{\csname vnsep@#1\endcsname}
\newcommand{\vnnewsep}[2]{\expandafter\newlength\csname vnsep@#1\endcsname
  \global\csname vnsep@#1\endcsname=#2\relax
  \typeout{VN-SEP #1 = \the\csname vnsep@#1\endcsname}%
  \ifdim\csname vnsep@#1\endcsname<1pt\PackageWarning{v10tables}{column separation of #1 below 1pt}\fi}
% \vnsd: zero-width strut, depth = array strut depth + \vnsubpad; used only on rows whose subscripts sit on a rule
\newlength{\vnsubpad}\setlength{\vnsubpad}{1.5pt}
\newcommand{\vnsd}{\rule[-\dimexpr\dp\@arstrutbox+\vnsubpad\relax]{0pt}{\dimexpr\dp\@arstrutbox+\vnsubpad\relax}}
\makeatother
\newcommand{\vncheckleft}[3]{\ifdim#1<#2\PackageWarning{v10tables}{baseline Method column too narrow in #3}\fi}
\newcommand{\vncaption}[1]{\par\vspace{3pt}\noindent#1\par}
"""

SAMPLE_HEAD = r"""\documentclass[10pt,letterpaper]{article}
% \textwidth of the paper (ieeeconf 10pt letter): 505.89pt; \textheight 696pt (measured from ieeeconf.cls).
\usepackage[paperwidth=8.5in,paperheight=11in,textwidth=505.89pt,textheight=696pt,centering,headheight=12pt,headsep=8pt,footskip=18pt]{geometry}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb}
\input{v10_tables_preamble}
\setlength{\parindent}{0pt}
\pagestyle{plain}
\begin{document}
\typeout{VN-TEXTWIDTH=\the\textwidth}\typeout{VN-TEXTHEIGHT=\the\textheight}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics", type=Path, default=DEF_METRICS)
    ap.add_argument("--candidates", type=Path, default=DEF_CANDIDATES)
    ap.add_argument("--out", type=Path, default=HERE)
    ap.add_argument("--sample-alpha0", type=float, default=0.15)
    ap.add_argument("--sample-alpha", default="0.05,0.10")
    ap.add_argument("--sample-delta", default="0.05,0.10")
    ap.add_argument("--sweep-full", default="", help="comma list of alpha0 values for the complete sweep tables")
    ap.add_argument("--page-lines", type=int, default=40, help="max value lines of a main table per page")
    ap.add_argument("--no-compile", action="store_true")
    a = ap.parse_args()
    out = a.out.resolve()
    frag = out / "fragments"
    frag.mkdir(parents=True, exist_ok=True)
    D = Data(a.metrics)

    # check: this Clopper-Pearson bound reproduces the stored upper95 fields
    chk = 0
    for (fam, *_), (vid, df, c) in list(D.idx.items())[:400]:
        e = c["realized_transition_events"]["q_CBF"]
        if e["assessed_n"]:
            assert abs(cp_upper(e["count"], e["assessed_n"]) - e["upper95"]) < 1e-9, (vid, e)
            chk += 1
    print(f"Clopper-Pearson check against stored upper95: {chk} cells OK")

    (out / "v10_tables_preamble.tex").write_text(PREAMBLE)
    files = []

    def emit(name, tex):
        (frag / name).write_text(tex)
        files.append(name)

    alphas = [float(x) for x in a.sample_alpha.split(",")]
    deltas = [float(x) for x in a.sample_delta.split(",")]
    note = rf"Sample rows: $\alpha\in\{{{','.join(num(x) for x in alphas)}\}}$, $\delta\in\{{{','.join(num(x) for x in deltas)}\}}$."
    tex, nom = None, None
    tex, _, nom = topfive_fragment(D, a.candidates, "topfive.tex", "tf")
    emit("topfive.tex", tex)
    for scen in ("corridor", "maze"):
        sc = scen[0]
        a0 = a.sample_alpha0
        tag = f"{int(round(a0 * 100)):03d}"
        tex, _ = sweep_fragment(D, scen, a0, alphas, deltas, KS, f"sweep_{scen}_a{tag}_sample.tex", f"sw{sc}",
                                a.page_lines, note)
        emit(f"sweep_{scen}_a{tag}_sample.tex", tex)
        for L_mu, lt in ((0.0005, "L0005"), (0.005, "L005")):
            tex, _ = checks_fragment(D, scen, L_mu, f"checks_{scen}_{lt}.tex", f"ck{sc}{lt[1:]}", 24)
            emit(f"checks_{scen}_{lt}.tex", tex)
        tex, _ = anchor_param_fragment(D, scen, "lmu", f"lmu_{scen}.tex", f"lm{sc}", 32)
        emit(f"lmu_{scen}.tex", tex)
        for L_mu, lt in ((0.0005, "L0005"), (0.005, "L005")):
            tex, _ = anchor_param_fragment(D, scen, "phi", f"phi_{scen}_{lt}.tex", f"ph{sc}{lt[1:]}", 40, L_mu)
            emit(f"phi_{scen}_{lt}.tex", tex)
        tex, _ = anchor_param_fragment(D, scen, "b", f"b_{scen}.tex", f"bn{sc}", 40)
        emit(f"b_{scen}.tex", tex)
    full_files = []
    for a0s in [x for x in a.sweep_full.split(",") if x.strip()]:
        a0 = float(a0s)
        tag = f"{int(round(a0 * 100)):03d}"
        al = [x / 100 for x in range(5, int(round(a0 * 100)), 5)]
        for scen in ("corridor", "maze"):
            name = f"sweep_{scen}_a{tag}_full.tex"
            tex, _ = sweep_fragment(D, scen, a0, al, DELTAS, KS, name, f"sf{scen[0]}{tag}", a.page_lines, "")
            (frag / name).write_text(tex)
            full_files.append(name)

    sample_files = files  # full sweeps are not part of the sample document
    entries = REC.entries
    values = dict(
        generator="build_tables.py", metrics=recorded_path(a.metrics), metrics_sha256=D.sha256,
        candidates=recorded_path(a.candidates), n30_nomination=nom,
        rounding=dict(normalized=2, percent=2, seconds=2, metres=3, counts_mean=2, upper_bound_pct="ceil 1 decimal",
                      rule="half-up on the stored double (as papertables.fmt_half_up)"),
        fragments=sample_files + full_files, cells=entries)
    (out / "tables_values.json").write_text(json.dumps(values, indent=1, default=str))
    print(f"{len(entries)} printed cells recorded")

    order = ["topfive.tex"]
    for scen in ("corridor", "maze"):
        order += [f for f in files if f.startswith(("sweep_" + scen,))]
        order += [f for f in files if f.startswith("checks_" + scen)]
        order += [f for f in files if f.startswith("lmu_" + scen)]
        order += [f for f in files if f.startswith("phi_" + scen)]
        order += [f for f in files if f.startswith("b_" + scen)]
    doc = SAMPLE_HEAD + "".join(rf"\input{{fragments/{f}}}" + "\n" for f in order) + r"\end{document}" + "\n"
    (out / "v10_tables_sample.tex").write_text(doc)
    if a.no_compile:
        return
    for _ in range(1):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "v10_tables_sample.tex"],
                           cwd=out, capture_output=True, text=True)
    log = (out / "v10_tables_sample.log").read_text(errors="replace")
    if p.returncode:
        print(log[-3000:])
        sys.exit("pdflatex failed")
    tw = re.search(r"VN-TEXTWIDTH=([\d.]+pt)", log).group(1)
    over = len(re.findall(r"Overfull \\[hv]box", log))
    under = len(re.findall(r"Underfull \\[hv]box", log))
    seps = re.findall(r"VN-SEP (\S+) = ([\d.\-]+)pt", log)
    pages = re.search(r"Output written on .*?\((\d+) pages", log).group(1)
    blocks = sum((frag / f).read_text().count("\\clearpage") for f in order)
    print(f"page blocks {blocks}; ", end="")
    print(f"textwidth {tw}; pages {pages}; Overfull {over}; Underfull {under}; "
          f"min column separation {min(float(s) for _, s in seps):.2f}pt")
    warn = len(re.findall(r"Package v10tables Warning", log))
    print(f"v10tables warnings {warn}")
    if tw != "505.89pt" or over or under or int(pages) != blocks or warn:
        sys.exit("layout check failed")


if __name__ == "__main__":
    main()
