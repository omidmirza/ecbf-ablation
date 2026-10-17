#!/usr/bin/env python3
"""Independent check of the N200 table fragments: re-reads the exports and re-derives every printed string
recorded in tables_values.json.gz, without importing build_tables.py.

For every recorded cell it
  * finds the view by (source, view_id) in final/common, new_anchors/common or final/temporal, checks the
    recorded configuration, n = 200 and the seeds 970001-970200;
  * recomputes the printed string: mean +- SD (half-up after rounding the double to d+9 decimals; decimals per
    metric, Section VI precision m_c 3 / q_use 2), event "a (u)" or "a/m (u)" with the conservative count
    (possible events; realized counts; trajectories with any excess), u = one-sided 95% Clopper-Pearson bound in %
    rounded up, \\checkmark iff the row is a certified variant, the event is q_cert or q_dep, and u <= 100 delta
    (finite K) or 100 delta/(2 b ln^2 2) (K = infinity); temporal 10^3 V / 100 J (3 decimals), assessed sizes;
    N-A when not applicable; "--" only when no matched baseline exists;
  * recomputes the bold marks (best mean among the four variants of one setting, ties all bold, arrow direction,
    never for m_cert, tau_use, p_max or conditional metrics);
  * checks that the strings appear in the fragment in the recorded order and that the fragment contains no other
    mean +- SD or event strings than the recorded ones;
  * recomputes the configuration counts of the design tables;
  * compares Section VI Table I (sec6_combined: default setting and baselines at alpha0 = 0.15) with the previous Table I (appendices/sim_tables/joint_main.tex):
    the same means and SDs, and counts equal to the upper end of a previous range.
Structure checks:
  * checks, detailed-event, delta, K, b and phi tables and the Section VI short table show only the paired
    variant of each setting (Default V-C, S1 V-C, S2 V-C, S3 V-L, S4 V-Delta, S5 V-C); paper parameter tables show
    Default and S3 only, report parameter tables all six settings; the full top-five table all four variants at S1-S5;
  * each recorded row label (alpha0, alpha, delta, K) equals the configuration of its cell and is printed;
  * key boxes only in table_key, sec6_short, sec6_combined and sec6_maze_arxiv (one each);
  * sec6_maze_arxiv (arXiv Table I, Maze only) holds S1-S5 with the four certified variants and the five
    baselines at alpha0 = 0.30; its chosen-variant rows and its baselines print the Maze cells of sec6_short in the
    shared columns (bold aside); its main table names each setting in a block row (as topfive_full);
  * fragments of earlier layouts are not generated; colour vl = 009E73 in appendices/n200_assets.tex.
arXiv checks and stress fragments: checks_{scen}_arxiv (S1-S5 x L_mu x {K, inf}, 20 rows) and
stress_{scen}_arxiv (S1-S5 x phi x L_mu, 50 rows, dagger iff phi < 1): every column named in the "% OMITTED" line is
re-derived as 0 or N-A in every row from the exports, no printed column is all-zero; every baseline box is centred
at \\linewidth, its Method column absorbs exactly the value columns it drops, and a first-column Method header cell
restates the left box rule (editions and catalog).
Earlier layout of the arXiv Tables VIII-X (arxiv_layout_earlier(), build_tables.Layout align/max_sep/inner_rules/stretch): in
the four arXiv fragments the value columns are right-aligned (w{r}) under centred headers (\\multicolumn{1}{c...}),
\\vnnock is defined before use and ends exactly the cells without a check mark in the columns that have check marks,
the column separation is capped at 6 pt and the settings row has the table's natural width (not \\linewidth), no rule
separates the rows of one innermost block, the stress tables label their first column "Method (alpha0, alpha, delta,
K)", and the checks tables set \\arraystretch 1.1.
Catalog: the same cell checks for the catalog fragments (catalog_values.json.gz, generated_n200/catalog/),
plus catalog_checks(): four variants per setting in every all-variant table, full main-grid coverage of the sweep
and compact tables, and the compact best variant and its count (k) re-derived from the exports.
Setting tuples: tuple_checks() - no setting tuple split over two lines or stacked, block rows are one-line l
columns; block rows name the setting tuple without the swept delta or K. The rows match it on every other field.
Uniform layout (Table X continued): uniform_layout_checks() on every edition fragment except the key, the design tables and the
catalog-only N30 record - no rotated or stacked labels, every column centred (w{c}), no capped column separation,
equal value-column widths when they fit (the equalization block before every separation), full-width settings and
block rows with the fixed padding, at least one block row per main table (Section VI tables: the setting in the
settings row or one row per setting instead), every baseline box centred at \\linewidth with a complete Method header
box and either the grid of its table (aligned) or its own balanced grid (independent), \\arraystretch within
[0.93, 1.4]; arxiv_layout_uniform() replaces the earlier layout checks of the arXiv checks and stress fragments (block
rows per setting or per (setting, L_mu), rows L_mu or phi, the setting's K and K = infinity as column groups,
\\vnnock iff a cell lacks the check mark of its column, the stress tables on two pages).
Corridor table and arXiv appendix: main_corridor_arxiv_v8 (arXiv Table IV) is checked cell by cell as
main_maze_arxiv_v8 (MAIN_V8_FRAGS) and for its structure (25 rows, 150 cells, layout, bold marks re-derived, the
Unfiltered T_mis cell at the horizon); the report fragments stress_{scen}_cases print the cells of
stress_{scen}_arxiv and pass the same stress checks; check_assumption_totals() re-derives the assumption-check totals and the main-grid event
counts of appendices/generated_n200/paper/assumption_checks_arxiv.tex from the exports; check_arxiv_parameter_table()
checks the arXiv branch of appendices/f_tab_setup.tex (Table III) against the exports, the extended rows and, when
ECBF_PREVIOUS_SETUP names a file, the earlier arXiv rows (otherwise that comparison is skipped with a NOTE).
Inputs: the exports named in the values files are resolved against ECBF_N200_NUMBERS (default
<repository>/N200_numbers; not included); the optional earlier parameter table against ECBF_PREVIOUS_SETUP.
Usage: python3 check_tables.py
"""
from __future__ import annotations

import json
import math
import os
import re
from collections import defaultdict
from decimal import ROUND_CEILING, ROUND_HALF_UP, Decimal
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1]
GEN = SRC / "appendices" / "generated_n200"
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_tables/<this file>)
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports (not included)


def resolve_export(p):
    """A path recorded in a values file: relative "N200_numbers/<rest>" resolves against NUM; absolute as is."""
    p = Path(p)
    return p if p.is_absolute() else NUM.joinpath(*p.parts[1:])


def gen_text(path):
    """Plain (uncoloured) variant of a generated fragment (colour switch, table_colour.py): the fragment is checked
    exactly as the editions typeset it by default; identity for a fragment without the switch."""
    import sys
    sys.path.insert(0, str(HERE))
    import table_colour
    return table_colour.plain_variant(Path(path).read_text())
PCT = {"q_CBF", "q_Y", "q_use", "q_unsafe", "q_nocert", "q_notrans", "p_max"}
DEC = {"m_c": 2, "m_Y": 2, "m_int": 2, "n_certificates": 2, "q_CBF": 2, "q_Y": 2, "q_use": 1, "q_unsafe": 2,
       "q_nocert": 1, "q_notrans": 1, "h_min": 2, "clearance_min": 2, "tau_use": 2, "p_max": 1, "T_mis": 2}
SEC6_DEC = {"m_c": 3, "q_use": 2}
PK_DEC = 3   # Table I column pk = statistics.p_max as a fraction, 3 decimals, plain
UP = {"q_use", "h_min", "clearance_min"}
NOBOLD = {"n_certificates", "tau_use", "p_max"}
VARIANTS = {"V-U", "V-D", "V-Cd", "V-C"}
EVSRC = {"qcert": ("trajectory_events", "qcert"), "qdep": ("trajectory_events", "qdep"),
         "rCBF": ("realized_transition_events", "q_CBF"), "runsafe": ("realized_transition_events", "q_unsafe"),
         "qA2": ("risk_increase_checks", "direct"), "qpw": ("risk_increase_checks", "pairwise")}
TMP = {"Vc": ("normalized_margin", "temporal_variance", 1000), "Jc": ("normalized_margin", "step_rms", 100),
       "Vp": ("estimated_applied_risk", "temporal_variance", 1000), "Jp": ("estimated_applied_risk", "step_rms", 100)}
SEEDS = list(range(970001, 970201))


def rnd(v, d):
    q = Decimal(("{:.%df}" % (d + 9)).format(v)).quantize(Decimal(10) ** -d, rounding=ROUND_HALF_UP)
    t = format(q, f".{d}f")
    return t[1:] if t.startswith("-") and q.is_zero() else t


def ceil_pct(u):
    return format((Decimal(repr(u)) * 100).quantize(Decimal("0.1"), rounding=ROUND_CEILING), ".1f")


def cp_upper(x, n):
    """Exact one-sided 95% Clopper-Pearson bound via exact binomial tails and bisection."""
    if x >= n:
        return 1.0
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(60):
        p = (lo + hi) / 2
        tail = sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(x + 1))
        if tail > Fraction(1, 20):
            lo = p
        else:
            hi = p
        lo, hi = Fraction(float(lo)), Fraction(float(hi))
    return float(hi)


STAR = HERE / "starred_setting.json"


def check_starred(e, fail):
    """A cell of row S* (sec6_combined) equals starred_setting.json, and the JSON equals the cell of
    the source table at the recorded line and column (one count = upper end of a printed range; bound as printed;
    \\checkmark iff bound <= 100 delta)."""
    S = json.loads(STAR.read_text())
    c = next((x for x in S["cells"] if x["spec"] == e["spec"]), None)
    if c is None or e["source_line"] != c["line"] or e["row"].get("label") != "(0.25, 0.225, 0.15, 50)":
        fail("starred cell not in JSON", e)
        return
    line = (SRC / S["source_file"]).read_text().splitlines()[c["line"] - 1]
    parts = [p.strip() for p in line.rstrip().rstrip("\\").split(" & ")]
    if parts[0] != S["method"] or parts[c["column"]] != c["source_text"]:
        fail(f"source cell {S['source_file']}:{c['line']} col {c['column']} = {parts[c['column']]!r}", e)
        return
    src = re.sub(r"\\reviewedit\{\\mathbf\{([^}]*)\}\}", r"\1", c["source_text"])
    if c["kind"] == "meanSD":
        m = re.fullmatch(r"\$(.+)\\pm (\S+)\$", src)
        ok = m and (m.group(1), m.group(2)) == (c["mean"], c["sd"])
        want = f"${c['mean']}\\vpm {c['sd']}$"
    else:
        m = re.fullmatch(r"(\d+)(?:--(\d+))?/200 \((\d+\.\d)\\%\)", src)
        # Consistent with the submitted values: a cell with count_rule
        # "confirmed" prints the lower end of the source range (confirmed events) with the source bound, which is
        # the bound of the upper end; every other cell prints the upper end with its own bound.
        hi = int(m.group(2) or m.group(1)) if m else None
        cnt = int(m.group(1)) if (m and c.get("count_rule") == "confirmed") else hi
        ok = m and cnt == c["count"] and m.group(3) == c["upper_pct"]
        ok = ok and ceil_pct(cp_upper(hi, 200)) == c["upper_pct"]
        want = f"{c['count']} ({c['upper_pct']})"
        if Decimal(c["upper_pct"]) <= Decimal(repr(100 * S["configuration"]["delta"])):
            want += r"\,\checkmark"
    if not ok:
        fail(f"JSON value differs from the source cell {c['source_text']!r}", e)
    if e["printed"] != want:
        fail(f"starred string, want {want}", e)


# Fragment star_maze (Appendix G-I.2, full edition): S* with the four certified variants in the full-metric layout
# of Tables XIII/XIV; spec -> (source block, source column, source unit -> MET unit).
STAR_VARIANT_NAME = {"V-U": "V-U", "V-D": r"V-$\Delta$", "V-Cd": "V-L", "V-C": "V-C"}
STAR_VARIANT_ROWS = {"V-U": (13, 27, 41, 55), "V-D": (14, 28, 42, 56), "V-Cd": (15, 29, 43, 57), "V-C": (16, 30, 44, 58)}
STAR_SPEC_BLOCK = {"m_c": (0, 2, 1), "q_use": (0, 3, 1), "p_max": (0, 4, 100), "tau_use": (0, 5, 1),
                   "T_mis": (0, 1, 1), "q_CBF": (1, 2, 1), "q_unsafe": (1, 3, 1), "h_min": (1, 5, Fraction(1, 1000)),
                   "qdep": (2, 1, None), "qcert": (2, 3, None), "n_certificates": (2, 5, 1),
                   "m_Y": (3, 1, 1), "q_Y": (3, 2, 1), "m_int": (3, 3, 1), "q_nocert": (3, 4, 1), "q_notrans": (3, 5, 1)}
STAR_NONE = {"clearance_min"}   # not in the source table: printed "--"


def star_expected(spec, delta):
    """Independent re-derivation of the four printed strings of one star_maze column from the source table.
    Means and SDs: converted exactly (fractions) to the unit of the column; printed with DEC[spec] decimals (half up)
    when the interval [x - u/2, x + u/2) of every source number of the column lies in one rounding cell, else with the
    source digits; bold: lowest/highest printed mean of the four variants (UP), never for NOBOLD; ties all bold."""
    lines = (SRC / json.loads(STAR.read_text())["source_file"]).read_text().splitlines()
    blk, col, f = STAR_SPEC_BLOCK[spec]
    cells = {}
    for v, rows in STAR_VARIANT_ROWS.items():
        parts = [p.strip() for p in lines[rows[blk] - 1].rstrip().rstrip("\\").split(" & ")]
        if parts[0] != STAR_VARIANT_NAME[v]:
            return None
        cells[v] = parts[col]
    out = {}
    if f is None:
        for v, cell in cells.items():
            txt = re.sub(r"^\\reviewedit\{\\textbf\{(.*?)\}\}", r"\1", cell)
            m = re.fullmatch(r"(\d+)(?:--(\d+))?/200 \((\d+\.\d)\\%\)", txt)
            if not m:
                return None
            hi = int(m.group(2) or m.group(1))
            if ceil_pct(cp_upper(hi, 200)) != m.group(3):
                return None
            cnt = int(m.group(1)) if spec == "qcert" else hi
            out[v] = f"{cnt} ({m.group(3)})" + (r"\,\checkmark" if Decimal(m.group(3)) <= Decimal(repr(100 * delta)) else "")
        return out
    d = DEC[spec]

    def parse(t):
        lt = t.startswith("{<}")
        body = t[3:] if lt else t
        nd = len(body.split(".")[1]) if "." in body else 0
        x, u = Fraction(body) * f, Fraction(1, 10 ** nd) * f
        sd_ = 0
        while u * 10 ** sd_ != int(u * 10 ** sd_):
            sd_ += 1
        src = format(Decimal(x.numerator) / Decimal(x.denominator), f".{sd_}f")
        if lt:
            return ("{<}" + src, ("{<}" + src) if sd_ <= d else None, x - u / 10)
        if sd_ <= d:
            return src, src, x
        lo_cell = math.floor((x - u / 2) * 10 ** d + Fraction(1, 2))       # half up of the lower end
        hi_val = (x + u / 2) * 10 ** d + Fraction(1, 2)                      # the upper end is excluded
        hi_cell = math.ceil(hi_val) - 1
        if lo_cell != hi_cell:
            return src, None, x
        r = format(Decimal(lo_cell).scaleb(-d), f".{d}f")
        return src, r, Fraction(lo_cell, 10 ** d)

    nums = {}
    for v, cell in cells.items():
        if cell == "---":
            out[v] = r"\na"
            continue
        plain = re.sub(r"\\reviewedit\{\\mathbf\{([^}]*)\}\}", r"\1", cell)
        m = re.fullmatch(r"\$(.+?)\\pm\s*(\S+)\$", plain)
        if not m:
            return None
        nums[v] = (parse(m.group(1)), parse(m.group(2)))
    rounded = all(p[1] is not None for a in nums.values() for p in a)
    pr = {v: ((a[0][1], a[1][1]) if rounded else (a[0][0], a[1][0]), a[0][2] if not rounded or a[0][1].startswith("{<}")
              else Fraction(Decimal(a[0][1]))) for v, a in nums.items()}
    best = None
    if len(pr) == 4 and spec not in NOBOLD:
        ks = [k for _, k in pr.values()]
        best = max(ks) if spec in UP else min(ks)
    for v, ((mean, sd), k) in pr.items():
        mm = rf"\mathbf{{{mean}}}" if best is not None and k == best else mean
        out[v] = f"${mm}\\vpm {sd}$"
    return out


def check_starred_variant(e, fail):
    """A cell of star_maze equals the re-derivation from the source table (star_expected); qcert prints the lower
    end of the source range (confirmed events), qdep the upper end, each with the source bound
    (= Clopper-Pearson bound of the upper end); \\checkmark iff bound <= 100 delta."""
    S = json.loads(STAR.read_text())
    if e["fragment"] != "star_maze" or e["source_file"] != S["source_file"]:
        fail("star_maze cell outside its fragment", e)
        return
    if e["spec"] in STAR_NONE:
        if e["kind"] != "starred_missing" or e["printed"] != "--":
            fail("star_maze: a column without source values must print --", e)
        return
    if e["spec"] not in STAR_SPEC_BLOCK:
        fail("star_maze spec without a source column", e)
        return
    blk, col, _ = STAR_SPEC_BLOCK[e["spec"]]
    if e["source_line"] != STAR_VARIANT_ROWS[e["method"]][blk] or e["source_column"] != col:
        fail(f"star_maze source line {e['source_line']} / column {e['source_column']}", e)
        return
    want = star_expected(e["spec"], S["configuration"]["delta"])
    if want is None:
        fail("star_maze source column unreadable", e)
    elif e["printed"] != want[e["method"]]:
        fail(f"star_maze string, want {want[e['method']]}", e)


VIEWS = {}


def load_views(V):
    """Views of the exports named in a values file (loaded once; the SHA-256 must equal the build's)."""
    import hashlib
    for src, path in V["inputs"].items():
        if any(k[0] == src for k in VIEWS):
            continue
        raw = resolve_export(path).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == V["sha256"][src], f"{path} changed since the build"
        M = json.loads(raw)
        del raw
        for v in M["views"]:
            VIEWS[(src, v["view_id"])] = v
        del M
    return VIEWS


# arXiv Table I (main_maze_arxiv): mean +- SD with two decimals,
# no bold; q_dep and q_cert as the percentage of the 200 trajectories with a possible event (one decimal) and the
# one-sided 95% bound (%) in parentheses, no check marks; '---' where not applicable (q_cert of the baselines).
# Every certified bound must meet its level (the caption states it).
MAIN_MAZE = "main_maze_arxiv"


def check_main_maze(e, c, fail):
    grp, _, sub = e["source"].partition(".")
    cert = e["method"] in VARIANTS and not e["row"].get("baseline")
    if grp == "statistics":
        s = c["statistics"][sub]
        if not s["applicable"] or s["observed_n"] == 0 or s["mean"] is None:
            want = "---"
        else:
            sc = 100 if sub in PCT else 1
            want = f"${rnd(s['mean'] * sc, 2)}\\vpm {rnd(s['sd'] * sc, 2)}$"
    elif grp == "trajectory_events":
        x = c[grp][sub]
        n, a = x["applicable_n"], x["possible"]
        if not n or (sub == "qcert" and not cert):
            want = "---"
        else:
            if n != 200:
                fail("main_maze_arxiv: event not assessed on all 200 trajectories", e)
            u = cp_upper(a, 200)
            if abs(u - x["upper95_possible"]) > 1e-9:
                fail("main_maze_arxiv: stored bound differs from Clopper-Pearson", e)
            up = ceil_pct(u)
            want = f"{rnd(100 * a / 200, 1)} ({up})"
            if cert:
                cf = c["configuration"]
                lev = cf["delta"] if cf["K"] != "infinity" else cf["delta"] / (2 * cf["budget_normalizer"] *
                                                                              math.log(2) ** 2)
                if Decimal(up) > Decimal(repr(100 * lev)):
                    fail("main_maze_arxiv: certified bound above its level (caption claims all meet it)", e)
    else:
        fail("main_maze_arxiv: unexpected source", e)
        return
    if e["printed"] != want:
        fail(f"main_maze_arxiv string, want {want}", e)


# main_maze_arxiv_v2: the rows of main_maze_arxiv with six metrics in two lines per method
# (line 1 m_c, q_CBF, qdep; line 2 T_mis, pk, qcert). m_c and q_CBF two decimals, T_mis two decimals, pk =
# statistics.p_max available mean and SD (fraction, three decimals); events as in main_maze_arxiv.
MAIN_MAZE_V2 = "main_maze_arxiv_v2"
MAIN_MAZE_V2_LINES = {1: ("m_c", "q_CBF", "qdep"), 2: ("T_mis", "pk", "qcert")}


def check_main_maze_v2(e, c, fail):
    grp, _, sub = e["source"].partition(".")
    spec = e["spec"]
    if spec not in MAIN_MAZE_V2_LINES.get(e["row"].get("line"), ()):
        fail("main_maze_arxiv_v2: metric on the wrong line", e)
    if grp == "statistics":
        want_sub = "p_max" if spec == "pk" else spec
        if sub != want_sub:
            fail("main_maze_arxiv_v2: wrong statistic", e)
            return
        s = c["statistics"][sub]
        mk, sk = ("available_mean", "available_sd") if spec == "pk" else ("mean", "sd")
        if not s["applicable"] or s["observed_n"] == 0 or s[mk] is None:
            want = "---"
        else:
            d = {"m_c": 2, "q_CBF": 2, "T_mis": 2, "pk": 3}[spec]
            sc = 100 if spec == "q_CBF" else 1
            want = f"${rnd(s[mk] * sc, d)}\\vpm {rnd(s[sk] * sc, d)}$"
        if e["printed"] != want:
            fail(f"main_maze_arxiv_v2 string, want {want}", e)
    elif grp == "trajectory_events" and sub == spec:
        check_main_maze(e, c, fail)   # same rule and level check as main_maze_arxiv
    else:
        fail("main_maze_arxiv_v2: unexpected source", e)


# main_maze_arxiv_v3: one line per method, six metrics, baselines then S1-S5 x four variants with a rule between
# the baselines and the variants. m_c mean +- SD 2 dec, T_mis mean 1 dec, q_CBF mean +- SD 1 dec (%), pk =
# statistics.p_max available mean 2 dec; events as main_maze_arxiv with '100' for a full count or bound.
MAIN_MAZE_V3 = "main_maze_arxiv_v3"
MAIN_MAZE_V3_SPECS = ("m_c", "T_mis", "q_CBF", "pk", "qdep", "qcert")
# main_maze_arxiv_v4 = main_maze_arxiv_v3 with the units on the header line and T_mis as mean +- SD (1 decimal).
MAIN_MAZE_V4 = "main_maze_arxiv_v4"
# main_maze_arxiv_v5 = the cells of main_maze_arxiv_v4, \footnotesize, tabcolsep 2.5 pt,
# one header line of symbols only (units in the caption).
MAIN_MAZE_V5 = "main_maze_arxiv_v5"
# main_maze_arxiv_v6 = the cells of main_maze_arxiv_v5 in the report-table style (bold
# setting rows with a rule below, variants unindented, direction arrows, bold lowest printed mean of the four variants
# of one setting in m_c/c_max, T_mis, q_CBF^+). Each cell is re-derived without its bold mark, then the marks.
MAIN_MAZE_V6 = "main_maze_arxiv_v6"
MAIN_MAZE_V6_BOLD = ("m_c", "T_mis", "q_CBF")
# main_maze_arxiv_v7 = the rows of main_maze_arxiv_v6 with m_{p^max} as mean +- SD (statistics.p_max available mean
# and SD, fraction, 2 decimals), units and arrows on the one header line, every column split into r@{}l halves
# (mean | {}\mathbin{\vpm}SD, count | \ (bound), '---' over both halves), bold on the lowest printed mean of the
# four variants of one setting in m_c/c_max, T_mis, q_CBF^+, m_{p^max}. Every cell is re-derived from the exports.
MAIN_MAZE_V7 = "main_maze_arxiv_v7"
MAIN_MAZE_V7_BOLD = ("m_c", "T_mis", "q_CBF", "pk")
MAIN_MAZE_V7_HEAD = (r"Method & \multicolumn{2}{c}{$\widehat m_c/c_{\max}\!\downarrow$} & "
                     r"\multicolumn{2}{c}{$\widehat T_{\rm mis}\!\downarrow$\,(s)} & "
                     r"\multicolumn{2}{c}{$\widehat q^{+}_{\rm CBF}\!\downarrow$\,(\%)} & "
                     r"\multicolumn{2}{c}{$\widehat m_{p^{\max}}\!\downarrow$} & "
                     r"\multicolumn{2}{c}{$\widehat q_{\rm dep}\!\downarrow$\,(\%)} & "
                     r"\multicolumn{2}{c}{$\widehat q_{\rm cert}\!\downarrow$\,(\%)}\\\midrule")


def v7_join(t):
    """A v7 cell (two r@{}l halves) as the one-cell string of v6: '$m\\vpm s$', 'p (u)' or '---'."""
    if t == r"\multicolumn{2}{c}{---}":
        return "---"
    m = re.fullmatch(r"\$([^$]*)\$&\$\{\}\\mathbin\{\\vpm\}([0-9.]+)\$", t)
    if m:
        return f"${m.group(1)}\\vpm {m.group(2)}$"
    m = re.fullmatch(r"(\d+(?:\.\d)?)&\\ \((\d+(?:\.\d)?)\)", t)
    if m:
        return f"{m.group(1)} ({m.group(2)})"
    return None


def check_main_maze_v7(e, c, fail):
    """One cell of main_maze_arxiv_v7, re-derived from the exports (bold marks are checked in arxiv_checks)."""
    one = v7_join(e["printed"])
    if one is None:
        fail(f"{MAIN_MAZE_V7}: not a split cell", e)
        return
    one = re.sub(r"\\mathbf\{([^{}]*)\}", r"\1", one)
    grp, _, sub = e["source"].partition(".")
    spec = e["spec"]
    if grp == "statistics" and spec in ("m_c", "T_mis", "q_CBF", "pk"):
        if sub != ("p_max" if spec == "pk" else spec):
            fail(f"{MAIN_MAZE_V7}: wrong statistic", e)
            return
        s = c["statistics"][sub]
        mk, sk = ("available_mean", "available_sd") if spec == "pk" else ("mean", "sd")
        if not s["applicable"] or s["observed_n"] == 0 or s[mk] is None:
            want = "---"
        else:
            d = {"m_c": 2, "T_mis": 1, "q_CBF": 1, "pk": 2}[spec]
            sc = 100 if spec == "q_CBF" else 1
            want = f"${rnd(s[mk] * sc, d)}\\vpm {rnd(s[sk] * sc, d)}$"
        if one != want:
            fail(f"{MAIN_MAZE_V7} string, want {want}", e)
    elif grp == "trajectory_events" and sub == spec:   # one decimal everywhere (100.0, no '100' shorthand)
        if one != "---" and re.search(r"(?<![\d.])100(?![\d.])", one):
            fail(f"{MAIN_MAZE_V7}: 100 instead of 100.0", e)
        check_main_maze(dict(e, printed=one), c, fail)
    else:
        fail(f"{MAIN_MAZE_V7}: unexpected source", e)


# main_maze_arxiv_v8 = the rows, views, statistics, rounding and bold marks of main_maze_arxiv_v7 as a one-column table
# (\footnotesize, tabcolsep >= 1.5 pt, no scaling, 245.7 pt): header arrows '{\downarrow}' without units, cells mean | {\scriptsize$\vpm SD$} with the SD
# leading zero dropped in m_c/c_max, T_mis, m_{p^max}, events count | \,(bound) with '100' for a full count or bound
# and one decimal otherwise, '---' over both halves. Every cell is re-derived from the exports.
MAIN_MAZE_V8 = "main_maze_arxiv_v8"
# main_corridor_arxiv_v8: arXiv Table IV, the same table for the Corridor; every cell rule of main_maze_arxiv_v8 applies
MAIN_CORR_V8 = "main_corridor_arxiv_v8"
MAIN_V8_FRAGS = (MAIN_MAZE_V8, MAIN_CORR_V8)
MAIN_MAZE_V8_NOLEAD = ("m_c", "T_mis", "pk")
MAIN_MAZE_V8_HEAD = (r"Method & \multicolumn{2}{c}{$\widehat m_c/c_{\max}{\downarrow}$} & "
                     r"\multicolumn{2}{c}{$\widehat T_{\rm mis}{\downarrow}$} & "
                     r"\multicolumn{2}{c}{$\widehat q^{+}_{\rm CBF}{\downarrow}$} & "
                     r"\multicolumn{2}{c}{$\widehat m_{p^{\max}}{\downarrow}$} & "
                     r"\multicolumn{2}{c}{$\widehat q_{\rm dep}{\downarrow}$} & "
                     r"\multicolumn{2}{c}{$\widehat q_{\rm cert}{\downarrow}$}\\\midrule")


def v8_join(t, spec):
    """A v8 cell as the one-cell string of v7 joined ('$m\\vpm s$' with the full SD, 'p (u)' with 100.0, '---'),
    or None when the cell breaks a v8 formatting rule."""
    if t == r"\multicolumn{2}{c}{---}":
        return "---"
    m = re.fullmatch(r"\$([^$]*)\$&\{\\scriptsize\$\\vpm ([0-9.]+)\$\}", t)
    if m:
        sd = m.group(2)
        if spec in MAIN_MAZE_V8_NOLEAD:
            if not re.fullmatch(r"\.\d+", sd):
                return None
            sd = "0" + sd
        elif sd.startswith("."):
            return None
        return f"${m.group(1)}\\vpm {sd}$"
    m = re.fullmatch(r"(\d+(?:\.\d)?)&\\,\((\d+(?:\.\d)?)\)", t)
    if m:
        p_, u_ = m.groups()
        if "100.0" in (p_, u_) or any("." not in x and x != "100" for x in (p_, u_)):
            return None   # '100' exactly for a full count or bound, one decimal otherwise
        full = lambda x: "100.0" if x == "100" else x
        return f"{full(p_)} ({full(u_)})"
    return None


def check_main_maze_v8(e, c, fail):
    """One cell of main_maze_arxiv_v8 or main_corridor_arxiv_v8, re-derived from the exports through the
    v7 rule (bold marks checked in arxiv_checks)."""
    frag = e["fragment"] if e["fragment"] in MAIN_V8_FRAGS else MAIN_MAZE_V8
    one = v8_join(e["printed"], e["spec"])
    if one is None:
        fail(f"{frag}: cell breaks the v8 format", e)
        return
    one = re.sub(r"\\mathbf\{([^{}]*)\}", r"\1", one)
    if one == "---":
        v7 = r"\multicolumn{2}{c}{---}"
    elif one.startswith("$"):
        mm = re.fullmatch(r"\$([^$\\]*)\\vpm ([0-9.]+)\$", one)
        v7 = f"${mm.group(1)}$&${{}}\\mathbin{{\\vpm}}{mm.group(2)}$"
    else:
        mm = re.fullmatch(r"(\S+) \((\S+)\)", one)
        v7 = f"{mm.group(1)}&\\ ({mm.group(2)})"
    check_main_maze_v7(dict(e, printed=v7), c, lambda msg, ee: fail(msg.replace(MAIN_MAZE_V7, frag), e))



def check_extended_compact(e, c, fail):
    if e["spec"] not in ("qdep", "qcert"):
        check_main_maze_v8(e, c, fail)
        return
    ev = c["trajectory_events"][e["spec"]]
    if not ev["applicable_n"]:
        want = r"\multicolumn{2}{c}{---}"
    else:
        count = str(ev["possible"])
        cfg = c["configuration"]
        if (cfg["benchmark"], cfg["method"], cfg["alpha0"], e["spec"]) == ("maze", "QT", 0.15, "qdep"):
            count = f"{ev['confirmed']}--{ev['possible']}"
        up = ceil_pct(cp_upper(ev["possible"], 200))
        up = "100" if up == "100.0" else up
        want = f"{count}&\\,({up})"
    if e["printed"] != want:
        fail(f"extended compact event: wanted {want}", e)


def check_starred_compact(e, fail):
    S = json.loads(STAR.read_text())
    c = next(x for x in S["cells"] if x["spec"] == e["json_spec"])
    original = f"${c['mean']}\\vpm {c['sd']}$" if c["kind"] == "meanSD" else f"{c['count']} ({c['upper_pct']})\\,\\checkmark"
    check_starred(dict(e, spec=e["json_spec"], printed=original), fail)
    if c["kind"] == "meanSD":
        d = {"m_c": 2, "T_mis": 1, "q_CBF": 1, "pk": 2}[e["spec"]]
        mean, sd = rnd(float(c["mean"]), d), rnd(float(c["sd"]), d)
        if e["spec"] in ("m_c", "T_mis", "pk"):
            sd = sd[1:]
        want = f"${mean}$&{{\\scriptsize$\\vpm {sd}$}}"
    else:
        want = f"{c['count']}&\\,({c['upper_pct']})"
    if e["printed"] != want:
        fail(f"starred compact: wanted {want}", e)


def check_percent_sweep(e, c, fail):
    grp, _, sub = e["source"].partition(".")
    ev = c[grp][sub]
    count = ev["possible"] if grp == "trajectory_events" else ev["any_excess"]
    up = ceil_pct(cp_upper(count, 200))
    cfg = c["configuration"]
    level = (cfg["delta"] if cfg["K"] != "infinity" else cfg["delta"] / (2 * cfg["budget_normalizer"] * math.log(2) ** 2)) if cfg["method"] in VARIANTS else None
    checked = level is not None and e["spec"] in ("qcert", "qdep") and Decimal(up) <= Decimal(repr(100 * level))
    show_check = checked and e["fragment"] != "stress_representative_arxiv_percent"
    want = f"{rnd(100 * count / 200, 1)} ({up})" + (r"\,\checkmark" if show_check else "")
    if e["printed"] != want or e["count"] != count:
        fail(f"percentage sweep: wanted {want}", e)


def extended_compact_checks(cells, gen, scen="maze"):
    frag = f"main_{scen}_extended_v8"
    es = [e for e in cells if e["target"] == "paper" and e["fragment"] == frag]
    if not es:
        return 0
    errors = []
    rows = {(e["row"].get("setting", "Baseline"), e["row"]["alpha0"], e["method"]) for e in es}
    # Read method IDs from the already independent-checked baseline cells: fixed baselines
    # have no certification setting and appear once at alpha0=0.15.
    baselines = {(e["row"]["alpha0"], e["method"]) for e in es if e["row"].get("baseline")}
    if len(baselines) != 7 or {(a, m) for a, m in baselines if m in ("ACI", "QT")} != {(a, m) for a in (0.15, 0.30) for m in ("ACI", "QT")}:
        errors.append("baseline groups")
    expected_rows = 29 if scen == "maze" else 28
    if len(es) != 6 * expected_rows or len(rows) != expected_rows:
        errors.append(f"{len(es)} cells/{len(rows)} methods, expected{6 * expected_rows}/{expected_rows}")
    for setting in ("S1", "S2", "S3", "S4", "S5"):
        if {e["method"] for e in es if e["row"].get("setting") == setting} != VARIANTS:
            errors.append(f"four variants missing at {setting}")
    for setting in (("Default", "S*") if scen == "maze" else ("Default",)):
        if {e["method"] for e in es if e["row"].get("setting") == setting} != {"V-C"}:
            errors.append(f"reference variant missing at {setting}")
    arx = {(e["row"].get("setting"), e["method"], e["spec"]): e["printed"] for e in cells
           if e["target"] == "paper" and e["fragment"] == f"main_{scen}_arxiv_v8" and e["row"].get("setting")}
    for e in es:
        key = e["row"].get("setting"), e["method"], e["spec"]
        if key in arx and e["spec"] not in ("qdep", "qcert") and e["printed"] != arx[key]:
            errors.append(f"mean/SD/bold mismatch {key}")
    txt = gen_text(gen / "paper" / f"{frag}.tex")
    if r"\footnotesize" not in txt or r"\tabcolsep}{1.2pt}" not in txt or "resizebox" in txt:
        errors.append("compact typography")
    for err in errors:
        print("FAIL extended compact", err)
    print(f"extended {scen} compact: {len(es)} cells, {len(rows)} methods; shared means/SD/bold agree with arXiv")
    return len(errors)

def check_main_maze_v3(e, c, fail, frag=MAIN_MAZE_V3, sd_specs=("m_c", "q_CBF")):
    grp, _, sub = e["source"].partition(".")
    spec = e["spec"]
    if spec not in MAIN_MAZE_V3_SPECS:
        fail(f"{frag}: unexpected column", e)
        return
    if grp == "statistics":
        if sub != ("p_max" if spec == "pk" else spec):
            fail(f"{frag}: wrong statistic", e)
            return
        s = c["statistics"][sub]
        mk = "available_mean" if spec == "pk" else "mean"
        if not s["applicable"] or s["observed_n"] == 0 or s[mk] is None:
            want = "---"
        else:
            d = {"m_c": 2, "T_mis": 1, "q_CBF": 1, "pk": 2}[spec]
            sc = 100 if spec == "q_CBF" else 1
            if spec in sd_specs:
                want = f"${rnd(s['mean'] * sc, d)}\\vpm {rnd(s['sd'] * sc, d)}$"
            else:
                want = f"${rnd(s[mk] * sc, d)}$"
        if e["printed"] != want:
            fail(f"{frag} string, want {want}", e)
    elif grp == "trajectory_events" and sub == spec:
        x = c[grp][sub]
        printed = e["printed"]
        if printed != "---":   # undo the '100' shorthand, then the main_maze_arxiv rule applies
            m = re.fullmatch(r"(\S+) \((\S+)\)", printed)
            if not m:
                fail(f"{frag}: event string", e)
                return
            p_, u_ = m.groups()
            if p_ == "100.0" or u_ == "100.0":
                fail(f"{frag}: 100.0 instead of 100", e)
            e = dict(e, printed=f"{'100.0' if p_ == '100' else p_} ({'100.0' if u_ == '100' else u_})")
        check_main_maze(e, c, fail)
    else:
        fail(f"{frag}: unexpected source", e)


# ------------------------------------------- assumption-check totals and arXiv parameter table (Appendix F-D)
AC_FLAGS = ("support_excess_observed", "fallback_infeasibility_observed", "filter_residual_failure_observed",
            "terminal_fallback_feasibility_unknown", "missing_tail")
AC_WANT_MACROS = ("acCfg", "acTraj", "acExec", "acChecked", "acNoChk", "acUnresDep", "acUnresCert", "acCert",
                  "gcDep", "gcMissFin", "gcMissInf", "gcMissInfLow", "gcMissInfCorr")


def check_assumption_totals(entries, views, gen):
    """Re-derives, from the export views (final and new_anchors, read by load_views), the
    assumption-check totals and the main-grid event counts of Appendix F-D of the arXiv, and compares them with the recorded values
    (kind a6_total) and with the printed macro file paper/assumption_checks_arxiv.tex. Scope: the certified variants of
    final main_grid_200 (all), of final detailed_200, budget_normalizer_200 and stress_drift_200 at the default anchor
    (anchor_id current_default), and of every new_anchors view; one configuration per configuration_id; phi = 1.
    Every failure count the text states as zero must be zero. Returns the number of failures."""
    nerr = 0

    def bad(msg):
        nonlocal nerr
        nerr += 1
        print("FAIL assumption checks", msg)

    cert = {"V-U", "V-D", "V-Cd", "V-C"}
    by_id = {}
    for (src, vid), v in views.items():
        if src not in ("final", "new"):
            continue
        d, c = v["definition"], v["cells"][0]
        if c["configuration"].get("method") not in cert:
            continue
        fam = d["family"]
        if src == "final" and not (fam == "main_grid_200" or (
                fam in ("detailed_200", "budget_normalizer_200", "stress_drift_200")
                and d.get("anchor_id") == "current_default")):
            continue
        cid = c["configuration_id"]
        if cid in by_id:
            o = by_id[cid]
            if any(o[f] != c[f] for f in ("risk_increase_checks", "numerical_ambiguity", "assumption_observed_flags")):
                bad(f"configuration {cid} carries different check fields in two views")
            continue
        by_id[cid] = c
    phi1 = [c for c in by_id.values() if float(c["configuration"]["bound_scale"]) == 1.0]
    want = defaultdict(int)
    zero = defaultdict(int)
    for c in phi1:
        r, na, out = c["risk_increase_checks"], c["numerical_ambiguity"], c["outcomes"]
        want["acCfg"] += 1
        want["acTraj"] += c["n"]
        want["acExec"] += c["assumption_assessed_population"]["executed_decisions"]
        want["acChecked"] += r["direct"]["assessed_decisions"]
        zero["checked direct != pairwise"] += r["direct"]["assessed_decisions"] != r["pairwise"]["assessed_decisions"]
        no = r["direct"]["assessed_decisions"] == 0
        want["acNoChk"] += no
        zero["V-C without a checked decision"] += no and c["configuration"]["method"] == "V-C"
        want["acUnresDep"] += na["deployment_intervals_unresolved"]
        want["acUnresCert"] += na["certificate_events_unresolved"]
        certificates = c["statistics"]["n_certificates"]
        count = certificates["mean"] * certificates["observed_n"]
        if certificates["observed_n"] != c["n"] or certificates["missing_n"] or abs(count - round(count)) > 1e-6:
            bad(f"certificate population or nonintegral total in {c['configuration_id']}")
        want["acCert"] += round(count)
        for side in ("direct", "pairwise"):
            zero[f"{side} excess decisions"] += r[side]["excess_decisions"]
            zero[f"{side} trajectories with an excess"] += r[side]["any_excess"]
        for f in AC_FLAGS:
            zero[f] += int(c["assumption_observed_flags"][f])
        zero["risk_numerical_errors"] += na["risk_numerical_errors"]
        zero["contact"] += out["contact"]["count"]
        zero["infeas"] += out["infeas"]["count"]
        zero["trajectories without the goal"] += c["n"] - out["goal"]["count"]
    for k, x in zero.items():
        if x:
            bad(f"{k} = {x} (the text states 0)")
    # the 3000 main-grid cells (certified, L_mu = 0.0005)
    bound = {}

    def misses(x, lev):
        k_ = (x["possible"], x["applicable_n"])
        if k_ not in bound:
            bound[k_] = cp_upper(*k_)
        if abs(bound[k_] - x["upper95_possible"]) > 1e-9:
            bad(f"stored upper95_possible {x['upper95_possible']} differs from Clopper-Pearson {bound[k_]}")
        return Decimal(ceil_pct(bound[k_])) > Decimal(repr(100 * lev))

    ncell = 0
    for (src, vid), v in views.items():
        if src != "final" or v["definition"]["family"] != "main_grid_200":
            continue
        c = v["cells"][0]
        cf = c["configuration"]
        if cf["method"] not in cert:
            continue
        ncell += 1
        if abs(cf["L_mu"] - 0.0005) > 1e-12 or float(cf["bound_scale"]) != 1.0:
            bad(f"main-grid cell {vid} not at L_mu = 0.0005, phi = 1")
        te = c["trajectory_events"]
        if te["qdep"]["possible"] > 0:
            want["gcDep"] += 1
            if cf["alpha0"] > 0.25 + 1e-12:
                bad(f"deployment event at alpha0 = {cf['alpha0']} (the text states alpha0 <= 0.25)")
        inf = cf["K"] == "infinity"
        lev = cf["delta"] / (2 * cf["budget_normalizer"] * math.log(2) ** 2) if inf else cf["delta"]
        mm = [ev for ev in ("qcert", "qdep") if te[ev]["applicable_n"] and misses(te[ev], lev)]
        if not mm:
            continue
        if not inf:
            want["gcMissFin"] += 1
            if mm != ["qcert"] or te["qcert"]["possible"] != 2 or te["qcert"]["applicable_n"] != 200 or \
                    cf["benchmark"] != "corridor" or (cf["alpha"], cf["delta"], cf["K"]) != (0.15, 0.025, 10):
                bad(f"finite-K miss {vid} is not a 2/200 q_cert count at (0.15, 0.025, 10) in the Corridor")
        else:
            want["gcMissInf"] += 1
            if abs(cf["delta"] - 0.025) < 1e-12:
                want["gcMissInfLow"] += 1
            elif cf["benchmark"] == "corridor" and abs(cf["alpha"] - 0.05) < 1e-12 and cf["delta"] in (0.05, 0.10):
                want["gcMissInfCorr"] += 1
            else:
                bad(f"K = infinity miss {vid} outside the two groups of the text")
    n_low = sum(1 for (src, _), v in views.items() if src == "final" and v["definition"]["family"] == "main_grid_200"
                and v["cells"][0]["configuration"]["method"] in cert and v["cells"][0]["configuration"]["K"] == "infinity"
                and abs(v["cells"][0]["configuration"]["delta"] - 0.025) < 1e-12)
    if ncell != 3000 or want["gcMissInfLow"] != n_low:
        bad(f"{ncell} main-grid cells (want 3000); {want['gcMissInfLow']} of {n_low} K = infinity cells with "
            f"delta = 0.025 miss")
    fmt = lambda n: f"{int(n):,}".replace(",", "{,}")
    rec = {e["name"]: e for e in entries}
    if sorted(rec) != sorted(AC_WANT_MACROS) or len(entries) != len(AC_WANT_MACROS):
        bad(f"recorded macros {sorted(rec)}")
    path = gen / "paper" / "assumption_checks_arxiv.tex"
    txt = path.read_text() if path.exists() else ""
    printed = dict(re.findall(r"^\\providecommand\{\\(\w+)\}\{\}\\renewcommand\{\\\1\}\{([^{}]*(?:\{,\}[^{}]*)*)\}$",
                              txt, re.M))
    body = [l for l in txt.splitlines() if l and not l.startswith("%")]
    if sorted(printed) != sorted(AC_WANT_MACROS) or len(body) != len(AC_WANT_MACROS):
        bad(f"macro file: {sorted(printed)}, {len(body)} non-comment lines")
    for m in AC_WANT_MACROS:
        w = fmt(want[m])
        if rec.get(m, {}).get("value") != want[m] or rec.get(m, {}).get("printed") != w or printed.get(m) != w:
            bad(f"\\{m}: re-derived {w}, recorded {rec.get(m, {}).get('printed')}, printed {printed.get(m)}")
    print("assumption-check totals: " + ", ".join(f"{m} {fmt(want[m])}" for m in AC_WANT_MACROS) +
          f" ({len(by_id) - len(phi1)} configurations with phi < 1 outside the totals)")
    return nerr


def strip_wrappers(s, names=("ownerfix", "releaseedit", "orderfix")):
    """Removes \\name{...} wrappers (balanced braces), keeping their content; nested wrappers too."""
    for name in names:
        tok = "\\" + name + "{"
        while tok in s:
            i = s.index(tok)
            j, depth = i + len(tok), 1
            while depth:
                depth += {"{": 1, "}": -1}.get(s[j], 0)
                j += 1
            s = s[:i] + s[i + len(tok):j - 1] + s[j:]
    return s


def brace_arg(s, i):
    """(content, end) of the brace group that opens at s[i] == '{'."""
    assert s[i] == "{", s[i:i + 20]
    j, depth = i + 1, 1
    while depth:
        depth += {"{": 1, "}": -1}.get(s[j], 0)
        j += 1
    return s[i + 1:j - 1], j


def expand_paper(s):
    """The extended paper text of a table branch: \\paperonly{X}, \\reponly{X} -> X (each row is compared alone),
    \\versionchoice{report}{A}{B} -> B, comment ends '%\\n' removed."""
    s = re.sub(r"%\n", "", s)
    for name in ("paperonly", "reponly"):
        tok = "\\" + name + "{"
        while tok in s:
            i = s.index(tok)
            a, j = brace_arg(s, i + len(tok) - 1)
            s = s[:i] + a + s[j:]
    tok = r"\versionchoice{report}{"
    while tok in s:
        i = s.index(tok)
        _, j = brace_arg(s, i + len(tok) - 1)
        b, k = brace_arg(s, j)
        s = s[:i] + b + s[k:]
    return s


def table_rows(branch):
    """Rows (list of cells, whitespace normalized, \\multicolumn unwrapped) of the tabular body of a branch."""
    i = branch.index(r"Quantity & Corridor & Maze\\\midrule") + len(r"Quantity & Corridor & Maze\\\midrule")
    body = branch[i:branch.index(r"\bottomrule")]
    rows = []
    for r in body.split("\\\\"):
        r = re.sub(r"\s+", " ", r).strip()
        r = re.sub(r"^\\midrule\s*", "", r).strip()
        if not r:
            continue
        cells = [x.strip() for x in r.split(" & ")]
        out = []
        for x in cells:
            if x.startswith(r"\multicolumn{"):
                _, j = brace_arg(x, len(r"\multicolumn"))
                _, j = brace_arg(x, j)
                x, _ = brace_arg(x, j)
            out.append(x)
        rows.append(out)
    return rows


PARAM_NEW_ROWS = ("Start; goal; goal tolerance", r"State bounds; braking $a_{\mathrm{br}}$",
                  r"Cone radius $\rho_j$; body radius $\rho_{\mathrm{body}}$", "Nose distance $l$", "Nominal input",
                  "Filter weight",
                  r"$\varrho_{\min}$ in \eqref{eq:unicycle-minimum-distance}; $K_x^z$; $K_u^z$ in "
                  r"\eqref{eq:unicycle-first-order-constants}, \eqref{eq:unicycle-second-order-constants}",
                  r"$L_{\mathrm{geom}}(c)$ in \eqref{eq:unicycle-clearance-sensitivity} over $\Gamma$",
                  r"$L_x^{\mathrm{geom}}$; $L_u^{\mathrm{geom}}$ in \eqref{eq:unicycle-risk-constants}",
                  r"ACI $k_{\mathrm{ACI}}$; $N^{\mathrm{win}}$; initial margin", r"QT $k_{\mathrm{QT}}$; initial margin")
# The arXiv names the unicycle constants L_x^geom, L_u^geom; the extended keeps L_x, L_u
# (cited by the submitted version). The values of the row must agree.
ARXIV_RENAMED = {r"$L_x^{\mathrm{geom}}$; $L_u^{\mathrm{geom}}$ in \eqref{eq:unicycle-risk-constants}":
                 r"$L_x$; $L_u$ in \eqref{eq:unicycle-risk-constants}"}
PARAM_GRID_ROWS = (r"Main grid: $\alpha_0$; $\alpha$", r"Main grid: $\delta$; $K$",
                   r"Default bound factor $\phi$ and budget normalizer $b$ (\appsecref{app:sim-assessments}); "
                   r"$L_\mu$ in \eqref{eq:distribution-drift-rate}", r"Selected settings $(\alpha_0,\alpha,\delta,K)$")
PARAM_SELECTED = {"S1": (0.30, 0.25, 0.20, 50), "S2": (0.30, 0.25, 0.15, 50), "S3": (0.30, 0.20, 0.20, 50),
                  "S4": (0.30, 0.25, 0.025, 100), "S5": (0.30, 0.25, 0.10, 10)}
PARAM_B = 2.11   # the printed budget normalizer (kappa rounded up)
# optional: an earlier appendices/f_tab_setup.tex without the eleven added rows in its arXiv branch (not included;
# check (c) is skipped when the file does not exist)
PREVIOUS_SETUP = Path(os.environ.get("ECBF_PREVIOUS_SETUP", SRC / "appendices" / "f_tab_setup_previous.tex"))


def check_arxiv_parameter_table(views):
    """arXiv Table III (the \\ifarxiv branch of appendices/f_tab_setup.tex, read only), after
    removing the \\ownerfix and \\releaseedit wrappers: (a) the main-grid, default and selected-setting rows equal the
    values re-derived from final main_grid_200; (b) the eleven added rows equal the extended rows of the same
    file (paper text) cell by cell, character for character; (c) the previous arXiv rows equal the arXiv branch of the
    earlier file (PREVIOUS_SETUP); no other row. Without that file, (c) and the row count are skipped with a NOTE.
    Before the grid rows exist: a NOTE, no check."""
    nerr = 0

    def bad(msg):
        nonlocal nerr
        nerr += 1
        print("FAIL Table III", msg)

    src_txt = (SRC / "appendices/f_tab_setup.tex").read_text()
    a0, a1 = src_txt.index("\\ifarxiv\n"), src_txt.index("\n\\else\n")
    arx = strip_wrappers(src_txt[a0:a1])
    if "Main grid" not in arx:
        print("NOTE Table III: the arXiv branch has no 'Main grid' row yet; check not run")
        return 0
    ext = expand_paper(strip_wrappers(src_txt[a1:src_txt.index("\\end{table*}", a1)])).replace("\\leavevmode", "")
    rows = table_rows(arx)
    byl = {r[0]: r for r in rows}
    if len(byl) != len(rows):
        bad("duplicate row labels")
    ext_rows = {r[0]: r for r in table_rows(ext)}
    if PREVIOUS_SETUP.exists():
        old = PREVIOUS_SETUP.read_text()
        old_rows = table_rows(strip_wrappers(old[old.index("\\ifarxiv\n"):old.index("\n\\else\n")]))
        old_keep = [r for r in old_rows if r[0] != "State bounds"]   # the old State bounds row became a (b) row
    else:
        old_keep = None
        print(f"NOTE Table III: {PREVIOUS_SETUP.name} not found (ECBF_PREVIOUS_SETUP); check (c) and the row count "
              "skipped")
    # (c) previous rows
    for r in old_keep or []:
        if byl.get(r[0]) != r:
            bad(f"(c) previous row {r[0]!r}: {byl.get(r[0])} != {r}")
    # (b) the eleven added rows
    for lab in PARAM_NEW_ROWS:
        if lab not in byl:
            bad(f"(b) row {lab!r} missing")
        elif byl[lab][1:] != (ext_rows.get(ARXIV_RENAMED.get(lab, lab)) or [None])[1:] or ARXIV_RENAMED.get(lab, lab) not in ext_rows:
            bad(f"(b) row {lab!r}: {byl[lab]} != extended {ext_rows.get(ARXIV_RENAMED.get(lab, lab))}")
    if old_keep is not None:
        extra = set(byl) - set(PARAM_NEW_ROWS) - set(PARAM_GRID_ROWS) - {r[0] for r in old_keep}
        if extra or len(rows) != len(old_keep) + len(PARAM_NEW_ROWS) + len(PARAM_GRID_ROWS):
            bad(f"rows not in (a)-(c): {sorted(extra)}; {len(rows)} rows")
    # (a) the grid rows from the exports
    grid = [v["cells"][0]["configuration"] for (src, _), v in views.items()
            if src == "final" and v["definition"]["family"] == "main_grid_200"]
    fv = lambda x: (lambda s_: s_ if len(s_.split(".")[1]) >= 2 else f"{float(x):.2f}")(f"{float(x):.4f}".rstrip("0"))
    fk = lambda x: r"\infty" if x == "infinity" else str(int(float(x)))
    a0s = sorted({c["alpha0"] for c in grid})
    pairs = sorted({(c["alpha0"], c["alpha"]) for c in grid})
    for x in a0s:   # alpha in {0.05, 0.10, ..., alpha0 - 0.05}
        al = sorted(a for a0_, a in pairs if a0_ == x)
        if [round(a, 6) for a in al] != [round(0.05 * k, 6) for k in range(1, round(x / 0.05))]:
            bad(f"(a) alpha set at alpha0 = {x}: {al} is not {{0.05, ..., alpha0 - 0.05}}")
    deltas = sorted({c["delta"] for c in grid})
    ks = sorted({c["K"] for c in grid}, key=lambda k: float("inf") if k == "infinity" else float(k))
    sett = {(c["alpha0"], c["alpha"], c["delta"], c["K"]) for c in grid}
    want_a = {
        PARAM_GRID_ROWS[0]: r"$\{" + ",".join(fv(x) for x in a0s) + r"\}$; $\{0.05,0.10,\ldots,\alpha_0-0.05\}$ ($" +
                            f"{len(pairs)}$ pairs)",
        PARAM_GRID_ROWS[1]: r"$\{" + ",".join(fv(x) for x in deltas) + r"\}$; $\{" + ",".join(fk(k) for k in ks) +
                            r"\}$ ($" + f"{len(sett)}$ settings)",
    }
    phis = {float(c["bound_scale"]) for c in grid}
    lms = {c["L_mu"] for c in grid}
    bs = {c["budget_normalizer"] for c in grid}
    if phis != {1.0} or lms != {0.0005} or bs != {PARAM_B}:
        bad(f"(a) main-grid defaults phi {phis}, L_mu {lms}, b {bs}")
    want_a[PARAM_GRID_ROWS[2]] = (f"${min(phis):g}$; $\\kappa$, rounded up to ${PARAM_B}$; ${fv(min(lms))}$"
                                  if len(phis) == len(lms) == 1 else None)
    for st, t in PARAM_SELECTED.items():
        if t not in {(a, b_, d, (int(float(k)) if k != "infinity" else k)) for a, b_, d, k in sett}:
            bad(f"(a) selected setting {st} {t} is not a main-grid setting")
    want_a[PARAM_GRID_ROWS[3]] = ", ".join(rf"\Tname{{{i}}}~$(" + ",".join(fv(x) for x in t[:3]) + f",{t[3]})$"
                                          for i, (st, t) in enumerate(PARAM_SELECTED.items(), 1))
    if len(pairs) != 15 or len(sett) != 375:
        bad(f"(a) {len(pairs)} pairs, {len(sett)} settings (want 15, 375)")
    for lab, w in want_a.items():
        r = byl.get(lab)
        if r is None or len(r) != 2 or r[1] != w:
            bad(f"(a) row {lab!r}: {r} (want [{lab!r}, {w!r}])")
    n_prev = "previous rows not compared" if old_keep is None else f"{len(old_keep)} previous"
    print(f"Table III (arXiv): {len(rows)} rows checked ({n_prev}, {len(PARAM_NEW_ROWS)} added, "
          f"{len(PARAM_GRID_ROWS)} grid and default rows)")
    return nerr


def check_file(values_path, gen, edition):
    """Checks one values file: tables_values.json.gz (editions, gen = generated_n200) or catalog_values.json.gz
    (catalog, gen = generated_n200/catalog). Returns the number of errors."""
    import gzip
    with gzip.open(values_path, "rt") as fh:
        V = json.load(fh)
    print(f"== {values_path.name} ({'editions' if edition else 'catalog'})")
    views = load_views(V)
    GEN = gen
    nerr = 0

    def fail(msg, e):
        nonlocal nerr
        nerr += 1
        if nerr <= 30:
            print("FAIL", msg, {k: e.get(k) for k in ("target", "fragment", "table", "method", "spec", "printed",
                                                     "view_id")})

    nstar = 0
    ac = [e for e in V["cells"] if e["kind"] == "a6_total"]   # checked by check_assumption_totals()
    cells = [e for e in V["cells"] if e["kind"] != "a6_total"]
    if edition:
        nerr += check_assumption_totals(ac, views, GEN)
        nerr += check_arxiv_parameter_table(views)
    elif ac:
        fail("assumption-check totals in the catalog values", ac[0])
    groups = defaultdict(list)
    msd = {}
    seen_views = set()
    for i, e in enumerate(cells):
        kind = e["kind"]
        if kind == "design_values":
            continue
        if kind == "starred_compact":
            check_starred_compact(e, fail)
            continue
        if kind == "starred":
            nstar += 1
            check_starred(e, fail)
            continue
        if kind in ("starred_variant", "starred_missing"):
            nstar += 1
            check_starred_variant(e, fail)
            continue
        if kind == "not_simulated":
            r = e["row"]
            hit = [v for (s, _), v in views.items() if s == "final" and v["definition"]["family"] in
                   ("matched_baselines_200", "matched_ACI_200") and v["cells"][0]["configuration"]["method"] == e["method"]
                   and abs(v["cells"][0]["configuration"]["alpha0"] - r["alpha0"]) < 1e-9
                   and abs(v["cells"][0]["configuration"]["L_mu"] - r["L_mu"]) < 1e-9
                   and v["cells"][0]["configuration"]["benchmark"] in r["scenario"].split("+")]
            if hit or e["printed"] != "--":
                fail("dash although simulated", e)
            continue
        v = views.get((e["src"], e["view_id"]))
        if v is None or len(v["cells"]) != 1:
            fail("view missing", e)
            continue
        c = v["cells"][0]
        if c["configuration"] != e["configuration"] or c["configuration"]["method"] != e["method"]:
            fail("configuration mismatch", e)
        if e["fragment"] in ("main_maze_extended_v8", "main_corridor_extended_v8"):
            check_extended_compact(e, c, fail)
            continue
        if kind == "event_pct_sweep":
            check_percent_sweep(e, c, fail)
            continue
        if e["fragment"] == MAIN_MAZE:
            check_main_maze(e, c, fail)
            continue
        if e["fragment"] == MAIN_MAZE_V2:
            check_main_maze_v2(e, c, fail)
            continue
        if e["fragment"] == MAIN_MAZE_V3:
            check_main_maze_v3(e, c, fail)
            continue
        if e["fragment"] == MAIN_MAZE_V4:
            check_main_maze_v3(e, c, fail, MAIN_MAZE_V4, ("m_c", "q_CBF", "T_mis"))
            continue
        if e["fragment"] == MAIN_MAZE_V5:
            check_main_maze_v3(e, c, fail, MAIN_MAZE_V5, ("m_c", "q_CBF", "T_mis"))
            continue
        if e["fragment"] in MAIN_V8_FRAGS:   # split cell re-derived (bold marks checked below)
            check_main_maze_v8(e, c, fail)
            continue
        if e["fragment"] == MAIN_MAZE_V7:   # split cell re-derived (bold marks checked below)
            check_main_maze_v7(e, c, fail)
            continue
        if e["fragment"] == MAIN_MAZE_V6:   # the cell without its bold mark (marks checked below)
            check_main_maze_v3(dict(e, printed=re.sub(r"\\mathbf\{([^{}]*)\}", r"\1", e["printed"])), c, fail,
                               MAIN_MAZE_V6, ("m_c", "q_CBF", "T_mis"))
            continue
        if (e["src"], e["view_id"]) not in seen_views:
            seen_views.add((e["src"], e["view_id"]))
            seeds = c.get("seeds") if e["src"] != "final_temporal" else c["statistics"].get("seeds")
            if c["n"] != 200 or seeds != SEEDS or not v["definition"]["family"].endswith("_200"):
                fail("n / seeds / family", e)
        grp, _, sub = e["source"].partition(".")
        cert = e["method"] in VARIANTS and not e["row"].get("baseline")
        if grp == "statistics":
            s = c["statistics"][sub]
            pk = e["spec"].split("#")[0].split("@")[0] == "pk"
            if pk and (sub != "p_max" or not e["fragment"].startswith("sec6_")):
                fail("peak-risk column pk outside Table I or not statistics.p_max", e)
            mk, sk = ("available_mean", "available_sd") if pk else ("mean", "sd")
            if not s["applicable"] or s["observed_n"] == 0 or (pk and s[mk] is None):
                if kind != "NA" or e["printed"] != r"\na":
                    fail("should be N-A", e)
                continue
            if pk:   # Table I's peak-risk column, a fraction with 3 decimals (available mean and SD)
                d, sc = PK_DEC, 1
            else:
                d = SEC6_DEC.get(sub, DEC[sub]) if e["fragment"].startswith("sec6_") else DEC[sub]
                sc = 100 if sub in PCT else 1
            m, sd = rnd(s[mk] * sc, d), rnd(s[sk] * sc, d)
            if sub in ("h_min", "clearance_min") and s["mean"] < 0 and Decimal(m) == 0:   # sign kept
                m = "-" + m
            msd[i] = (m, sd)
            if cert and sub not in NOBOLD and not s["conditional"]:
                cfg = {k: x for k, x in c["configuration"].items() if k != "method"}
                gk = (e["target"], e["fragment"], e["table"], json.dumps(cfg, sort_keys=True), e["spec"])
                groups[gk].append((i, Decimal(m), sub))
        elif grp == "temporal":
            t = c["statistics"]
            if sub == "nK":
                want = f"{t['normalized_margin']['assessed_n']}, {t['estimated_applied_risk']['assessed_n']}"
            else:
                q, stat, sc = TMP[sub]
                st = t[q][stat]
                if not st or st.get("n", 0) < 2 or st.get("mean") is None:
                    want = r"\na"
                else:
                    want = f"${rnd(st['mean'] * sc, 3)}\\vpm {rnd(st['sd'] * sc, 3)}$"
            if e["printed"] != want:
                fail(f"temporal string, want {want}", e)
        else:
            gname, key = EVSRC[e["spec"].split("@")[0]]
            assert (grp, sub) == (gname, key), (grp, sub, gname, key)
            x = c[grp][sub]
            if grp == "trajectory_events":
                n, a, u = x["applicable_n"], x["possible"], x["upper95_possible"]
                assert not n or x["possible"] >= x["confirmed"]
            elif grp == "realized_transition_events":
                n, a, u = x["assessed_n"], x["count"], None
            else:
                n, a = x["assessed_trajectory"], x["any_excess"]
                u = cp_upper(a, 200) if n else None
            if not n:
                if e["printed"] != r"\na":
                    fail("zero assessed must be N-A", e)
                continue
            up = ceil_pct(u) if u is not None else None
            want = (f"{a}" if n == 200 else f"{a}/{n}") + (f" ({up})" if up is not None else "")
            if e.get("count_rule") == "range":   # the submission's range cell (build_tables.SUBMISSION_RANGES)
                cf = c["configuration"]
                allowed = (e["target"] in ("paper", "report")
                           and e["fragment"] in ("sec6_combined", "default_maze", "checks_maze_L0005",
                                                 "checks_all_maze_L0005", "topfive_events")
                           and grp == "trajectory_events" and e["spec"].split("@")[0] == "qdep"
                           and cf["method"] == "QT"
                           and cf["benchmark"] == "maze" and abs(cf["alpha0"] - 0.15) < 1e-9
                           and abs(cf["L_mu"] - 0.0005) < 1e-12 and n == 200)
                if not allowed or x["confirmed"] == a:
                    fail("range cell outside the submission's Table I cell", e)
                want = f"{x['confirmed']}--{a} ({up})"
            if cert and e["spec"].split("@")[0] in ("qcert", "qdep"):
                cfg = c["configuration"]
                lev = cfg["delta"] if cfg["K"] != "infinity" else cfg["delta"] / (2 * cfg["budget_normalizer"] *
                                                                                  math.log(2) ** 2)
                if Decimal(up) <= Decimal(repr(100 * lev)):
                    want += r"\,\checkmark"
            if e["printed"] != want:
                fail(f"event string, want {want}", e)
    expected_bold = {i: False for i in msd}
    for gk, lst in groups.items():
        if len(lst) != 4:
            continue
        vals = [val for _, val, _ in lst]
        best = max(vals) if lst[0][2] in UP else min(vals)
        for (i, val, _) in lst:
            expected_bold[i] = val == best
    for i, (m, sd) in msd.items():
        mm = rf"\mathbf{{{m}}}" if expected_bold[i] else m
        want = f"${mm}\\vpm {sd}$"
        if cells[i]["printed"] != want:
            fail(f"mean/SD string, want {want}", cells[i])
    # design tables: no counts; the printed value sets are re-derived from the exports
    def cfgs(src, fams):
        return [v["cells"][0]["configuration"] for (s_, _), v in views.items() if s_ == src
                and v["definition"]["family"] in fams]

    def fmt_v(x):
        if x == "infinity":
            return r"\infty"
        s_ = f"{float(x):.4f}".rstrip("0")
        return s_ if len(s_.split(".")[1]) >= 2 else f"{float(x):.2f}"

    def fmt_k(x):
        return r"\infty" if x == "infinity" else str(int(float(x)))

    def vset(cs, k, f=fmt_v):
        vs = sorted({c[k] for c in cs if c.get(k) is not None}, key=lambda x: (isinstance(x, str), float(x)
                                                                              if not isinstance(x, str) else 0))
        return r"\{" + ",".join(f(v) for v in vs) + r"\}"

    gsh = lambda x: f"{float(x):g}"
    want_rows = {False: ["Main grid", "Baselines", "Default setting", "Selected settings", "Budget normalizer",
                         "Bound factor", "Test period"],
                 True: ["Main grid", "Baselines", "Selected settings", "Budget normalizer", "Bound factor",
                        "Test period"]}
    dsg = defaultdict(list)
    for e in cells:
        if e["kind"] == "design_values":
            dsg[(e["target"], e["fragment"])].append(e)
    for (target, frag), es in dsg.items():
        arx = frag.startswith("arxiv")
        if [e["row"]["study"] for e in es] != want_rows[arx]:
            fail(f"design rows {[e['row']['study'] for e in es]}", es[0])
        txt = gen_text(GEN / target / f"{frag}.tex")
        if "Count" in txt or "triples" in txt or "Comparison & Varied values & Fixed values" not in txt:
            fail("design table header (Count column or key line left, or wrong header)", es[0])
        for e in es:
            st, pr = e["row"]["study"], e["printed"]
            if st == "Main grid":
                m = cfgs("final", {"main_grid_200"})
                need = [rf"$\delta\in{vset(m, 'delta')}$", rf"$K\in{vset(m, 'K', fmt_k)}$"]
                for a0 in sorted({c["alpha0"] for c in m}):
                    al = sorted({c["alpha"] for c in m if c["alpha"] is not None and c["alpha0"] == a0})
                    rhs = rf"$\alpha={fmt_v(al[0])}$" if len(al) == 1 else \
                        r"$\alpha\in\{" + ",".join(fmt_v(x) for x in al) + r"\}$"
                    need.append(rf"$\alpha_0={fmt_v(a0)}$: & {rhs}")
            elif st == "Baselines":
                b_ = cfgs("final", {"matched_baselines_200", "matched_ACI_200"})
                need = [rf"$L_\mu\in{vset(b_, 'L_mu')}$",
                        rf"$\alpha_0\in{vset([c for c in b_ if c['method'] == 'QT'], 'alpha0')}$"]
                aci = {c["alpha0"] for c in b_ if c["method"] == "ACI"} - {c["alpha0"] for c in b_ if c["method"] == "QT"}
                need.append(r"ACI also at $\alpha_0\in\{" + ",".join(fmt_v(x) for x in sorted(aci)) + r"\}$")
            elif st in ("Budget normalizer", "Bound factor", "Test period"):
                fam = {"Budget normalizer": "budget_normalizer_200", "Bound factor": "stress_drift_200",
                       "Test period": "detailed_200"}[st]
                n_ = cfgs("new", {fam})
                # the default b = 2.11 prints as \kappa (b = kappa by default, F-D)
                gsh_b = lambda x: r"\kappa" if abs(float(x) - 2.11) < 1e-9 else gsh(x)
                need = {"Budget normalizer": [rf"$b\in{vset(n_, 'budget_normalizer', gsh_b)}$"],
                        "Bound factor": [rf"$\phi\in{vset(n_, 'bound_scale', gsh)}$", rf"$L_\mu\in{vset(n_, 'L_mu')}$"],
                        "Test period": [rf"$L_\mu\in{vset(n_, 'L_mu')}$"]}[st]
                if st == "Test period" and "infinity" not in {c["K"] for c in n_}:
                    fail("detailed_200 has no K = infinity", e)
                tops = {(0.3, 0.25, 0.2), (0.3, 0.25, 0.15), (0.3, 0.2, 0.2), (0.3, 0.25, 0.025), (0.3, 0.25, 0.1)}
                if {(c["alpha0"], c["alpha"], c["delta"]) for c in n_} != tops:
                    fail(f"{fam}: new settings are not S1-S5", e)
            elif st == "Selected settings":
                t5 = [(0.30, 0.25, 0.20, 50), (0.30, 0.25, 0.15, 50), (0.30, 0.20, 0.20, 50), (0.30, 0.25, 0.025, 100),
                      (0.30, 0.25, 0.10, 10)]
                need = [rf"\Tname{{{i}}}: & $(" + ",".join(fmt_v(x) for x in t[:3]) + f",{t[3]})$"
                        for i, t in enumerate(t5, 1)]
            else:
                need = []
            for x in need:
                if x not in pr:
                    fail(f"design values: missing {x}", e)
    # order and completeness inside each fragment
    byfrag = defaultdict(list)
    for e in cells:
        byfrag[(e["target"], e["fragment"])].append(e)
    # an event string stands alone or as the argument of \multirow{L}{*}{...} (one-member column of a box)
    pat_ev = re.compile(r"(?:(?<![\w{])|(?<=\{\*\}\{))\d+(?:/\d+)? \(\d+\.\d\)")
    for (target, frag), es in byfrag.items():
        txt = gen_text(GEN / target / f"{frag}.tex")
        body = "\n".join(l for l in txt.splitlines() if not l.startswith((r"\vnfit", r"\vnnew", "%")))
        pos = 0
        for e in es:
            j = body.find(e["printed"], pos)
            if j < 0:
                nerr += 1
                print("FAIL order", target, frag, e["printed"])
                break
            pos = j + len(e["printed"])
        n_pm = body.count(r"\vpm")
        n_rec = sum(1 for e in es if r"\vpm" in e["printed"])
        n_ev = len(pat_ev.findall(body))
        if frag in (MAIN_MAZE_V3, MAIN_MAZE_V4, MAIN_MAZE_V5, MAIN_MAZE_V6):   # '100' for a full count or bound
            n_ev = len(re.findall(r"(?<![\w{])(?:\d+\.\d|100) \((?:\d+\.\d|100)\)", body))
        if frag in MAIN_V8_FRAGS:   # count and bound in the two halves, a thin space before '('
            n_ev = len(re.findall(r"(?<![\w{])(?:\d+\.\d|100)&\\,\((?:\d+\.\d|100)\)", body))
        if frag in ("main_maze_extended_v8", "main_corridor_extended_v8"):
            n_ev = len(re.findall(r"(?<![\w{])\d+(?:--\d+)?&\\,\((?:\d+\.\d|100)\)", body))
        if frag.endswith("_arxiv_percent"):
            n_ev = len(re.findall(r"(?:(?<![\w{])|(?<=\{\*\}\{))\d+\.\d \(\d+\.\d\)", body))
        if frag == MAIN_MAZE_V7:   # count and bound in the two halves of a split column
            n_ev = len(re.findall(r"(?<![\w{])(?:\d+\.\d|100)&\\ \((?:\d+\.\d|100)\)", body))
        n_rec_ev = sum(1 for e in es if (e["kind"] in ("event", "event_pct", "event_count_compact", "event_pct_sweep") or (e["kind"] in ("starred", "starred_variant", "starred_compact") and r"\vpm" not in
                                                                   e["printed"])) and "(" in e["printed"])   # bound-free counts (CBF, unsafe) are checked cell by cell above
        if n_pm != n_rec or n_ev != n_rec_ev:
            nerr += 1
            print("FAIL completeness", target, frag, n_pm, n_rec, n_ev, n_rec_ev)
    nerr += baseline_box_checks(gen)
    nerr += tuple_checks(gen)
    if edition:
        nerr += uniform_layout_checks(gen)
        nerr += arxiv_checks(cells, byfrag, views, gen)
        nerr += secondary_arxiv_checks(cells, views, gen)
        nerr += sum(extended_compact_checks(cells, gen, scen) for scen in ("maze", "corridor"))
    if not edition:
        nerr += catalog_checks(cells, byfrag, views, gen)
        print(f"{len(cells)} recorded cells, {len(seen_views)} views, {sum(expected_bold.values())} bold marks, "
              f"{sum(1 for e in cells if r'checkmark' in e['printed'])} check marks")
        print("catalog errors:", nerr)
        return nerr
    # ---------------------------------------------------------------- structure checks
    PAIRED = {"Default": "V-C", "S1": "V-C", "S2": "V-C", "S3": "V-Cd", "S4": "V-D", "S5": "V-C"}
    SETTINGS = {"Default": (0.15, 0.10, 0.05, 50), "S1": (0.30, 0.25, 0.20, 50), "S2": (0.30, 0.25, 0.15, 50),
                "S3": (0.30, 0.20, 0.20, 50), "S4": (0.30, 0.25, 0.025, 100), "S5": (0.30, 0.25, 0.10, 10)}

    def fmt_setting(x):
        s_ = f"{float(x):.4f}".rstrip("0")
        return s_ if len(s_.split(".")[1]) >= 2 else f"{float(x):.2f}"

    def fmt_tuple(c):
        k = r"$\infty$" if c["K"] == "infinity" else str(int(float(c["K"])))
        return "(" + ", ".join([fmt_setting(c["alpha0"]), fmt_setting(c["alpha"]), fmt_setting(c["delta"]), k]) + ")"

    def fmt_fields(c, fields):   # fmt_tuple restricted to `fields`
        return "(" + ", ".join((r"$\infty$" if c[f] == "infinity" else str(int(float(c[f])))) if f == "K"
                               else fmt_setting(c[f]) for f in fields) + ")"

    def struct_fail(msg):
        nonlocal nerr
        nerr += 1
        print("FAIL structure", msg)

    paired_frag = re.compile(r"^(checks_\w+|stress_\w+_(?:arxiv|cases)|detailed_events_\w+|delta_\w+|period_\w+|b_\w+|phi_\w+|sec6_short|sec6_combined)$")
    seen_settings = defaultdict(set)
    nlab = 0
    for (target, frag), es in byfrag.items():
        txt = gen_text(GEN / target / f"{frag}.tex")
        for e in es:
            row = e.get("row") or {}
            if row.get("baseline") or e["kind"] == "design_values":
                continue
            if row.get("starred"):
                # star_maze (G-I.2) holds the S* setting with the four certified variants
                ok = (frag == "main_maze_extended_v8" and e["kind"] == "starred_compact") or (frag == "sec6_combined" and e["kind"] == "starred") or (frag == "star_maze" and e["kind"] in ("starred_variant", "starred_missing"))
                if not ok or row.get("label") not in txt:
                    struct_fail(f"{target}/{frag}: starred row misplaced or label not printed")
                continue
            sets_ = [x for x in str(row.get("setting", "")).split("+") if x]
            if paired_frag.match(frag):
                for st in sets_:
                    seen_settings[(target, frag)].add(st)
                    if PAIRED[st] != e["method"]:
                        struct_fail(f"{target}/{frag}: {e['method']} at {st} is not the chosen variant {PAIRED[st]}")
                        break
                if not sets_:
                    struct_fail(f"{target}/{frag}: row without setting")
            if frag.startswith(("topfive_full_", "sec6_short", "sec6_combined", "sec6_maze_arxiv")):
                for st in sets_:
                    c = e["configuration"]
                    if (c["alpha0"], c["alpha"], c["delta"], c["K"]) != SETTINGS[st]:
                        struct_fail(f"{target}/{frag}: configuration of {st} differs")
            if "label" in row:
                nlab += 1
                c = dict(e["configuration"])
                if e["spec"].endswith("@inf"):   # K = infinity column of a row labelled with the setting's K
                    if c["K"] != "infinity":
                        struct_fail(f"{target}/{frag}: K = infinity column holds K = {c['K']}")
                    c["K"] = SETTINGS[sets_[0]][3] if sets_ else c["K"]
                want = fmt_tuple(c)
                if row["label"] != want:
                    struct_fail(f"{target}/{frag}: label {row['label']} != configuration {want}")
                if row.get("block"):   # the block row names the setting tuple without the swept field
                    bf = tuple(row.get("block_fields", ("alpha0", "alpha", "delta", "K")))
                    swf = {"delta": "delta", "period": "K"}.get(row.get("swept"))
                    if swf and swf in bf:
                        struct_fail(f"{target}/{frag}: swept {swf} printed in the block row")
                    blk = [fmt_fields(dict(zip(("alpha0", "alpha", "delta", "K"), SETTINGS[st])), bf) for st in sets_]
                    if row["block"] != blk:
                        struct_fail(f"{target}/{frag}: block {row['block']} != settings {blk}")
                    if any(b_ not in txt for b_ in blk):
                        struct_fail(f"{target}/{frag}: block tuple {blk} not printed")
                    sw = {"delta": 2, "period": 3}.get(row.get("swept"))
                    st_t = list(SETTINGS[sets_[0]])
                    if sw is not None:
                        st_t[sw] = c["K"] if sw == 3 else c["delta"]
                    if fmt_tuple(dict(zip(("alpha0", "alpha", "delta", "K"), st_t))) != want:
                        struct_fail(f"{target}/{frag}: row {want} differs from its block beyond the swept value")
                elif want not in txt:
                    struct_fail(f"{target}/{frag}: label {want} not printed")
        nkey = gen_text(GEN / target / f"{frag}.tex").count("%KEYBOX")
        if nkey != (1 if frag in ("sec6_short", "sec6_combined", "sec6_maze_arxiv") else 0):
            struct_fail(f"{target}/{frag}: {nkey} key boxes")
    for target in ("paper", "report"):
        for f_ in GEN.glob(f"{target}/*.tex"):
            if re.match(r"^(sweep_|topfive_(corridor|maze|events)|sec6_topfive|phi_\w+_L00|lmu_)", f_.stem):
                struct_fail(f"fragment of an earlier layout still generated: {target}/{f_.name}")
        tk = gen_text(GEN / target / "table_key.tex")
        if tk.count("%KEYBOX") != 1:
            struct_fail(f"{target}/table_key: needs exactly one key box")
    for (target, frag), sts in seen_settings.items():
        if re.match(r"^(delta|period|b|phi)_", frag):
            want = {"Default", "S3"} if target == "paper" and not frag.startswith("phi_") else set(SETTINGS)
            if sts != want:
                struct_fail(f"{target}/{frag}: settings {sorted(sts)} != {sorted(want)}")
    for target in ("paper", "report"):
        for scen in ("corridor", "maze"):
            es = byfrag.get((target, f"topfive_full_{scen}"), [])
            got = {(e["row"].get("setting"), e["method"]) for e in es if not e["row"].get("baseline")}
            specs = {e["spec"] for e in es}
            if got != {(t, v) for t in ("S1", "S2", "S3", "S4", "S5") for v in VARIANTS} or not {"qdep", "qcert"} <= specs:
                struct_fail(f"{target}/topfive_full_{scen}: rows or q_dep/q_cert columns incomplete")
    if not re.search(r"\\definecolor\{vl\}\{HTML\}\{009E73\}", (SRC / "appendices/n200_assets.tex").read_text()):
        struct_fail("n200_assets.tex: V-L colour vl is not 009E73")
    print(f"structure: {nlab} row labels checked")
    # the starred row exists exactly once, in Maze, and sec6_short holds no starred cell
    es_star = [e for e in cells if e["kind"] == "starred"]
    S = json.loads(STAR.read_text())
    if len(es_star) != len(S["cells"]) or {e["fragment"] for e in es_star} != {"sec6_combined"}:
        struct_fail(f"starred cells: {len(es_star)} recorded, fragments {sorted({e['fragment'] for e in es_star})}")
    # the rows S1-S5 and the ACI/QT rows at alpha0 = 0.30 of sec6_combined print exactly the cells of sec6_short
    # (arXiv Table I); its rows at the default setting and its baselines at alpha0 = 0.15 print the cells of the
    # previous full Table I (compared with joint_main.tex below); the alpha0-free baselines appear once per scenario.
    def cell_key(e):
        r = e["row"]
        return (e["configuration"]["benchmark"], str(r.get("setting") or r.get("alpha0")), e["method"], e["spec"])
    a = sorted((cell_key(e), e["printed"]) for e in cells if e["fragment"] == "sec6_short" and e["target"] == "paper"
               and (not e["row"].get("baseline") or e["method"] in ("ACI", "QT")))
    b = sorted((cell_key(e), e["printed"]) for e in cells if e["fragment"] == "sec6_combined" and e["kind"] != "starred"
               and (e["row"].get("setting") in ("S1", "S2", "S3", "S4", "S5")
                    or (e["row"].get("baseline") and e["method"] in ("ACI", "QT") and abs(e["row"]["alpha0"] - 0.30) < 1e-9)))
    if a != b:
        struct_fail("sec6_combined differs from sec6_short in the rows S1-S5 or the ACI/QT rows at alpha0 = 0.30")
    # arXiv Table I (sec6_maze_arxiv) = Maze, S1-S5 x four certified variants, baselines at
    # alpha0 = 0.30; its chosen-variant rows and its baselines equal the Maze cells of sec6_short (bold aside)
    unb = lambda t: re.sub(r"\\mathbf\{([^}]*)\}", r"\1", t)
    sm = [e for e in cells if e["fragment"] == "sec6_maze_arxiv" and e["target"] == "paper"]
    got = {(e["row"].get("setting"), e["method"]) for e in sm if not e["row"].get("baseline")}
    if got != {(t, v) for t in ("S1", "S2", "S3", "S4", "S5") for v in VARIANTS}:
        struct_fail(f"sec6_maze_arxiv: certified rows {sorted(got)}")
    gotb = sorted({(e["method"], e["row"]["alpha0"]) for e in sm if e["row"].get("baseline")})
    if gotb != sorted((b_, 0.30) for b_ in ("Unfiltered", "Zero margin", "Conservative", "ACI", "QT")):
        struct_fail(f"sec6_maze_arxiv: baseline rows {gotb}")
    if {e["configuration"]["benchmark"] for e in sm} != {"maze"}:
        struct_fail("sec6_maze_arxiv: a row outside the Maze")
    if sorted({e["spec"] for e in sm}) != sorted(["m_c#3", "q_CBF", "pk", "qdep", "qcert", "T_mis"]):
        struct_fail(f"sec6_maze_arxiv: columns {sorted({e['spec'] for e in sm})}")
    ss = {cell_key(e): e["printed"] for e in cells if e["fragment"] == "sec6_short" and e["target"] == "paper"
          and e["configuration"]["benchmark"] == "maze"}
    ncross = 0
    for e in sm:
        if not e["row"].get("baseline") and PAIRED[e["row"]["setting"]] != e["method"]:
            continue
        k_ = cell_key(e)
        if k_ not in ss or unb(ss[k_]) != unb(e["printed"]):
            struct_fail(f"sec6_maze_arxiv {k_}: {e['printed']} differs from sec6_short {ss.get(k_)}")
        ncross += 1
    print(f"sec6_maze_arxiv: {len(sm)} cells, {ncross} compared with sec6_short (Maze)")
    # main_maze_arxiv = Maze, baselines at alpha0 = 0.30, S3 and S5 x four certified variants, four columns;
    # its S3 and S5 rows and its baselines hold the same cells (views) as sec6_maze_arxiv
    mm = [e for e in cells if e["fragment"] == MAIN_MAZE and e["target"] == "paper"]
    got = {(e["row"].get("setting"), e["method"]) for e in mm if not e["row"].get("baseline")}
    if got != {(t, v) for t in ("S3", "S5") for v in VARIANTS}:
        struct_fail(f"main_maze_arxiv: certified rows {sorted(got)}")
    gotb = sorted({(e["method"], e["row"]["alpha0"]) for e in mm if e["row"].get("baseline")})
    if gotb != sorted((b_, 0.30) for b_ in ("Unfiltered", "Zero margin", "Conservative", "ACI", "QT")):
        struct_fail(f"main_maze_arxiv: baseline rows {gotb}")
    if {e["configuration"]["benchmark"] for e in mm} != {"maze"} or \
            sorted({e["spec"] for e in mm}) != sorted(["m_c", "q_CBF", "qdep", "qcert"]):
        struct_fail("main_maze_arxiv: scenario or columns")
    vid = {(e["method"], str(e["row"].get("setting"))): e["view_id"] for e in sm}
    for e in mm:
        if vid.get((e["method"], str(e["row"].get("setting")))) != e["view_id"]:
            struct_fail(f"main_maze_arxiv: {e['method']} {e['row'].get('setting')} not the cell of sec6_maze_arxiv")
            break
    txt_mm = gen_text(GEN / "paper" / f"{MAIN_MAZE}.tex")
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark)|\|", txt_mm) or \
            r"\begin{tabular}{@{}lcccc@{}}" not in txt_mm:
        struct_fail("main_maze_arxiv: not a plain booktabs table (rules, boxes, colours or check marks)")
    print(f"main_maze_arxiv: {len(mm)} cells")
    # main_maze_arxiv_v2 = the rows and views of main_maze_arxiv, six metrics on two lines per method
    m2 = [e for e in cells if e["fragment"] == MAIN_MAZE_V2 and e["target"] == "paper"]
    rows_mm = {(e["method"], str(e["row"].get("setting")), e["view_id"]) for e in mm}
    rows_m2 = {(e["method"], str(e["row"].get("setting")), e["view_id"]) for e in m2}
    if rows_m2 != rows_mm or len(m2) != 6 * len(rows_mm):
        struct_fail(f"main_maze_arxiv_v2: rows differ from main_maze_arxiv or not six cells per row ({len(m2)})")
    for e in m2:
        if e["spec"] not in MAIN_MAZE_V2_LINES.get(e["row"].get("line"), ()):
            struct_fail(f"main_maze_arxiv_v2: {e['spec']} on line {e['row'].get('line')}")
            break
    same = {(e["method"], str(e["row"].get("setting")), e["spec"]): e["printed"] for e in mm}
    for e in m2:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        if k_ in same and same[k_] != e["printed"]:
            struct_fail(f"main_maze_arxiv_v2 {k_}: {e['printed']} differs from main_maze_arxiv {same[k_]}")
    pk7 = {(e["method"], str(e["row"].get("setting")), e["spec"]): e["printed"] for e in sm
           if e["spec"] in ("pk", "T_mis")}
    ncmp = 0
    for e in m2:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        if k_ in pk7:
            ncmp += 1
            if unb(pk7[k_]).replace(r"\na", "---") != e["printed"]:
                struct_fail(f"main_maze_arxiv_v2 {k_}: {e['printed']} differs from sec6_maze_arxiv {pk7[k_]}")
    txt_m2 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V2}.tex")
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox)|\|", txt_m2) \
            or r"\begin{tabular}{@{}lccc@{}}" not in txt_m2 or r"{\small" not in txt_m2:
        struct_fail("main_maze_arxiv_v2: not a plain booktabs table at \\small (rules, boxes, colours, scaling)")
    print(f"main_maze_arxiv_v2: {len(m2)} cells, {ncmp} T_mis/pk cells compared with sec6_maze_arxiv")
    # main_maze_arxiv_v3 = baselines (alpha0 = 0.30) and S1-S5 x four variants, the views of sec6_maze_arxiv,
    # six cells per row; values shared with main_maze_arxiv(_v2) at the same precision must agree
    m3 = [e for e in cells if e["fragment"] == MAIN_MAZE_V3 and e["target"] == "paper"]
    got3 = {(e["row"].get("setting"), e["method"]) for e in m3 if not e["row"].get("baseline")}
    if got3 != {(f"S{i}", v) for i in range(1, 6) for v in VARIANTS}:
        struct_fail(f"main_maze_arxiv_v3: certified rows {sorted(got3)}")
    gotb3 = sorted({(e["method"], e["row"]["alpha0"]) for e in m3 if e["row"].get("baseline")})
    if gotb3 != gotb:
        struct_fail(f"main_maze_arxiv_v3: baseline rows {gotb3}")
    rows3 = {(e["method"], str(e["row"].get("setting"))) for e in m3}
    if len(m3) != 6 * len(rows3) or len(rows3) != 25 or \
            {e["configuration"]["benchmark"] for e in m3} != {"maze"}:
        struct_fail(f"main_maze_arxiv_v3: not six cells per row of 25 Maze rows ({len(m3)} cells)")
    for k_ in rows3:
        if sorted(e["spec"] for e in m3 if (e["method"], str(e["row"].get("setting"))) == k_) != \
                sorted(MAIN_MAZE_V3_SPECS):
            struct_fail(f"main_maze_arxiv_v3: columns of {k_}")
            break
    for e in m3:
        if vid.get((e["method"], str(e["row"].get("setting")))) != e["view_id"]:
            struct_fail(f"main_maze_arxiv_v3: {e['method']} {e['row'].get('setting')} not the cell of sec6_maze_arxiv")
            break
    same1 = {(e["method"], str(e["row"].get("setting")), e["spec"]): e["printed"] for e in mm if e["spec"] == "m_c"}
    for e in m3:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        if k_ in same1 and same1[k_] != e["printed"]:
            struct_fail(f"main_maze_arxiv_v3 {k_}: {e['printed']} differs from main_maze_arxiv {same1[k_]}")
    txt_m3 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V3}.tex")
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox)|\|", txt_m3) \
            or r"\begin{tabular}{@{}lcccccc@{}}" not in txt_m3 or r"{\small" not in txt_m3 \
            or txt_m3.count(r"\midrule") != 2 or txt_m3.count(r"\cmidrule{1-7}") != 4:
        struct_fail("main_maze_arxiv_v3: not the plain booktabs layout (rules, boxes, colours, scaling)")
    print(f"main_maze_arxiv_v3: {len(m3)} cells in {len(rows3)} rows")
    # main_maze_arxiv_v4 = the rows, rules and columns of main_maze_arxiv_v3, the views of sec6_maze_arxiv, one
    # header line with the units, T_mis mean +- SD; every cell except T_mis printed as in main_maze_arxiv_v3
    m4 = [e for e in cells if e["fragment"] == MAIN_MAZE_V4 and e["target"] == "paper"]
    rows4 = {(e["method"], str(e["row"].get("setting"))) for e in m4}
    if rows4 != rows3 or len(m4) != 6 * len(rows4):
        struct_fail(f"main_maze_arxiv_v4: rows differ from main_maze_arxiv_v3 ({len(m4)} cells)")
    for k_ in rows4:
        if sorted(e["spec"] for e in m4 if (e["method"], str(e["row"].get("setting"))) == k_) != \
                sorted(MAIN_MAZE_V3_SPECS):
            struct_fail(f"main_maze_arxiv_v4: columns of {k_}")
            break
    for e in m4:
        if vid.get((e["method"], str(e["row"].get("setting")))) != e["view_id"]:
            struct_fail(f"main_maze_arxiv_v4: {e['method']} {e['row'].get('setting')} not the cell of sec6_maze_arxiv")
            break
    same3 = {(e["method"], str(e["row"].get("setting")), e["spec"]): e["printed"] for e in m3}
    for e in m4:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        if e["spec"] == "T_mis":
            if not (e["printed"] == same3.get(k_) == "---" or e["printed"].startswith(same3.get(k_, "?")[:-1] + r"\vpm ")):
                struct_fail(f"main_maze_arxiv_v4 {k_}: T_mis mean {e['printed']} differs from main_maze_arxiv_v3 {same3.get(k_)}")
        elif same3.get(k_) != e["printed"]:
            struct_fail(f"main_maze_arxiv_v4 {k_}: {e['printed']} differs from main_maze_arxiv_v3 {same3.get(k_)}")
    txt_m4 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V4}.tex")
    head4 = [l for l in txt_m4.splitlines() if l.startswith("Method & ")]
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox|makebox)|\|", txt_m4) \
            or r"\begin{tabular}{@{}lcccccc@{}}" not in txt_m4 \
            or not re.search(r"\{\\(small|footnotesize)\\setlength", txt_m4) \
            or txt_m4.count(r"\midrule") != 2 or txt_m4.count(r"\cmidrule{1-7}") != 4 \
            or len(head4) != 1 or not head4[0].endswith(r"\\\midrule") \
            or [h.split("$")[-1].strip() for h in head4[0][:-len(r"\\\midrule")].split(" & ")[1:]] != \
            ["", "(s)", r"(\%)", "", r"(\%)", r"(\%)"]:
        struct_fail("main_maze_arxiv_v4: not the one-header-line booktabs layout (units, rules, font, scaling)")
    print(f"main_maze_arxiv_v4: {len(m4)} cells in {len(rows4)} rows")
    # main_maze_arxiv_v5 = every row and cell of main_maze_arxiv_v4 (same views, same strings), \footnotesize,
    # tabcolsep 2.5 pt, one header line of the six symbols without units, the m_{p^max} header in a zero-width box
    m5 = [e for e in cells if e["fragment"] == MAIN_MAZE_V5 and e["target"] == "paper"]
    rows5 = {(e["method"], str(e["row"].get("setting"))) for e in m5}
    if rows5 != rows4 or len(m5) != 6 * len(rows5):
        struct_fail(f"main_maze_arxiv_v5: rows differ from main_maze_arxiv_v4 ({len(m5)} cells)")
    same4 = {(e["method"], str(e["row"].get("setting")), e["spec"]): (e["printed"], e["view_id"]) for e in m4}
    for e in m5:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        if same4.get(k_) != (e["printed"], e["view_id"]):
            struct_fail(f"main_maze_arxiv_v5 {k_}: {e['printed']} differs from main_maze_arxiv_v4 {same4.get(k_)}")
    txt_m5 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V5}.tex")
    head5 = [l for l in txt_m5.splitlines() if l.startswith("Method & ")]
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox)|\|", txt_m5) \
            or r"\begin{tabular}{@{}lcccccc@{}}" not in txt_m5 \
            or r"{\footnotesize\setlength{\tabcolsep}{2.5pt}" not in txt_m5 \
            or txt_m5.count(r"\midrule") != 2 or txt_m5.count(r"\cmidrule{1-7}") != 4 \
            or txt_m5.count(r"\makebox") != 1 or len(head5) != 1 or not head5[0].endswith(r"\\\midrule") \
            or any(h.split("$")[-1].strip() not in ("", "}") for h in head5[0][:-len(r"\\\midrule")].split(" & ")[1:]):
        struct_fail("main_maze_arxiv_v5: not the one-header-line booktabs layout (units, rules, font, spacing, scaling)")
    print(f"main_maze_arxiv_v5: {len(m5)} cells in {len(rows5)} rows")
    # main_maze_arxiv_v6 = every row and cell of main_maze_arxiv_v5 (same views, same strings once the bold
    # marks are removed); bold exactly on the lowest printed mean of the four variants of one setting in m_c/c_max,
    # T_mis, q_CBF^+ (ties all bold); one header line with the report's arrows; setting rows bold, a \cmidrule
    # directly below each, variants unindented
    m6 = [e for e in cells if e["fragment"] == MAIN_MAZE_V6 and e["target"] == "paper"]
    rows6 = {(e["method"], str(e["row"].get("setting"))) for e in m6}
    if rows6 != rows5 or len(m6) != 6 * len(rows6):
        struct_fail(f"main_maze_arxiv_v6: rows differ from main_maze_arxiv_v5 ({len(m6)} cells)")
    same5 = {(e["method"], str(e["row"].get("setting")), e["spec"]): (e["printed"], e["view_id"]) for e in m5}
    unbold = lambda t: re.sub(r"\\mathbf\{([^{}]*)\}", r"\1", t)
    for e in m6:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        if same5.get(k_) != (unbold(e["printed"]), e["view_id"]):
            struct_fail(f"main_maze_arxiv_v6 {k_}: {e['printed']} differs from main_maze_arxiv_v5 {same5.get(k_)}")
    nb6 = 0
    for st in sorted({str(e["row"].get("setting")) for e in m6 if not e["row"].get("baseline")}):
        for sp in MAIN_MAZE_V3_SPECS:
            grp6 = [e for e in m6 if str(e["row"].get("setting")) == st and e["spec"] == sp]
            if len(grp6) != 4:
                struct_fail(f"main_maze_arxiv_v6 {st} {sp}: {len(grp6)} variants")
                continue
            means = [re.match(r"\$(?:\\mathbf\{)?([0-9.]+)", e["printed"]) for e in grp6]
            want = set()
            if sp in MAIN_MAZE_V6_BOLD and all(means):
                vals = [Decimal(m.group(1)) for m in means]
                want = {i for i, v in enumerate(vals) if v == min(vals)}
            got = {i for i, e in enumerate(grp6) if r"\mathbf" in e["printed"]}
            nb6 += len(got)
            if got != want:
                struct_fail(f"main_maze_arxiv_v6 {st} {sp}: bold on {sorted(got)}, want {sorted(want)}")
    if any(r"\mathbf" in e["printed"] for e in m6 if e["row"].get("baseline")):
        struct_fail("main_maze_arxiv_v6: bold baseline cell")
    txt_m6 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V6}.tex")
    head6 = [l for l in txt_m6.splitlines() if l.startswith("Method & ")]
    want_head6 = (r"Method & $\widehat m_c/c_{\max}\!\downarrow$ & $\widehat T_{\rm mis}\!\downarrow$ & "
                  r"$\widehat q^{+}_{\rm CBF}\!\downarrow$ & \makebox[0pt]{\hspace{4.8pt}$\widehat m_{p^{\max}}$} & $\widehat q_{\rm dep}\!\downarrow$ & "
                  r"$\widehat q_{\rm cert}\!\downarrow$\\\midrule")
    setrows6 = re.findall(r"^\\multicolumn\{7\}\{@\{\}l\}\{\\textbf\{(S[1-5])\}: .*\}\\\\\\cmidrule\{1-7\}$", txt_m6, re.M)
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox|quad)|\|", txt_m6) or txt_m6.count(r"\makebox") != 1 \
            or r"\begin{tabular}{@{}lcccccc@{}}" not in txt_m6 \
            or not re.search(r"\{\\footnotesize\\setlength\{\\tabcolsep\}\{(2|2\.5)pt\}", txt_m6) \
            or txt_m6.count(r"\midrule") != 2 or txt_m6.count(r"\cmidrule{1-7}") != 9 \
            or head6 != [want_head6] or setrows6 != ["S1", "S2", "S3", "S4", "S5"]:
        struct_fail("main_maze_arxiv_v6: not the report-style one-header-line booktabs layout")
    print(f"main_maze_arxiv_v6: {len(m6)} cells in {len(rows6)} rows, {nb6} bold marks")
    # main_maze_arxiv_v7 = every row of main_maze_arxiv_v6 (same views); joined and unbolded, every cell equals
    # its v6 cell except m_{p^max}, whose mean equals the v6 mean and which carries the SD; bold exactly on the lowest
    # printed mean of the four variants of one setting in m_c/c_max, T_mis, q_CBF^+, m_{p^max} (ties all bold); the one
    # header line with symbols, arrows and units; booktabs, no scaling, \small or \footnotesize, tabcolsep >= 2 pt
    m7 = [e for e in cells if e["fragment"] == MAIN_MAZE_V7 and e["target"] == "paper"]
    rows7 = {(e["method"], str(e["row"].get("setting"))) for e in m7}
    if rows7 != rows6 or len(m7) != 6 * len(rows7):
        struct_fail(f"main_maze_arxiv_v7: rows differ from main_maze_arxiv_v6 ({len(m7)} cells)")
    same6 = {(e["method"], str(e["row"].get("setting")), e["spec"]): (unbold(e["printed"]), e["view_id"]) for e in m6}
    for e in m7:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        one = v7_join(e["printed"])
        one = None if one is None else unbold(one)
        w6, vid6 = same6.get(k_, (None, None))
        if w6 and e["spec"] in ("qdep", "qcert"):   # v7 prints 100.0 where v6 prints the shorthand 100
            w6 = re.sub(r"(?<![\d.])100(?![\d.])", "100.0", w6)
        if vid6 != e["view_id"]:
            struct_fail(f"main_maze_arxiv_v7 {k_}: view differs from main_maze_arxiv_v6")
        elif e["spec"] == "pk":
            if not (one == w6 == "---" or (one or "").startswith(w6[:-1] + r"\vpm ")):
                struct_fail(f"main_maze_arxiv_v7 {k_}: pk {one} differs from main_maze_arxiv_v6 {w6}")
        elif one != w6:
            struct_fail(f"main_maze_arxiv_v7 {k_}: {one} differs from main_maze_arxiv_v6 {w6}")
    nb7 = 0
    for st in sorted({str(e["row"].get("setting")) for e in m7 if not e["row"].get("baseline")}):
        for sp in MAIN_MAZE_V3_SPECS:
            grp7 = [e for e in m7 if str(e["row"].get("setting")) == st and e["spec"] == sp]
            if len(grp7) != 4:
                struct_fail(f"main_maze_arxiv_v7 {st} {sp}: {len(grp7)} variants")
                continue
            means = [re.match(r"\$(?:\\mathbf\{)?([0-9.]+)", e["printed"]) for e in grp7]
            want = set()
            if sp in MAIN_MAZE_V7_BOLD and all(means):
                vals = [Decimal(m.group(1)) for m in means]
                want = {i for i, v in enumerate(vals) if v == min(vals)}
            got = {i for i, e in enumerate(grp7) if r"\mathbf" in e["printed"]}
            nb7 += len(got)
            if got != want:
                struct_fail(f"main_maze_arxiv_v7 {st} {sp}: bold on {sorted(got)}, want {sorted(want)}")
    if any(r"\mathbf" in e["printed"] for e in m7 if e["row"].get("baseline")):
        struct_fail("main_maze_arxiv_v7: bold baseline cell")
    txt_m7 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V7}.tex")
    head7 = [l for l in txt_m7.splitlines() if l.startswith("Method & ")]
    setrows7 = re.findall(r"^\\multicolumn\{13\}\{@\{\}l\}\{\\textbf\{(S[1-5])\}: .*\}\\\\\\cmidrule\{1-13\}$", txt_m7, re.M)
    sep7 = re.search(r"\{\\(small|footnotesize)\\setlength\{\\tabcolsep\}\{([0-9.]+)pt\}", txt_m7)
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox|makebox|quad)|\|", txt_m7) \
            or r"\begin{tabular}{@{}l" + "r@{}l" * 6 + "@{}}" not in txt_m7 \
            or not sep7 or float(sep7.group(2)) < 2 \
            or txt_m7.count(r"\midrule") != 2 or txt_m7.count(r"\cmidrule{1-13}") != 9 \
            or head7 != [MAIN_MAZE_V7_HEAD] or setrows7 != ["S1", "S2", "S3", "S4", "S5"]:
        struct_fail("main_maze_arxiv_v7: not the split-column one-header-line booktabs layout")
    print(f"main_maze_arxiv_v7: {len(m7)} cells in {len(rows7)} rows, {nb7} bold marks")
    # main_maze_arxiv_v8 = every row of main_maze_arxiv_v7 (same views); joined, every cell equals its v7 cell
    # (same bold marks); the one header line with tight arrows and no units; booktabs, no scaling, \footnotesize,
    # tabcolsep >= 1.5 pt, every SD in \scriptsize
    m8 = [e for e in cells if e["fragment"] == MAIN_MAZE_V8 and e["target"] == "paper"]
    rows8 = {(e["method"], str(e["row"].get("setting"))) for e in m8}
    if rows8 != rows7 or len(m8) != 6 * len(rows8):
        struct_fail(f"main_maze_arxiv_v8: rows differ from main_maze_arxiv_v7 ({len(m8)} cells)")
    same7 = {(e["method"], str(e["row"].get("setting")), e["spec"]): (v7_join(e["printed"]), e["view_id"]) for e in m7}
    nb8 = 0
    for e in m8:
        k_ = (e["method"], str(e["row"].get("setting")), e["spec"])
        one = v8_join(e["printed"], e["spec"])
        w7, vid7 = same7.get(k_, (None, None))
        nb8 += r"\mathbf" in e["printed"]
        if vid7 != e["view_id"]:
            struct_fail(f"main_maze_arxiv_v8 {k_}: view differs from main_maze_arxiv_v7")
        elif one is None or one != w7:
            struct_fail(f"main_maze_arxiv_v8 {k_}: {e['printed']} differs from main_maze_arxiv_v7 {w7}")
    if nb8 != nb7:
        struct_fail(f"main_maze_arxiv_v8: {nb8} bold marks, main_maze_arxiv_v7 has {nb7}")
    txt_m8 = gen_text(GEN / "paper" / f"{MAIN_MAZE_V8}.tex")
    head8 = [l for l in txt_m8.splitlines() if l.startswith("Method & ")]
    setrows8 = re.findall(r"^\\multicolumn\{13\}\{@\{\}l\}\{\\textbf\{(S[1-5])\}: .*\}\\\\\\cmidrule\{1-13\}$", txt_m8, re.M)
    sep8 = re.search(r"\{\\footnotesize\\setlength\{\\tabcolsep\}\{([0-9.]+)pt\}", txt_m8)
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox|makebox|quad)|\|", txt_m8) \
            or r"\begin{tabular}{@{}l" + "r@{}l" * 6 + "@{}}" not in txt_m8 \
            or not sep8 or float(sep8.group(1)) < 1.5 \
            or txt_m8.count(r"\midrule") != 2 or txt_m8.count(r"\cmidrule{1-13}") != 9 \
            or txt_m8.count(r"\vpm") != txt_m8.count(r"{\scriptsize$\vpm ") \
            or head8 != [MAIN_MAZE_V8_HEAD] or setrows8 != ["S1", "S2", "S3", "S4", "S5"]:
        struct_fail("main_maze_arxiv_v8: not the one-column split-column one-header-line booktabs layout")
    print(f"main_maze_arxiv_v8: {len(m8)} cells in {len(rows8)} rows, {nb8} bold marks")
    # main_corridor_arxiv_v8 (arXiv Table IV) = the rows of main_maze_arxiv_v8 for the Corridor: the five
    # baselines at alpha0 = 0.30 (L_mu = 0.0005) and S1-S5 with the four certified variants (main_grid_200), 25 rows
    # of six cells; bold exactly on the lowest printed mean of the four variants of one setting in m_c/c_max, T_mis,
    # q_CBF^+, m_{p^max} (ties all bold), never on a baseline; the layout and header of main_maze_arxiv_v8; the
    # Unfiltered T_mis cell prints the horizon 70.0 +- .0 s (contact in every trajectory, Table IV caption)
    mc = [e for e in cells if e["fragment"] == MAIN_CORR_V8 and e["target"] == "paper"]
    rowsc = {(e["method"], str(e["row"].get("setting"))) for e in mc}
    want_rc = {(b_, "None") for b_ in ("Unfiltered", "Zero margin", "Conservative", "ACI", "QT")} | \
        {(v_, st) for st in ("S1", "S2", "S3", "S4", "S5") for v_ in VARIANTS}
    if rowsc != want_rc or len(mc) != 150:
        struct_fail(f"{MAIN_CORR_V8}: {len(rowsc)} rows, {len(mc)} cells (want 25 rows, 150 cells)")
    for e in mc:
        cf = e["configuration"]
        if cf["benchmark"] != "corridor" or abs(cf["L_mu"] - 0.0005) > 1e-12:
            struct_fail(f"{MAIN_CORR_V8}: {e['method']} not a Corridor cell at L_mu = 0.0005")
        if e["row"].get("baseline"):
            if abs(cf["alpha0"] - 0.30) > 1e-9 and e["method"] in ("ACI", "QT"):
                struct_fail(f"{MAIN_CORR_V8}: baseline {e['method']} not at alpha0 = 0.30")
            if r"\mathbf" in e["printed"]:
                struct_fail(f"{MAIN_CORR_V8}: bold baseline cell")
        else:
            st = e["row"].get("setting")
            if (cf["alpha0"], cf["alpha"], cf["delta"], cf["K"]) != SETTINGS.get(st) or \
                    views[(e["src"], e["view_id"])]["definition"]["family"] != "main_grid_200":
                struct_fail(f"{MAIN_CORR_V8}: {st} {e['method']} not the main-grid cell of the setting")
    nbc = 0
    for st in ("S1", "S2", "S3", "S4", "S5"):
        for sp in MAIN_MAZE_V3_SPECS:
            grp_ = [e for e in mc if str(e["row"].get("setting")) == st and e["spec"] == sp]
            if len(grp_) != 4:
                struct_fail(f"{MAIN_CORR_V8} {st} {sp}: {len(grp_)} variants")
                continue
            means = [re.match(r"\$(?:\\mathbf\{)?([0-9.]+)", e["printed"]) for e in grp_]
            want = set()
            if sp in MAIN_MAZE_V7_BOLD and all(means):
                vals = [Decimal(m.group(1)) for m in means]
                want = {i for i, v in enumerate(vals) if v == min(vals)}
            got = {i for i, e in enumerate(grp_) if r"\mathbf" in e["printed"]}
            nbc += len(got)
            if got != want:
                struct_fail(f"{MAIN_CORR_V8} {st} {sp}: bold on {sorted(got)}, want {sorted(want)}")
    unf = [e["printed"] for e in mc if e["method"] == "Unfiltered" and e["spec"] == "T_mis"]
    if unf != [r"$70.0$&{\scriptsize$\vpm .0$}"]:
        struct_fail(f"{MAIN_CORR_V8}: Unfiltered T_mis cell {unf} (want the horizon 70.0 +- .0)")
    txt_mc = gen_text(GEN / "paper" / f"{MAIN_CORR_V8}.tex")
    headc = [l for l in txt_mc.splitlines() if l.startswith("Method & ")]
    setrowsc = re.findall(r"^\\multicolumn\{13\}\{@\{\}l\}\{\\textbf\{(S[1-5])\}: .*\}\\\\\\cmidrule\{1-13\}$", txt_mc, re.M)
    sepc = re.search(r"\{\\footnotesize\\setlength\{\\tabcolsep\}\{([0-9.]+)pt\}", txt_mc)
    if re.search(r"\\(Xhline|hline|vrule|cline|rowcolor|cellcolor|textcolor|checkmark|resizebox|scalebox|makebox|quad)|\|", txt_mc) \
            or r"\begin{tabular}{@{}l" + "r@{}l" * 6 + "@{}}" not in txt_mc \
            or not sepc or float(sepc.group(1)) < 1.5 \
            or txt_mc.count(r"\midrule") != 2 or txt_mc.count(r"\cmidrule{1-13}") != 9 \
            or txt_mc.count(r"\vpm") != txt_mc.count(r"{\scriptsize$\vpm ") \
            or headc != [MAIN_MAZE_V8_HEAD] or setrowsc != ["S1", "S2", "S3", "S4", "S5"]:
        struct_fail(f"{MAIN_CORR_V8}: not the layout of main_maze_arxiv_v8")
    print(f"{MAIN_CORR_V8}: {len(mc)} cells in {len(rowsc)} rows, {nbc} bold marks")
    for scen in ("corridor", "maze"):
        got = [(e["row"].get("setting"), e["method"], e["row"].get("alpha0")) for e in cells
               if e["fragment"] == "sec6_combined" and e["configuration"]["benchmark"] == scen and e["spec"] == "m_c#3"]
        want = ([(None, b_, 0.15) for b_ in ("Unfiltered", "Zero margin", "Conservative")]
                + [(None, b_, a0) for b_ in ("ACI", "QT") for a0 in (0.15, 0.30)])
        got_b = sorted((m, a0) for s_, m, a0 in got if s_ is None)
        if got_b != sorted((m, a0) for _, m, a0 in want):
            struct_fail(f"sec6_combined {scen}: baseline rows {got_b}")
        got_v = [(s_, m) for s_, m, _ in got if s_ is not None]
        want_v = [(st, PAIRED[st]) for st in ("Default", "S1", "S2", "S3", "S4", "S5")] + (
            [("S*", "V-C")] if scen == S["scenario"] else [])
        if got_v != want_v:
            struct_fail(f"sec6_combined {scen}: certified rows {got_v}")
    # Table I of both editions has the peak-risk column pk in every row (both scenarios,
    # baselines and certified rows), N-A for Unfiltered and c_max, and its header has no arrow
    for frag in ("sec6_short", "sec6_combined", "sec6_maze_arxiv"):
        rows_all = {(e["configuration"]["benchmark"] if e["kind"] != "starred" else "maze", e["method"],
                     str(e["row"].get("setting") or e["row"].get("alpha0"))) for e in cells
                    if e["fragment"] == frag and e["target"] == "paper"}
        pks = {(e["configuration"]["benchmark"] if e["kind"] != "starred" else "maze", e["method"],
                str(e["row"].get("setting") or e["row"].get("alpha0"))): e for e in cells
               if e["fragment"] == frag and e["target"] == "paper" and e["spec"] == "pk"}
        if set(pks) != rows_all or not pks:
            struct_fail(f"{frag}: pk column missing in rows {sorted(rows_all - set(pks))}")
        for k_, e in pks.items():
            if k_[1] in ("Unfiltered", "Conservative") and e["printed"] != r"\na":
                struct_fail(f"{frag}: pk of {k_} must be N-A")
        txt = gen_text(GEN / "paper" / f"{frag}.tex")
        if r"\widehat m_{p^{\max}}$" not in txt or re.search(r"m_\{p\^\{\\max\}\}\\!\\(down|up)arrow", txt):
            struct_fail(f"{frag}: pk header missing or with an arrow")
        print(f"pk column {frag}: {len(pks)} cells")
    print(f"starred row: {nstar} cells checked against {STAR.name} and {S['source_file']}")
    # the V-C row of star_maze (G-I.2) prints the cells of row S* of sec6_combined (bold aside; a mean
    # or SD of row S* may print with more digits there, and then rounds half up to the star_maze string), the
    # fragment has the four certified variants and every column of Tables XIII/XIV
    sv = {e["spec"]: re.sub(r"\\mathbf\{([^}]*)\}", r"\1", e["printed"]) for e in cells
          if e["fragment"] == "star_maze" and e["method"] == "V-C"}
    s6 = {e["spec"].partition("#")[0]: e["printed"] for e in cells if e["fragment"] == "sec6_combined" and e["kind"] == "starred"}
    for sp, p6 in s6.items():
        if sp == "pk":   # S* peak risk as a fraction (Table I) equals the star_maze p_max cell in % (x 100)
            m6 = re.fullmatch(r"\$(.+)\\vpm (.+)\$", p6)
            m7 = re.fullmatch(r"\$(.+)\\vpm (.+)\$", sv.get("p_max") or "")
            if not (m6 and m7 and all(Decimal(a) * 100 == Decimal(b) for a, b in zip(m6.groups(), m7.groups()))):
                struct_fail(f"row S* pk {p6} differs from star_maze V-C p_max {sv.get('p_max')}")
            continue
        p7 = sv.get(sp)
        m6, m7 = re.fullmatch(r"\$(.+)\\vpm (.+)\$", p6), re.fullmatch(r"\$(.+)\\vpm (.+)\$", p7 or "")
        if m6 and m7:
            same = all(b == a or (not a.startswith("{<}") and "." in b and
                                  rnd(float(a), len(b.split(".")[1])) == b and len(b) < len(a))
                       for a, b in zip(m6.groups(), m7.groups()))
        else:
            same = p6 == p7
        if not same:
            struct_fail(f"star_maze V-C {sp} {p7} differs from row S* of sec6_combined {p6}")
    if set(sv) != set(STAR_SPEC_BLOCK) | STAR_NONE:
        struct_fail(f"star_maze columns {sorted(sv)}")
    if sorted({e["method"] for e in cells if e["fragment"] == "star_maze"}) != sorted(VARIANTS):
        struct_fail("star_maze: not the four certified variants")
    # Section VI Table I versus the previous Table I
    old = (SRC / "appendices/sim_tables/joint_main.tex").read_text()
    names = {"V-U": "V-U", "V-$\\Delta$": "V-D", "V-L": "V-Cd", "V-C": "V-C", "Unfiltered": "Unfiltered",
             "$c{=}0$": "Zero margin", "$c_{\\max}$": "Conservative", "ACI": "ACI", "QT": "QT"}
    order = ["m_c", "q_use", "q_CBF", "q_unsafe", "qdep", "qcert"]
    oldvals = {}
    for scen, blk in zip(("corridor", "maze"), old.split(r"\textit{Maze: all nine methods}")):
        for line in blk.splitlines():
            parts = [p.strip() for p in line.rstrip("\\").split(" & ")]
            if parts[0] in names and len(parts) == 7:
                for met, p in zip(order, parts[1:]):
                    oldvals[(scen, names[parts[0]], met)] = p
    ncmp = 0
    for e in cells:
        if e["target"] != "paper" or e["fragment"] != "sec6_combined" or e["kind"] == "starred":
            continue
        if e["row"].get("setting") not in (None, "Default") or (e["row"].get("baseline") and abs(e["row"]["alpha0"] - 0.15) > 1e-9):
            continue
        met = e["spec"].split("#")[0]
        o = oldvals.get((e["configuration"]["benchmark"], e["method"], met))
        if o is None:
            continue
        ncmp += 1
        p = e["printed"].replace(r"\mathbf{", "").replace("}", "").replace(r"\vpm ", r"\pm").replace("$", "")
        p = re.sub(r"\\,\\checkmark", "", p)
        oo = o.replace(r"\mathbf{", "").replace("}", "").replace("$", "").replace(r"\mbox{", "").replace(" ", "")
        oo = oo.replace("\\pm", "\\pm")
        if met in ("qdep", "qcert"):
            if o.strip() == "---":
                if e["printed"] != r"\na":
                    nerr += 1
                    print("FAIL Table I", e["method"], met, o, e["printed"])
                continue
            mo = re.match(r"(\d+)(?:--(\d+))?/200\((\d+\.\d)\)", oo)
            want = f"{mo.group(2) or mo.group(1)} ({mo.group(3)})"
            if e.get("count_rule") == "range":   # printed as in the submission's Table I: the full source range
                want = f"{mo.group(1)}--{mo.group(2)} ({mo.group(3)})" if mo.group(2) else want
            if not p.startswith(want):
                nerr += 1
                print("FAIL Table I", e["method"], met, o, e["printed"])
        else:
            if o.strip() == "---":
                if e["printed"] != r"\na":
                    nerr += 1
                    print("FAIL Table I", e["method"], met, o, e["printed"])
                continue
            if p.replace(" ", "") != oo:
                nerr += 1
                print("FAIL Table I", e["configuration"]["benchmark"], e["method"], met, o, e["printed"])
    print(f"{len(cells)} recorded cells, {len(seen_views)} views, {sum(expected_bold.values())} bold marks, "
          f"{sum(1 for e in cells if r'checkmark' in e['printed'])} check marks, {ncmp} Table I values compared")
    print("edition errors:", nerr)
    return nerr


ARX_SPECS = {"checks": ["qcert", "qdep", "rCBF", "runsafe", "qA2", "qpw"],
             "stress": ["m_c", "q_CBF", "qdep", "qA2", "qpw"]}


def zero_from_view(c, spec):
    """Re-derivation of an omitted column from the export cell: 0 or N-A (event: zero count or nothing assessed;
    mean: N-A or 0.00 +- 0.00 at the printed decimals)."""
    if spec in EVSRC:
        grp, sub = EVSRC[spec]
        x = c[grp][sub]
        if grp == "trajectory_events":
            return not x["applicable_n"] or x["possible"] == 0
        if grp == "realized_transition_events":
            return not x["assessed_n"] or x["count"] == 0
        return not x["assessed_trajectory"] or x["any_excess"] == 0
    st = c["statistics"][spec]
    if not st["applicable"] or st["observed_n"] == 0 or st["mean"] is None:
        return True
    sc = 100 if spec in PCT else 1
    return Decimal(rnd(st["mean"] * sc, DEC[spec])) == 0 and Decimal(rnd(st["sd"] * sc, DEC[spec])) == 0


def arxiv_checks(cells, byfrag, views, gen):
    """The arXiv checks and stress-test fragments (omitted columns, rows, layout)."""
    nerr = 0

    def bad(msg):
        nonlocal nerr
        nerr += 1
        print("FAIL arxiv", msg)

    SET = {(0.3, 0.25, 0.2, 50): "S1", (0.3, 0.25, 0.15, 50): "S2", (0.3, 0.2, 0.2, 50): "S3",
           (0.3, 0.25, 0.025, 100): "S4", (0.3, 0.25, 0.1, 10): "S5"}
    jobs = [(scen, kind, "paper", f"{kind}_{scen}_arxiv") for scen in ("corridor", "maze") for kind in ("checks", "stress")]
    jobs += [(scen, "stress", "report", f"stress_{scen}_cases") for scen in ("corridor", "maze")]   # report stress tables
    for scen, kind, target, frag in jobs:
            es = byfrag.get((target, frag), [])
            if target == "report":   # the report prints exactly the cells of the arXiv stress table
                pa = [(e["view_id"], e["spec"], e["printed"]) for e in byfrag.get(("paper", f"stress_{scen}_arxiv"), [])]
                if [(e["view_id"], e["spec"], e["printed"]) for e in es] != pa:
                    bad(f"report/{frag}: cells differ from paper/stress_{scen}_arxiv")
            path = gen / target / f"{frag}.tex"
            if not es or not path.exists():
                bad(f"{frag}: not generated")
                continue
            txt = gen_text(path)
            m = re.search(r"^% OMITTED \(0 or N-A in every row\): (.*)$", txt, re.M)
            if not m:
                bad(f"{frag}: no OMITTED line")
                continue
            gone = [] if m.group(1).strip() == "none" else [x.strip() for x in m.group(1).split(",")]
            printed = [x for x in ARX_SPECS[kind] if any(e["spec"].split("@")[0] == x for e in es)]   # @K/@inf columns
            if set(printed) | set(gone) != set(ARX_SPECS[kind]) or set(printed) & set(gone):
                bad(f"{frag}: printed {printed} + omitted {gone} != {ARX_SPECS[kind]}")
            rows = {}
            for e in es:
                rows.setdefault((e["src"], e["view_id"]), views[(e["src"], e["view_id"])]["cells"][0])
            want_n = 20 if kind == "checks" else 50
            fam = "detailed_200" if kind == "checks" else "stress_drift_200"
            if len(rows) != want_n:
                bad(f"{frag}: {len(rows)} rows, want {want_n}")
            got = set()
            for (src, vid), c in rows.items():
                cf = c["configuration"]
                if views[(src, vid)]["definition"]["family"] != fam or cf["benchmark"] != scen:
                    bad(f"{frag}: {vid} family/scenario")
                K = cf["K"] if kind == "stress" else None
                t = (cf["alpha0"], cf["alpha"], cf["delta"],
                     cf["K"] if kind == "stress" or cf["K"] != "infinity" else None)
                name = next((n for k, n in SET.items() if k[:3] == t[:3] and (t[3] is None or k[3] == t[3])), None)
                if name is None:
                    bad(f"{frag}: {vid} not at S1-S5")
                    continue
                got.add((name, cf["L_mu"], cf["K"] == "infinity" if kind == "checks" else cf["bound_scale"]))
                for x in gone:
                    if not zero_from_view(c, x):
                        bad(f"{frag}: omitted column {x} is not 0/N-A at {vid}")
            if kind == "checks":
                want = {(n, lm, inf) for n in SET.values() for lm in (0.0005, 0.005) for inf in (False, True)}
            else:
                want = {(n, lm, p) for n in SET.values() for lm in (0.0005, 0.005) for p in (0.0, 0.25, 0.5, 0.75, 1.0)}
            if got != want:
                bad(f"{frag}: row set differs ({len(got)} of {len(want)})")
            for x in printed:
                if all(zero_from_view(c, x) for c in rows.values()):
                    bad(f"{frag}: printed column {x} is 0/N-A in every row")
            if kind == "stress":
                body = "\n".join(l for l in txt.splitlines() if not l.startswith((r"\vnfit", r"\vnnew", "%")))
                nd = len(re.findall(r"\$0(?:\.\d+)?\^\\dagger\$", body))
                n1 = len(re.findall(r"(?m)^\$1\$ &", body))   # phi is the first column
                if nd != 40 or n1 != 10:
                    bad(f"{frag}: {nd} dagger rows (want 40), {n1} phi = 1 rows (want 10)")
            nerr += arxiv_layout_uniform(frag, kind, txt, bad)
    print("arXiv checks and stress fragments checked")
    return nerr


def arxiv_layout_earlier(frag, kind, txt, bad):
    """Earlier layout of the arXiv checks and stress fragments (see the module docstring). Calls bad() per failure;
    returns 0 (bad() counts)."""
    lines = txt.splitlines()
    nock_def = next((i for i, l in enumerate(lines) if l == r"\providecommand{\vnnock}{\phantom{\,\checkmark}}"), None)
    first_use = next((i for i, l in enumerate(lines) if r"\vnnock" in l and not l.startswith(r"\providecommand")), None)
    if nock_def is None or (first_use is not None and first_use < nock_def):
        bad(f"{frag}: \\vnnock not defined before its first use")
    if not re.search(r"\\ifdim\\vnsep\{(\w+)\}>6pt\\vnnewsep\{\1\}\{6pt\}\\fi", txt):
        bad(f"{frag}: column separation not capped at 6pt")
    if (kind == "checks") != (r"\renewcommand{\arraystretch}{1.1}" in txt):
        bad(f"{frag}: arraystretch 1.1 in the checks tables only")
    k0 = next(i for i, l in enumerate(lines) if l.startswith(r"{\setlength{\tabcolsep}") and r"\begin{tabular}" in l)
    spec = re.search(r"\\begin\{tabular\}\{(.*?)\}\\Xhline", lines[k0]).group(1)
    cols = re.findall(r"w\{([crl])\}\{\\vnw\{\w+?(t|l|meth|v\d+)\}\}", spec)
    if [a for a, _ in cols[:3]] != ["c", "c", "c"] or any(a != "r" for a, _ in cols[3:]) or len(cols) < 4:
        bad(f"{frag}: level and Method columns centred, value columns right-aligned: {cols}")
    nv = len(cols) - 3
    if r"\linewidth-2\vnbox" in txt or "+ " + str(2 * len(cols) - 2) + r"\tabcolsep + " + str(nv) + r"\arrayrulewidth" not in txt:
        bad(f"{frag}: settings row not at the table's natural width")
    head_i = next(i for i, l in enumerate(lines) if i > k0 and l.endswith(r"\vnsd\\"))
    head = lines[head_i][:-len(r"\vnsd\\")].split(" & ")
    want0 = r"\begin{tabular}{@{}c@{}}Method\\$(\alpha_0,\alpha,\delta,K)$\end{tabular}" if kind == "stress" else ""
    if head[0] != want0:
        bad(f"{frag}: first-column header {head[0]!r}")
    for n, c in enumerate(head[3:]):
        right = "|" if n < nv - 1 else r"!{\vrule width\vnbox}"
        if not c.startswith(r"\multicolumn{1}{c" + right + "}{"):
            bad(f"{frag}: value header {n} not centred: {c[:40]}")
    body = [l for l in lines[head_i + 1:] if l.endswith("\\\\") and " & " in l]
    rules = [l for l in lines[head_i + 1:] if l.startswith((r"\Xcline", r"\Xhline"))]
    if any(l.startswith(r"\Xcline{3-") for l in rules):
        bad(f"{frag}: a rule inside an innermost block")
    want_rows = 20 if kind == "checks" else 50
    # rules after the header: the header rule, 9 between the 10 (setting, L_mu) blocks, the closing box rule
    if len(body) != want_rows or len(rules) != 11:
        bad(f"{frag}: {len(body)} rows (want {want_rows}), {len(rules)} rules (want 11)")
    cells = [l[:-2].split(" & ") for l in body]
    for j in range(3, 3 + nv):
        col = [r[j] for r in cells]
        has = any(r"\checkmark" in c for c in col)
        for c in col:
            want_pad = has and r"\checkmark" not in c
            if c.endswith(r"\vnnock") != want_pad or c.count(r"\vnnock") > 1:
                bad(f"{frag}: column {j}: {c!r} (\\vnnock iff no check mark in a column with check marks)")
    return 0


def arxiv_layout_uniform(frag, kind, txt, bad):
    """Layout of the arXiv checks and stress fragments under the uniform layout of Table X
    (uniform_layout_checks() covers the common rules), here the structure: checks - first column L_mu, column groups 'K of the
    setting' and 'K = infinity' with the same event columns, 5 block rows (one per setting), 10 rows (two L_mu per
    setting), one page; stress - first column phi, 10 block rows (one per (setting, L_mu)), 50 rows, two pages of five
    blocks; \\vnnock exactly on the cells without a check mark in a column with check marks. Returns 0 (bad() counts)."""
    lines = txt.splitlines()
    nock_def = next((i for i, l in enumerate(lines) if l == r"\providecommand{\vnnock}{\phantom{\,\checkmark}}"), None)
    first_use = next((i for i, l in enumerate(lines) if r"\vnnock" in l and not l.startswith(r"\providecommand")), None)
    if first_use is not None and (nock_def is None or first_use < nock_def):
        bad(f"{frag}: \\vnnock not defined before its first use")
    heads = [i for i, l in enumerate(lines) if l.endswith(r"\vnsd\\")]
    if not heads:
        bad(f"{frag}: no header line")
        return 0
    head = lines[heads[0]][:-len(r"\vnsd\\")].split(" & ")
    want0 = r"$\phi$" if kind == "stress" else ""
    if head[0] != want0:
        bad(f"{frag}: first-column header {head[0]!r} (want {want0!r})")
    if kind == "checks":
        grp = lines[heads[0] - 1]
        if not grp.startswith(r"\multirow{2}{*}{$L_\mu$} & \multicolumn{") or "$K$ of the setting" not in grp \
                or r"$K=\infty$" not in grp:
            bad(f"{frag}: column groups 'K of the setting' / 'K = infinity' missing")
        nv = len(head) - 1
        if nv % 2 or head[1:1 + nv // 2] != head[1 + nv // 2:]:
            bad(f"{frag}: the two K groups do not show the same columns")
    blocks = [l for l in lines if l.startswith(r"\multicolumn{") and r"@{\hskip 4pt}l!{\vrule width\vnbox}}" in l]
    rows = [l for l in lines if l.endswith("\\\\") and " & " in l and not l.endswith(r"\vnsd\\")
            and not l.startswith(r"\multirow{2}{*}{$L_\mu$}")]
    nb, nr, npg = (5, 10, 1) if kind == "checks" else (10, 50, 2)
    if len(blocks) != nb or len(rows) != nr or txt.count(r"\ntablebreak") != npg - 1:
        bad(f"{frag}: {len(blocks)} block rows (want {nb}), {len(rows)} rows (want {nr}), "
            f"{txt.count(chr(92) + 'ntablebreak') + 1} pages (want {npg})")
    if kind == "checks":
        for l in rows:
            if not l.startswith(("$0.0005$ & ", "$0.005$ & ")):
                bad(f"{frag}: row label {l[:20]!r} is not an L_mu value")
                break
    cells = [l[:-2].split(" & ") for l in rows]
    for j in range(1, len(head)):
        col = [r[j] for r in cells if len(r) > j]
        has = any(r"\checkmark" in c for c in col)
        for c in col:
            want_pad = has and r"\checkmark" not in c
            if c.endswith(r"\vnnock") != want_pad or c.count(r"\vnnock") > 1:
                bad(f"{frag}: column {j}: {c!r} (\\vnnock iff no check mark in a column with check marks)")
                break
    return 0


UNIFORM_EXEMPT = re.compile(r"^(main_(?:maze|corridor)_extended_v8|table_key(?:_arxiv)?|design_overview|arxiv_design_overview|n30_record_\w+|main_maze_arxiv|main_maze_arxiv_v2|main_maze_arxiv_v3|main_maze_arxiv_v4|main_maze_arxiv_v5|main_maze_arxiv_v6|main_maze_arxiv_v7|main_maze_arxiv_v8|main_corridor_arxiv_v8|assumption_checks_arxiv)$")   # includes Table IV and the assumption-check macro file
BLOCK_ROW = re.compile(r"^\\multicolumn\{\d+\}\{!\{\\vrule width\\vnbox\}(?:@\{\\hskip [\d.]+pt\})?l!\{\\vrule width\\vnbox\}\}")


def uniform_layout_checks(gen):
    """The uniform layout of Table X (continued) in every edition fragment (paper and
    report) except the key, the design tables and the catalog-only N30 record. Returns the number of failures."""
    nerr, nfrag, nbox = 0, 0, 0

    def bad(path, msg):
        nonlocal nerr
        nerr += 1
        print("FAIL layout", path, msg)

    for path in sorted(gen.glob("paper/*.tex")) + sorted(gen.glob("report/*.tex")):
        if UNIFORM_EXEMPT.match(path.stem):
            continue
        nfrag += 1
        rel = f"{path.parent.name}/{path.name}"
        txt = gen_text(path)
        lines = txt.splitlines()
        if r"\rotatebox" in txt:
            bad(rel, "rotated label")
        if re.search(r"\\ifdim\\vnsep\{\w+\}>", txt):
            bad(rel, "capped column separation (table narrower than the line)")
        m = re.search(r"\\renewcommand\{\\arraystretch\}\{([\d.]+)\}", txt)
        if m and not 0.93 <= float(m.group(1)) <= 1.4:
            bad(rel, f"arraystretch {m.group(1)} outside [0.93, 1.4]")
        # every layout or independent box: equalization block right before its separation
        for i, l in enumerate(lines):
            if l.startswith(r"\vnnewsep{"):
                tag = re.match(r"\\vnnewsep\{(\w+)\}", l).group(1)
                if not (i >= 1 and lines[i - 1].startswith(r"\ifdim\dimexpr\linewidth-") and
                        rf"\vnw{{{tag}eqv}}" in lines[i - 1]):
                    bad(rel, f"layout {tag}: no equalization of the value columns before its separation")
                if r"\dimexpr(\linewidth - " not in l:
                    bad(rel, f"layout {tag}: separation not computed from \\linewidth")
        tabs = [i for i, l in enumerate(lines) if r"\begin{tabular}{!{\vrule width\vnbox}" in l
                and r"\begin{tabular}{!{\vrule width\vnbox}p{" not in l]   # not the key box
        for n_, i in enumerate(tabs):
            l = lines[i]
            spec = re.search(r"\\begin\{tabular\}\{(.*?)\}\\Xhline", l).group(1)
            if re.search(r"w\{[lr]\}", spec):
                bad(rel, f"table {n_ + 1}: a column not centred")
            end = next(j for j in range(i, len(lines)) if lines[j].startswith(r"\Xhline{\vnbox}\end{tabular}"))
            body = lines[i + 1:end + 1]
            if not path.stem.startswith(("sec6_short", "sec6_combined")) or r"!{\hskip\arrayrulewidth}" in l:   # sec6_short, sec6_combined: line above
                if not body or not body[0].startswith(r"\multicolumn{") or r"@{\hskip 4pt}p{\dimexpr" not in body[0]:
                    bad(rel, f"table {n_ + 1}: no full-width settings row with the fixed padding")
            is_box = r"!{\hskip\arrayrulewidth}" in spec
            if is_box:
                nbox += 1
                if r"\makebox[\linewidth][c]{" not in l:
                    bad(rel, f"box {n_ + 1}: not centred at \\linewidth")
                hd = next((x for x in body if x.endswith(r"\vnsd\\")), "")
                if not hd.startswith(r"\multicolumn{1}{!{\vrule width\vnbox}w{c}{"):
                    bad(rel, f"box {n_ + 1}: Method header without its left box rule")
                sep = re.search(r"\\setlength\{\\tabcolsep\}\{\\vnsep\{(\w+)\}\}", l).group(1)
                first = re.match(r"!\{\\vrule width\\vnbox\}w\{c\}\{\\vnw\{(\w+)\}\}", spec)
                if not first:
                    bad(rel, f"box {n_ + 1}: first column")
                elif first.group(1).endswith("meth"):   # aligned: the grid of the main table
                    main = next((lines[j] for j in range(end, len(lines)) if r"\begin{tabular}{" in lines[j]
                                 and r"!{\hskip\arrayrulewidth}" not in lines[j]), "")
                    mspec = re.search(r"\\begin\{tabular\}\{(.*?)\}\\Xhline", main)
                    if not mspec or spec.replace(r"!{\hskip\arrayrulewidth}", "|") != mspec.group(1) or \
                            rf"\vnsep{{{sep}}}" not in main:
                        bad(rel, f"box {n_ + 1}: aligned box with another grid than its table")
                elif not re.search(rf"\\vnnewsep\{{{sep}\}}", txt) or not first.group(1).startswith(sep):
                    bad(rel, f"box {n_ + 1}: independent box without its own separation")
            else:
                nblk = sum(1 for x in body if BLOCK_ROW.match(x))
                # Section VI: one setting per table or one row per setting; the arXiv Maze table names each
                # setting in a block row (four variants per setting, as topfive_full)
                sec6 = path.stem.startswith("sec6_") and path.stem != "sec6_maze_arxiv"
                if not nblk and not sec6:
                    bad(rel, f"table {n_ + 1}: no block row naming the setting")
                if sec6 and nblk:
                    bad(rel, "Section VI table: the setting is named in the settings row or its own column")
    print(f"uniform layout: {nfrag} edition fragments, {nbox} baseline boxes checked")
    return nerr


SPLIT_TUPLE = re.compile(r"\(\d+\.\d+, \d+\.\d+,(?=\s*(?:\\\\|&|\}|$))", re.M)
STACKED_TUPLE = re.compile(r"\\begin\{tabular\}\{@\{\}c@\{\}\}[^&\n]*\$\\alpha_0=")


def tuple_checks(gen):
    """No setting tuple is split over two lines (no '(a, b,' ending a cell or line, no
    stacked 'alpha_0=' lines); every block row is a one-line l column. Editions and catalog."""
    nerr, nblk = 0, 0
    for path in sorted(gen.rglob("*.tex")):
        if "_superseded" in path.parts:
            continue
        txt = gen_text(path)
        for pat, what in ((SPLIT_TUPLE, "split tuple"), (STACKED_TUPLE, "stacked tuple")):
            m = pat.search(txt)
            if m:
                nerr += 1
                print("FAIL tuple", what, path.relative_to(gen), txt[max(0, m.start() - 40):m.end() + 20].replace("\n", " "))
        # the swept field is never printed in the block-row tuple (delta sweeps: no delta, K sweeps and the
        # detailed checks with rows or column groups at the setting's K and K = infinity: no K)
        swept_sym = (r"(\alpha_0,\alpha,\delta" if path.stem.startswith("delta_") else
                     ",K)$" if path.stem.startswith(("period_", "checks_", "detailed_events_")) else None)
        for l in txt.splitlines():
            if BLOCK_ROW.match(l):   # also with the fixed left padding @{\hskip 4pt}
                nblk += 1
                if swept_sym and swept_sym in l:
                    nerr += 1
                    print("FAIL tuple: swept field in block row", path.relative_to(gen), l[:120])
                if r"\\" not in l or r"\par" in l:
                    nerr += 1
                    print("FAIL tuple: block row", path.relative_to(gen), l[:80])
    print(f"setting tuples: {nblk} block rows checked, no split tuples" if not nerr else f"setting tuples: {nerr} failures")
    return nerr


def baseline_box_checks(gen):
    """Baseline boxes centred at full \\linewidth, Method column = left part of the main table plus
    every dropped value column, complete Method header box."""
    nerr, nbox = 0, 0
    for path in sorted(gen.glob("paper/*.tex")) + sorted(gen.glob("report/*.tex")):
        txt = gen_text(path)
        if r"\vncheckleft" not in txt:
            continue
        if r"\makebox[\linewidth][l]" in txt:
            nerr += 1
            print("FAIL box: left-aligned box", path.relative_to(gen))
        lines = txt.splitlines()
        for i, l in enumerate(lines):
            m = re.search(r"\\vncheckleft\{(.*?)\}\{\\vnw\{(\w+?)bname\}\}", l)
            if not m:
                continue
            nbox += 1
            gid = m.group(2)
            spec = re.search(r"\\begin\{tabular\}\{(.*?)\}\\Xhline", l).group(1)
            kept = set(re.findall(rf"w\{{c\}}\{{\\vnw\{{{gid}v(\d+)\}}\}}", spec))
            filled = set(re.findall(rf"\+ \\vnw\{{{gid}v(\d+)\}} \+ 2\\tabcolsep \+ \\arrayrulewidth", m.group(1)))
            main = next((x for x in lines[i + 1:] if re.match(r"\{\\setlength\{\\tabcolsep\}.*\\begin\{tabular\}", x)
                         and r"\vncheckleft" not in x), None)
            allv = set(re.findall(rf"w\{{c\}}\{{\\vnw\{{{gid}v(\d+)\}}\}}", main)) if main else None
            if allv is None or kept | filled != allv or kept & filled:
                nerr += 1
                print("FAIL box width", path.relative_to(gen), gid, sorted(kept), sorted(filled), allv and sorted(allv))
            if re.match(r"!\{\\vrule width\\vnbox\}w\{c\}\{\\dimexpr", spec):   # Method is the first column
                j = i + 1
                while j < len(lines) and not lines[j].startswith(r"\Xhline{0.6pt}"):
                    if r"\multicolumn{1}{w{c}{\dimexpr" in lines[j]:
                        nerr += 1
                        print("FAIL box: Method header without left rule", path.relative_to(gen), gid)
                        break
                    j += 1
    print(f"baseline boxes: {nbox} checked")
    return nerr


def compact_winner(cfgs):
    """cfgs: method -> cell of the four variants at one main-grid setting. Re-derives the compact rule of
    catalog_tables.py: most best means (bold rule; ties all best) among m_c, q_CBF, q_unsafe, T_mis (all lower is
    better), ties by printed m_c, q_unsafe, q_CBF, T_mis, then V-U, V-D, V-Cd, V-C."""
    order = ["V-U", "V-D", "V-Cd", "V-C"]
    mets = ["m_c", "q_unsafe", "q_CBF", "T_mis"]
    vals = {}
    for v in order:
        st = cfgs[v]["statistics"]
        row = {}
        for m in mets:
            s_ = st[m]
            ok = s_["applicable"] and s_["observed_n"] > 0 and s_["mean"] is not None
            row[m] = (Decimal(rnd(s_["mean"] * (100 if m in PCT else 1), DEC[m])) if ok else None,
                      ok and not s_["conditional"])
        vals[v] = row
    wins = {v: 0 for v in order}
    for m in mets:
        if all(vals[v][m][1] for v in order):
            best = min(vals[v][m][0] for v in order)
            for v in order:
                wins[v] += vals[v][m][0] == best
    key = lambda v: (-wins[v], tuple(vals[v][m][0] if vals[v][m][0] is not None else Decimal("Infinity")
                                     for m in mets), order.index(v))
    w = min(order, key=key)
    return w, wins[w]


def catalog_checks(cells, byfrag, views, gen):
    """Catalog structure: all-variant tables hold the four certified variants for every setting and column;
    sweep and compact tables cover the main grid at their alpha0; the compact variant and its (k) re-derived from
    the exports and printed in table order."""
    nerr = 0

    def bad(msg):
        nonlocal nerr
        nerr += 1
        if nerr <= 30:
            print("FAIL catalog", msg)

    by_cfg = {}
    for (s, vid), v in views.items():
        if s == "final" and v["definition"]["family"] == "main_grid_200":
            c = v["cells"][0]
            k = json.dumps({a: b for a, b in c["configuration"].items() if a != "method"}, sort_keys=True)
            by_cfg.setdefault(k, {})[c["configuration"]["method"]] = c
    VNAME = {"V-U": r"\textcolor{vu}{V-U}", "V-D": r"\textcolor{vd}{V-$\Delta$}", "V-Cd": r"\textcolor{vl}{V-L}",
             "V-C": r"\textcolor{vc}{V-C}"}
    ncomp = 0
    for (target, frag), es in byfrag.items():
        var = [e for e in es if not (e.get("row") or {}).get("baseline")]
        a0 = int(frag[-3:]) / 100 if frag.startswith(("sweep_", "compact_")) else None
        if a0 is not None:
            want = len(range(5, int(round(a0 * 100)), 5)) * 25
            cfgs = {json.dumps({k: x for k, x in e["configuration"].items() if k != "method"}, sort_keys=True)
                    for e in var}
            if len(cfgs) != want or any(abs(json.loads(c)["alpha0"] - a0) > 1e-9 for c in cfgs):
                bad(f"{target}/{frag}: {len(cfgs)} settings, want {want} at alpha0 = {a0}")
        if frag.startswith("compact_"):
            rows = []
            for e in var:
                if e["spec"] != "m_c":
                    continue
                k = json.dumps({a: b for a, b in e["configuration"].items() if a != "method"}, sort_keys=True)
                w, n = compact_winner(by_cfg[k])
                ncomp += 1
                if e["method"] != w or e["row"].get("wins") != n:
                    bad(f"{target}/{frag}: {k} variant {e['method']} ({e['row'].get('wins')}), want {w} ({n})")
                rows.append(f"{VNAME[w]} ({n})")
            txt = "\n".join(l for l in gen_text(gen / target / f"{frag}.tex").splitlines()
                            if not l.startswith((r"\vnfit", r"\vnnew", "%")))
            printed = re.findall(r"(\\textcolor\{v[udlc]\}\{V-[^}]*\} \(\d\))", txt)
            if printed != rows:
                bad(f"{target}/{frag}: printed variant column differs from the re-derived one "
                    f"({len(printed)} printed, {len(rows)} rows)")
            continue
        groups = defaultdict(set)
        for e in var:
            c = {k: x for k, x in e["configuration"].items() if k != "method"}
            if e["spec"].endswith("@inf"):
                c["K"] = "K-of-setting"
            groups[(e["table"], json.dumps(c, sort_keys=True), e["spec"].split("@")[0])].add(e["method"])
        for g, ms in groups.items():
            if ms != VARIANTS:
                bad(f"{target}/{frag}: {g[1]} {g[2]} has variants {sorted(ms)}")
                break
    print(f"catalog structure: {ncomp} compact rows re-derived")
    return nerr


def secondary_arxiv_checks(cells, views, gen):
    """Check arXiv-only coverage and percentage/SD semantics independently of the generator.

    check_file already re-derives every displayed mean, SD, scale and bold mark from the export. Here require
    all five settings/all four variants, five matched baselines, exact metric coverage, 200 temporal samples,
    the two-page layout, and a frequency key free of the historical count-cell notation.
    """
    nerr = 0
    def bad(msg):
        nonlocal nerr
        nerr += 1
        print("FAIL arxiv secondary", msg)
    specs = {"m_Y", "m_int", "n_certificates", "q_Y", "q_unsafe", "q_use", "q_nocert", "q_notrans", "tau_use",
             "h_min", "clearance_min", "Vc@temporal", "Jc@temporal", "Vp@temporal", "Jp@temporal"}
    settings = {"S1": (.30, .25, .20, 50), "S2": (.30, .25, .15, 50), "S3": (.30, .20, .20, 50),
                "S4": (.30, .25, .025, 100), "S5": (.30, .25, .10, 10)}
    expected_rows = {(s, v) for s in settings for v in VARIANTS} | {
        (None, v) for v in ("Unfiltered", "Zero margin", "Conservative", "ACI", "QT")}
    for scen in ("corridor", "maze"):
        frag = f"secondary_{scen}_arxiv"
        es = [e for e in cells if e["fragment"] == frag]
        rows = defaultdict(set)
        for e in es:
            cf, st = e["configuration"], e["row"].get("setting")
            rows[(st, e["method"])].add(e["spec"])
            if e["target"] != "paper" or cf["benchmark"] != scen or cf["L_mu"] != .0005:
                bad(f"{frag}: wrong edition/scenario/drift for {st}/{e['method']}")
            if st and tuple(cf[x] for x in ("alpha0", "alpha", "delta", "K")) != settings[st]:
                bad(f"{frag}: wrong setting {st}")
            if not st and e["method"] in ("ACI", "QT") and cf["alpha0"] != .30:
                bad(f"{frag}: unmatched baseline {e['method']}")
            if e["kind"] in ("event", "count"):
                bad(f"{frag}: count-form event cell {e['spec']}")
            if e["source"].startswith("temporal."):
                v = views[(e["src"], e["view_id"])]["cells"][0]
                metric = e["spec"].split("@")[0]
                q = "normalized_margin" if metric in ("Vc", "Jc") else "estimated_applied_risk"
                n = v["statistics"][q]["assessed_n"]
                if n != (0 if e["method"] == "Unfiltered" else 200):
                    bad(f"{frag}: {q} population {n} for {e['method']}")
            if e["spec"] in PCT and e["kind"] == "meanSD" and e["scale"] != 100:
                bad(f"{frag}: decision fraction without percent scale")
        if set(rows) != expected_rows or any(sp != (specs if st else specs - {"n_certificates", "q_nocert", "q_notrans"})
                                               for (st, _), sp in rows.items()) or len(es) != 360:
            bad(f"{frag}: coverage is {len(rows)} rows and {len(es)} cells, expected 25 rows / 360 applicable-column cells")
        txt = gen_text(gen / "paper" / f"{frag}.tex")
        if txt.count(r"\ntablebreak") != 1:
            bad(f"{frag}: expected two page blocks")
    key = gen_text(gen / "paper" / "table_key_arxiv.tex")
    if key.count("%KEYBOX") != 1 or "in percent" not in key or "confidence bound" not in key \
            or "sample SD" not in key or "$a/m$" in key or "counts:" in key:
        bad("arXiv key lacks percentage/uncertainty definitions or retains count-cell notation")
    print("arXiv secondary: 720 cells, 50 rows, percentages and temporal sample sizes checked")
    return nerr


def colour_switch_checks(gen):
    """Colour switch: every fragment of table_colour.COLOUR_TABLES carries exactly one switch (coloured variant,
    then the plain variant), no other fragment carries one, and the plain variant has no colour."""
    import sys
    sys.path.insert(0, str(HERE))
    import table_colour as tc
    n, nerr = 0, 0
    for target in ("paper", "report"):
        for path in sorted(gen.glob(f"{target}/*.tex")):
            txt = path.read_text()
            want = path.stem in tc.COLOUR_TABLES
            got = [txt.count(m) for m in (tc.MARK_DEF, tc.MARK_IF.replace("NAME", path.stem), tc.MARK_ELSE,
                                          tc.MARK_END)]
            plain = tc.plain_variant(txt)
            if got != ([1, 1, 1, 1] if want else [0, 0, 0, 0]) or "\\rgcell" in plain or "COLOUR-SWITCH" in plain \
                    or "vncolourlegend" in plain:
                print(f"FAIL colour switch {target}/{path.name}: markers {got}, wanted {'one each' if want else 'none'}")
                nerr += 1
            n += want
    print(f"colour switch: {n} fragments with coloured and plain variants")
    return nerr


def main():
    nerr = check_file(HERE / "tables_values.json.gz", GEN, edition=True)
    nerr += colour_switch_checks(GEN)
    cat = HERE / "catalog_values.json.gz"
    if cat.exists():
        nerr += check_file(cat, GEN / "catalog", edition=False)
    print("ERRORS:", nerr)
    raise SystemExit(1 if nerr else 0)


if __name__ == "__main__":
    main()
