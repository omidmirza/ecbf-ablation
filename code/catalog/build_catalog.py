#!/usr/bin/env python3
"""Build the N200 display catalog and the three sweep documents.

Inputs (read only):
  appendices/generated_n200/catalog/CATALOG_TABLES.json   catalog table fragments (tools/n200_tables, --catalog)
  appendices/generated_n200/{paper,report}/*.tex           edition table fragments (tools/n200_tables)
  tools/n200_tables/tables_values.json.gz                  data of the edition fragments (views per fragment)
  figures/n200/catalog/CATALOG_FIGURES.json                catalog figures (paper and report styles)
  figures/n200/{paper,report}/*.pdf                        figures as used in the editions
  sections/, appendices/, report/ (*.tex)                  where each edition table and figure is restored today

Writes catalog/generated/<doc>_body_<style>.tex (generated bodies) and compiles, from the paper source root, the
wrappers catalog/<doc>_<style>.tex with latexmk (pdflatex) into pdf/catalog/build/, then copies the PDFs to
pdf/catalog/:
  display_catalog_paper.pdf, display_catalog_report.pdf   every table and figure in that style
  sweeps_full_{paper,report}.pdf, sweeps_compact_{paper,report}.pdf, sweeps_a030_{paper,report}.pdf
Each document has a table of contents and a list of tables (the display catalog also a list of figures); each
caption gives the file, the data source and the exact restore line. The log check fails on any LaTeX error, undefined
reference, missing asset (N200-PENDING-ASSET), Overfull box, table-layout warning (of v10tables, the table package
of the paper source) or table with an unexpected page count, apart from the known layout results in KNOWN_LAYOUT.

Usage (from anywhere): python3 catalog/build_catalog.py [--only display_catalog_paper ...] [--no-compile]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

CAT = Path(__file__).resolve().parent
SRC = CAT.parent
GEN = SRC / "appendices" / "generated_n200"
FIG = SRC / "figures" / "n200"
OUT = SRC / "pdf" / "catalog"
BUILD = OUT / "build"
BODY = CAT / "generated"
STYLES = ("paper", "report")
# The script runs from catalog/ of the paper source tree (SRC is the tree root). REPO is the repository root for
# code/catalog/build_catalog.py; in the paper source tree it is the parent of the tree, so set ECBF_N200_NUMBERS.
REPO = CAT.parents[1]
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports (not included)
sys.path.insert(0, str(SRC / "tools" / "n200_tables"))
import table_colour  # noqa: E402  (a fragment with the colour switch holds two variants; count the plain one)
SCN = {"corridor": "Corridor", "maze": "Maze"}
WINNER_ASSET = "New3"   # \\winnerasset of the paper source: the asset label (New3) of setting S3
# Captions name the settings S1-S5 (asset labels New1-New5 in file names, T1-T5 in CATALOG_FIGURES.json)
SET = {"Default": "Default", "New1": "S1", "New2": "S2", "New3": "S3", "New4": "S4", "New5": "S5"}
TS = {"T1": "S1", "T2": "S2", "T3": "S3", "T4": "S4", "T5": "S5"}
SETT = {"Default": "(0.15, 0.10, 0.05, 50)", "S1": "(0.30, 0.25, 0.20, 50)", "S2": "(0.30, 0.25, 0.15, 50)",
        "S3": "(0.30, 0.20, 0.20, 50)", "S4": "(0.30, 0.25, 0.025, 100)", "S5": "(0.30, 0.25, 0.10, 10)"}
METS = {"m_c": r"$\widehat m_c/c_{\max}$", "q_CBF": r"$\widehat q_{\rm CBF}^+$", "q_unsafe": r"$\widehat q_{\rm unsafe}^+$",
        "T_mis": r"$\widehat T_{\rm mis}$"}


# ----------------------------------------------------------------------------------------------- text helpers
def tt(s: str) -> str:
    """Literal text in \\texttt with break points after / and _ (exact characters when printed)."""
    rep = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "_": r"\_\allowbreak{}", "$": r"\$", "#": r"\#",
           "%": r"\%", "&": r"\&", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "/": r"/\allowbreak{}"}
    return r"\texttt{" + "".join(rep.get(c, c) for c in s) + "}"


def rel(p: str) -> str:
    """A recorded input path as printed: relative paths as recorded, an absolute export path as N200_numbers/<rest>."""
    q = Path(p)
    if q.is_absolute() and NUM in q.parents:
        return "N200_numbers/" + q.relative_to(NUM).as_posix()
    return p


def a0tex(a0) -> str:
    return f"{a0:.2f}"


# ----------------------------------------------------------------------------------------------- usage in editions
def usage():
    """name -> list of 'file:line' where \\ntable{name} or \\nfigure{...name} appears in the edition sources."""
    use = defaultdict(list)
    files = [p for d in ("sections", "appendices", "report") for p in (SRC / d).rglob("*.tex")
             if "_superseded" not in p.parts and "generated_n200" not in p.parts and p.name != "n200_assets.tex"]
    for p in sorted(files):
        for i, line in enumerate(p.read_text(errors="replace").splitlines(), 1):
            if line.lstrip().startswith("%"):
                continue
            # a size group around the table, {\footnotesize\ntable{name}}, belongs to its restore line
            for m in re.finditer(r"(\{\\(?:footnotesize|small|scriptsize))?\\ntable\{([^}]*)\}(?(1)\})", line):
                use[("table", m.group(2))].append((f"{p.relative_to(SRC)}:{i}", m.group(0)))
            for m in re.finditer(r"\\nfigure(\[[^\]]*\])?\{([^}]*)\}", line):
                name = m.group(2).replace(r"\winnerasset", WINNER_ASSET)
                use[("figure", name)].append((f"{p.relative_to(SRC)}:{i}", m.group(0)))
    return use


# ----------------------------------------------------------------------------------------------- edition tables
ED_SKIP = {"assumption_checks_arxiv"}   # edition fragments without a table (macros of the assumption-check totals)
ED_TITLES = [   # first match wins: the specific names come before the general patterns
    (r"^table_key$", "Key to all tables", "key"),
    (r"^table_key_arxiv$", "Key to all tables (arXiv)", "key"),
    (r"^sec6_default$", "Section VI: default setting, four variants and baselines (Table I of the extended version)", "sec6"),
    (r"^sec6_short$", r"Section VI: chosen variant at S1--S5, baselines at $\alpha_0=0.30$ (reference of the extended version's Table I)", "sec6"),
    (r"^sec6_maze_arxiv$", r"Section VI: Maze, S1--S5 with the four variants, baselines at $\alpha_0=0.30$ (arXiv Table I)", "sec6"),
    (r"^main_maze_arxiv_v8$", r"Section VI: Maze, S1--S5, four variants, baselines at $\alpha_0=0.30$ (arXiv Table I)", "sec6"),
    (r"^main_corridor_arxiv_v8$", r"Corridor, S1--S5, four variants, baselines at $\alpha_0=0.30$ (arXiv Table IV)", "sec6"),
    (r"^main_(corridor|maze)_extended_v8$", r"Section VI: {0}, S1--S5, four variants, baselines at $\alpha_0=0.30$ (extended version)", "sec6"),
    (r"^secondary_(corridor|maze)_arxiv$", r"Secondary and temporal metrics, S1--S5, all variants and baselines, {0} (arXiv)", "default_topfive"),
    (r"^stress_representative_arxiv_percent$", r"Bound factor $\phi$ and distribution-drift bound $L_\mu$, S3, V-L, both scenarios (arXiv Table V)", "parameter"),
    (r"^stress_(corridor|maze)_arxiv_percent$", r"Bound factor $\phi$ and distribution-drift bound $L_\mu$, S1--S5, {0} (arXiv, percentages)", "parameter"),
    (r"^stress_(corridor|maze)_cases$", r"Bound factor $\phi$ and distribution-drift bound $L_\mu$, S1--S5, {0} (report)", "parameter"),
    (r"^topfive_full_(corridor|maze)_arxiv_percent$", r"S1--S5, four variants, all metrics with $\widehat q_{{\rm dep}}$, $\widehat q_{{\rm cert}}$, {0} (arXiv, percentages)", "default_topfive"),
    (r"^topfive_full_(corridor|maze)_onepage$", r"S1--S5, four variants, all metrics with $\widehat q_{{\rm dep}}$, $\widehat q_{{\rm cert}}$, {0} (one page)", "default_topfive"),
    (r"^variation_(corridor|maze)_arxiv$", r"Temporal variation $10^3\widehat V$, $100\widehat J$, four variants, {0} (arXiv)", "variation"),
    (r"^sec6_short_star$", r"Section VI: chosen variant at S1--S5 and row S$^\ast$ ($\alpha=0.225$, not on the main grid), extended version", "sec6"),
    (r"^sec6_combined$", r"Section VI: baselines; default setting, S1--S5 and S$^\ast$, each with one variant (extended version)", "sec6"),
    (r"^star_maze$", r"S$^\ast$ ($\alpha=0.225$), four variants, Maze", "default_topfive"),
    (r"^checks_(\w+)_arxiv$", r"Checks, S1--S5, {0} (arXiv)", "checks"),
    (r"^stress_(\w+)_arxiv$", r"Bound factor $\phi$ and distribution-drift bound $L_\mu$, S1--S5, {0} (arXiv)", "parameter"),
    (r"^design_overview$", "Design of the comparisons (extended version and report)", "design"),
    (r"^arxiv_design_overview$", "Design of the comparisons (arXiv)", "design"),
    (r"^default_(\w+)$", "Default setting, four variants, all metrics, {0}", "default_topfive"),
    (r"^topfive_full_(\w+)$", r"S1--S5, four variants, all metrics with $\widehat q_{{\rm dep}}$, $\widehat q_{{\rm cert}}$, {0}", "default_topfive"),
    (r"^delta_(\w+)$", r"Parameter $\delta$, chosen variant, {0}", "parameter"),
    (r"^period_(\w+)$", r"Test period $K$, chosen variant, {0}", "parameter"),
    (r"^b_(\w+)$", r"Budget normalizer $b$, chosen variant, {0}", "parameter"),
    (r"^phi_(\w+)$", r"Bound factor $\phi\times L_\mu$, chosen variant, {0}", "parameter"),
    (r"^detailed_events_(\w+)$", r"Detailed events at $K$ and $K=\infty$, chosen variant, {0}", "checks"),
    (r"^checks_(\w+)_L0005_nodef$", r"Checks, chosen variant, S1--S5, $L_\mu=0.0005$, {0} (arXiv)", "checks"),
    (r"^checks_(\w+)_L0005$", r"Checks, chosen variant, $L_\mu=0.0005$, {0}", "checks"),
    (r"^checks_(\w+)_L005$", r"Checks, chosen variant, $L_\mu=0.005$, {0}", "checks"),
    (r"^variation_(\w+)$", r"Temporal variation $10^3\widehat V$, $100\widehat J$, four variants, {0}", "variation"),
    (r"^n30_record_(\w+)$", "Settings selected on the first 30 trajectories, {0}", "n30"),
]
ED_GROUPS = [("key", "Key"), ("sec6", "Section VI tables"), ("design", "Design of the comparisons"),
             ("default_topfive", "Default and selected settings"), ("parameter", "Parameter tables (chosen variant)"),
             ("checks", "Checks and detailed events (chosen variant)"), ("variation", "Temporal variation"),
             ("n30", "Settings selected on the first 30 trajectories")]


def ed_title(name):
    for pat, title, grp in ED_TITLES:
        m = re.match(pat, name)
        if m:
            return title.format(*(SCN.get(g, g) for g in m.groups())), grp
    raise KeyError(name)


def edition_data():
    """(target, fragment) -> data summary from tables_values.json.gz."""
    with gzip.open(SRC / "tools/n200_tables/tables_values.json.gz", "rt") as fh:
        V = json.load(fh)
    d = defaultdict(lambda: dict(srcs=set(), fams=set(), views=set(), star=False, design=False))
    for e in V["cells"]:
        x = d[(e["target"], e["fragment"])]
        if e.get("view_id"):
            x["srcs"].add(e["src"]), x["fams"].add(e["family"]), x["views"].add(e["view_id"])
        if e["kind"] == "starred":
            x["star"] = True
        if e["kind"] == "design_count":
            x["design"] = True
    return d, {k: rel(v) for k, v in V["inputs"].items()}


ABBR = {"final": "[F]", "new": "[A]", "final_temporal": "[T]"}   # defined in \\catreading (catalog_macros.tex)


def abbr(path):
    return "[A]" if "new_anchors" in path else "[T]" if "temporal" in path else "[F]"


def data_text(srcs, fams, views, inputs, where):
    files = ", ".join(ABBR[s] for s in sorted(srcs))
    ids = (": " + ", ".join(tt(v[:12]) for v in sorted(views))) if 0 < len(views) <= 4 else \
        f" (ids in {where})"
    return f"{files}; families {', '.join(tt(f) for f in sorted(fams))}; {len(views)} views{ids}"


# ----------------------------------------------------------------------------------------------- captions
BREAKS = {}   # label of a table block -> expected pages (number of \\ntablebreak + 1 in its fragment)
# Known layout results of the current table fragments, reported but not counted as failures: two column separations
# below 1 pt in the paper display catalog, and the report-style topfive_full tables, which take one page more than
# their breaks give.
KNOWN_LAYOUT = {"layout_warn": {"display_catalog_paper": 2},
                "pages": {("display_catalog_report", "-ntable-topfive-full-corridor-"): 3,
                          ("display_catalog_report", "-ntable-topfive-full-maze-"): 3}}


def table_block(title, path, data, restore, used):
    where = ("; used in " + ", ".join(tt(u[0]) for u in used)) if used else "; not used in the current editions"
    style = "report" if "/report/" in path else "paper"
    BREAKS[(style, re.sub(r"[^A-Za-z0-9]", "-", restore))] = table_colour.plain_variant((SRC / path).read_text()).count(r"\ntablebreak") + 1
    long = (rf"{title}. \emph{{File:}} {tt(path)}. \emph{{Data:}} {data}. \emph{{Restore:}} {tt(restore)}{where}.")
    lab = re.sub(r"[^A-Za-z0-9]", "-", restore)
    return (r"\begin{simbreaktable}\caption[" + title + "]{" + long + "}\\label{cs:" + lab + "}\n" + restore +
            r"\par\label{ce:" + lab + "}\n\\end{simbreaktable}\n" r"\clearpage" + "\n")


def figure_block(title, path, data, show, restore, used):
    where = ("; used in " + ", ".join(tt(u[0]) for u in used)) if used else ""
    long = rf"{title}. \emph{{File:}} {tt(path)}. \emph{{Data:}} {data}. \emph{{Restore:}} {tt(restore)}{where}."
    return (r"\catfigure{" + show + "}{" + title + "}{" + long + "}\n")


# ----------------------------------------------------------------------------------------------- sections
def tables_catalog(T, style, kinds, sec_level="section"):
    """Catalog table fragments of one style, filtered by type; organised by type, scenario and parameter."""
    out = []
    order = {"default_topfive": 0, "parameter": 1, "checks": 2, "events": 3, "sweep_compact": 4, "sweep_full": 5}
    rows = [t for t in T["tables"] if t["target"] == style and t["type"] in kinds]
    rows.sort(key=lambda t: (order[t["type"]], t["parameters"].get("scenario", ""), t["parameters"].get("alpha0", 0),
                             t["name"]))
    cur = None
    for t in rows:
        sub = None
        if t["type"].startswith("sweep"):
            sub = f"{SCN[t['parameters']['scenario']]}, $\\alpha_0={a0tex(t['parameters']['alpha0'])}$"
        elif t["type"] == "parameter":
            sub = SCN[t["parameters"]["scenario"]]
        if sub and sub != cur:
            out.append(rf"\subsection{{{sub}}}")
            cur = sub
        d = t["data"]
        files = ", ".join(abbr(f) for f in d["files"])
        ids = (": " + ", ".join(tt(v[:12]) for v in d["view_ids"])) if 0 < d["n_views"] <= 4 else " (ids in [C])"
        data = f"{files}; families {', '.join(tt(f) for f in d['families'])}; {d['n_views']} views{ids}"
        out.append(table_block(t["title"], t["file"], data, t["restore"], None))
    return "\n".join(out)


def edition_tables(style, use, ED, inputs):
    names = sorted(p.stem for p in (GEN / style).glob("*.tex") if p.stem not in ED_SKIP)
    out = []
    for gid, gtitle in ED_GROUPS:
        sel = [n for n in names if ed_title(n)[1] == gid]
        if not sel:
            continue
        out.append(rf"\subsection{{{gtitle}}}" + (r"\label{app:gallery-n30}" if gid == "n30" else ""))
        for n in sel:
            title, _ = ed_title(n)
            x = ED[(style, n)]
            path = f"appendices/generated_n200/{style}/{n}.tex"
            if x["design"]:
                data = "configuration counts of [F] and [A] (recomputed by check\\_tables.py)"
            elif x["views"]:
                data = data_text(x["srcs"], x["fams"], x["views"], inputs, "[V]")
                if x["star"]:
                    data += (f"; row S$^\\ast$: {tt('tools/n200_tables/starred_setting.json')} (values and source lines)")
            else:
                data = "notation only"
            used = use.get(("table", n))
            restore = used[0][1] if used and used[0][1].startswith("{") else rf"\ntable{{{n}}}"
            out.append(table_block(title, path, data, restore, used))
    return "\n".join(out)


def setting(e):
    """Setting of a catalog figure as printed in captions (S1-S5 for the recorded T1-T5)."""
    return TS.get(e.get("setting"), e.get("setting"))


def fig_title(e):
    k = e["kind"]
    sc = SCN.get(e.get("scenario") or "", "")
    st = setting(e)
    if k.startswith("distribution"):
        return f"Per-trajectory distributions, {sc}, {st} {SETT.get(st, '')}"
    if k.startswith("sweep plot (compact"):
        return f"Sweep plot, compact, {st} {SETT.get(st, '')}, Corridor and Maze"
    if k.startswith("sweep plot"):
        return f"Sweep plot, {sc}, {st} {SETT.get(st, '')}"
    if k.startswith("trajectory map"):
        form = k[k.index("(") + 1:-1]
        return f"Trajectory map ({form}), {sc}, {st} {SETT.get(st, '')}"
    if k.startswith("heatmap"):
        ds = ", ".join(f"{x:g}" for x in e["deltas"])
        return (f"Heatmap {METS.get(e['metric'], e['metric'])}, {sc}, $\\alpha_0={a0tex(e['alpha0'])}$, "
                f"$\\delta\\in\\{{{ds}\\}}$")
    raise KeyError(k)


def fig_data(e):
    s = re.sub(r"^.*/maps_top5/", "maps_top5/", e["source_script"])   # recorded path of the map scripts, as maps_top5/<script>
    s = re.sub(r"^(?:.*/)?(?:source|code)/(?=tools/)", "", s)          # scripts of tools/ by their path in the tree
    k = e["kind"]
    if k.startswith("distribution"):
        d = (f"{tt('tools/n200_catalog/setting_assessments.json')} (saved per-trajectory assessments, hash-checked "
             f"against {tt('N200_data/N200_INDEX.sqlite')})")
    elif k.startswith("sweep"):
        d = (f"{tt('tools/n200_catalog/sweep_cells.json')} from [F] and [A]")
    elif k.startswith("heatmap"):
        d = "[F] (main\\_grid\\_200)"
    else:
        d = (f"{tt('N200_numbers/final/common/population.jsonl')} and the saved "
             f"trajectory payloads (seed 970001 highlighted)")
    return d + f" via {tt(s)}"


def fig_show(path_from_kind, e_width, style):
    """\\nfigure line used in the catalog: natural size (the figure is drawn at its final width)."""
    return rf"\nfigure[{e_width:.2f}pt]{{{path_from_kind}}}"


def catalog_figures(F, style, use):
    figs = [e for e in F["figures"] if e["edition"] == style]
    out = []
    order = [("distribution", "Per-trajectory distributions"), ("sweep plot", "Sweep plots"),
             ("trajectory map", "Trajectory maps"), ("heatmap (alpha x K, five", "Heatmaps, five $\\delta$ slices"),
             ("heatmap (alpha x K, two", "Heatmaps, two $\\delta$ slices (arXiv style)")]
    mo = {"m_c": 0, "q_CBF": 1, "q_unsafe": 2, "T_mis": 3}
    for pre, title in order:
        sel = [e for e in figs if e["kind"].startswith(pre)]
        if not sel:
            continue
        out.append(r"\needspace{0.45\textheight}" + rf"\subsection{{{title}}}")
        sel.sort(key=lambda e: (e.get("scenario") or "", mo.get(e.get("metric"), 9), e.get("alpha0") or 0,
                                tuple(e.get("deltas") or ()), list(SETT).index(setting(e)) if setting(e) in SETT else 9,
                                e["kind"]))
        cur = None
        for e in sel:
            if pre.startswith("heatmap"):
                sub = f"{SCN[e['scenario']]}, {METS[e['metric']]}"
            else:
                sub = SCN.get(e.get("scenario") or "", "Corridor and Maze")
            if sub != cur:
                out.append(r"\needspace{0.45\textheight}" + rf"\catsubsub{{{sub}}}")
                cur = sub
            rest = e["file"][len(style) + 1:]
            arg = f"../catalog/{style}/{rest}"
            width = r"[\columnwidth]" if e["width_texpt"] < 300 else ""
            restore = rf"\nfigure{width}{{{arg}}}"
            out.append(figure_block(fig_title(e), f"figures/n200/catalog/{e['file']}", fig_data(e),
                                    fig_show(arg, e["width_texpt"], style), restore, None))
    return "\n".join(out)


EDFIG = [
    (r"^map_(corridor|maze)_(\w+?)_(text|column|compact|pair)(_zoom|_fit|_main|_main_abcd)?\.pdf$", "map"),
    (r"^distributions_(corridor|maze)_(Default|New\d)(_sep)?\.pdf$", "dist"),
    (r"^distributions_(Default|New\d)_half\.pdf$", "disthalf"),
    (r"^distributions_default_(corridor|maze)_row(\d)\.pdf$", "defdist"),
    (r"^default_margin_distributions\.pdf$", "defmargin"),
    (r"^decisions_default_(corridor|maze)_K50\.pdf$", "decisions"),
    (r"^theorem_checks_(corridor|maze)_(New\d|Default)\.pdf$", "theorem"),
    (r"^sweepplot_T3_compact()\.pdf$", "sweepc"),
    (r"^sweepplot_T3_compact(_twin)\.pdf$", "sweepc"),
    (r"^sweepbox_T3_(compact|full)\.pdf$", "sweepc"),
    (r"^sweepplot_(corridor|maze)_(Default|New\d)\.pdf$", "sweep"),
    (r"^sweepbox_(corridor|maze)_(Default|New\d)\.pdf$", "sweepb"),
    (r"^heatmap_(corridor|maze)_(m_c|q_CBF|q_unsafe|T_mis)_a0_(0p\d+)_d_(all|[\dp.-]+)\.pdf$", "heatmap"),
]
EDFIG_GROUPS = [("map", "Trajectory maps"), ("dist", "Per-trajectory distributions"), ("disthalf", None),
                ("defdist", "Default setting: distributions, decisions and margins"),
                ("defmargin", None), ("decisions", None), ("theorem", "Theorem checks"),
                ("sweepc", "Sweep plots"), ("sweep", None), ("sweepb", None), ("heatmap", "Heatmaps")]


def edfig_info(name, style):
    for pat, kind in EDFIG:
        m = re.match(pat, name)
        if m:
            g = m.groups()
            if kind == "map":
                st = SET.get(g[1], g[1])
                form = {"_main_abcd": "main, panels (a)-(d)"}.get(g[3], (g[3] or "")[1:])
                t = f"Trajectory map ({g[2]}{', ' + form if form else ''}), {SCN[g[0]]}, {st} {SETT.get(st, '')}"
                script = "maps_top5/build_compact_main.py" if g[3] == "_main_abcd" else "maps_top5/build_maps_top5.py"
                d = (f"{tt('N200_numbers/final/common/population.jsonl')} and saved "
                     f"trajectory payloads via "
                     + (tt(script)
                        if style == "paper" else tt("tools/n200_maps/build_maps_report" + ("_pair" if g[2] == "pair" else "") + ".py")))
            elif kind == "dist":
                st = SET[g[1]]
                t = f"Per-trajectory distributions, {SCN[g[0]]}, {st} {SETT[st]}" + (" (separate form)" if g[2] else "")
                d = f"{tt('tools/n200_figures/setting_assessments.json')} via {tt('tools/n200_figures/build_distribution_figures.py')}"
            elif kind == "disthalf":
                st = SET[g[0]]
                t = f"Per-trajectory distributions, Corridor and Maze, {st} {SETT[st]}"
                d = f"{tt('tools/n200_figures/setting_assessments.json')} via {tt('tools/n200_figures/build_distribution_figures.py')} (form half)"
            elif kind in ("defdist", "defmargin", "decisions"):
                t = {"defdist": lambda: f"Default setting, per-trajectory distributions (row {g[1]}), {SCN[g[0]]}",
                     "defmargin": lambda: "Default setting, per-trajectory margins $m_{c,k}/c_{\\max}$, Corridor and Maze",
                     "decisions": lambda: f"Default setting, decision fractions, {SCN[g[0]]}, $K=50$"}[kind]()
                d = f"{tt('tools/n200_figures/default_assessments.json')} via {tt('tools/n200_figures/build_default_figures.py')}"
            elif kind == "theorem":
                st = SET[g[1]]
                t = f"Theorem checks, {SCN[g[0]]}, {st} {SETT[st]} " + ("V-C" if st == "Default" else "V-L")
                d = (("[F]" if st == "Default" else "[A]") + " (detailed\\_200) via "
                     f"{tt('tools/n200_figures/build_theorem_checks_box.py')}")
            elif kind in ("sweep", "sweepc"):
                if kind == "sweep":
                    st = SET[g[1]]
                    t = f"Sweep plot, {SCN[g[0]]}, {st} {SETT[st]}"
                    d = f"{tt('tools/n200_figures/sweep_cells.json')} via {tt('tools/n200_figures/build_sweeps.py')}"
                elif name.startswith("sweepbox_"):
                    t = f"Sweep box plot, {g[0]}, S3 {SETT['S3']}, Corridor and Maze"
                    d = f"{tt('tools/n200_figures/sweepbox_T3_values.json')} via {tt('tools/n200_figures/build_sweeps_box.py')}"
                elif g[0]:
                    t = f"Sweep plot, compact, two axes, S3 {SETT['S3']}, Corridor and Maze"
                    d = f"{tt('tools/n200_figures/sweep_cells.json')} via {tt('tools/n200_figures/build_sweeps_twin.py')}"
                else:
                    t = f"Sweep plot, compact, S3 {SETT['S3']}, Corridor and Maze"
                    d = f"{tt('tools/n200_figures/sweep_cells.json')} via {tt('tools/n200_figures/build_sweeps.py')}"
            elif kind == "sweepb":
                st = SET[g[1]]
                t = f"Sweep box plot, {SCN[g[0]]}, {st} {SETT[st]}"
                d = f"{tt('tools/n200_figures/sweepbox_ext_values.json')} via {tt('tools/n200_figures/build_sweeps_box_ext.py')}"
            else:
                a0 = float(g[2].replace("p", "."))
                ds = "all five" if g[3] == "all" else g[3].replace("p", ".").replace("-", ", ")
                t = f"Heatmap {METS[g[1]]}, {SCN[g[0]]}, $\\alpha_0={a0tex(a0)}$, $\\delta$: {ds}"
                d = (f"[F] (main\\_grid\\_200) via "
                     f"{tt('tools/n200_heatmaps/build_heatmaps.py')}")
            return kind, t, d
    raise KeyError(name)


def edition_figures(style, use):
    files = sorted(p.name for p in (FIG / style).glob("*.pdf"))
    out = []
    for gk, gtitle in EDFIG_GROUPS:
        sel = [f for f in files if edfig_info(f, style)[0] == gk]
        if not sel:
            continue
        if gtitle:
            out.append(rf"\subsection{{{gtitle}}}")
        for f in sel:
            _, title, data = edfig_info(f, style)
            w = pdf_width(FIG / style / f)
            used = use.get(("figure", f))
            if used:
                restore = used[0][1].replace(r"\winnerasset", WINNER_ASSET)
            else:
                restore = (r"\nfigure[\columnwidth]{" if w < 300 else r"\nfigure{") + f + "}"
            out.append(figure_block(title, f"figures/n200/{style}/{f}", data, fig_show(f, w, style), restore, used))
    return "\n".join(out)


def pdf_width(p):
    """Page width in TeX pt (pdfinfo reports PostScript points)."""
    out = subprocess.run(["pdfinfo", str(p)], capture_output=True, text=True).stdout
    return float(re.search(r"Page size:\s+([\d.]+) x", out).group(1)) * 72.27 / 72


# ----------------------------------------------------------------------------------------------- documents
def bodies(T, F, use, ED, inputs):
    for style in STYLES:
        sty = "paper (ieeeconf, one column, \\textwidth 505.89\\,pt, 10\\,pt)" if style == "paper" else \
            "report (article, 11\\,pt Times, \\textwidth 505.89\\,pt)"
        # display catalog
        parts = [r"\section{Tables in the editions}",
                 r"Fragments of " + tt(f"appendices/generated_n200/{style}/") + " exactly as the editions input them.",
                 edition_tables(style, use, ED, inputs),
                 r"\section{Catalog tables: default and selected settings}",
                 tables_catalog(T, style, {"default_topfive"}),
                 r"\section{Catalog tables: parameters, all four variants}",
                 tables_catalog(T, style, {"parameter"}),
                 r"\section{Catalog tables: checks and events, all four variants}",
                 tables_catalog(T, style, {"checks", "events"}),
                 r"\section{Sweep tables, compact (best variant)}",
                 tables_catalog(T, style, {"sweep_compact"}),
                 r"\section{Sweep tables, full}",
                 tables_catalog(T, style, {"sweep_full"}),
                 r"\section{Figures in the editions}",
                 "Figures of " + tt(f"figures/n200/{style}/") + ".",
                 edition_figures(style, use),
                 r"\section{Catalog figures}",
                 "Figures of " + tt(f"figures/n200/catalog/{style}/") + " (" + tt("CATALOG_FIGURES.json") + ").",
                 catalog_figures(F, style, use)]
        (BODY / f"display_catalog_body_{style}.tex").write_text(
            "% Generated by catalog/build_catalog.py; do not edit by hand.\n" + "\n".join(parts) + "\n")
        for doc, kinds, filt in (("sweeps_full", {"sweep_full"}, None), ("sweeps_compact", {"sweep_compact"}, None),
                                 ("sweeps_a030", {"sweep_compact", "sweep_full"}, 0.30)):
            TT = dict(T, tables=[t for t in T["tables"] if filt is None or abs(t["parameters"].get("alpha0", 0) - filt) < 1e-9])
            if doc == "sweeps_a030":
                body = "\n".join([r"\section{Compact sweep tables at $\alpha_0=0.30$}",
                                  tables_catalog(TT, style, {"sweep_compact"}),
                                  r"\section{Full sweep tables at $\alpha_0=0.30$}",
                                  tables_catalog(TT, style, {"sweep_full"})])
            else:
                body = "\n".join([rf"\section{{{'Full' if doc == 'sweeps_full' else 'Compact'} sweep tables}}",
                                  tables_catalog(TT, style, kinds)])
            (BODY / f"{doc}_body_{style}.tex").write_text(
                "% Generated by catalog/build_catalog.py; do not edit by hand.\n" + body + "\n")


WRAP = r"""%% Wrapper generated by catalog/build_catalog.py: %(doc)s in the %(style)s style.
%% Compile from the paper source root, as build_catalog.py does.
\def\catstyle{%(style)s}
\input{catalog/preamble_%(style)s}
\hypersetup{pdftitle={Anytime-Valid Certification of Robustness Margins for Safe Control via E-Values: %(title)s},
  pdfauthor={Panagiotis S. Trakas, Omid Mirzaeedodangeh, Lars Lindemann}}
\begin{document}
\input{catalog/%(doc)s}
\end{document}
"""


TITLES = {"display_catalog": "N200 Display Catalog", "sweeps_full": "Full Sweep Tables",
          "sweeps_compact": "Compact Sweep Tables", "sweeps_a030": "Sweep Tables at alpha0 = 0.30"}


def wrappers():
    out = []
    for doc in ("display_catalog", "sweeps_full", "sweeps_compact", "sweeps_a030"):
        for style in STYLES:
            p = CAT / f"{doc}_{style}.tex"
            p.write_text(WRAP % dict(doc=doc, style=style, title=f"{TITLES[doc]} ({style} style)"))
            out.append(p)
    return out


# ----------------------------------------------------------------------------------------------- compile and check
def compile_doc(p):
    BUILD.mkdir(parents=True, exist_ok=True)
    cmd = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "-f-",
           f"-outdir={BUILD.relative_to(SRC)}", str(p.relative_to(SRC))]
    r = subprocess.run(cmd, cwd=SRC, capture_output=True, text=True)
    log = (BUILD / (p.stem + ".log")).read_text(errors="replace")
    res = dict(doc=p.stem, ok=r.returncode == 0)
    res["errors"] = len(re.findall(r"^! ", log, re.M))
    res["undefined"] = len(re.findall(r"(Reference|Citation) `[^']*' on page \d+ undefined", log))
    res["pending"] = len(re.findall(r"N200-PENDING-ASSET", log))
    res["overfull"] = len(re.findall(r"^Overfull \\[hv]box", log, re.M))
    res["layout_warn"] = len(re.findall(r"Package v10tables Warning", log))   # table-layout warnings
    res["width"] = len(re.findall(r"N200-ASSET-WIDTH", log))
    # every table fragment spans exactly (number of \\ntablebreak + 1) pages from its caption to its end
    aux = (BUILD / (p.stem + ".aux")).read_text(errors="replace")
    pg = {(m.group(1), m.group(2)): int(m.group(3))
          for m in re.finditer(r"\\newlabel\{(cs|ce):([^}]*)\}\{\{[^}]*\}\{(\d+)\}", aux)}
    res["page_mismatch"] = []
    for (k, lab), a in pg.items():
        if k != "cs":
            continue
        want = BREAKS.get((p.stem.rsplit("_", 1)[1], lab))
        got = pg[("ce", lab)] - a + 1
        if want is not None and got != want:
            res["page_mismatch"].append((lab, got, want))
    m = re.search(r"Output written on .*?\((\d+) pages", log)
    res["pages"] = int(m.group(1)) if m else None
    res["tables"] = len(re.findall(r"\\contentsline \{table\}", (BUILD / (p.stem + ".lot")).read_text())) \
        if (BUILD / (p.stem + ".lot")).exists() else 0
    res["figures"] = len(re.findall(r"\\contentsline \{figure\}", (BUILD / (p.stem + ".lof")).read_text())) \
        if (BUILD / (p.stem + ".lof")).exists() else 0
    if res["ok"]:
        shutil.copy2(BUILD / (p.stem + ".pdf"), OUT / (p.stem + ".pdf"))
    else:
        print(log[-2500:])
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", help="wrapper stems to compile (default: all)")
    ap.add_argument("--no-compile", action="store_true")
    a = ap.parse_args()
    T = json.loads((GEN / "catalog/CATALOG_TABLES.json").read_text())
    F = json.loads((FIG / "catalog/CATALOG_FIGURES.json").read_text())
    for e in F["figures"]:   # every catalog figure must exist with the recorded hash
        p = FIG / "catalog" / e["file"]
        assert hashlib.sha256(p.read_bytes()).hexdigest() == e["sha256"], p
    for t in T["tables"]:
        assert hashlib.sha256((SRC / t["file"]).read_bytes()).hexdigest() == t["sha256"], t["file"]
    BODY.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    ED, inputs = edition_data()
    bodies(T, F, usage(), ED, inputs)
    docs = wrappers()
    if a.no_compile:
        return
    if a.only:
        docs = [d for d in docs if d.stem in a.only]
    res = [compile_doc(d) for d in docs]
    bad = False
    for r in res:
        print(f"{r['doc']}: ok={r['ok']} pages={r['pages']} tables={r['tables']} figures={r['figures']} "
              f"errors={r['errors']} undefined={r['undefined']} pending={r['pending']} overfull={r['overfull']} "
              f"layout-warnings={r['layout_warn']} asset-width-notes={r['width']}")
        known = [x for x in r["page_mismatch"] if KNOWN_LAYOUT["pages"].get((r["doc"], x[0])) == x[1]]
        if r["page_mismatch"]:
            print("  page mismatch (label, pages, expected):", r["page_mismatch"][:20],
                  f"({len(known)} known)" if known else "")
        bad |= (not r["ok"] or r["errors"] or r["undefined"] or r["pending"] or r["overfull"]
                or r["layout_warn"] > KNOWN_LAYOUT["layout_warn"].get(r["doc"], 0)
                or len(known) < len(r["page_mismatch"]))
    (OUT / "BUILD_RESULTS.json").write_text(json.dumps(res, indent=1))
    print("CATALOG BUILD:", "FAILED" if bad else "PASSED")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
