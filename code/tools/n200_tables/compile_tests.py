#!/usr/bin/env python3
"""Compile every generated N200 table fragment at its target width and check the layout.

Paper test: the real paper preamble (ieeeconf 10 pt, paper_helper, release_helper with n200_assets);
Section VI tables as table* floats in two-column mode, Appendix G fragments in one-column mode inside
simbreaktable (as in g_*.tex), the arXiv checks pair inside one simdisplay (as in arxiv_gallery.tex).
Report test: the ablation-report preamble (article 11 pt, geometry letterpaper margin 0.75in = 505.89 pt text width
times, report_helper).
Each fragment starts on a fresh page after a two-line stand-in caption; the check requires
pages(fragment) = number of \\ntablebreak + 1, zero Overfull/Underfull boxes, no "Float too large",
no v10tables warning, every column separation >= 1 pt, and no N200-FRAGMENT-CAPTION.
Every page block fits its page (PAGEFILL <= page goal: no table runs off the page);
the fill of every page block is printed; the arXiv pair is checks_{corridor,maze}_arxiv in one simdisplay (as in
arxiv_gallery.tex, at \\normalsize); the arXiv stress tables are breakable (simbreaktable, two pages each).
Compiles from the source root (cwd) exactly like build.py; outputs go to tools/n200_tables/test/out.
Colour switch: the default run typesets the plain variants (as every edition by default); --colour defines
\\releasecolortables (as editions/*_colored.tex; test_*_colour.tex) and applies the same checks to the coloured
variants.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1]
GEN = SRC / "appendices" / "generated_n200"
TEST = HERE / "test"
OUT = TEST / "out"
sys.path.insert(0, str(HERE))
import table_colour  # noqa: E402  (colour switch: plain variant of a fragment)
COLOUR = ""   # "\\def\\releasecolortables{all}" with --colour (set in main)

CAP = (r"\caption{Stand-in caption of two lines for the layout test: the writers own the captions; "
       r"this text only reserves the vertical space that a caption of two lines takes on the first page.}")

SEC6 = ["sec6_short", "sec6_combined", "sec6_maze_arxiv", "main_maze_arxiv_v7"]   # sec6_combined: one table in place of sec6_default and sec6_short_star; sec6_maze_arxiv: arXiv Table I (Maze); main_maze_arxiv_v7: arXiv Table I as a table* (one header line with units does not fit 245.7 pt)
SEC6_COL = ["main_maze_extended_v8", "main_corridor_extended_v8", "main_maze_arxiv", "main_maze_arxiv_v2", "main_maze_arxiv_v3", "main_maze_arxiv_v4", "main_maze_arxiv_v5", "main_maze_arxiv_v6", "main_maze_arxiv_v8", "main_corridor_arxiv_v8"]   # main_corridor_arxiv_v8: arXiv Table IV (Corridor) in one column as Table I; main_maze_arxiv_v8: arXiv Table I in one column; each a one-column table float at the arXiv \\columnwidth (245.7 pt)
NO_TABLE = ("assumption_checks_arxiv",)   # macro files among the paper fragments
ARXIV_PAIR = ["checks_corridor_arxiv", "checks_maze_arxiv"]   # the pair of arXiv Table VIII (arxiv_gallery.tex)
FOOTNOTE = ("topfive_full_corridor_onepage", "topfive_full_maze_onepage")   # set in \footnotesize (arxiv_restored_figures.tex)


def paper_doc(names):
    body = [r"\typeout{TEST-TEXTWIDTH=\the\textwidth}\typeout{TEST-TEXTHEIGHT=\the\textheight}"]
    body.append(r"\section{Two-column text}\lipsumlike")
    for n in SEC6:
        body.append(rf"\begin{{table*}}[p]\centering\caption{{Stand-in caption ({n.replace('_', ' ')}) of two lines for "
                    r"the layout test of the Section VI tables at the text width of the paper; the writers own the "
                    r"captions.}" + rf"\typeout{{SEC6-FLOAT {n}}}\ntable{{{n}}}\end{{table*}}")
    for n in SEC6_COL:   # natural width; an Overfull box in the 245.7 pt minipage fails the test
        body.append(rf"\begin{{table}}[p]\centering\begin{{minipage}}{{245.7pt}}\centering\caption{{Stand-in caption "
                    r"of the one-column Section VI table.}" + rf"\typeout{{SEC6-FLOAT {n}}}\ntable{{{n}}}"
                    r"\end{minipage}\end{table}")
    body.append(r"\clearpage\onecolumn\raggedbottom\section{One-column appendix}")
    body.append(r"\typeout{TEST-G-LINEWIDTH=\the\linewidth}")
    body.append(r"\clearpage\typeout{FRAG-START arxiv_pair \thepage}\begin{simdisplay}{table}" + CAP +
                r"\normalsize" + r"\par\medskip".join(rf"\ntable{{{n}}}" for n in ARXIV_PAIR) +
                r"\end{simdisplay}\par\typeout{PAGEFILL \the\pagetotal\space of \the\pagegoal}"
                r"\clearpage\typeout{FRAG-END arxiv_pair \thepage}")
    for n in names:
        if n in SEC6 or n in SEC6_COL or n in ARXIV_PAIR:
            continue
        tab = rf"{{\footnotesize\ntable{{{n}}}}}" if n in FOOTNOTE else rf"\ntable{{{n}}}"
        body.append(rf"\clearpage\typeout{{FRAG-START {n} \thepage}}\begin{{simbreaktable}}{CAP}{tab}"
                    rf"\end{{simbreaktable}}\par\typeout{{PAGEFILL \the\pagetotal\space of \the\pagegoal}}\clearpage\typeout{{FRAG-END {n} \thepage}}")
    return (r"""\documentclass[letterpaper,10pt,journal,twoside]{ieeeconf}
\def\paperentryclassloaded{1}
\def\paperversion{2}
\def\releaseedition{full}
\def\releasechanges{0}
""" + COLOUR + r"""
\input{paper_helper}
\input{release_helper}
\newcommand{\lipsumlike}{Body text at the paper's 10 pt size for comparison with the table text.
The Section VI tables are table* floats at the text width; Appendix G is one column.\par}
\begin{document}
\makeatletter\let\ntablebreak@orig\ntablebreak
\renewcommand{\ntablebreak}{\par\typeout{PAGEFILL \the\pagetotal\space of \the\pagegoal}\ntablebreak@orig}
\makeatother
""" + "\n".join(body) + "\n\\end{document}\n")


def report_doc(names):
    body = [r"\typeout{TEST-TEXTWIDTH=\the\textwidth}\typeout{TEST-TEXTHEIGHT=\the\textheight}"]
    body.append(r"\section{Report}Body text at the report's 11 pt Times size for comparison with the table text.\par")
    for n in names:
        body.append(rf"\clearpage\typeout{{FRAG-START {n} \thepage}}\begin{{simbreaktable}}{CAP}\ntable{{{n}}}"
                    rf"\end{{simbreaktable}}\par\typeout{{PAGEFILL \the\pagetotal\space of \the\pagegoal}}\clearpage\typeout{{FRAG-END {n} \thepage}}")
    return (r"""\documentclass[11pt]{article}
\pdfoutput=1
\usepackage[T1]{fontenc}
\usepackage[font=small,labelfont=bf,labelsep=colon,justification=justified,singlelinecheck=false,aboveskip=4pt,belowskip=2pt]{caption}
\usepackage{subcaption}
\usepackage{booktabs,array}
\usepackage{placeins}
\usepackage{setspace}
\usepackage[letterpaper,margin=0.75in]{geometry}\usepackage{etoolbox}
\usepackage{times}
""" + COLOUR + r"""
\input{report/report_helper}
\input{report/template/compact_envs}
\begin{document}
\renewcommand{\baselinestretch}{0.97}
\makeatletter\let\ntablebreak@orig\ntablebreak
\renewcommand{\ntablebreak}{\par\typeout{PAGEFILL \the\pagetotal\space of \the\pagegoal}\ntablebreak@orig}
\makeatother
""" + "\n".join(body) + "\n\\end{document}\n")


def run(texname, names):
    OUT.mkdir(parents=True, exist_ok=True)
    rel = (TEST / texname).relative_to(SRC)
    for _ in range(2):
        p = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            f"-output-directory={OUT.relative_to(SRC)}", str(rel)],
                           cwd=SRC, capture_output=True, text=True)
    log = (OUT / texname.replace(".tex", ".log")).read_text(errors="replace")
    ok = p.returncode == 0
    if not ok:
        print(log[-3000:])
    res = dict(file=texname, ok=ok)
    res["textwidth"] = re.search(r"TEST-TEXTWIDTH=([\d.]+pt)", log).group(1)
    res["textheight"] = re.search(r"TEST-TEXTHEIGHT=([\d.]+pt)", log).group(1)
    cur, ov, un = "preamble", [], []
    for line in log.splitlines():
        m = re.match(r"(FRAG-START|SEC6-FLOAT) (\S+)", line)
        if m:
            cur = m.group(2)
        if line.startswith("Overfull \\"):
            ov.append(f"[{cur}] {line}")
        if line.startswith("Underfull \\"):
            un.append(f"[{cur}] {line}")
    res["overfull"], res["underfull"] = ov, un
    res["fill"] = re.findall(r"PAGEFILL ([\d.]+)pt of ([\d.]+)pt", log)
    # page fill per fragment; a block that ends below the page goal ran off the page
    fills, cur2 = [], None
    for line in log.splitlines():
        m = re.match(r"FRAG-START (\S+)", line)
        if m:
            cur2 = m.group(1)
        m = re.match(r"PAGEFILL ([\d.]+)pt of ([\d.]+)pt", line)
        if m:
            fills.append((cur2, float(m.group(1)), float(m.group(2))))
    res["fills"] = fills
    res["overflow"] = [(n, a, b) for n, a, b in fills if a > b + 0.5]
    res["float_too_large"] = len(re.findall(r"Float too large", log))
    res["v10warn"] = re.findall(r"Package v10tables Warning: ([^\n]*)", log)
    res["fragment_caption"] = len(re.findall(r"N200-FRAGMENT-CAPTION", log))
    res["pending"] = re.findall(r"N200-PENDING-ASSET: ([^\n]*)", log)
    seps = [(a, float(b)) for a, b in re.findall(r"VN-SEP (\S+) = ([\d.\-]+)pt", log)]
    res["min_sep"] = min(seps, key=lambda x: x[1]) if seps else None
    res["seps"] = seps
    starts = dict(re.findall(r"FRAG-START (\S+) (\d+)", log))
    ends = dict(re.findall(r"FRAG-END (\S+) (\d+)", log))
    bad = []
    for n in starts:
        pages = int(ends[n]) - int(starts[n])
        if n == "arxiv_pair":
            want = 1
        else:
            kind = "paper" if texname.startswith("test_paper") else "report"
            want = table_colour.plain_variant((GEN / kind / f"{n}.tex").read_text()).count(r"\ntablebreak") + 1
        if pages != want:
            bad.append((n, pages, want))
    res["page_mismatch"] = bad
    res["fragments"] = len(starts)
    res["pages"] = int(re.search(r"Output written on .*?\((\d+) pages", log).group(1)) if ok else None
    return res


def main():
    global COLOUR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--colour", action="store_true", help="test the coloured variants (\\releasecolortables)")
    a = ap.parse_args()
    COLOUR = "\\def\\releasecolortables{all}" if a.colour else ""
    sfx = "_colour" if a.colour else ""
    # assumption_checks_arxiv is a macro file (prints nothing), not a table
    paper = sorted(p.stem for p in (GEN / "paper").glob("*.tex") if p.stem not in NO_TABLE)
    report = sorted(p.stem for p in (GEN / "report").glob("*.tex"))
    TEST.mkdir(exist_ok=True)
    (TEST / f"test_paper{sfx}.tex").write_text(paper_doc(paper))
    (TEST / f"test_report{sfx}.tex").write_text(report_doc(report))
    fail = False
    for tex, names in ((f"test_paper{sfx}.tex", paper), (f"test_report{sfx}.tex", report)):
        r = run(tex, names)
        print(f"{tex}: ok={r['ok']} pages={r['pages']} fragments={r['fragments']} textwidth={r['textwidth']} "
              f"textheight={r['textheight']}")
        print(f"  Overfull {len(r['overfull'])}; Underfull {len(r['underfull'])}; float too large "
              f"{r['float_too_large']}; v10tables warnings {len(r['v10warn'])}; fragment captions "
              f"{r['fragment_caption']}; pending {len(r['pending'])}; min sep {r['min_sep']}")
        for x in r["overfull"][:10] + r["underfull"][:5] + r["v10warn"][:10]:
            print("   ", x)
        if r["page_mismatch"]:
            print("  page mismatch (fragment, pages, expected):", r["page_mismatch"])
        fl = {}
        for n, a, b in r["fills"]:
            fl.setdefault(n, []).append(f"{100 * a / b:.0f}%")
        print("  page fill per block:", "; ".join(f"{n} {' '.join(v)}" for n, v in fl.items()))
        if r["overflow"]:
            print("  OVERFLOW (fragment, fill, goal):", r["overflow"])
        if (not r["ok"] or r["overfull"] or r["underfull"] or r["float_too_large"] or r["v10warn"]
                or r["fragment_caption"] or r["page_mismatch"] or r["overflow"]
                or (r["min_sep"] and r["min_sep"][1] < 1)):
            fail = True
    print(f"LAYOUT CHECK{' (coloured variants)' if sfx else ''}:", "FAILED" if fail else "PASSED")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
