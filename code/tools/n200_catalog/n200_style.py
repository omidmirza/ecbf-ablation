"""Shared method styles, widths and layout helpers of the N200 data figures.

Colours: V-U #56B4E9, V-Delta #0072B2, V-L #009E73, V-C #E69F00, ACI #CC79A7, QT #F0E442, Unfiltered #D55E00,
c=0 #7F7F7F, c_max #B0B0B0. Markers and line styles as in the paper's method legend. Figures are saved through the
figure helper ieee_figure.py (strict save: vector PDF, zero margin, text at least the body size, line widths).
Widths: paper \\textwidth 505.89 pt at the 10 pt body; report \\textwidth 505.89 pt (geometry letterpaper, margin
0.75in) at the 11 pt body. The figure width is set in TeX pt (width / 72.27 in) and the saved
PDF must span it within 1 TeX pt (\\nfigure includes N200 assets at natural size).
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

if sys.platform == "darwin":   # MacTeX and Homebrew (latex, pdfinfo) when not on the PATH
    os.environ["PATH"] = "/Library/TeX/texbin:/opt/homebrew/bin:/usr/local/bin:" + os.environ.get("PATH", "")
sys.dont_write_bytecode = True
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/<dir>/n200_style.py)
# The figure helper ieee_figure.py and the N200 summary exports are not part of this repository.
FIGURE_STYLE = Path(os.environ.get("ECBF_FIGURE_STYLE", REPO / "figure_style"))
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))
sys.path.insert(0, str(FIGURE_STYLE))


def recorded_path(p):
    """A path as recorded in manifests and caches: '<default folder>/<rest>' under an input folder named by an ECBF_*
    variable, else relative to the repository root, else the file name (no absolute local paths)."""
    p = Path(p).resolve()
    for var, name in (("ECBF_N200_NUMBERS", "N200_numbers"), ("ECBF_N200_DATA", "N200_data"),
                      ("ECBF_FIGURE_STYLE", "figure_style"), ("ECBF_DISPLAYS", "data")):
        if os.environ.get(var):
            try:
                return f"{name}/{p.relative_to(Path(os.environ[var]).resolve()).as_posix()}"
            except ValueError:
                pass
    try:
        return p.relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return p.name


import ieee_figure as fs  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TEX_PT_PER_BP = 72.27 / 72.0
KINDS = {"paper": (505.89, 10.0), "report": (505.89, 11.0)}

CERTIFIED = ("V-U", "V-D", "V-Cd", "V-C")
BASELINES = ("Zero margin", "Conservative", "ACI", "QT")            # Unfiltered control omitted from plots
LABEL = {"Zero margin": r"$c{=}0$", "Conservative": r"$c_{\max}$", "ACI": "ACI", "QT": "QT", "V-U": "V-U",
         "V-D": r"V-$\Delta$", "V-Cd": "V-L", "V-C": "V-C", "Unfiltered": "Unfiltered"}
#           colour,    marker, line style
STYLE = {"Unfiltered": ("#D55E00", "X", ":"),
         "Zero margin": ("#7F7F7F", "v", "--"),
         "Conservative": ("#B0B0B0", "s", "-."),
         "ACI": ("#CC79A7", "D", "-."),
         "QT": ("#F0E442", "P", ":"),
         "V-U": ("#56B4E9", "o", "--"),
         "V-D": ("#0072B2", "h", "--"),
         "V-Cd": ("#009E73", "p", "-."),
         "V-C": ("#E69F00", "^", "-")}
REF_COLOUR = "#333333"
MS, MEW = 3.0, 0.5            # markers on curves about 3 pt, edges 0.5 pt
LW_DATA, LW_REF, LW_BOX = 2.0, 1.8, 1.2
SETTINGS = {  # label: (alpha0, alpha, delta, K, paired variant code, paper name)
    "Default": (0.15, 0.10, 0.05, 50, "V-C"),
    "New1": (0.30, 0.25, 0.20, 50, "V-C"),
    "New2": (0.30, 0.25, 0.15, 50, "V-C"),
    "New3": (0.30, 0.20, 0.20, 50, "V-Cd"),
    "New4": (0.30, 0.25, 0.025, 100, "V-D"),
    "New5": (0.30, 0.25, 0.10, 10, "V-C"),
}
TLABEL = {"Default": "Default", "New1": "S1", "New2": "S2", "New3": "S3", "New4": "S4", "New5": "S5"}
SCEN = {"corridor": "Corridor", "maze": "Maze"}


def edge(colour):
    """Marker edge: black on light fills (yellow, light grey), white on dark fills (figs._edge)."""
    r, g, b = matplotlib.colors.to_rgb(colour)
    return "black" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.8 else "white"


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def setup(kind: str):
    """Register the width, set every text size to the body size; returns (width_in, width_texpt, txt)."""
    width_texpt, txt = KINDS[kind]
    width_in = width_texpt / 72.27
    if not any(abs(width_in - w) < 1e-6 for w in fs.TEMPLATE_WIDTHS_IN):
        fs.TEMPLATE_WIDTHS_IN = tuple(fs.TEMPLATE_WIDTHS_IN) + (width_in,)
    fs.MIN_TEXT_PT = 10.95 if txt > 10 else 10.0
    fs.setup()
    for key in ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize",
                "legend.fontsize", "legend.title_fontsize"):
        plt.rcParams[key] = txt
    plt.rcParams["figure.constrained_layout.use"] = False
    plt.rcParams["text.latex.preamble"] = (r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{bm}"
                                           rf"\DeclareMathSizes{{{txt:g}}}{{{txt:g}}}{{{txt:g}}}{{{txt:g}}}")
    return width_in, width_texpt, txt


def fill_box(fig, left=0.08, right=0.99, bottom=0.08, top=0.99, iters=6):
    """Adjust the subplot box so that the tight bounding box spans exactly the figure width and height."""
    W, H = fig.get_size_inches()
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    for _ in range(iters):
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        gl, gr = bb.x0 / W, (W - bb.x1) / W
        gb, gt = bb.y0 / H, (H - bb.y1) / H
        if max(abs(gl), abs(gr), abs(gb), abs(gt)) < 2e-4:
            break
        sp = fig.subplotpars
        fig.subplots_adjust(left=sp.left - gl, right=sp.right + gr, bottom=sp.bottom - gb, top=sp.top + gt)


def check_width(pdf: Path, target_texpt: float):
    out = subprocess.run([shutil.which("pdfinfo") or "pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    w, h = [float(x) for x in out.split("Page size:")[1].split("pts")[0].split(" x ")]
    w_tex, h_tex = w * TEX_PT_PER_BP, h * TEX_PT_PER_BP
    if abs(w_tex - target_texpt) > 1.0:
        raise ValueError(f"{pdf.name}: width {w_tex:.2f} TeX pt, target {target_texpt}")
    return round(w_tex, 2), round(h_tex, 2)


def _display_points(ax, line, n_interp=24):
    xy = line.get_xydata()
    if len(xy) == 0:
        return np.zeros((0, 2))
    pts = line.get_transform().transform(np.asarray(xy, dtype=float))
    pts = pts[np.all(np.isfinite(pts), axis=1)]
    if line.get_linestyle() in ("None", "none", "") or len(pts) < 2:
        return pts
    seg = [pts[:1]]
    for a, b in zip(pts[:-1], pts[1:]):
        t = np.linspace(0, 1, n_interp)[1:, None]
        seg.append(a + t * (b - a))
    return np.vstack(seg)


def legend_clear(ax, leg, pad_px=1.5, ignore=()):
    """True when the legend box contains no point of any line of the axes (data, markers, reference lines) and
    lies inside the axes."""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    bb = leg.get_window_extent(r)
    abb = ax.get_window_extent(r)
    if bb.x0 < abb.x0 - 0.5 or bb.x1 > abb.x1 + 0.5 or bb.y0 < abb.y0 - 0.5 or bb.y1 > abb.y1 + 0.5:
        return False
    for ln in ax.get_lines():
        if ln in ignore or not ln.get_visible():
            continue
        p = _display_points(ax, ln)
        ms_px = ln.get_markersize() * fig.dpi / 72.0 / 2 if ln.get_marker() not in (None, "None", "", " ") else 0
        m = pad_px + ms_px + ln.get_linewidth() * fig.dpi / 72.0 / 2
        inside = ((p[:, 0] > bb.x0 - m) & (p[:, 0] < bb.x1 + m) & (p[:, 1] > bb.y0 - m) & (p[:, 1] < bb.y1 + m))
        if inside.any():
            return False
    for other in ax.texts:
        if other.get_visible() and other.get_text().strip():
            ob = other.get_window_extent(r)
            if ob.overlaps(bb):
                return False
    return True


LOCS = ("upper right", "upper left", "lower right", "lower left", "center right", "center left", "upper center",
        "lower center", "center")


def place_legend(candidates, handles, txt, ncols=(2, 1), ignore=(), **kw):
    """Try each axes in `candidates`, each column count and each location; keep the first legend that covers no
    line. Returns (axes, legend) or (None, None)."""
    for ax in candidates:
        for nc in ncols:
            for loc in LOCS:
                leg = ax.legend(handles=handles, loc=loc, ncol=nc, fontsize=txt, **kw)
                if legend_clear(ax, leg, ignore=ignore):
                    return ax, leg
                leg.remove()
            # finer vertical positions along the left and right edges
            for y in np.arange(0.0, 1.0001, 0.025):
                for loc, x in (("lower left", 0.0), ("lower right", 1.0), ("lower center", 0.5)):
                    leg = ax.legend(handles=handles, loc=loc, bbox_to_anchor=(x, y), ncol=nc, fontsize=txt, **kw)
                    if legend_clear(ax, leg, ignore=ignore):
                        return ax, leg
                    leg.remove()
    return None, None


def text_overlaps(fig, gap_pt=2.0):
    """Pairs of visible texts (tick labels, labels, legend texts, titles) closer than gap_pt, and texts that reach
    into another panel."""
    fig.canvas.draw()
    g = gap_pt * fig.dpi / 72.0
    r = fig.canvas.get_renderer()
    items = []
    tick_all, tick_shown = set(), []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                tick_all.update({id(tk.label1), id(tk.label2)})
            lo, hi = sorted(axis.get_view_interval())
            eps = 1e-9 * max(1.0, abs(hi - lo))
            for tk, loc in zip(axis.get_major_ticks(), axis.get_majorticklocs()):
                if lo - eps <= loc <= hi + eps:
                    tick_shown += [lab for lab in (tk.label1, tk.label2) if lab.get_visible()]
    texts = [t for t in fig.findobj(matplotlib.text.Text) if id(t) not in tick_all] + tick_shown
    for t in texts:
        if t.get_visible() and t.get_text().strip():
            try:
                bb = t.get_window_extent(r)
            except Exception:
                continue
            if bb.width > 0 and bb.height > 0:
                items.append((t, bb))
    bad = []
    for t, bb in items:                      # a text may not reach into another panel
        for ax in fig.axes:
            if ax is t.axes or not ax.get_visible():
                continue
            ab = ax.get_window_extent(r)
            if bb.x0 < ab.x1 - 0.5 and ab.x0 < bb.x1 - 0.5 and bb.y0 < ab.y1 - 0.5 and ab.y0 < bb.y1 - 0.5:
                bad.append((t.get_text()[:25], "axes"))
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            ti, tj = items[i][0], items[j][0]
            if ti.axes is None or tj.axes is None or ti.axes is tj.axes:
                hit = a.x0 < b.x1 - 0.5 and b.x0 < a.x1 - 0.5 and a.y0 < b.y1 - 0.5 and b.y0 < a.y1 - 0.5
            else:                                   # texts of different panels keep gap_pt apart
                hit = a.x0 < b.x1 + g and b.x0 < a.x1 + g and a.y0 < b.y1 + g and b.y0 < a.y1 + g
            if hit:
                bad.append((ti.get_text()[:25], tj.get_text()[:25]))
    return bad
