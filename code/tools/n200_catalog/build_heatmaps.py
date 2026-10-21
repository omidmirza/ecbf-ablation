#!/usr/bin/env python3
"""alpha-by-K heatmaps of the N200 main grid.

Data: N200_numbers/final/common/metrics.json (folder ECBF_N200_NUMBERS); --width-texpt and --text-pt set the
target width and body size (paper 505.89 pt / 10 pt, report 505.89 pt / 11 pt); the layout scales with the
text size; --rotate-delta prints the delta labels vertically to give the cells more width.

One figure = one scenario, one alpha0, one metric and a list of delta slices. Panel rows are delta
slices, panel columns are the four variants V-U, V-Delta, V-L, V-C; each panel is an alpha-by-K grid
(alpha rows 0.05.. < alpha0, K columns 10, 50, 100, 200, infinity). Every cell prints the mean with the
sample SD below it. One colour bar per figure; the colour range is fixed per metric per scenario (the
minimum and maximum over the whole main grid of that scenario, every alpha0). Red-to-blue map: blue is the better
end, red the worse end, per metric direction;
ink black or white, whichever contrasts more with the cell (at least 4.58:1 for any colour).

Vector PDF at the paper \\textwidth (505.89 TeX pt), LaTeX Computer Modern, all text at
the 10 pt body size including sub/superscripts (\\DeclareMathSizes), saved by the figure helper
ieee_figure.save() in strict mode. Read-only on the data; writes only into --out.

Default call builds the sample set (Corridor and Maze, alpha0 = 0.30, delta in {0.05, 0.10},
metrics q_unsafe and m_c). Full set example:
    python3 build_heatmaps.py --alpha0 0.10 0.15 0.20 0.25 0.30 \\
        --deltas 0.025 0.05 0.10 0.15 0.20 --metrics q_unsafe m_c q_CBF T_mis --deltas-per-figure 2
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path

if sys.platform == "darwin":   # MacTeX and Homebrew (latex, pdfinfo) when not on the PATH
    os.environ["PATH"] = "/Library/TeX/texbin:/opt/homebrew/bin:/usr/local/bin:" + os.environ.get("PATH", "")
sys.dont_write_bytecode = True

REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/<dir>/build_heatmaps.py)
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
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = NUM / "final/common/metrics.json"
POPULATION = DATA.with_name("population.jsonl")
FAMILY = "main_grid_200"
N_MC = 200

CMAP_NAME, CMAP_LO, CMAP_HI = "RdBu", 0.06, 0.94      # red (worse) to blue (better); ColorBrewer RdBu (colour-blind safe), ends trimmed
NA_FACE = "#d0d0d0"
TEXT_PT = 10.0                                       # body size of the paper, in bp (= 10.04 TeX pt)
TEX_PT_PER_BP = 72.27 / 72.0
TEXTWIDTH_TEX_PT = 505.89
W = TEXTWIDTH_TEX_PT / TEX_PT_PER_BP                 # drawn width in bp, so \includegraphics scale ~ 1

METHODS = ("V-U", "V-D", "V-Cd", "V-C")
METHOD_LABEL = {"V-U": r"V-U", "V-D": r"V-$\Delta$", "V-Cd": r"V-L", "V-C": r"V-C"}
# the paper's method colours and markers
METHOD_STYLE = {"V-U": ("#56B4E9", "o"), "V-D": ("#0072B2", "h"), "V-Cd": ("#009E73", "p"), "V-C": ("#E69F00", "^")}
K_LEVELS = (10, 50, 100, 200, "infinity")
K_LABEL = {10: "10", 50: "50", 100: "100", 200: "200", "infinity": r"$\infty$"}
ALPHAS = (0.05, 0.10, 0.15, 0.20, 0.25)

# key -> (symbol, unit kind, direction); direction drives arrows, colour sense and bold
METRICS = {
    "T_mis": (r"\widehat T_{\mathrm{mis}}", "s", "down"),
    "m_c": (r"\widehat m_c/c_{\max}", None, "down"),
    "m_Y": (r"\widehat m_Y/c_{\max}", None, "down"),
    "m_int": (r"\widehat m_{\|u-u^{\mathrm{nom}}\|^2}", None, "down"),
    "q_CBF": (r"\widehat q_{\mathrm{CBF}}^+", "fraction", "down"),
    "q_Y": (r"\widehat q_Y^+", "fraction", "down"),
    "q_use": (r"\widehat q_{\mathrm{use}}^+", "fraction", "up"),
    "q_unsafe": (r"\widehat q_{\mathrm{unsafe}}^+", "fraction", "down"),
    "q_nocert": (r"\widehat q_{\mathrm{nocert}}^+", "fraction", "down"),
    "q_notrans": (r"\widehat q_{\mathrm{nodep}}^+", "fraction", "down"),
    "tau_use": (r"\widehat m_{\tau^{\mathrm{use}}}", "s", "down"),
    "p_max": (r"\widehat m_{p^{\max}}", "fraction", "down"),
    "h_min": (r"\widehat m_{h^{\min}}", "m", "up"),
    "clearance_min": (r"\widehat m_{\psi^{\min}}", "m", "up"),
    "n_certificates": (r"\widehat m_{\mathrm{cert}}", None, None),
}
NO_BOLD = {"n_certificates", "tau_use", "p_max"}     # m_cert and the metrics averaged over a subset of trajectories

# layout in bp (paper values; configure() rescales them for another body size and width)
CELL_H = 22.0
LEFT_DELTA = 24.5          # delta column of the left margin
LEFT_ALPHA = 20.5          # alpha column of the left margin
LEFT = LEFT_DELTA + LEFT_ALPHA
SPINE = 1.2
CB_TICKS = 11.0            # colour-bar tick labels above the bar
CB_BAR = 6.0
CB_GAP = 6.0
TITLE_ROW = 13.0           # variant names, delta and alpha headers
BOTTOM = 12.0              # K tick labels
SD_OFFSET = 5.2            # half distance between the mean and SD lines
ROTATE_DELTA = False
BOLD_NARROW = False
FIG_WIDTH_IN = fs.OVERRIDE_TEXT_WIDTH_IN


def configure(width_texpt: float, text_pt: float, rotate_delta: bool) -> None:
    """Rescale the layout for a target width (TeX pt) and body size (bp)."""
    global TEXT_PT, TEXTWIDTH_TEX_PT, W, CELL_H, LEFT_DELTA, LEFT_ALPHA, LEFT, CB_TICKS, CB_BAR, CB_GAP
    global TITLE_ROW, BOTTOM, SD_OFFSET, ROTATE_DELTA, FIG_WIDTH_IN
    k = text_pt / 10.0
    TEXT_PT, TEXTWIDTH_TEX_PT = text_pt, width_texpt
    W = width_texpt / TEX_PT_PER_BP
    CELL_H, LEFT_ALPHA = 22.0 * k, 20.5 * k
    LEFT_DELTA = 13.0 * k if rotate_delta else 24.5 * k
    LEFT = LEFT_DELTA + LEFT_ALPHA
    CB_TICKS, CB_BAR, CB_GAP, TITLE_ROW, BOTTOM = 11.0 * k, 6.0 * k, 6.0 * k, 13.0 * k, 12.0 * k
    SD_OFFSET = 5.2 * k
    ROTATE_DELTA = rotate_delta
    FIG_WIDTH_IN = width_texpt / 72.0
    if not any(abs(FIG_WIDTH_IN - w) < 1e-6 for w in fs.TEMPLATE_WIDTHS_IN):
        fs.TEMPLATE_WIDTHS_IN = tuple(fs.TEMPLATE_WIDTHS_IN) + (FIG_WIDTH_IN,)   # measured report width
    fs.MIN_TEXT_PT = min(text_pt, 10.95) if text_pt > 10 else fs.MIN_TEXT_PT


def rel_lum(rgb) -> float:
    """WCAG 2 relative luminance of an sRGB colour (components 0..1)."""
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb[:3]]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def ink_for(face) -> str:
    """Black or white text, whichever has the higher WCAG contrast ratio with the cell colour."""
    L = rel_lum(face)
    return "black" if (L + 0.05) / 0.05 >= 1.05 / (L + 0.05) else "white"


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_grid(path: Path) -> dict:
    data = json.loads(path.read_text())
    index = {}
    for view in data["views"]:
        definition = view["definition"]
        if definition["family"] != FAMILY:
            continue
        if definition["n"] != N_MC or len(view["cells"]) != 1:
            raise ValueError(f"unexpected view {view['view_id']}")
        cell = view["cells"][0]
        c = cell["configuration"]
        if cell["n"] != N_MC or len(set(cell["seeds"])) != N_MC:
            raise ValueError(f"view {view['view_id']} is not a complete N={N_MC} cell")
        if (c["L_mu"], c["bound_scale"], c["budget_normalizer"], c["rho_buf"]) != (0.0005, 1.0, 2.11, 0.0):
            raise ValueError(f"main-grid view {view['view_id']} off the declared fixed settings")
        key = (c["benchmark"], c["method"], round(c["alpha0"], 4), round(c["alpha"], 4), round(c["delta"], 4), c["K"])
        if key in index:
            raise ValueError(f"duplicate main-grid cell {key}")
        index[key] = (view["view_id"], cell)
    return index


def check_population(path: Path, view_ids: set) -> dict:
    """Seeds per used view in the per-trajectory map (must be exactly 200 distinct)."""
    seeds = {v: set() for v in view_ids}
    with open(path) as f:
        for line in f:
            i = line.find('"n200_view_id":"')          # population key of the N200 export
            vid = line[i + 16:i + 80] if i >= 0 else line[line.find('"view_id":"') + 11:line.find('"view_id":"') + 75]
            if vid in seeds:
                seeds[vid].add(json.loads(line)["seed"])
    bad = {v: len(s) for v, s in seeds.items() if len(s) != N_MC}
    if bad:
        raise ValueError(f"population map does not hold {N_MC} seeds for {len(bad)} views, e.g. {list(bad.items())[:3]}")
    return {v: N_MC for v in seeds}


def usable(stat) -> bool:
    return bool(stat["applicable"]) and stat["observed_n"] > 0 and stat["mean"] is not None


def fmt(x: float) -> str:
    """Compact number: 3 significant digits, at most 2 decimals, at most 4 characters below 1000."""
    if x == 0:
        return "0"
    d = max(0, min(2, 2 - math.floor(math.log10(abs(x)))))
    s = f"{x:.{d}f}"
    while len(s.lstrip("-")) > 4 and d > 0:
        d -= 1
        s = f"{x:.{d}f}"
    return s


def choose_scale(metric: str, values: list) -> int:
    """Display scale per metric per scenario: the smallest nonzero mean keeps two significant digits."""
    unit = METRICS[metric][1]
    nonzero = [abs(v) for v in values if v]
    candidates = {"fraction": (100, 10**3, 10**4, 10**5, 10**6), "m": (1, 100, 1000)}.get(unit, (1, 10, 100, 1000, 10**4))
    if not nonzero:
        return candidates[0]
    for scale in candidates:
        if min(nonzero) * scale >= 0.0995 and max(nonzero) * scale < 999.5:
            return scale
    raise ValueError(f"no compact display scale for {metric}: range {min(nonzero)}..{max(nonzero)}")


def unit_text(metric: str, scale: int) -> str:
    unit = METRICS[metric][1]
    if unit == "fraction":
        return r"(\%)" if scale == 100 else rf"($10^{{-{round(math.log10(scale))}}}$)"
    if unit == "m":
        return {1: "(m)", 100: "(cm)", 1000: "(mm)"}[scale]
    if unit == "s":
        return "(s)" if scale == 1 else rf"($10^{{-{round(math.log10(scale))}}}$\,s)"
    return "" if scale == 1 else rf"($10^{{-{round(math.log10(scale))}}}$)"


def heading(metric: str, scale: int) -> str:
    symbol, _, direction = METRICS[metric]
    arrow = {"down": r"\,\downarrow", "up": r"\,\uparrow", None: ""}[direction]
    u = unit_text(metric, scale)
    return "$" + symbol + arrow + "$" + (" " + u if u else "")


def scenario_ranges(index: dict, scenario: str, metric: str):
    """Fixed colour range (display units) and scale over the scenario's whole main grid."""
    stats = [cell["statistics"][metric] for (b, *_), (_, cell) in index.items() if b == scenario]
    means = [s["mean"] for s in stats if usable(s)]
    if not means:
        raise ValueError(f"{metric} has no applicable values in {scenario}")
    scale = choose_scale(metric, means)
    lo, hi = min(means) * scale, max(means) * scale
    if hi <= lo:
        hi = lo + 1.0
    return scale, lo, hi


def build_figure(index, scenario, alpha0, deltas, metric, out_dir: Path, ranges, dtag_override=None):
    alphas = [a for a in ALPHAS if a < alpha0 - 1e-9]
    if not alphas:
        raise ValueError(f"no alpha below alpha0={alpha0}")
    direction = METRICS[metric][2]
    scale, lo, hi = ranges
    na, nk, nd = len(alphas), len(K_LEVELS), len(deltas)
    panel_w = (W - LEFT - SPINE / 2) / len(METHODS)
    cell_w = panel_w / nk
    panel_h = na * CELL_H
    top = CB_TICKS + CB_BAR + CB_GAP + TITLE_ROW
    H = top + nd * panel_h + BOTTOM

    # collect cells
    grid = {}
    for di, delta in enumerate(deltas):
        for ai, alpha in enumerate(alphas):
            for ki, K in enumerate(K_LEVELS):
                row = []
                for method in METHODS:
                    key = (scenario, method, round(alpha0, 4), round(alpha, 4), round(delta, 4), K)
                    if key not in index:
                        raise KeyError(f"missing main-grid cell {key}")
                    view_id, cell = index[key]
                    stat = cell["statistics"][metric]
                    ok = usable(stat)
                    mean = stat["mean"] * scale if ok else None
                    sd = stat["sd"] * scale if ok and stat["sd"] is not None else None
                    if ok and stat["mean"] != 0 and float(fmt(mean)) == 0:
                        raise ValueError(f"nonzero mean prints as zero: {key}")
                    row.append(dict(view_id=view_id, configuration_id=cell["configuration_id"],
                                    population_sha256=cell["population_sha256"], mean_raw=stat["mean"] if ok else None,
                                    sd_raw=stat["sd"] if ok else None, mean=mean, sd=sd, usable=ok,
                                    applicable=stat["applicable"], observed_n=stat["observed_n"], total_n=stat["total_n"],
                                    conditional=stat["conditional"],
                                    mean_text=fmt(mean) if ok else "N-A",
                                    sd_text=(fmt(sd) if sd is not None else "") if ok else ""))
                # bold: best displayed mean among the variants, ties all bold
                bold = [False] * len(METHODS)
                if metric not in NO_BOLD and direction is not None:
                    shown = [float(e["mean_text"]) if e["usable"] else None for e in row]
                    avail = [v for v in shown if v is not None]
                    if avail:
                        best = min(avail) if direction == "down" else max(avail)
                        bold = [v is not None and v == best for v in shown]
                for e, b in zip(row, bold):
                    e["bold"] = b
                for mi, e in enumerate(row):
                    grid[(di, ai, ki, mi)] = e

    fig, ax0 = fs.figure(width=FIG_WIDTH_IN, height=H / 72.0, layout="none")
    ax0.remove()
    plt.rcParams["text.latex.preamble"] += rf"\DeclareMathSizes{{{TEXT_PT:g}}}{{{TEXT_PT:g}}}{{{TEXT_PT:g}}}{{{TEXT_PT:g}}}"
    _base = matplotlib.colormaps[CMAP_NAME]
    cmap = LinearSegmentedColormap.from_list("RdBu_trim", _base(np.linspace(CMAP_LO, CMAP_HI, 256)))
    if direction == "down":   # RdBu runs red (low) to blue (high); blue = better: low values of "down" metrics
        cmap = cmap.reversed()
    norm = Normalize(lo, hi)
    fw, fh = FIG_WIDTH_IN * 72, H

    def frac(x0, y0, w, h):  # bp from the top-left corner -> figure fraction
        return [x0 / fw, 1 - (y0 + h) / fh, w / fw, h / fh]

    texts = []          # (artist, bounds in bp: x0, y0top, x1, y1top) for the fit check

    def ftext(x, y, s, bounds=None, **kw):
        t = fig.text(x / fw, 1 - y / fh, s, fontsize=TEXT_PT, **kw)
        texts.append((t, bounds))
        return t

    # colour bar row: heading on the left, bar to the right
    head = ftext(0, CB_TICKS + CB_BAR, heading(metric, scale), ha="left", va="baseline")
    cb_x0 = LEFT + panel_w * 1.35
    cb_w = W - cb_x0 - 12
    cax = fig.add_axes(frac(cb_x0, CB_TICKS, cb_w, CB_BAR))
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.solids.set_rasterized(False)
    cb.outline.set_linewidth(0.8)
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.tick_params(axis="x", labelsize=TEXT_PT, length=2.5, width=0.8, pad=1.0, direction="out")
    cax.grid(False)
    ticks = cb.get_ticks()
    ticks = [t for t in ticks if lo - 1e-9 <= t <= hi + 1e-9]
    cb.set_ticks(ticks, labels=[f"{t:g}" for t in ticks])

    # left-margin headers
    y_title = CB_TICKS + CB_BAR + CB_GAP + TITLE_ROW / 2
    ftext(LEFT_DELTA / 2, y_title, r"$\delta$", ha="center", va="center")
    ftext(LEFT_DELTA + LEFT_ALPHA / 2, y_title, r"$\alpha$", ha="center", va="center")
    ftext(LEFT_DELTA + LEFT_ALPHA / 2, top + nd * panel_h + BOTTOM / 2 + 0.5, r"$K$", ha="center", va="center")

    for mi, method in enumerate(METHODS):
        x0 = LEFT + mi * panel_w
        colour, marker = METHOD_STYLE[method]
        t = ftext(x0 + panel_w / 2 + 4, y_title, METHOD_LABEL[method], ha="center", va="center")
        fig.add_artist(Line2D([(x0 + panel_w / 2 - 20) / fw], [1 - y_title / fh], marker=marker, markersize=4.5,
                              markerfacecolor=colour, markeredgecolor="black", markeredgewidth=0.5,
                              linestyle="none", transform=fig.transFigure))

    for di, delta in enumerate(deltas):
        y_panel = top + di * panel_h
        ftext(LEFT_DELTA / 2, y_panel + panel_h / 2, f"{delta:.2f}" if round(delta, 2) == delta else f"{delta:.3f}",
              ha="center", va="center", rotation=90 if ROTATE_DELTA else 0)
        for ai, alpha in enumerate(alphas):
            ftext(LEFT_DELTA + LEFT_ALPHA / 2, y_panel + (ai + 0.5) * CELL_H, f"{alpha:.2f}", ha="center", va="center")
        for mi, method in enumerate(METHODS):
            x0 = LEFT + mi * panel_w
            ax = fig.add_axes(frac(x0, y_panel, panel_w, panel_h))
            ax.set_xlim(0, nk)
            ax.set_ylim(na, 0)
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(r"$K$")
            ax.set_ylabel(r"$\alpha$")
            ax.xaxis.label.set_visible(False)
            ax.yaxis.label.set_visible(False)
            for s in ax.spines.values():
                s.set_linewidth(SPINE)
            for ai in range(na):
                for ki in range(nk):
                    e = grid[(di, ai, ki, mi)]
                    face = cmap(norm(e["mean"])) if e["usable"] else matplotlib.colors.to_rgba(NA_FACE)
                    ax.add_patch(matplotlib.patches.Rectangle((ki, ai), 1, 1, facecolor=face, edgecolor="white" if e["usable"] else "#9a9a9a",
                                                              linewidth=0.5,
                                                              hatch=None if e["usable"] else "////"))
                    ink = ink_for(face)
                    cx = x0 + (ki + 0.5) * cell_w
                    cy = y_panel + (ai + 0.5) * CELL_H
                    bounds = (x0 + ki * cell_w, y_panel + ai * CELL_H, x0 + (ki + 1) * cell_w, y_panel + (ai + 1) * CELL_H)
                    m = e["mean_text"]
                    if e["bold"]:   # cmb (bold, regular widths) when the cells are narrow, else cmbx
                        m = (r"{\fontseries{b}\selectfont " + m + "}") if BOLD_NARROW else (r"\textbf{" + m + "}")
                    if e["usable"]:
                        ftext(cx, cy - SD_OFFSET, m, bounds, ha="center", va="center", color=ink)
                        ftext(cx, cy + SD_OFFSET, e["sd_text"], bounds, ha="center", va="center", color=ink)
                    else:
                        ftext(cx, cy, "N-A", bounds, ha="center", va="center", color="black")
            if di == nd - 1:
                for ki, K in enumerate(K_LEVELS):
                    ftext(x0 + (ki + 0.5) * cell_w, top + nd * panel_h + BOTTOM / 2 + 0.5, K_LABEL[K],
                          ha="center", va="center")

    # fit and overlap checks (display units -> bp)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    px_per_bp = fig.dpi / 72.0
    boxes = []
    for t, bounds in texts:
        if not t.get_text().strip():
            continue
        bb = t.get_window_extent(r)
        x0, x1 = bb.x0 / px_per_bp, bb.x1 / px_per_bp
        y0, y1 = fh - bb.y1 / px_per_bp, fh - bb.y0 / px_per_bp
        if x0 < -0.01 or x1 > W + 0.01:
            raise ValueError(f"text {t.get_text()!r} leaves the template width ({x0:.2f}..{x1:.2f} bp)")
        if bounds is not None:
            bx0, by0, bx1, by1 = bounds
            if x0 < bx0 + 0.6 or x1 > bx1 - 0.6 or y0 < by0 + 0.4 or y1 > by1 - 0.4:
                raise ValueError(f"text {t.get_text()!r} does not fit its cell at {TEXT_PT} pt: text {x0:.2f},{y0:.2f},{x1:.2f},{y1:.2f} cell {bx0:.2f},{by0:.2f},{bx1:.2f},{by1:.2f}")
        boxes.append((t.get_text(), x0, y0, x1, y1))
    for lab in cax.get_xticklabels():
        if lab.get_text().strip():
            bb = lab.get_window_extent(r)
            boxes.append((lab.get_text(), bb.x0 / px_per_bp, fh - bb.y1 / px_per_bp, bb.x1 / px_per_bp, fh - bb.y0 / px_per_bp))
            if bb.x1 / px_per_bp > W + 0.01:
                raise ValueError("colour-bar tick label leaves the template width")
    hb = head.get_window_extent(r)
    if hb.x1 / px_per_bp > cb_x0 - 4:
        raise ValueError("heading runs into the colour bar")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            if a[1] < b[3] and b[1] < a[3] and a[2] < b[4] and b[2] < a[4]:
                raise ValueError(f"texts overlap: {a[0]!r} and {b[0]!r}")

    tag = lambda x: f"{x:.3f}".rstrip("0").replace(".", "p")
    dtag = dtag_override or "-".join(tag(d) for d in deltas)
    stem = out_dir / f"heatmap_{scenario}_{metric}_a0_{tag(alpha0)}_d_{dtag}"
    pdf = fs.save(fig, stem, png=False, strict=True)
    plt.close(fig)
    min_text = min(t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip())
    cells = []
    for (di, ai, ki, mi), e in sorted(grid.items()):
        cells.append(dict(variant=METHODS[mi], variant_label=METHOD_LABEL[METHODS[mi]], delta=deltas[di],
                          alpha=alphas[ai], K=K_LEVELS[ki], statistic=f"cells[0].statistics.{metric}",
                          **{k: e[k] for k in ("view_id", "configuration_id", "population_sha256", "mean_raw", "sd_raw",
                                               "applicable", "observed_n", "total_n", "conditional", "mean_text",
                                               "sd_text", "bold")}))
    return dict(pdf=recorded_path(pdf), pdf_sha256=sha_file(pdf), scenario=scenario, alpha0=alpha0, deltas=deltas, metric=metric,
                heading=heading(metric, scale), display_scale=scale, colour_limits_display=[lo, hi],
                colour_good_end="blue", direction=direction, height_bp=H, width_bp=W,
                min_text_size_bp=min_text, cells=cells)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=DATA)
    ap.add_argument("--population", type=Path, default=POPULATION)
    ap.add_argument("--skip-population-check", action="store_true")
    ap.add_argument("--out", type=Path, default=HERE)
    ap.add_argument("--scenarios", nargs="+", default=["corridor", "maze"])
    ap.add_argument("--alpha0", nargs="+", type=float, default=[0.30])
    ap.add_argument("--deltas", nargs="+", type=float, default=[0.05, 0.10])
    ap.add_argument("--deltas-per-figure", type=int, default=2, help="2 = 8 panels (half page), 4 = 16 panels (page)")
    ap.add_argument("--metrics", nargs="+", default=["q_unsafe", "m_c"], choices=sorted(METRICS))
    ap.add_argument("--scale", nargs="*", default=[], help="override display scale, e.g. q_CBF=100")
    ap.add_argument("--width-texpt", type=float, default=505.89, help="target width in TeX pt")
    ap.add_argument("--text-pt", type=float, default=10.0, help="figure text size in bp (>= body size)")
    ap.add_argument("--rotate-delta", action="store_true")
    ap.add_argument("--cell-h", type=float, default=None, help="cell height in bp (default 22 x text/10)")
    ap.add_argument("--sd-offset", type=float, default=None, help="mean/SD line offset in bp (default 5.2 x text/10)")
    ap.add_argument("--bold-narrow", action="store_true", help="bold best means in cmb (same digit widths as cmr)")
    ap.add_argument("--tag", default=None, help="delta tag in the file name (default: the delta values)")
    args = ap.parse_args(argv)
    configure(args.width_texpt, args.text_pt, args.rotate_delta)
    global BOLD_NARROW
    BOLD_NARROW = args.bold_narrow
    global CELL_H, SD_OFFSET
    if args.cell_h:
        CELL_H = args.cell_h
    if args.sd_offset:
        SD_OFFSET = args.sd_offset
    overrides = {k: int(v) for k, v in (s.split("=") for s in args.scale)}

    index = load_grid(args.data)
    out_dir = args.out.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    figures = []
    used = set()
    for scenario in args.scenarios:
        for metric in args.metrics:
            scale, lo, hi = scenario_ranges(index, scenario, metric)
            if metric in overrides:
                f = overrides[metric] / scale
                scale, lo, hi = overrides[metric], lo * f, hi * f
            for alpha0 in args.alpha0:
                for i in range(0, len(args.deltas), args.deltas_per_figure):
                    fig = build_figure(index, scenario, alpha0, args.deltas[i:i + args.deltas_per_figure], metric,
                                       out_dir, (scale, lo, hi), args.tag)
                    used.update(c["view_id"] for c in fig["cells"])
                    figures.append(fig)
                    print(fig["pdf"], flush=True)
    population = None
    if not args.skip_population_check:
        check_population(args.population, used)
        population = dict(path=recorded_path(args.population), views_checked=len(used), seeds_per_view=N_MC)
    manifest = dict(generator=recorded_path(Path(__file__)), generator_sha256=sha_file(Path(__file__).resolve()),
                    figure_helper_sha256=sha_file(FIGURE_STYLE / "ieee_figure.py"),
                    data=recorded_path(args.data), data_sha256=sha_file(args.data), family=FAMILY, n_mc=N_MC,
                    width_texpt=TEXTWIDTH_TEX_PT, text_pt=TEXT_PT, rotate_delta=ROTATE_DELTA, cell_h_bp=CELL_H, bold_font="cmb" if BOLD_NARROW else "cmbx",
                    fixed_settings=dict(L_mu=0.0005, b=2.11, phi=1.0, rho_buf=0.0),
                    cell_text="mean (upper line, bold = best displayed mean among V-U, V-Delta, V-L, V-C at the same "
                              "(alpha, delta, K) by the arrow direction, ties all bold) and sample SD (lower line), "
                              "in display units = raw x display_scale; N-A when applicable=false or observed_n=0",
                    colour_range="min..max of the means over the scenario's whole main grid (all alpha0), per metric",
                    population_check=population, figures=figures)
    (out_dir / "heatmap_values.json").write_text(json.dumps(manifest, indent=1) + "\n")


if __name__ == "__main__":
    main()
