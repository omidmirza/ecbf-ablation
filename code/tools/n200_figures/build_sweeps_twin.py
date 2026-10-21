#!/usr/bin/env python3
"""Twin-axis compact sweep figure (not printed; the arXiv paper uses sweepbox_T3_compact.pdf of build_sweeps_box.py).

sweepplot_T3_compact_twin.pdf: setting T3 = (alpha0, alpha, delta, K) = (0.30, 0.20, 0.20, 50); rows Corridor and
Maze; columns the swept alpha, delta, K, b (all other parameters at T3, phi = 1, L_mu = 0.0005). Each panel has two
y-axes: left m_c/c_max (solid lines, filled markers), right q+_CBF in % (dashed lines, open markers, symmetric
logarithmic scale). Two variants only: V-L (paired, #009E73) and V-C (#E69F00). Dotted vertical line: the value at
T3. Thin horizontal reference lines on the left axis: the c_max and QT baselines at alpha0 = 0.30
(matched_baselines_200). One legend inside one panel serves all panels. Data and style: build_sweeps.Data and
n200_style (strict helper save, width 505.89 TeX pt, text 10 pt). Height cap 260 pt.

Usage: python3 build_sweeps_twin.py --out out_paper
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs
import build_sweeps as B

from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
LABEL = "New3"
VARIANTS = ("V-Cd", "V-C")
SCENS = ("corridor", "maze")
COLS = ("alpha", "delta", "K", "b")
LEFT, RIGHT = "m_c", "q_CBF"
REF_BASE = ("Conservative", "QT")
DASH = (0, (3, 1.5))
LW_RIGHT = 1.6
HEIGHT_PT = 215.0
NAME = "sweepplot_T3_compact_twin"


def series(data, scen, col, m, row):
    pts = data.sweep(LABEL, scen, col, m)
    keep = [(x, c) for x, c in pts if B.value(c, row) is not None]
    return [x for x, _ in keep], [B.value(c, row) for _, c in keep], [c["view_id"] for _, c in keep]


def style_x(ax, col):
    B.style_axis(ax, col, LEFT, 10, LABEL)
    ax.grid(True, which="major", axis="both")


def build(out: Path):
    width_in, width_tex, txt = S.setup("paper")
    data = B.Data()
    a0 = S.SETTINGS[LABEL][0]
    fig, axes = plt.subplots(len(SCENS), len(COLS), figsize=(width_in, HEIGHT_PT / 72.27), sharex="col",
                             sharey="row", squeeze=False)
    twins = np.empty_like(axes)
    info = dict(figure=NAME, setting=LABEL, tuple=S.SETTINGS[LABEL][:4], variants=list(VARIANTS),
                rows=list(SCENS), cols=list(COLS), left=LEFT, right=RIGHT, series=[], baselines={})
    qmax = 0.0
    for i, scen in enumerate(SCENS):
        for j, col in enumerate(COLS):
            ax = axes[i, j]
            tx = ax.twinx()
            twins[i, j] = tx
            ax.axvline(B.setting_x(LABEL, col), color=S.REF_COLOUR, ls=":", lw=S.LW_REF, zorder=1)
            for m in REF_BASE:
                v = B.value(data.baseline(scen, a0, m), LEFT)
                c, _, ls = S.STYLE[m]
                ax.axhline(v, color=c, ls=ls, lw=S.LW_REF, zorder=1.5)
                info["baselines"].setdefault(scen, {})[m] = v
            for m in VARIANTS:
                c, mk, _ = S.STYLE[m]
                xs, ys, vids = series(data, scen, col, m, LEFT)
                ax.plot(xs, ys, color=c, ls="-", lw=S.LW_DATA, marker=mk, ms=S.MS, mfc=c, mec=S.edge(c),
                        mew=S.MEW, zorder=B.ZORDER[m])
                info["series"].append(dict(scen=scen, col=col, method=m, metric=LEFT, x=xs, y=ys, view_ids=vids))
                xs, ys, vids = series(data, scen, col, m, RIGHT)
                tx.plot(xs, ys, color=c, ls=DASH, lw=LW_RIGHT, marker=mk, ms=S.MS, mfc="white", mec=c, mew=0.8,
                        zorder=B.ZORDER[m])
                qmax = max([qmax] + ys)
                info["series"].append(dict(scen=scen, col=col, method=m, metric=RIGHT, x=xs, y=ys, view_ids=vids))
            style_x(ax, col)
            ax.set_ylim(-0.06, 1.1)
            ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
            tx.set_yscale("symlog", linthresh=0.1, linscale=1.0)
            tx.set_ylim(-0.02, 60)
            tx.set_yticks([0, 0.1, 1, 10], ["0", "0.1", "1", "10"])
            tx.yaxis.set_minor_locator(plt.NullLocator())
            tx.grid(False)
            ax.tick_params(labelsize=txt)
            tx.tick_params(labelsize=txt)
            # every axes carries its labels (strict save); only the outer ones are shown
            ax.set_ylabel(S.SCEN[scen] + "\n" + B.YLAB[LEFT])
            tx.set_ylabel(B.YLAB[RIGHT], rotation=270, va="bottom")
            ax.set_xlabel(B.XLAB[col])
            tx.set_xlabel(B.XLAB[col])
            tx.xaxis.label.set_visible(False)
            if j > 0:
                ax.tick_params(labelleft=False)
                ax.yaxis.label.set_visible(False)
            if j < len(COLS) - 1:
                tx.tick_params(labelright=False)
                tx.yaxis.label.set_visible(False)
            if i < len(SCENS) - 1:
                ax.xaxis.label.set_visible(False)
    if qmax >= 50:
        raise ValueError("q+_CBF above 50 %: widen the right axis")
    fig.align_ylabels(axes[:, 0])
    fig.align_ylabels(twins[:, -1])
    S.fill_box(fig, left=0.09, right=0.93, bottom=0.12, top=0.995)
    fig.subplots_adjust(wspace=0.06, hspace=0.08)
    S.fill_box(fig, left=fig.subplotpars.left, right=fig.subplotpars.right, bottom=fig.subplotpars.bottom,
               top=fig.subplotpars.top)

    handles = [Line2D([], [], color=S.STYLE[m][0], ls="", marker=S.STYLE[m][1], ms=S.MS + 1, mfc=S.STYLE[m][0],
                      mec=S.edge(S.STYLE[m][0]), mew=S.MEW, label=S.LABEL[m]) for m in VARIANTS]
    handles += [Line2D([], [], color="black", ls="-", lw=S.LW_DATA, label="left"),
                Line2D([], [], color="black", ls=DASH, lw=LW_RIGHT, label="right")]
    handles += [Line2D([], [], color=S.STYLE[m][0], ls=S.STYLE[m][2], lw=S.LW_REF, label=S.LABEL[m])
                for m in REF_BASE]
    common = dict(handlelength=1.6, handletextpad=0.3, columnspacing=0.6, borderaxespad=0.2, borderpad=0.25,
                  labelspacing=0.15, framealpha=1.0)
    order = [(axes[i, j], twins[i, j]) for i in range(len(SCENS)) for j in range(len(COLS))]
    found = None
    for ncol in (3, 2):
        for ax, tx in order:
            for loc in S.LOCS:
                leg = tx.legend(handles=handles, loc=loc, ncol=ncol, fontsize=txt, **common)
                fig.canvas.draw()
                if S.legend_clear(tx, leg) and S.legend_clear(ax, leg):
                    found = (ax, tx, loc, ncol)
                    break
                leg.remove()
            if found:
                break
        if found:
            break
    if not found:
        raise ValueError("no free place for the legend")
    bad = S.text_overlaps(fig)
    # a twin axes shares its panel with its host: texts of the host reaching into its own twin are fine
    bad = [b for b in bad if b[1] != "axes"]
    if bad:
        raise ValueError(f"overlapping texts {bad[:6]}")
    pdf = fs.save(fig, out / NAME, png=False, strict=True)
    plt.close(fig)
    w, h = S.check_width(pdf, width_tex)
    if h > 260.5:
        raise ValueError(f"height {h} pt above 260")
    k = [p[0] for p in order].index(found[0])
    info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=h, text_pt=txt,
                legend_panel=f"{SCENS[k // len(COLS)]}/{COLS[k % len(COLS)]}", legend_loc=found[2],
                legend_ncol=found[3])
    print(NAME, f"{w} x {h} TeX pt", "legend:", info["legend_panel"], found[2], "ncol", found[3], flush=True)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=HERE / "out_paper")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    info = build(a.out)
    (a.out / "sweeps_manifest_paper_compact_twin.json").write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        style_sha256=S.sha_file(HERE / "n200_style.py"), inputs={S.recorded_path(B.CELLS): S.sha_file(B.CELLS)},
        figures=[info]), indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
