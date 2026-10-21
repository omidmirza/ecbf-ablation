#!/usr/bin/env python3
"""Per-trajectory distribution figures.

distributions_<scen>_<label>.pdf (label Default, New1..New5 = T1..T5): six panels in the layout of the paper's
distribution figure: mission time T_mis,k (s); mean robustness margin
m_c,k/c_max; CBF violations q+_CBF,k (%); safe-set violations q+_unsafe,k (%); estimated peak exceedance risk
p-hat_k^max with alpha0 dotted; minimum barrier value h_k^min (m). Per method (c=0, c_max, ACI, QT at the setting's
alpha0; V-U, V-Delta, V-L, V-C at the setting): jittered per-trajectory points (left) and a box plot (right; box =
quartiles, whiskers = 1.5 IQR, no fliers), N_MC = 200 trajectories each. Unfiltered control omitted (caption).
p-hat^max is undefined where no smaller-margin decision is executed (c_max): N-A.
Data: setting_assessments.json (load_setting_assessments.py: saved assessments, hash-checked against
N200_INDEX.sqlite, means equal to the N200 export). Saved through the figure helper (strict save), with the paper's
method colours.

Forms (--form): "separate" = the figures above (one per scenario, 3 x 2 panels, centred panel
titles; --suffix _sep writes distributions_<scen>_<label>_sep.pdf); "half" = distributions_<label>_half.pdf, one
figure per setting with both scenarios, default --layout 4x3 (rows 1-2 Corridor, rows 3-4 Maze, 3 metric columns in
the order of the separate figures, bold scenario heading set vertically at the left, centred on its block; 3x4 = the two
blocks side by side), text width, at most 0.46 \\textheight tall, so that two fit one page with their captions
(24 panels per page). Centred panel titles; method names under the bottom row only (all panels share the method axis).
The paper and the report print the half form: --form half --kind paper|report.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs

from matplotlib.lines import Line2D  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402

HERE = Path(__file__).resolve().parent
CACHE = HERE / "setting_assessments.json"
METHODS = S.BASELINES + S.CERTIFIED
PANELS = (("T_mis", "Mission time", r"$T_{\mathrm{mis},k}$ (s)", 1.0),
          ("m_c", "Mean robustness margin", r"$m_{c,k}/c_{\max}$", 1.0),
          ("q_CBF", "CBF violations", r"$q^+_{\mathrm{CBF},k}$ (\%)", 100.0),
          ("q_unsafe", "Safe-set violations", r"$q^+_{\mathrm{unsafe},k}$ (\%)", 100.0),
          ("p_max", "Estimated peak exceedance risk", r"$\widehat p_k^{\max}$", 1.0),
          ("h_min", "Minimum barrier value", r"$h_k^{\min}$ (m)", 1.0))


def light(c, f=0.55):
    r, g, b = mcolors.to_rgb(c)
    return (1 - f * (1 - r), 1 - f * (1 - g), 1 - f * (1 - b))


def strip_box(ax, data, scale, txt, na_y):
    rng = np.random.default_rng(0)                  # jitter positions only
    stats = {}
    for i, m in enumerate(METHODS):
        vals = np.array([v * scale for v in data[m] if v is not None], dtype=float)
        colour, marker, _ = S.STYLE[m]
        if len(vals) == 0:
            ax.text(i, na_y, "N-A", ha="center", va="center", fontsize=txt)
            stats[m] = None
            continue
        x = i - 0.17 + rng.uniform(-0.12, 0.12, len(vals))
        ax.plot(x, vals, linestyle="none", marker=marker, markersize=S.MS, markerfacecolor=colour,
                markeredgecolor=S.edge(colour), markeredgewidth=S.MEW, alpha=0.8, label="_nolegend_", zorder=3)
        bp = ax.boxplot(vals, positions=[i + 0.23], widths=0.24, whis=1.5, showfliers=False, patch_artist=True,
                        manage_ticks=False, zorder=4)
        for part in ("whiskers", "caps", "medians"):
            for ln in bp[part]:
                ln.set_linewidth(S.LW_BOX)
                ln.set_color("black")
        for b in bp["boxes"]:
            b.set_facecolor(light(colour))
            b.set_edgecolor("black")
            b.set_linewidth(S.LW_BOX)
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        stats[m] = dict(n=int(len(vals)), mean=float(vals.mean()), median=float(med), q1=float(q1), q3=float(q3),
                        min=float(vals.min()), max=float(vals.max()))
    ax.set_xlim(-0.6, len(METHODS) - 0.4)
    ax.grid(axis="x", visible=False)
    return stats


CLEAR_PT = 3.0       # labels and legends keep this distance from data, lines and frames


def text_clear(ax, t, pad_pt=CLEAR_PT):
    """True when text t keeps pad_pt from every line point and box of ax and from the frames."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb, abb = t.get_window_extent(r), ax.get_window_extent(r)
    pad = pad_pt * fig.dpi / 72.0
    if bb.x0 < abb.x0 + pad or bb.x1 > abb.x1 - pad or bb.y0 < abb.y0 + pad or bb.y1 > abb.y1 - pad:
        return False
    for ln in ax.get_lines():
        if not ln.get_visible():
            continue
        pts = S._display_points(ax, ln)
        if len(pts) == 0:
            continue
        ms = (ln.get_markersize() / 2 + ln.get_markeredgewidth() / 2) * fig.dpi / 72.0 \
            if ln.get_marker() not in (None, "None", "", " ") else 0.0
        lw = ln.get_linewidth() * fig.dpi / 72.0 / 2 if ln.get_linestyle() not in ("None", "none", "") else 0.0
        m = pad + ms + lw                 # points without a line: marker radius only
        if ((pts[:, 0] > bb.x0 - m) & (pts[:, 0] < bb.x1 + m) & (pts[:, 1] > bb.y0 - m) & (pts[:, 1] < bb.y1 + m)).any():
            return False
    for p in ax.patches:
        if p.get_visible() and p.get_zorder() > 0 and p.get_window_extent(r).padded(pad).overlaps(bb):
            return False
    return True


def place_na(axes):
    """Each "N-A" label (a method without values in a panel) at the first height, from the middle of the panel
    outwards, where it keeps CLEAR_PT from every point, line, box and frame."""
    out = []
    for ax in axes:
        nas = [t for t in ax.texts if t.get_text() == "N-A"]
        # upright first; set along the method slot (all labels of the panel) where neighbouring labels would touch
        rots = (0, 90)
        if len(nas) > 1:
            fig = ax.figure
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            wd = max(t.get_window_extent(r).width for t in nas) + 2.0 * fig.dpi / 72.0
            xs = sorted(ax.transData.transform((t.get_position()[0], 0))[0] for t in nas)
            if min(b - a for a, b in zip(xs[:-1], xs[1:])) < wd:
                rots = (90,)
        for t in nas:
            x = t.get_position()[0]
            t.set_transform(ax.get_xaxis_transform())
            for rot, f in [(rot, f) for rot in rots for f in
                           [0.5] + [v for d in np.arange(0.04, 0.5, 0.04) for v in (0.5 + d, 0.5 - d)]]:
                t.set_rotation(rot)
                t.set_position((x, f))
                if text_clear(ax, t):
                    out.append(dict(panel=ax.get_title(), x=float(x), axes_fraction=round(float(f), 2),
                                    rotation=rot))
                    break
            else:
                raise ValueError(f"{ax.get_title()}: no clear place for N-A")
    return out


def style_panel(ax, key, a0, compact=False):
    """Axis scales, ticks and reference lines of one panel (identical in both forms)."""
    if key in ("q_CBF", "q_unsafe"):
        ax.set_yscale("symlog", linthresh=0.1, linscale=0.8)
        ax.set_ylim(-0.03, 160)
        ax.set_yticks([0, 0.1, 1, 10, 100], ["0", "0.1", "1", "10", "100"])
        ax.yaxis.set_minor_locator(plt.NullLocator())
    if key in ("m_c", "p_max"):
        ax.set_ylim(-0.06, 1.1)
        ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
    if key == "p_max":
        ax.axhline(a0, color=S.REF_COLOUR, ls=":", lw=S.LW_REF, zorder=2, label="_nolegend_")
    if key == "h_min":
        ax.axhline(0.0, color="#555555", ls="-", lw=0.8, zorder=1, label="_nolegend_")
    if compact and key in ("T_mis", "h_min"):
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=3, min_n_ticks=3))


def build(kind, out, labels, scens, suffix=""):
    width_in, width_tex, txt = S.setup(kind)
    cache = json.loads(CACHE.read_text())["cells"]
    manifest = []
    # paper 316 pt, so that two figures and their captions fit one float page (no overflow)
    height_pt = 316.0 if kind == "paper" else 360.0
    for label in labels:
        a0 = S.SETTINGS[label][0]
        for scen in scens:
            name = f"distributions_{scen}_{label}{suffix}"
            fig, axes = plt.subplots(3, 2, figsize=(width_in, height_pt / 72.27), sharex=True, squeeze=False)
            fig.subplots_adjust(wspace=0.25, hspace=0.30)
            info = dict(figure=name, setting=label, tuple=S.SETTINGS[label][:4], paired=S.SETTINGS[label][4],
                        alpha0=a0, methods=list(METHODS), n_mc=200, panels={})
            for k, (key, title, ylab, scale) in enumerate(PANELS):
                ax = axes.flat[k]
                data = {m: [r[key] for r in cache[f"{label}|{scen}|{m}"]["values"]] for m in METHODS}
                for m in METHODS:
                    if len(data[m]) != 200:
                        raise ValueError(f"{name} {m}: {len(data[m])} trajectories")
                info["panels"][key] = strip_box(ax, data, scale, txt, na_y=0.5)
                ax.set_title(title, loc="center", pad=2, fontsize=txt)
                ax.set_ylabel(ylab)
                ax.set_xlabel("Method")
                ax.set_xticks(range(len(METHODS)), [S.LABEL[m] for m in METHODS])
                ax.tick_params(axis="x", length=0, pad=3)
                if k < 4:
                    ax.xaxis.label.set_visible(False)
                style_panel(ax, key, a0)
            S.fill_box(fig, left=0.08, right=0.995, bottom=0.08, top=0.97)
            # points, whiskers and reference lines keep CLEAR_PT from the frames
            info["frame_widening_pt"] = {key: [0.0, 0.0] for key, *_ in PANELS}
            for _ in range(2):
                for k, ax in enumerate(axes.flat):
                    w = S.frame_clearance(ax, CLEAR_PT)
                    info["frame_widening_pt"][PANELS[k][0]] = [round(a + b, 2) for a, b in
                                                               zip(info["frame_widening_pt"][PANELS[k][0]], w)]
                S.fill_box(fig, left=0.08, right=0.995, bottom=0.08, top=0.97)
            info["na_label"] = place_na(axes.flat)
            h = [Line2D([], [], color=S.REF_COLOUR, ls=":", lw=S.LW_REF, label=r"$\alpha_0$")]
            ax_p = axes.flat[4]
            axl, leg = S.place_legend([ax_p], h, txt, ncols=(1,), pad_pt=CLEAR_PT, handlelength=1.8,
                                      borderaxespad=0.25, framealpha=1.0)
            if leg is None:
                raise ValueError(f"{name}: no free place for the alpha0 legend")
            bad = S.text_overlaps(fig)
            if bad:
                raise ValueError(f"{name}: overlapping texts {bad[:6]}")
            pdf = fs.save(fig, out / name, png=False, strict=True)
            plt.close(fig)
            w, hh = S.check_width(pdf, width_tex)
            info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=hh, text_pt=txt,
                        legend_loc=leg._loc)
            manifest.append(info)
            print(name, kind, f"{w} x {hh} TeX pt", flush=True)
    return manifest


# Half-page form (the form printed in the paper and the report): short centred titles, as printed; the y labels carry
# the symbols of the separate figures.
HALF_TITLES = {"T_mis": "Mission time", "m_c": "Robustness margin", "q_CBF": "CBF violations",
               "q_unsafe": "Safe-set violations", "p_max": "Peak exceedance risk", "h_min": "Minimum barrier"}
HALF_HEIGHT_PT = {"paper": 316.0, "report": 312.0}     # <= 0.46 \textheight (paper 696 pt, report 686.6 pt)
HALF_TEXTHEIGHT_PT = {"paper": 696.0, "report": 686.565}
HEAD_GAP_PT = 2.0


LAYOUTS = {   # rows, columns, panel k of scenario block b -> (row, column), cells of block b
    "4x3": (4, 3, lambda b, k: (2 * b + k // 3, k % 3)),
    "3x4": (3, 4, lambda b, k: (k // 2, 2 * b + k % 2)),
}


def _fit_half(fig, axes, heads, blocks, sep_pt=4.0, block_pt=8.0, head_gap_pt=HEAD_GAP_PT, iters=6):
    """Place the panels by hand so that the figure has no margin. Each column gets the room its y labels need
    (plus sep_pt), each row the room of its titles and method labels (plus sep_pt); a scenario heading sits centred
    over its block, head_gap_pt above the titles of its first row; blocks side by side are block_pt further apart.
    All panels have the same size. blocks[b] = (rows, columns) of scenario block b."""
    W, H = fig.get_size_inches() * fig.dpi
    pt = fig.dpi / 72.0
    nr, nc = len(axes), len(axes[0])
    side = all(len(cols) == nc for _, cols in blocks)      # row blocks: headings at the left, set upright
    row_starts = set() if side else {min(rows) for rows, _ in blocks}
    col_breaks = {min(cols) for _, cols in blocks} - {0}
    for _ in range(iters):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        box = [[ax.get_window_extent(r) for ax in row] for row in axes]
        tbb = [[ax.get_tightbbox(r) for ax in row] for row in axes]
        ttl = [[ax.title.get_window_extent(r) for ax in row] for row in axes]
        dl = [max(box[i][c].x0 - tbb[i][c].x0 for i in range(nr)) for c in range(nc)]
        dr = [max(tbb[i][c].x1 - box[i][c].x1 for i in range(nr)) for c in range(nc)]
        tl = [max(box[i][c].x0 - ttl[i][c].x0 for i in range(nr)) for c in range(nc)]
        tr = [max(ttl[i][c].x1 - box[i][c].x1 for i in range(nr)) for c in range(nc)]
        up = [max(tbb[i][c].y1 - box[i][c].y1 for c in range(nc)) for i in range(nr)]
        dn = [max(box[i][c].y0 - tbb[i][c].y0 for c in range(nc)) for i in range(nr)]
        hh = max(h.get_window_extent(r).height for h in heads)
        hw = max(h.get_window_extent(r).width for h in heads)
        head = hh + head_gap_pt * pt
        x_left = hw + head_gap_pt * pt if side else 0.0
        gaps = [max(dl[c + 1], tr[c] + tl[c + 1]) + sep_pt * pt + (block_pt * pt if c + 1 in col_breaks else 0.0)
                for c in range(nc - 1)]
        w = (W - x_left - dl[0] - dr[-1] - sum(gaps)) / nc
        vg = [max(dn[i], 0.0) + up[i + 1] + sep_pt * pt + (head if i + 1 in row_starts else 0.0)
              + (block_pt * pt if side and any(i + 1 == min(rows) for rows, _ in blocks) else 0.0)
              for i in range(nr - 1)]
        top = H - (0.0 if side else head) - up[0]
        h = (top - dn[-1] - sum(vg)) / nr
        xs, x = [], x_left + dl[0]
        for c in range(nc):
            xs.append(x)
            x += w + (gaps[c] if c < nc - 1 else 0.0)
        ys, y = [], top
        for i in range(nr):
            ys.append(y)                                   # top edge of row i
            for c in range(nc):
                axes[i][c].set_position([xs[c] / W, (y - h) / H, w / W, h / H])
            y -= h + (vg[i] if i < nr - 1 else 0.0)
        for hd, (rows, cols) in zip(heads, blocks):
            i0, i1 = min(rows), max(rows)
            if side:                                       # centred on the panels of the block
                hd.set_position((0.0, 0.5 * (ys[i0] + ys[i1] - h) / H))
                continue
            x0, x1 = xs[min(cols)], xs[max(cols)] + w
            hd.set_position((0.5 * (x0 + x1) / W, (ys[i0] + up[i0] + head_gap_pt * pt) / H))


def build_half(kind, out, labels, scens, layout="4x3", names="bottom"):
    """names: 'block' prints the method names under the last row of each scenario block, 'bottom' under the
    last row of the figure only."""
    width_in, width_tex, txt = S.setup(kind)
    cache = json.loads(CACHE.read_text())["cells"]
    manifest = []
    height_pt = HALF_HEIGHT_PT[kind]
    if height_pt > 0.46 * HALF_TEXTHEIGHT_PT[kind] + 1e-6:
        raise ValueError(f"half height {height_pt} pt above 0.46 textheight")
    nr, nc, cell = LAYOUTS[layout]
    blocks = []
    for b in range(len(scens)):
        cells = [cell(b, k) for k in range(len(PANELS))]
        blocks.append((sorted({i for i, _ in cells}), sorted({c for _, c in cells})))
    name_rows = {max(rows) for rows, _ in blocks} if names == "block" else {nr - 1}
    for label in labels:
        a0 = S.SETTINGS[label][0]
        name = f"distributions_{label}_half"
        fig = plt.figure(figsize=(width_in, height_pt / 72.27))
        axes = [[None] * nc for _ in range(nr)]
        first = None
        for i in range(nr):
            for c in range(nc):
                ax = fig.add_axes([0.05 + 0.9 * c / nc, 0.05 + 0.9 * (nr - 1 - i) / nr, 0.5 / nc, 0.5 / nr],
                                  sharex=first)
                first = first or ax
                axes[i][c] = ax
        info = dict(figure=name, setting=label, tuple=S.SETTINGS[label][:4], paired=S.SETTINGS[label][4],
                    alpha0=a0, methods=list(METHODS), n_mc=200, layout=layout, method_names=names,
                    blocks={sc: blk for sc, blk in zip(scens, blocks)}, panels={})
        heads, flat = [], []
        for b, scen in enumerate(scens):
            if len(blocks[b][1]) == nc:      # row block: heading at the left edge, centred on the block
                heads.append(fig.text(0.0, 0.5, rf"\textbf{{{S.SCEN[scen]}}}", ha="left", va="center",
                                      rotation=90, fontsize=txt))
            else:
                heads.append(fig.text(0.5, 1.0, rf"\textbf{{{S.SCEN[scen]}}}", ha="center", va="bottom",
                                      fontsize=txt))
            info["panels"][scen] = {}
            for k, (key, _title, ylab, scale) in enumerate(PANELS):
                i, c = cell(b, k)
                ax = axes[i][c]
                flat.append((scen, key, ax))
                data = {m: [r[key] for r in cache[f"{label}|{scen}|{m}"]["values"]] for m in METHODS}
                for m in METHODS:
                    if len(data[m]) != 200:
                        raise ValueError(f"{name} {scen} {m}: {len(data[m])} trajectories")
                info["panels"][scen][key] = strip_box(ax, data, scale, txt, na_y=0.5)
                ax.set_title(HALF_TITLES[key], loc="center", pad=2, fontsize=txt)
                ax.set_ylabel(ylab, labelpad=2)
                ax.set_xticks(range(len(METHODS)), [S.LABEL[m] for m in METHODS])
                ax.set_xlabel("Method", labelpad=2)
                if i in name_rows:
                    ax.tick_params(axis="x", length=0, pad=2, labelrotation=90)
                else:                     # shared method axis: names under the last row of the block only
                    ax.tick_params(axis="x", length=0, labelbottom=False)
                    ax.xaxis.label.set_visible(False)
                style_panel(ax, key, a0, compact=True)
        _fit_half(fig, axes, heads, blocks)
        info["frame_widening_pt"] = {f"{sc}|{key}": [0.0, 0.0] for sc, key, _ in flat}
        for _ in range(2):
            for sc, key, ax in flat:
                w = S.frame_clearance(ax, CLEAR_PT)
                info["frame_widening_pt"][f"{sc}|{key}"] = [
                    round(a + b, 2) for a, b in zip(info["frame_widening_pt"][f"{sc}|{key}"], w)]
            _fit_half(fig, axes, heads, blocks)
        info["symlog_labels"] = {}
        for sc, key, ax in flat:          # 0 and 0.1 too close in a short panel: keep the 0.1 tick, drop its label
            if key in ("q_CBF", "q_unsafe"):
                fig.canvas.draw()
                r = fig.canvas.get_renderer()
                ys = {lab.get_text(): lab.get_window_extent(r) for lab in ax.get_yticklabels()}
                if ys["0.1"].y0 < ys["0"].y1 + 1.0 * fig.dpi / 72.0:
                    ax.set_yticks([0, 0.1, 1, 10, 100], ["0", "", "1", "10", "100"])
                    info["symlog_labels"][f"{sc}|{key}"] = "0, 1, 10, 100"
        info["na_label"] = place_na([ax for *_, ax in flat])
        h = [Line2D([], [], color=S.REF_COLOUR, ls=":", lw=S.LW_REF, label=r"$\alpha_0$")]
        info["legend_loc"] = {}
        for sc, key, ax in flat:
            if key != "p_max":
                continue
            _, leg = S.place_legend([ax], h, txt, ncols=(1,), pad_pt=CLEAR_PT, handlelength=1.4,
                                    handletextpad=0.4, borderpad=0.3, borderaxespad=0.25, framealpha=1.0)
            info["legend_loc"][sc] = None if leg is None else leg._loc
        if all(v is None for v in info["legend_loc"].values()):
            raise ValueError(f"{name}: no free place for the alpha0 legend")
        bad = S.text_overlaps(fig)
        if bad:
            raise ValueError(f"{name}: overlapping texts {bad[:6]}")
        pdf = fs.save(fig, out / name, png=False, strict=True)
        plt.close(fig)
        w, hh = S.check_width(pdf, width_tex)
        if hh > 0.46 * HALF_TEXTHEIGHT_PT[kind] + 0.5:
            raise ValueError(f"{name}: height {hh} TeX pt above 0.46 textheight")
        info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=hh, text_pt=txt)
        manifest.append(info)
        print(name, kind, layout, f"{w} x {hh} TeX pt", flush=True)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(S.KINDS), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--labels", nargs="+", default=list(S.SETTINGS))
    ap.add_argument("--scens", nargs="+", default=list(S.SCEN))
    ap.add_argument("--form", choices=("separate", "half"), default="separate",
                    help="separate: one 3 x 2 figure per scenario; half: one figure per setting (both scenarios)")
    ap.add_argument("--layout", choices=sorted(LAYOUTS), default="4x3", help="half form: rows x columns")
    ap.add_argument("--suffix", default="",
                    help="separate form: file-name suffix, e.g. _sep (keeps the files without it)")
    ap.add_argument("--names", choices=("block", "bottom"), default="bottom",
                    help="half form: method names under each scenario block or under the last row only")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.form == "half":
        if list(a.scens) != list(S.SCEN):
            raise SystemExit("--form half draws both scenarios")
        man = build_half(a.kind, a.out, a.labels, a.scens, a.layout, a.names)
    else:
        man = build(a.kind, a.out, a.labels, a.scens, a.suffix)
    tag = (a.suffix.strip("_") + "_" if a.suffix else "") if a.form == "separate" else "half_"
    (a.out / f"distributions_manifest_{tag}{a.kind}_{'_'.join(a.labels)}.json").write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        style_sha256=S.sha_file(HERE / "n200_style.py"), inputs={S.recorded_path(CACHE): S.sha_file(CACHE)},
        figures=man), indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
