#!/usr/bin/env python3
"""Theorem-check figures at S3 (they replace theorem_checks_<scen>_New3.pdf of build_default_figures.py).

Design (labels by symbol only, values kept clear of lines and edges):
  - y = percentage of the N_MC = 200 trajectories with the event (the trajectory frequency (80), in %);
  - marker at that percentage, bar to the one-sided 95 % Clopper-Pearson upper bound (81), in %;
  - the level as one dashed horizontal segment over the groups it applies to, labelled at its right end by its
    symbol only: "delta" at the setting's K = 50, "delta_1" at K = infinity ((100) in the paper); no arithmetic;
  - every label at least CLEAR_PT from any line, bar, marker and frame; every marker at least CLEAR_PT from the
    frame (y range set from the axes height in pt); no grid (grid lines at 10 % and 20 % would sit on the levels).

theorem_checks_<scen>_New3.pdf, one figure per scenario: S3 = (0.30, 0.20, 0.20, 50) V-L, b = 2.11, phi = 1, at
K in {50, inf} x L_mu in {0.0005, 0.005}. Two panels with a shared percentage axis: false certifications q_cert
(left) and deployment events q_dep (right). Marker = V-L (pentagon, #009E73) at 100 x possible / n (confirmed +
numerically possible events, as in the tables), bar to 100 x the Clopper-Pearson bound. The risk-increase checks q_A2
(direct) and q_pw (pairwise) are drawn only if some count is non-zero; at S3 every count is 0 (checked here), so the
figure omits them.
Style as build_sweeps_box*.py: text 11 pt (paper: sub/superscripts 10 pt),
marker 5 pt, bars 1.3 pt, level lines 1.8 pt, legend strip above the panels; strict helper save.

Data (read only; folder ECBF_N200_NUMBERS): N200_numbers/new_anchors/common/metrics.json, family detailed_200, role 3.

--setting Default: the same figure at the default setting (0.15, 0.10, 0.05, 50) with V-C
(triangle, #E69F00), b = 2.11, phi = 1; data N200_numbers/final/common/metrics.json, family
detailed_200, anchor current_default; theorem_checks_<scen>_Default.pdf; y ticks 0-5 % (levels 5 % and 2.47 %).
Usage: python3 build_theorem_checks_box.py --kind paper|report --out DIR [--setting S3|Default]
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs
import build_sweeps_box as X

from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
NEW_ANCHORS = S.NUM / "new_anchors/common/metrics.json"
FINAL = S.NUM / "final/common/metrics.json"
# setting -> (data file, anchor tag, method, (alpha0, alpha, delta), file label, legend text, y ticks)
SETTINGS = {
    "S3": (NEW_ANCHORS, "role3", "V-Cd", (0.3, 0.2, 0.2), "New3", r"V-L at S3", [0, 5, 10, 15, 20]),
    "Default": (FINAL, "current_default", "V-C", (0.15, 0.1, 0.05), "Default", r"V-C at the default setting",
                [0, 1, 2, 3, 4, 5]),
}
VL_C, VL_M, _ = S.STYLE["V-Cd"]
LEVEL_COLOUR = "#D55E00"
TXT = 11.0
LW = 1.3
MS = 5.0
MEW = 0.6
HEIGHT_PT = 168.0                     # two figures and their captions share one page
CLEAR_PT = 3.0                        # minimum gap between a label and any line, marker, bar or frame
LABEL_GAP_PT = 4.5                    # gap between a level segment's right end and its label
BOTTOM_PT = 7.0                       # 0 % to the bottom frame: marker radius 2.5 + edge + CLEAR_PT + frame
TOP_PT = 10.0                         # highest level to the top frame (its label is centred on it)
GROUPS = [("50", 0.0005), ("50", 0.005), ("infinity", 0.0005), ("infinity", 0.005)]
GLABEL = {g: (r"$K{=}50$" if g[0] == "50" else r"$K{=}\infty$") + "\n" + rf"$L_\mu{{=}}{g[1]:g}$" for g in GROUPS}
FAMS = (("qcert", r"False certifications $\widehat q_{\mathrm{cert}}$"),
        ("qdep", r"Deployment events $\widehat q_{\mathrm{dep}}$"))
SEGMENT = {"50": (-0.38, 1.38), "infinity": (1.62, 3.38)}     # level segment of each K block (groups 0-1, 2-3)
LEVEL_TEXT = {"50": r"$\delta$", "infinity": r"$\delta_1$"}
XMIN = -0.55
YTICKS = [0, 5, 10, 15, 20]


def rows_of(anchors, setting="S3"):
    _, tag, method, tup, _, _, _ = SETTINGS[setting]
    rows = {}
    for v in anchors["views"]:
        d = v["definition"]
        if d["family"] != "detailed_200" or tag not in (d.get("anchor_id") or ""):
            continue
        c = v["cells"][0]
        cf = c["configuration"]
        if cf["method"] != method:
            continue
        if (cf["alpha0"], cf["alpha"], cf["delta"], cf["budget_normalizer"], cf["bound_scale"]) != (*tup, 2.11, 1.0):
            raise ValueError(f"unexpected {setting} configuration {cf}")
        key = (cf["benchmark"], str(cf["K"]), cf["L_mu"])
        if key[1] not in ("50", "infinity"):   # GROUPS, GLABEL and the manifest name K = 50 (the setting's K)
            raise ValueError(f"{setting}: K = {cf['K']} is not the K = 50 of the group labels; set GROUPS from the setting")
        if c["n"] != 200 or key in rows:
            raise ValueError(f"{setting} cell {key}")
        rows[key] = dict(view_id=v["view_id"], configuration_id=c["configuration_id"], cell=c)
    return rows


def pt_px(fig, pt):
    return pt * fig.dpi / 72.0


def clear_of_lines(ax, bb, pad_pt):
    """True when the display box bb keeps pad_pt from every line point (bars, caps, markers, level segments) of ax and
    lies inside the axes with pad_pt to each frame (the frame is 1.2 pt wide: half of it counts)."""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    pad = pt_px(fig, pad_pt)
    abb = ax.get_window_extent(r)
    frame = pt_px(fig, ax.spines["left"].get_linewidth() / 2)
    if (bb.x0 < abb.x0 + frame + pad or bb.x1 > abb.x1 - frame - pad or bb.y0 < abb.y0 + frame + pad
            or bb.y1 > abb.y1 - frame - pad):
        return False
    for ln in ax.get_lines():
        if not ln.get_visible():
            continue
        p = S._display_points(ax, ln, n_interp=60)
        ms = pt_px(fig, ln.get_markersize() / 2 + ln.get_markeredgewidth() / 2) \
            if ln.get_marker() not in (None, "None", "", " ") else 0.0
        m = pad + ms + pt_px(fig, ln.get_linewidth() / 2)
        if ((p[:, 0] > bb.x0 - m) & (p[:, 0] < bb.x1 + m) & (p[:, 1] > bb.y0 - m) & (p[:, 1] < bb.y1 + m)).any():
            return False
    return True


def set_ranges(fig, axes, top_value):
    """y: 0 % sits BOTTOM_PT above the bottom frame and the highest level TOP_PT below the top frame; x: the label of
    the right-hand level segment fits with CLEAR_PT to the right frame. Uses the current axes size in pt."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = axes[0].get_window_extent(r)
    h_pt, w_pt = bb.height * 72.0 / fig.dpi, bb.width * 72.0 / fig.dpi
    span = top_value / (1.0 - (BOTTOM_PT + TOP_PT) / h_pt)
    ppu = h_pt / span
    axes[0].set_ylim(-BOTTOM_PT / ppu, top_value + TOP_PT / ppu)
    label_w = max(t.get_window_extent(r).width for ax in axes for t in ax.texts) * 72.0 / fig.dpi
    need_pt = LABEL_GAP_PT + label_w + CLEAR_PT + 1.0
    x1 = SEGMENT["infinity"][1]
    xspan = (x1 - XMIN) / (1.0 - need_pt / w_pt)
    for ax in axes:
        ax.set_xlim(XMIN, XMIN + xspan)


def build(kind, out: Path, anchors, setting="S3"):
    width_in, width_tex, body = S.setup(kind)
    for key in ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize",
                "legend.fontsize"):
        plt.rcParams[key] = TXT
    sub = 10 if body <= 10 else 11
    plt.rcParams["text.latex.preamble"] = (r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{bm}"
                                           rf"\DeclareMathSizes{{{TXT:g}}}{{{TXT:g}}}{{{sub}}}{{{sub}}}")
    rows = rows_of(anchors, setting)
    _, _, method, tup, flabel, legtext, yticks = SETTINGS[setting]
    vc, vm, _ = S.STYLE[method]
    delta, b = tup[2], 2.11
    delta1 = delta / (2 * b * math.log(2) ** 2)
    level = {"50": delta, "infinity": delta1}
    infos = []
    for scen in S.SCEN:
        # risk-increase checks: drawn only if some count is non-zero
        ri_any = sum(rows[(scen, *g)]["cell"]["risk_increase_checks"][k]["any_excess"]
                     for g in GROUPS for k in ("direct", "pairwise"))
        if ri_any:
            raise ValueError(f"{scen}: risk-increase excess counts are non-zero; add the q_A2/q_pw panel")
        fig, axes = plt.subplots(1, 2, figsize=(width_in, HEIGHT_PT / 72.27), sharey=True)
        fig.subplots_adjust(wspace=0.04)
        info = dict(figure=f"theorem_checks_{scen}_{flabel}", kind=kind,
                    setting=f"{setting} ({tup[0]:.2f},{tup[1]:.2f},{tup[2]:.2f},50) {S.LABEL[method]}", b=b,
                    phi=1.0, delta=delta, delta1=delta1, unit="percent of the N_MC = 200 trajectories", groups=[],
                    risk_increase_all_zero=True, risk_increase={}, level_labels={})
        for g in GROUPS:
            c = rows[(scen, *g)]["cell"]
            info["risk_increase"][f"K={g[0]},L_mu={g[1]}"] = {
                k: dict(any_excess=c["risk_increase_checks"][k]["any_excess"],
                        total_trajectories=c["risk_increase_checks"][k]["total_trajectories"],
                        assessed_decisions=c["risk_increase_checks"][k].get("assessed_decisions"))
                for k in ("direct", "pairwise")}
        for ai, (fam, title) in enumerate(FAMS):
            ax = axes[ai]
            for gi, g in enumerate(GROUPS):
                c = rows[(scen, *g)]["cell"]
                e = c["trajectory_events"][fam]
                n = e["applicable_n"]
                if n != 200:
                    raise ValueError(f"{scen} {g} {fam}: n = {n}")
                lev = level[g[0]]
                if abs(e["upper95_possible"] - X.cp_upper(e["possible"], n)) > 1e-9:
                    raise ValueError("Clopper-Pearson bound differs from the export")
                y, ub = 100.0 * e["possible"] / n, 100.0 * e["upper95_possible"]
                ax.plot([gi, gi], [y, ub], color=vc, lw=LW, solid_capstyle="butt", zorder=3)
                ax.plot([gi - 0.08, gi + 0.08], [ub, ub], color=vc, lw=LW, solid_capstyle="butt", zorder=3)
                ax.plot([gi], [y], ls="", marker=vm, ms=MS, mfc=vc, mec="black", mew=MEW, zorder=4)
                if ub > 100.0 * lev + 1e-9:
                    raise ValueError(f"{scen} {g} {fam}: bound above the level")
                if ai == 0:
                    info["groups"].append(dict(K=g[0], L_mu=g[1], view_id=rows[(scen, *g)]["view_id"],
                                               level=lev, level_percent=100.0 * lev))
                info["groups"][gi][fam] = dict(confirmed=e["confirmed"], possible=e["possible"], n=n,
                                               percent=y, upper95=e["upper95_possible"], upper95_percent=ub,
                                               meets_level=ub <= 100.0 * lev + 1e-12)
            for K, (x0, x1) in SEGMENT.items():
                v = 100.0 * level[K]
                ax.plot([x0, x1], [v, v], color=LEVEL_COLOUR, ls="--", lw=S.LW_REF, zorder=2)
                ax.annotate(LEVEL_TEXT[K], xy=(x1, v), xytext=(LABEL_GAP_PT, 0), textcoords="offset points",
                            ha="left", va="center", fontsize=TXT, zorder=5)
            ax.set_title(title, fontsize=TXT, pad=4)
            ax.set_xticks(range(4), [GLABEL[g] for g in GROUPS])
            ax.tick_params(axis="x", length=0, pad=3)
            ax.grid(False)
            ax.set_xlabel("Test period $K$ and drift bound $L_\\mu$")
            ax.xaxis.label.set_visible(False)       # the tick labels name K and L_mu
        axes[0].set_yticks(yticks, [f"{t:g}" for t in yticks])
        axes[0].set_ylabel(r"Trajectories (\%)")
        axes[1].set_ylabel(r"Trajectories (\%)")
        axes[1].yaxis.label.set_visible(False)
        axes[1].tick_params(labelleft=False)
        handles = [(Line2D([], [], color=vc, lw=LW), Line2D([], [], ls="", marker=vm, ms=MS, mfc=vc,
                                                            mec="black", mew=MEW))]
        labels = [legtext + r", bar to the 95\% upper confidence bound"]
        leg = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False,
                         fontsize=TXT, handlelength=1.6, handletextpad=0.4, borderaxespad=0.0, borderpad=0.1,
                         handler_map={tuple: HandlerTuple(ndivide=1, pad=0.0)})
        top_value = 100.0 * max(level.values())
        for _ in range(3):                     # the ranges depend on the axes size and the size on the tick labels
            layout(fig, axes, leg)
            set_ranges(fig, axes, top_value)
        layout(fig, axes, leg)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        for ax in axes:
            for t in ax.texts:
                if not clear_of_lines(ax, t.get_window_extent(r), CLEAR_PT):
                    raise ValueError(f"{scen}: label {t.get_text()} closer than {CLEAR_PT} pt to a line or frame")
                info["level_labels"].setdefault(ax.get_title(), []).append(t.get_text())
        ppu = axes[0].get_window_extent(r).height * 72.0 / fig.dpi / np.ptp(axes[0].get_ylim())
        low_gap = (0.0 - axes[0].get_ylim()[0]) * ppu - MS / 2 - MEW / 2 - axes[0].spines["bottom"].get_linewidth() / 2
        if low_gap < CLEAR_PT - 1e-6:
            raise ValueError(f"{scen}: markers at 0 only {low_gap:.2f} pt from the bottom frame")
        bad = S.text_overlaps(fig)
        if bad:
            raise ValueError(f"overlapping texts {bad[:6]}")
        pdf = fs.save(fig, out / info["figure"], png=False, strict=True)
        info.update(ylim=list(axes[0].get_ylim()), xlim=list(axes[0].get_xlim()), marker_frame_gap_pt=low_gap)
        plt.close(fig)
        w, h = S.check_width(pdf, width_tex)
        info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=h, text_pt=TXT)
        infos.append(info)
        print(info["figure"], kind, w, "x", h, flush=True)
    return infos


def layout(fig, axes, leg, gap_pt=4.0, iters=14):
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.18, top=0.80)
    W, H = fig.get_size_inches()
    for _ in range(iters):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        bb = fig.get_tightbbox(r)
        sp = fig.subplotpars
        dl, dr, db = -bb.x0 / W, (W - bb.x1) / W, -bb.y0 / H
        title_top = max(ax.title.get_window_extent(r).y1 for ax in axes) / (H * fig.dpi)
        dt = (leg.get_window_extent(r).y0 / (H * fig.dpi) - gap_pt / 72.0 / H) - title_top
        if max(abs(dl), abs(dr), abs(db), abs(dt)) < 1e-4:
            break
        fig.subplots_adjust(left=sp.left + dl, right=sp.right + dr, bottom=sp.bottom + db, top=sp.top + dt)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(S.KINDS), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--setting", choices=sorted(SETTINGS), default="S3")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    src = SETTINGS[a.setting][0]
    infos = build(a.kind, a.out, json.loads(src.read_text()), a.setting)
    tag = "" if a.setting == "S3" else f"_{a.setting}"
    (a.out / f"figures_manifest_{a.kind}_checks_box{tag}.json").write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        inputs={S.recorded_path(src): S.sha_file(src)}, figures=infos), indent=1) + "\n")


if __name__ == "__main__":
    main()
